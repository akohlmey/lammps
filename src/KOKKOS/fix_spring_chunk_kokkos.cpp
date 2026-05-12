// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_spring_chunk_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "compute_chunk_atom.h"
#include "compute_com_chunk.h"
#include "error.h"
#include "memory.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr double SMALL = 1.0e-10;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSpringChunkKokkos<DeviceType>::FixSpringChunkKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixSpringChunk(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringChunkKokkos<DeviceType>::init()
{
  FixSpringChunk::init();

  if (utils::strmatch(update->integrate_style,"^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Cannot (yet) use respa with fix spring/chunk/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSpringChunkKokkos<DeviceType>::post_force(int /*vflag*/)
{
  double dx, dy, dz, r;

  // chunk/atom and com/chunk computes require CPU data

  atomKK->sync(Host, X_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  // lock chunk compute on first use

  if (com0 == nullptr) cchunk->lock(this, update->ntimestep, -1);

  // compute current centers of mass for each chunk

  ccom->compute_array();

  nchunk = cchunk->nchunk;
  int *ichunk = cchunk->ichunk;
  double *masstotal = ccom->masstotal;
  double **com = ccom->array;

  // on first call, allocate arrays and store initial COM

  if (com0 == nullptr) {
    memory->create(com0, nchunk, 3, "spring/chunk:com0");
    memory->create(fcom, nchunk, 3, "spring/chunk:fcom");

    for (int m = 0; m < nchunk; m++) {
      com0[m][0] = com[m][0];
      com0[m][1] = com[m][1];
      com0[m][2] = com[m][2];
    }
  }

  // compute fcom = force on each COM divided by masstotal

  esprings = 0.0;
  for (int m = 0; m < nchunk; m++) {
    dx = com[m][0] - com0[m][0];
    dy = com[m][1] - com0[m][1];
    dz = com[m][2] - com0[m][2];
    r = sqrt(dx*dx + dy*dy + dz*dz);
    r = MAX(r, SMALL);

    if (masstotal[m] != 0.0) {
      fcom[m][0] = k_spring*dx/r / masstotal[m];
      fcom[m][1] = k_spring*dy/r / masstotal[m];
      fcom[m][2] = k_spring*dz/r / masstotal[m];
      esprings += 0.5*k_spring*r*r;
    } else {
      fcom[m][0] = fcom[m][1] = fcom[m][2] = 0.0;
    }
  }

  // copy ichunk and fcom to device views

  int nlocal = atom->nlocal;
  l_nchunk = nchunk;

  // allocate device views if too small
  if ((int)d_ichunk.extent(0) < nlocal)
    d_ichunk = Kokkos::View<int*, DeviceType>("spring_chunk:ichunk", nlocal);
  if ((int)d_fcom.extent(0) < nchunk)
    d_fcom = Kokkos::View<double*[3], Kokkos::LayoutRight, DeviceType>(
             "spring_chunk:fcom", nchunk);

  // copy via subviews to avoid extent mismatch
  auto h_ichunk_view =
    Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(
                 ichunk, nlocal);
  auto d_ichunk_sub = Kokkos::subview(d_ichunk, Kokkos::make_pair(0, nlocal));
  Kokkos::deep_copy(d_ichunk_sub, h_ichunk_view);

  auto h_fcom_view =
    Kokkos::View<double*[3], Kokkos::LayoutRight, Kokkos::HostSpace,
                 Kokkos::MemoryTraits<Kokkos::Unmanaged>>(&fcom[0][0], nchunk);
  auto d_fcom_sub = Kokkos::subview(d_fcom, Kokkos::make_pair(0, nchunk), Kokkos::ALL);
  Kokkos::deep_copy(d_fcom_sub, h_fcom_view);

  // apply restoring forces on device

  atomKK->sync(execution_space, F_MASK | MASK_MASK | TYPE_MASK | RMASS_MASK);

  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  copymode = 1;
  if (atom->rmass) {
    rmass = atomKK->k_rmass.view<DeviceType>();
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringChunkRmass>(0,nlocal),*this);
  } else {
    mass = atomKK->k_mass.view<DeviceType>();
    type = atomKK->k_type.view<DeviceType>();
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType,TagFixSpringChunk>(0,nlocal),*this);
  }
  copymode = 0;

  atomKK->modified(execution_space, F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringChunkKokkos<DeviceType>::operator()(TagFixSpringChunk, const int &i) const
{
  const int m = d_ichunk[i] - 1;
  if (m < 0) return;
  if (m >= l_nchunk) return;
  const double massone = mass[type[i]];
  f(i,0) -= d_fcom(m,0) * massone;
  f(i,1) -= d_fcom(m,1) * massone;
  f(i,2) -= d_fcom(m,2) * massone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixSpringChunkKokkos<DeviceType>::operator()(TagFixSpringChunkRmass, const int &i) const
{
  const int m = d_ichunk[i] - 1;
  if (m < 0) return;
  if (m >= l_nchunk) return;
  const double massone = rmass[i];
  f(i,0) -= d_fcom(m,0) * massone;
  f(i,1) -= d_fcom(m,1) * massone;
  f(i,2) -= d_fcom(m,2) * massone;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixSpringChunkKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixSpringChunkKokkos<LMPHostType>;
#endif
}
