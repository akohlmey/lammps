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

#include "fix_move_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixMoveKokkos<DeviceType>::FixMoveKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixMove(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMoveKokkos<DeviceType>::initial_integrate(int vflag)
{
  atomKK->sync(Host, ALL_MASK);
  FixMove::initial_integrate(vflag);
  atomKK->modified(Host, ALL_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixMoveKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(Host, ALL_MASK);
  FixMove::final_integrate();
  atomKK->modified(Host, ALL_MASK);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixMoveKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixMoveKokkos<LMPHostType>;
#endif
}
