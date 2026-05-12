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

#include "fix_gle_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixGLEKokkos<DeviceType>::FixGLEKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixGLE(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGLEKokkos<DeviceType>::initial_integrate(int vflag)
{
  atomKK->sync(Host, X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);
  FixGLE::initial_integrate(vflag);
  atomKK->modified(Host, X_MASK | V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGLEKokkos<DeviceType>::final_integrate()
{
  atomKK->sync(Host, V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);
  FixGLE::final_integrate();
  atomKK->modified(Host, V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGLEKokkos<DeviceType>::initial_integrate_respa(int vflag, int ilevel, int iloop)
{
  atomKK->sync(Host, X_MASK | V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);
  FixGLE::initial_integrate_respa(vflag, ilevel, iloop);
  atomKK->modified(Host, X_MASK | V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixGLEKokkos<DeviceType>::final_integrate_respa(int ilevel, int iloop)
{
  atomKK->sync(Host, V_MASK | F_MASK | MASK_MASK | RMASS_MASK | TYPE_MASK);
  FixGLE::final_integrate_respa(ilevel, iloop);
  atomKK->modified(Host, V_MASK);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixGLEKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixGLEKokkos<LMPHostType>;
#endif
}
