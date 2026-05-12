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

#include "fix_press_langevin_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixPressLangevinKokkos<DeviceType>::FixPressLangevinKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixPressLangevin(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPressLangevinKokkos<DeviceType>::initial_integrate(int vflag)
{
  atomKK->sync(Host, ALL_MASK);
  FixPressLangevin::initial_integrate(vflag);
  atomKK->modified(Host, ALL_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPressLangevinKokkos<DeviceType>::post_integrate()
{
  atomKK->sync(Host, X_MASK | V_MASK | MASK_MASK);
  FixPressLangevin::post_integrate();
  atomKK->modified(Host, X_MASK | V_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPressLangevinKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(Host, F_MASK | MASK_MASK);
  FixPressLangevin::post_force(vflag);
  atomKK->modified(Host, F_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPressLangevinKokkos<DeviceType>::end_of_step()
{
  atomKK->sync(Host, ALL_MASK);
  FixPressLangevin::end_of_step();
  atomKK->modified(Host, ALL_MASK);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixPressLangevinKokkos<DeviceType>::pre_exchange()
{
  atomKK->sync(Host, X_MASK | IMAGE_MASK);
  FixPressLangevin::pre_exchange();
  atomKK->modified(Host, X_MASK);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixPressLangevinKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixPressLangevinKokkos<LMPHostType>;
#endif
}
