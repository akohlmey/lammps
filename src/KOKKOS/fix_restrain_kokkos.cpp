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

#include "fix_restrain_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixRestrainKokkos<DeviceType>::FixRestrainKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixRestrain(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixRestrainKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(Host, X_MASK | F_MASK);
  FixRestrain::post_force(vflag);
  atomKK->modified(Host, F_MASK);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class FixRestrainKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixRestrainKokkos<LMPHostType>;
#endif
}
