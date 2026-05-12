/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(aveforce/kk,FixAveForceKokkos<LMPDeviceType>);
FixStyle(aveforce/kk/device,FixAveForceKokkos<LMPDeviceType>);
FixStyle(aveforce/kk/host,FixAveForceKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_AVEFORCE_KOKKOS_H
#define LMP_FIX_AVEFORCE_KOKKOS_H

#include "fix_aveforce.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixAveForceKokkos : public FixAveForce {
 public:
  typedef DeviceType device_type;

  FixAveForceKokkos(class LAMMPS *, int, char **);
  void post_force(int) override;
  void post_force_respa(int, int, int) override;

 private:
  class AtomKokkos *atomKK;
  ExecutionSpace execution_space;
};

}    // namespace LAMMPS_NS

#endif
#endif
