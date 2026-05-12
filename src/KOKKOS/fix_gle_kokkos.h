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
FixStyle(gle/kk,FixGLEKokkos<LMPDeviceType>);
FixStyle(gle/kk/device,FixGLEKokkos<LMPDeviceType>);
FixStyle(gle/kk/host,FixGLEKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_GLE_KOKKOS_H
#define LMP_FIX_GLE_KOKKOS_H

#include "fix_gle.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixGLEKokkos : public FixGLE {
 public:
  typedef DeviceType device_type;

  FixGLEKokkos(class LAMMPS *, int, char **);
  ~FixGLEKokkos() override {}
  void initial_integrate(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;

 private:
  class AtomKokkos *atomKK;
  ExecutionSpace execution_space;
};

}    // namespace LAMMPS_NS

#endif
#endif
