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
FixStyle(MINIMIZE/kk,FixMinimizeKokkos);
FixStyle(MINIMIZE/kk/device,FixMinimizeKokkos);
FixStyle(MINIMIZE/kk/host,FixMinimizeKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_MINIMIZE_KOKKOS_H
#define LMP_FIX_MINIMIZE_KOKKOS_H

#include "fix_minimize.h"
#include "kokkos_type.h"
#include "kokkos_base.h"

namespace LAMMPS_NS {

class FixMinimizeKokkos : public FixMinimize, public KokkosBase {
  friend class MinLineSearchKokkos;

 public:
  FixMinimizeKokkos(class LAMMPS *, int, char **);
  ~FixMinimizeKokkos() override;
  void init() override {}

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  void add_vector_kokkos();
  DAT::t_kkfloat_1d request_vector_kokkos(int);
  void reset_coords();

  // LMP_KOKKOS_DUALVIEW_64BIT_WORKAROUND
  DAT::ttransform_kkfloat_2d_lr_big k_vectors;
  // LMP_KOKKOS_DUALVIEW_64BIT_WORKAROUND
  DAT::t_kkfloat_2d_lr_big d_vectors;
  // LMP_KOKKOS_DUALVIEW_64BIT_WORKAROUND
  HAT::t_double_2d_lr_big h_vectors;
};

}

#endif
#endif
