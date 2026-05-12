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
FixStyle(spring/chunk/kk,FixSpringChunkKokkos<LMPDeviceType>);
FixStyle(spring/chunk/kk/device,FixSpringChunkKokkos<LMPDeviceType>);
FixStyle(spring/chunk/kk/host,FixSpringChunkKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_SPRING_CHUNK_KOKKOS_H
#define LMP_FIX_SPRING_CHUNK_KOKKOS_H

#include "fix_spring_chunk.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixSpringChunk{};
struct TagFixSpringChunkRmass{};

template<class DeviceType>
class FixSpringChunkKokkos : public FixSpringChunk {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixSpringChunkKokkos(class LAMMPS *, int, char **);
  ~FixSpringChunkKokkos() override {}
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringChunk, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringChunkRmass, const int &) const;

 private:
  class AtomKokkos *atomKK;
  ExecutionSpace execution_space;

  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;

  // device views for chunk data
  Kokkos::View<int*, DeviceType> d_ichunk;
  Kokkos::View<double*[3], Kokkos::LayoutRight, DeviceType> d_fcom;
  int l_nchunk;
};

}    // namespace LAMMPS_NS

#endif
#endif
