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
FixStyle(rigid/small/kk,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/kk/device,FixRigidSmallKokkos<LMPDeviceType>);
FixStyle(rigid/small/kk/host,FixRigidSmallKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_RIGID_SMALL_KOKKOS_H
#define LMP_FIX_RIGID_SMALL_KOKKOS_H

#include "fix_rigid_small.h"
#include "kokkos_base.h"
#include "kokkos_type.h"
#include "comm_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixRigidSmallKokkos : public FixRigidSmall, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixRigidSmallKokkos(class LAMMPS *, int, char **);
  ~FixRigidSmallKokkos() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void setup_pre_neighbor() override;
  void pre_exchange() override;
  void pre_neighbor() override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void final_integrate() override;
  void write_restart_file(const char *) override;

  // host-fallback methods that read/modify body[] while it is resident on the
  // device; they sync the body data to/from the host as needed
  double compute_scalar() override;
  void deform(int) override;
  void zero_momentum() override;
  void zero_rotation() override;
  void *extract(const char *, int &) override;

  void grow_arrays(int) override;

  // device communication of body data (fixed stride per atom)
  int pack_forward_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &,
                               int, int *) override;
  void unpack_forward_comm_kokkos(int, int, DAT::tdual_double_1d &) override;
  int pack_reverse_comm_kokkos(int, int, DAT::tdual_double_1d &) override;
  void unpack_reverse_comm_kokkos(int, DAT::tdual_int_1d, DAT::tdual_double_1d &) override;

 protected:

  // device per-body data (a device-resident copy of the base-class body[] list)
  Kokkos::View<Body*, DeviceType> d_body;

  // true once the device copy d_body holds the canonical body data (after
  // setup / pre_neighbor); false while the host body[] is being (re)built
  bool body_resident_device = false;

  class CommKokkos *commKK;

  void set_xv_kokkos(int setx);
  void compute_forces_and_torques_kokkos();
  void copy_body_to_device();
  void copy_body_to_host();
  void sync_peratom_to_device();


  // per-atom rigid-body arrays, stored as DualViews so the host pointers in
  // the FixRigidSmall base class alias the host side of each view and the
  // existing CPU comm/exchange/sort paths keep working unchanged, while the
  // device kernels (added in later build stages) operate on the device side.

  DAT::tdual_int_1d k_bodyown;          // index of body owned by atom, else -1
  DAT::tdual_tagint_1d k_bodytag;       // ID of body the atom belongs to, else 0
  DAT::tdual_int_1d k_atom2body;        // index of body the atom is in, else -1
  DAT::tdual_imageint_1d k_xcmimage;    // internal image flags of body atoms
  DAT::tdual_double_2d_lr k_displace;   // displacement of atom in body coords

  typename AT::t_int_1d d_bodyown;
  typename AT::t_tagint_1d d_bodytag;
  typename AT::t_int_1d d_atom2body;
  typename AT::t_imageint_1d d_xcmimage;
  typename AT::t_double_2d_lr d_displace;
};

}    // namespace LAMMPS_NS

#endif
#endif
