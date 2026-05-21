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
FixStyle(ilves,FixIlves);
// clang-format on
#else

#ifndef LMP_FIX_ILVES_H
#define LMP_FIX_ILVES_H

#include "fix.h"

namespace LAMMPS_NS {

class FixIlves : public Fix {

 public:
  FixIlves(class LAMMPS *, int, char **);
  ~FixIlves() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void pre_neighbor() override;
  void post_neighbor() override;
  void post_force(int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  void reset_dt() override;
  bigint dof(int) override;

 protected:
  double tolerance;      // convergence tolerance
  int max_iter;          // max number of Newton iterations
  int output_every;      // print stats every this many steps (0 = never)
  bigint next_output;    // timestep of next output

  // settings from input command (bond types to constrain)
  int *bond_flag;           // bond_flag[i] = 1 if bond type i is constrained
  double *bond_distance;    // equilibrium distance for each bond type

  // settings from input command (angle types to constrain)
  int *angle_flag;           // angle_flag[i] = 1 if angle type i is constrained
  double *angle_distance;    // virtual end-to-end distance for each angle type
  bool has_angle;            // true if any angle types are constrained

  // atom-based arrays (allocated per atom, survive exchange)
  int *ilves_flag;    // 1 if atom participates in any ILVES constraint
  double **xshake;    // unconstrained position (working array)

  // constraint list (flat, rebuilt each neighbor list update)
  int n_constr;        // number of constraints owned by this proc
  int max_constr;      // allocated size
  int *c_atom1;        // local index of first atom (always local, < nlocal)
  int *c_atom2;        // local index of second atom (local or ghost)
  double *c_dist2;     // squared constraint distance
  double *c_lambda;    // accumulated Lagrange multiplier (reset each step)

  // timestep info (set in setup/reset_dt)
  double dtv, dtfsq;

  // local copies of commonly used atom-class pointers
  double **x, **v, **f;
  double *mass, *rmass;
  int *type;
  int nlocal;

  // statistics (per bond type)
  bigint *b_count, *b_count_all;
  double *b_ave, *b_ave_all;
  double *b_max, *b_max_all;
  double *b_min, *b_min_all;

  void add_constraint(int, int, double);
  void stats();
};

}    // namespace LAMMPS_NS

#endif
#endif
