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
FixStyle(gcmc,FixGCMC);
// clang-format on
#else

#ifndef LMP_FIX_GCMC_H
#define LMP_FIX_GCMC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGCMC : public Fix {
 public:
  FixGCMC(class LAMMPS *, int, char **);
  ~FixGCMC() override;
  int setmask() override;
  void init() override;
  void pre_exchange() override;
  void post_run() override;
  double compute_vector(int) override;
  double memory_usage() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  void *extract(const char *, int &) override;

 private:
  int molecule_group, molecule_group_bit;
  int molecule_group_inversebit;
  int exclusion_group, exclusion_group_bit;
  int ngcmc_type, nevery, seed;
  int ncycles, nexchanges, nmcmoves;
  double patomtrans, pmoltrans, pmolrotate, pmctot;
  int ngas;                // # of gas atoms on all procs
  int ngas_local;          // # of gas atoms on this proc
  int ngas_before;         // # of gas atoms on procs < this proc
  int exchmode;            // exchange ATOM or MOLECULE
  int movemode;            // move ATOM or MOLECULE
  class Region *region;    // gcmc region
  char *idregion;          // gcmc region id
  bool pressure_flag;      // true if user specified reservoir pressure
  bool charge_flag;        // true if user specified atomic charge
  bool full_flag;          // true if doing full system energy calculations

  int natoms_per_molecule;    // number of atoms in each inserted molecule
  int nmaxmolatoms;           // number of atoms allocated for molecule arrays

  int groupbitall;            // group bitmask for inserted atoms
  int ngroups;                // number of group-ids for inserted atoms
  char **groupstrings;        // list of group-ids for inserted atoms
  int ngrouptypes;            // number of type-based group-ids for inserted atoms
  char **grouptypestrings;    // list of type-based group-ids for inserted atoms
  int *grouptypebits;         // list of type-based group bitmasks
  int *grouptypes;            // list of type-based group types
  double ntranslation_attempts;
  double ntranslation_successes;
  double nrotation_attempts;
  double nrotation_successes;
  double ndeletion_attempts;
  double ndeletion_successes;
  double ninsertion_attempts;
  double ninsertion_successes;

  int mc_active;              // 1 during MC trials, otherwise 0

  int gcmc_nmax;
  int max_region_attempts;
  double gas_mass;
  double reservoir_temperature;
  double tfac_insert;
  double chemical_potential;
  double displace;
  double max_rotation_angle;
  double beta, zz, sigma, volume;
  double pressure, fugacity_coeff, charge;
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double region_xlo, region_xhi, region_ylo, region_yhi, region_zlo, region_zhi;
  double region_volume;
  double energy_stored;    // full energy of old/current configuration

  // fast local-energy path: evaluate a trial's energy change locally
  // (no per-trial neighbor-list rebuild / domain re-decomposition) instead
  // of recomputing the whole-system energy via energy_full().
  int fast_mode;           // 0=off, 1=fast, 2=validate (run both, compare)
  int fast_supported;      // 1 if the fast path is usable for this config
  double kspace_stored;    // cached reciprocal-space energy of current config
  double etail_stored;     // cached pair->etail of current config
  double mol_const;        // calibrated position-independent insertion constant
                           // (molecule internal bonded + intramol reciprocal)
  int mol_const_set;       // 1 once mol_const has been calibrated
  double fast_maxerr;      // max |E_fast - E_full| observed in validate mode
  bigint fast_nvalid;      // number of validate-mode comparisons
  // per-segment wall-time accumulators, reported in post_run()
  double time_full, time_molenergy, time_kspace, time_borders;

  double *sublo, *subhi;
  int *local_gas_list;
  double **cutsq;
  double **molcoords;
  double *molq;
  imageint *molimage;
  imageint imagezero;
  double overlap_cutoffsq;    // square distance cutoff for overlap
  int overlap_flag;
  int max_ngas;
  int min_ngas;

  double energy_intra;

  class Pair *pair;

  class RanPark *random_equal;
  class RanPark *random_unequal;

  class Atom *model_atom;

  class Molecule **onemols;
  int imol, nmol;
  class Fix *fixrigid, *fixshake;
  int rigidflag, shakeflag;
  char *idrigid, *idshake;
  int triclinic;    // 0 = orthog box, 1 = triclinic

  class Compute *c_pe;

  // private methods

  void options(int, char **);

  void attempt_atomic_translation();
  void attempt_atomic_deletion();
  void attempt_atomic_insertion();
  void attempt_molecule_translation();
  void attempt_molecule_rotation();
  void attempt_molecule_deletion();
  void attempt_molecule_insertion();
  void attempt_atomic_translation_full();
  void attempt_atomic_deletion_full();
  void attempt_atomic_insertion_full();
  void attempt_molecule_translation_full();
  void attempt_molecule_rotation_full();
  void attempt_molecule_deletion_full();
  void attempt_molecule_insertion_full();

  double energy(int, int, tagint, double *);
  double energy_full();
  double molecule_energy(tagint);

  // fast local-energy path helpers
  double kspace_energy_only();    // ENERGY_ONLY reciprocal recompute -> energy
  double long_range_delta();      // (kspace + tail) energy change vs the cache
  void sync_long_cache();         // refresh kspace_stored/etail_stored on accept
  void refresh_ghosts();          // rebuild ghosts without migration (no exchange)
  double validate_fast(double energy_before, double energy_after_fast);

  int pick_random_gas_atom();
  tagint pick_random_gas_molecule();
  void toggle_intramolecular(int);
  void update_gas_atoms_list();

  void grow_molecule_arrays(int);
};

}    // namespace LAMMPS_NS

#endif
#endif
