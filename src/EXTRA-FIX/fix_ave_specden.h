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
FixStyle(ave/specden,FixAveSpecden);
// clang-format on
#else

#ifndef LMP_FIX_AVE_SPECDEN_H
#define LMP_FIX_AVE_SPECDEN_H

#include "fix.h"
#include "safe_pointers.h"

namespace LAMMPS_NS {

class FixAveSpecden : public Fix {
 public:
  FixAveSpecden(class LAMMPS *, int, char **);
  ~FixAveSpecden() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  double compute_array(int, int) override;

 private:
  struct value_t {
    int which;         // type of data: COMPUTE, FIX, VARIABLE
    int argindex;      // 1-based index if data is vector, else 0
    int iarg;          // argument index in original argument list
    std::string id;    // compute/fix/variable ID
    union {
      class Compute *c;
      class Fix *f;
      int v;
    } val;
  };
  std::vector<value_t> values;

  int nvalues, nrepeat, nfreq;
  bigint nvalid, nvalid_last;
  SafeFilePtr fp;

  int ave, startstep, overwrite;
  double prefactor;
  bigint filepos;

  int firstindex;    // index in cvalues ring of earliest time sample
  int lastindex;     // index in cvalues ring of latest time sample
  int nsample;       // number of time samples collected in ring buffer
  int noutput;       // number of DFT windows accumulated (for ave running)
  int nfreqs;        // number of positive frequency bins = nrepeat/2 + 1

  double **cvalues;      // ring buffer of time samples [nrepeat][nvalues]
  double **specden;      // accumulated spectral density [nfreqs][nvalues]
  double **save_specden; // saved PSD values for compute_array() [nfreqs][nvalues]

  void compute_specden();
  bigint nextvalid();
};

}    // namespace LAMMPS_NS

#endif
#endif
