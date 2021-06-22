/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_specden_scalar.h"

#include "compute.h"
#include "error.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixSpecdenScalar::FixSpecdenScalar(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), data(nullptr), specden(nullptr)
{
  time_depend = 1;
}

/* ---------------------------------------------------------------------- */

int FixSpecdenScalar::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpecdenScalar::init() {
  if (modify->find_fix_by_style("dt/reset") >= 0)
    error->all(FLERR,"Cannot use fix specden/scalar with fix dt/reset",style);
}

/* ---------------------------------------------------------------------- */

void FixSpecdenScalar::end_of_step() {}
