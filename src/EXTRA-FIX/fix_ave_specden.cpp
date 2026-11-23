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

#include "fix_ave_specden.h"

#include "arg_info.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum { ONE, RUNNING };

/* ---------------------------------------------------------------------- */

FixAveSpecden::FixAveSpecden(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), cvalues(nullptr), specden(nullptr), save_specden(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix ave/specden", error);

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat = utils::inumeric(FLERR, arg[4], false, lmp);
  nfreq = utils::inumeric(FLERR, arg[5], false, lmp);

  time_depend = 1;
  global_freq = nfreq;

  // expand args if any have wildcard character "*"

  int ioffset = 6;
  int expand = 0;
  char **earg;
  char **oarg = arg;
  int *amap = nullptr;
  int nargnew = utils::expand_args(FLERR, narg - ioffset, &arg[ioffset], 0, earg, lmp, &amap);

  if (earg != &arg[ioffset]) expand = 1;
  arg = earg;

  // parse values (scalars from computes, fixes, or variables)

  int iarg = 0;
  while (iarg < nargnew) {
    ArgInfo argi(arg[iarg]);
    value_t val;

    if (expand) val.iarg = amap[iarg] + ioffset;
    else val.iarg = iarg + ioffset;

    if (argi.get_type() == ArgInfo::NONE) break;
    if ((argi.get_type() == ArgInfo::UNKNOWN) || (argi.get_dim() > 1))
      error->all(FLERR, val.iarg, "Unknown fix ave/specden data type: {}", arg[iarg]);

    val.which = argi.get_type();
    val.argindex = argi.get_index1();
    val.id = argi.get_name();
    val.val.c = nullptr;

    values.push_back(val);
    iarg++;
  }
  nvalues = values.size();

  // get argument offset if optional arguments are present

  if (iarg < nargnew) {
    for (int i = 0; i < narg; ++i) {
      if (strcmp(oarg[i], arg[iarg]) == 0)
        ioffset = i - iarg;
    }
  }

  // optional args

  ave = ONE;
  startstep = 0;
  prefactor = 1.0;
  overwrite = 0;
  char *title1 = nullptr;
  char *title2 = nullptr;
  char *title3 = nullptr;

  while (iarg < nargnew) {
    int errptr = iarg + ioffset;
    if (strcmp(arg[iarg],"ave") == 0) {
      if (iarg+2 > nargnew) utils::missing_cmd_args(FLERR, "fix ave/specden ave", error);
      if (strcmp(arg[iarg+1],"one") == 0) ave = ONE;
      else if (strcmp(arg[iarg+1],"running") == 0) ave = RUNNING;
      else error->all(FLERR, errptr+1, "Unknown fix ave/specden ave mode: {}", arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > nargnew) utils::missing_cmd_args(FLERR, "fix ave/specden start", error);
      startstep = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"prefactor") == 0) {
      if (iarg+2 > nargnew) utils::missing_cmd_args(FLERR, "fix ave/specden prefactor", error);
      prefactor = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if ((strcmp(arg[iarg],"file") == 0) || (strcmp(arg[iarg],"append") == 0)) {
      if (iarg+2 > nargnew)
        utils::missing_cmd_args(FLERR, fmt::format("fix ave/specden {}", arg[iarg]), error);
      if (comm->me == 0) {
        if (strcmp(arg[iarg],"file") == 0) fp = fopen(arg[iarg+1],"w");
        else fp = fopen(arg[iarg+1],"a");
        if (fp == nullptr)
          error->one(FLERR, errptr+1, "Cannot open fix ave/specden file {}: {}",
                     arg[iarg+1], utils::getsyserror());
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"overwrite") == 0) {
      overwrite = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"title1") == 0) {
      if (iarg+2 > nargnew) utils::missing_cmd_args(FLERR, "fix ave/specden title1", error);
      delete[] title1;
      title1 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title2") == 0) {
      if (iarg+2 > nargnew) utils::missing_cmd_args(FLERR, "fix ave/specden title2", error);
      delete[] title2;
      title2 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"title3") == 0) {
      if (iarg+2 > nargnew) utils::missing_cmd_args(FLERR, "fix ave/specden title3", error);
      delete[] title3;
      title3 = utils::strdup(arg[iarg+1]);
      iarg += 2;
    } else error->all(FLERR, errptr, "Unknown fix ave/specden keyword: {}", arg[iarg]);
  }

  // setup and error check

  if (nevery <= 0)
    error->all(FLERR, 3, "Illegal fix ave/specden nevery value: {}", nevery);
  if (nrepeat <= 0)
    error->all(FLERR, 4, "Illegal fix ave/specden nrepeat value: {}", nrepeat);
  if (nfreq <= 0)
    error->all(FLERR, 5, "Illegal fix ave/specden nfreq value: {}", nfreq);
  if (nfreq % nevery || nrepeat*nevery > nfreq)
    error->all(FLERR, Error::NOPOINTER, "Inconsistent fix ave/specden nevery/nrepeat/nfreq "
               "values: nfreq must be a multiple of nevery, and nrepeat*nevery <= nfreq");
  if (ave != RUNNING && overwrite)
    error->all(FLERR, Error::NOPOINTER,
               "Fix ave/specden overwrite keyword requires ave running");
  if (nvalues == 0)
    error->all(FLERR, Error::NOPOINTER, "No values specified for fix ave/specden");

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR, val.iarg, "Compute ID {} for fix ave/specden does not exist", val.id);
      if (val.argindex == 0 && val.val.c->scalar_flag == 0)
        error->all(FLERR, val.iarg,
                   "Fix ave/specden compute {} does not calculate a scalar", val.id);
      if (val.argindex && val.val.c->vector_flag == 0)
        error->all(FLERR, val.iarg,
                   "Fix ave/specden compute {} does not calculate a vector", val.id);
      if (val.argindex && val.argindex > val.val.c->size_vector)
        error->all(FLERR, val.iarg,
                   "Fix ave/specden compute {} vector is accessed out-of-range{}",
                   val.id, utils::errorurl(20));

    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR, val.iarg, "Fix ID {} for fix ave/specden does not exist", val.id);
      if (val.argindex == 0 && val.val.f->scalar_flag == 0)
        error->all(FLERR, val.iarg,
                   "Fix ave/specden fix {} does not calculate a scalar", val.id);
      if (val.argindex && val.val.f->vector_flag == 0)
        error->all(FLERR, val.iarg,
                   "Fix ave/specden fix {} does not calculate a vector", val.id);
      if (val.argindex && val.argindex > val.val.f->size_vector)
        error->all(FLERR, val.iarg,
                   "Fix ave/specden fix {} vector is accessed out-of-range{}",
                   val.id, utils::errorurl(20));
      if (nevery % val.val.f->global_freq)
        error->all(FLERR, val.iarg,
                   "Fix {} for fix ave/specden not computed at compatible time{}",
                   val.id, utils::errorurl(7));

    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR, val.iarg,
                   "Variable name {} for fix ave/specden does not exist", val.id);
      if (val.argindex == 0 && input->variable->equalstyle(val.val.v) == 0)
        error->all(FLERR, val.iarg,
                   "Fix ave/specden variable {} is not equal-style variable", val.id);
      if (val.argindex && input->variable->vectorstyle(val.val.v) == 0)
        error->all(FLERR, val.iarg,
                   "Fix ave/specden variable {} is not vector-style variable", val.id);
    }
  }

  // number of positive frequency bins (DC through Nyquist)
  nfreqs = nrepeat / 2 + 1;

  // print file comment lines

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (title1) fprintf(fp,"%s\n",title1);
    else fprintf(fp,"# Spectral density data for fix %s\n",id);
    if (title2) fprintf(fp,"%s\n",title2);
    else fprintf(fp,"# Timestep Number-of-freq-bins\n");
    if (title3) fprintf(fp,"%s\n",title3);
    else {
      fprintf(fp,"# Index Frequency Ncount");
      for (int i = 0; i < nvalues; i++)
        fprintf(fp," specden(%s)",arg[i]);
      fprintf(fp,"\n");
    }
    if (ferror(fp))
      error->one(FLERR, Error::NOPOINTER,
                 "Error writing fix ave/specden header: {}", utils::getsyserror());
    filepos = platform::ftell(fp);
  }

  delete[] title1;
  delete[] title2;
  delete[] title3;

  // if wildcard expansion occurred, free earg memory from expand_args()
  // wait to do this until after file comment lines are printed

  if (expand) {
    for (int i = 0; i < nargnew; i++) delete[] earg[i];
    memory->sfree(earg);
    memory->sfree(amap);
  }

  // allocate and initialize memory

  memory->create(cvalues, nrepeat, nvalues, "ave/specden:cvalues");
  memory->create(specden, nfreqs, nvalues, "ave/specden:specden");
  memory->create(save_specden, nfreqs, nvalues, "ave/specden:save_specden");

  for (int i = 0; i < nrepeat; i++)
    for (int j = 0; j < nvalues; j++)
      cvalues[i][j] = 0.0;

  for (int k = 0; k < nfreqs; k++)
    for (int j = 0; j < nvalues; j++)
      save_specden[k][j] = specden[k][j] = 0.0;

  // this fix produces a global array
  // columns: frequency, ncount, PSD for each input value

  array_flag = 1;
  size_array_rows = nfreqs;
  size_array_cols = nvalues + 2;
  extarray = 0;

  lastindex = -1;
  firstindex = 0;
  nsample = 0;
  noutput = 0;
  nvalid_last = -1;
  nvalid = nextvalid();
  modify->addstep_compute_all(nvalid);
}

/* ---------------------------------------------------------------------- */

FixAveSpecden::~FixAveSpecden()
{
  memory->destroy(cvalues);
  memory->destroy(specden);
  memory->destroy(save_specden);
}

/* ---------------------------------------------------------------------- */

int FixAveSpecden::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAveSpecden::init()
{
  // set current indices for all computes, fixes, variables

  for (auto &val : values) {
    if (val.which == ArgInfo::COMPUTE) {
      val.val.c = modify->get_compute_by_id(val.id);
      if (!val.val.c)
        error->all(FLERR, Error::NOLASTLINE,
                   "Compute ID {} for fix ave/specden does not exist", val.id);
    } else if (val.which == ArgInfo::FIX) {
      val.val.f = modify->get_fix_by_id(val.id);
      if (!val.val.f)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix ID {} for fix ave/specden does not exist", val.id);
    } else if (val.which == ArgInfo::VARIABLE) {
      val.val.v = input->variable->find(val.id.c_str());
      if (val.val.v < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Variable name {} for fix ave/specden does not exist", val.id);
    }
  }

  // need to reset nvalid if nvalid < ntimestep b/c minimize was performed

  if (nvalid < update->ntimestep) {
    lastindex = -1;
    firstindex = 0;
    nsample = 0;
    nvalid = nextvalid();
    modify->addstep_compute_all(nvalid);
  }
}

/* ----------------------------------------------------------------------
   only does something if nvalid = current timestep
------------------------------------------------------------------------- */

void FixAveSpecden::setup(int /*vflag*/)
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixAveSpecden::end_of_step()
{
  int i;

  // skip if not step which requires doing something

  bigint ntimestep = update->ntimestep;
  if (ntimestep != nvalid) return;
  nvalid_last = nvalid;

  // accumulate results of computes, fixes, variables to ring buffer
  // compute/fix/variable may invoke computes so wrap with clear/add

  modify->clearstep_compute();

  // lastindex = index in cvalues ring of latest time sample

  lastindex++;
  if (lastindex == nrepeat) lastindex = 0;

  i = 0;
  for (auto &val : values) {
    double scalar = 0.0;

    // invoke compute if not previously invoked

    if (val.which == ArgInfo::COMPUTE) {
      if (val.argindex == 0) {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_SCALAR)) {
          val.val.c->compute_scalar();
          val.val.c->invoked_flag |= Compute::INVOKED_SCALAR;
        }
        scalar = val.val.c->scalar;
      } else {
        if (!(val.val.c->invoked_flag & Compute::INVOKED_VECTOR)) {
          val.val.c->compute_vector();
          val.val.c->invoked_flag |= Compute::INVOKED_VECTOR;
        }
        scalar = val.val.c->vector[val.argindex-1];
      }

    // access fix fields, guaranteed to be ready

    } else if (val.which == ArgInfo::FIX) {
      if (val.argindex == 0)
        scalar = val.val.f->compute_scalar();
      else
        scalar = val.val.f->compute_vector(val.argindex-1);

    // evaluate equal-style or vector-style variable
    // if index exceeds vector length, use a zero value

    } else if (val.which == ArgInfo::VARIABLE) {
      if (val.argindex == 0)
        scalar = input->variable->compute_equal(val.val.v);
      else {
        double *varvec;
        int nvec = input->variable->compute_vector(val.val.v,&varvec);
        int index = val.argindex;
        if (index > nvec) scalar = 0.0;
        else scalar = varvec[index-1];
      }
    }

    cvalues[lastindex][i++] = scalar;
  }

  // firstindex = index in cvalues ring of earliest time sample
  // nsample = number of time samples in cvalues ring

  if (nsample < nrepeat) nsample++;
  else {
    firstindex++;
    if (firstindex == nrepeat) firstindex = 0;
  }

  nvalid += nevery;
  modify->addstep_compute(nvalid);

  // only compute and output at nfreq steps

  if (ntimestep % nfreq) return;

  // compute spectral density from current ring buffer contents

  compute_specden();

  // save results for compute_array() access

  for (int k = 0; k < nfreqs; k++)
    for (int j = 0; j < nvalues; j++)
      save_specden[k][j] = (noutput > 0) ? specden[k][j] / noutput : 0.0;

  // output result to file

  if (fp && comm->me == 0) {
    clearerr(fp);
    if (overwrite) platform::fseek(fp,filepos);
    utils::print(fp,"{} {}\n",ntimestep,nfreqs);
    double dt = nevery * update->dt;
    for (int k = 0; k < nfreqs; k++) {
      double freq = (dt > 0.0) ? k / (nrepeat * dt) : 0.0;
      fprintf(fp,"%d %g %d",k+1,freq,noutput);
      for (int j = 0; j < nvalues; j++)
        fprintf(fp," %g",save_specden[k][j]);
      fprintf(fp,"\n");
    }
    if (ferror(fp))
      error->one(FLERR, Error::NOLASTLINE,
                 "Error writing fix ave/specden data: {}", utils::getsyserror());
    fflush(fp);
    if (overwrite) {
      bigint fileend = platform::ftell(fp);
      if ((fileend > 0) && (platform::ftruncate(fp,fileend)))
        error->warning(FLERR, "Error while truncating output: {}", utils::getsyserror());
    }
  }

  // zero accumulation if ave = one

  if (ave == ONE) {
    for (int k = 0; k < nfreqs; k++)
      for (int j = 0; j < nvalues; j++)
        specden[k][j] = 0.0;
    noutput = 0;
  }
}

/* ----------------------------------------------------------------------
   compute spectral density via DFT of the time series in cvalues ring
   accumulate into specden[k][i] for ave = running
------------------------------------------------------------------------- */

void FixAveSpecden::compute_specden()
{
  // use nrepeat as the nominal DFT length for consistent frequency bins;
  // samples before the ring buffer is full are effectively zero-padded

  double inv_nrepeat = 1.0 / nrepeat;

  for (int k = 0; k < nfreqs; k++) {
    for (int j = 0; j < nvalues; j++) {
      double re = 0.0, im = 0.0;
      for (int n = 0; n < nsample; n++) {
        int idx = (firstindex + n) % nrepeat;
        double x = cvalues[idx][j];
        double angle = MY_2PI * k * n * inv_nrepeat;
        re += x * cos(angle);
        im -= x * sin(angle);
      }
      // squared magnitude of DFT coefficient, normalized by window length
      specden[k][j] += prefactor * (re*re + im*im) * inv_nrepeat;
    }
  }
  noutput++;
}

/* ----------------------------------------------------------------------
   return I,J array value
   column 0: frequency (cycles per time unit)
   column 1: number of accumulated DFT windows (noutput)
   column 2+: average PSD for each input value
------------------------------------------------------------------------- */

double FixAveSpecden::compute_array(int i, int j)
{
  if (j == 0) {
    double dt = nevery * update->dt;
    return (dt > 0.0) ? i / (nrepeat * dt) : 0.0;
  } else if (j == 1) return 1.0 * noutput;
  else if (noutput) return save_specden[i][j-2];
  return 0.0;
}

/* ----------------------------------------------------------------------
   nvalid = next step on which end_of_step does something
   this step if multiple of nevery, else next multiple
   startstep is lower bound
------------------------------------------------------------------------- */

bigint FixAveSpecden::nextvalid()
{
  bigint nvalid = update->ntimestep;
  if (startstep > nvalid) nvalid = startstep;
  if (nvalid % nevery) nvalid = (nvalid/nevery)*nevery + nevery;
  return nvalid;
}
