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

#include "element_map.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
//#include "elements.h"
#include "error.h"
#include "force.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ElementMap::ElementMap(LAMMPS *_lmp) : Pointers(_lmp)
{
  element.clear();
  element.resize(1);
  element[0] = "Unk";
  print();
}

/* ---------------------------------------------------------------------- */

void ElementMap::init()
{
  // resize element array if possible and needed

  if ((domain->box_exist == 0) || (element.size() > 1)) return;

  element.resize(atom->ntypes + 1);
}

/* ----------------------------------------------------------------------
   element command in input script
------------------------------------------------------------------------- */

void ElementMap::modify_emap(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Incorrect number of arguments for element command");

  init();
  const int ntypes = atom->ntypes;

  char *typestr = utils::expand_type(FLERR, arg[0], Atom::ATOM, lmp);
  const std::string str = typestr ? typestr : arg[0];
  delete[] typestr;

  int lo, hi;
  utils::bounds(FLERR, str, 1, ntypes, lo, hi, error);
  if ((lo < 1) || (hi > ntypes)) error->all(FLERR, "Invalid atom type {} for element", str);

  for (int itype = lo; itype <= hi; itype++) element[itype] = arg[1];

  print();
}

/* ----------------------------------------------------------------------
   get element for given atom type. guess from mass if not set.
------------------------------------------------------------------------- */

const std::string &ElementMap::find(int itype)
{
  init();

  // protect against out-of-bounds access
  if ((itype < 0) || (itype > atom->ntypes)) itype = 0;

  // guess element from mass, if not set, return "Unk" if not detected
  if (!element[itype].size()) { return element[0]; }

  return element[itype];
}

/* ----------------------------------------------------------------------
   set element for given atom type. skip if type is out-of-bounds.
------------------------------------------------------------------------- */

void ElementMap::set(int itype, const std::string &str)
{
  // protect against out-of-bounds access
  if ((itype < 1) || (itype >= (int) element.size())) return;

  element[itype] = str;
}

/* ---------------------------------------------------------------------- */

void ElementMap::print()
{
  utils::logmesg(lmp, "Current element map:\n");
  for (std::size_t i = 0; i < element.size(); ++i)
    utils::logmesg(lmp, " {:4}  {}\n", i, element[i]);
}

#if 0
/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void ElementMap::read_restart(FILE *fp)
{
  char *charlabel;

  for (int i = 0; i < natomtypes; i++) {
    charlabel = read_string(fp);
    typelabel[i] = charlabel;
    typelabel_map[charlabel] = i + 1;
    delete[] charlabel;
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void ElementMap::write_restart(FILE *fp)
{
  for (int i = 0; i < natomtypes; i++) write_string(element[i], fp);
}

/* ----------------------------------------------------------------------
   read a char string (including nullptr) and bcast it
   str is allocated here, ptr is returned, caller must deallocate
------------------------------------------------------------------------- */

char *ElementMap::read_string(FILE *fp)
{
  int n = read_int(fp);
  if (n < 0) error->all(FLERR, "Illegal size string or corrupt restart");
  char *value = new char[n];
  if (comm->me == 0) utils::sfread(FLERR, value, sizeof(char), n, fp, nullptr, error);
  MPI_Bcast(value, n, MPI_CHAR, 0, world);
  return value;
}

/* ----------------------------------------------------------------------
   write a flag and a C-style char string (including the terminating null
   byte) into the restart file
------------------------------------------------------------------------- */

void ElementMap::write_string(const std::string &str, FILE *fp)
{
  const char *cstr = str.c_str();
  int n = strlen(cstr) + 1;
  fwrite(&n, sizeof(int), 1, fp);
  fwrite(cstr, sizeof(char), n, fp);
}
#endif
