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

#ifndef LMP_ELEMENT_MAP_H
#define LMP_ELEMENT_MAP_H

#include "pointers.h"    // IWYU pragma: export

namespace LAMMPS_NS {

class ElementMap : protected Pointers {

 public:
  ElementMap(LAMMPS *lmp);

  void init();                           // initialize map
  void modify_emap(int, char **);        // element command in the input script
  const std::string &find(int);    // find element of given type
  void set(int, const std::string &);    // set element for given type

#if 0
  void read_restart(FILE *fp);
  void write_restart(FILE *);
#endif

 protected:
  std::vector<std::string> element;
  void print();
};
}    // namespace LAMMPS_NS
#endif
