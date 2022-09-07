/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(tersoff/hg,PairTERSOFFHG);
// clang-format on
#else

#ifndef LMP_PAIR_TERSOFF_HG_H
#define LMP_PAIR_TERSOFF_HG_H

#include "HGvector.h"
#include "pair_tersoff.h"
#include <map>
#include <vector>

#define X1_NGRIDPOINTS 5
#define X2_NGRIDPOINTS 5
#define X1_NGRIDSQUARES (X1_NGRIDPOINTS - 1)
#define X2_NGRIDSQUARES (X2_NGRIDPOINTS - 1)

// See "Humbird and Graves, JCP 120, 2405 (2004) for equations

namespace LAMMPS_NS {

class PairTERSOFFHG : public PairTersoff {

  static constexpr int NPARAMS_PER_LINE = 28;

 private:
  HGvector Fij;
  double bbar_ij;
  double bsp_ij, bsp_ji, Nconj_ij;
  int cnt, scrcount;
  double VA_ij, VR_ij, dVA_ij, dVR_ij;
  double iFv_ij, jFv_ij, force_ik, force_jk;
  std::vector<double> iFv_ik, iFv_jk;    //, iFv_ik2;
  std::vector<double> jFv_ik, jFv_jk;    //, jFv_jk2;
  std::map<int, std::map<int, double>> Nmap;

 public:
  PairTERSOFFHG(class LAMMPS *);
  ~PairTERSOFFHG() override {}
  void compute(int, int) override;
  void init_style() override;

 protected:
  int maxlocal;             // size of numneigh, firstneigh arrays
  int pgsize;               // size of neighbor page
  int oneatom;              // max # of neighbors for one atom
  MyPage<int> *ipage;       // neighbor list pages
  int *REBO_numneigh;       // # of pair neighbors for each atom
  int **REBO_firstneigh;    // ptr to 1st neighbor of each atom

  // Library file containing the spline coefficient
  void read_lib(Param *);
  void read_file(char *);

  void *xmalloc(size_t);

  // Coordination terms, Hij of Eq. A12
  void count_neigh();
  int nmax;
  double coordenergy[5][402], coordforce[5][402], coordnumber[5][402];
  double *NCl, *NSi;

  // communication functions
  int pack_flag;
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  void bicubic_genCoef(double y[X1_NGRIDPOINTS][X2_NGRIDPOINTS], Param *);
  void bcucof(double y[], double y1[], double y2[], double y12[], double d1, double d2,
              double c[4][4]);
  double BondOrder(int, int, double, double, int, int);
  void bicubicint(double x1, double x2, double *y, double *y1, double *y2, Param *);
  void bcuint(double x1l, double x1u, double x2l, double x2u, double x1, double x2, double *ansy,
              double *ansy1, double *ansy2, double c[4][4]);
};

}    // namespace LAMMPS_NS

#endif
#endif
