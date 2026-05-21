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

/* ----------------------------------------------------------------------
   ILVES: Iterative Linear-Velocity-Equation Solver for bond constraints

   Implements the ILVES algorithm from:
     P. Garcia-Risueno et al., J. Chem. Theory Comput. (2025)
     DOI: 10.1021/acs.jctc.5c01376

   Key advantages over SHAKE:
     1. Handles arbitrarily large constraint networks (chains, rings, etc.)
        unlike SHAKE which only handles small clusters (up to 4 atoms)
     2. Simultaneous (Jacobi-style) Newton update: all corrections computed
        from the same positions and applied together, giving quadratic
        convergence for tightly coupled constraint networks
     3. Parallel via MPI ghost communication between Newton iterations

   Angle constraints are treated as three bond constraints following
   the same approach used by fix shake: the angle A-B-C (B = central)
   is enforced by constraining:
     - Bond B-A to its equilibrium length
     - Bond B-C to its equilibrium length
     - A virtual bond A-C to sqrt(b1^2+b2^2-2*b1*b2*cos(theta0))
   All three become regular entries in the flat constraint list and
   are treated identically by the Newton iteration.

   Algorithm (each timestep):
     1. Compute unconstrained positions (same as SHAKE)
     2. Forward communicate ghost xshake positions
     3. Newton iteration loop:
        a. Compute residuals g_k = 0.5*(|s_k|^2 - d_k^2)
           where s_k = xshake[a] - xshake[b] (current direction)
        b. Compute corrections dl_k = -g_k / A_kk  (note: negative sign)
           where A_kk = (1/m_a + 1/m_b)*(r_k . s_k)
           and r_k = x[a] - x[b] (OLD positions)
        c. Apply ALL corrections simultaneously:
           xshake[a] += dl_k/m_a * r_k
           xshake[b] -= dl_k/m_b * r_k
        d. Forward communicate ghost xshake
        e. Check global convergence (MPI max of |g_k|/d_k^2)
     4. Apply constraint forces from accumulated Lagrange multipliers

   Note on parallel constraints: For constraint chains spanning multiple
   MPI ranks, convergence requires multiple Newton iterations (proportional
   to the number of inter-rank constraint boundaries). The full Schur-based
   parallel solver from the original GROMACS implementation can be added
   in future work to eliminate this scaling.
------------------------------------------------------------------------- */

#include "fix_ilves.h"

#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "label_map.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr int DELTA_CONSTR = 256;
static constexpr double SMALL_DENOM = 1.0e-10;
static constexpr double BIG = 1.0e20;

/* ---------------------------------------------------------------------- */

FixIlves::FixIlves(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), bond_flag(nullptr), bond_distance(nullptr), angle_flag(nullptr),
    angle_distance(nullptr), has_angle(false), ilves_flag(nullptr), xshake(nullptr),
    c_atom1(nullptr), c_atom2(nullptr), c_dist2(nullptr), c_lambda(nullptr), x(nullptr), v(nullptr),
    f(nullptr), mass(nullptr), rmass(nullptr), type(nullptr), b_count(nullptr),
    b_count_all(nullptr), b_ave(nullptr), b_ave_all(nullptr), b_max(nullptr), b_max_all(nullptr),
    b_min(nullptr), b_min_all(nullptr), a_count(nullptr), a_count_all(nullptr), a_ave(nullptr),
    a_max(nullptr), a_min(nullptr), a_ave_all(nullptr), a_max_all(nullptr), a_min_all(nullptr)
{
  if (narg < 8) utils::missing_cmd_args(FLERR, "fix ilves", error);

  // error checks

  if (atom->molecular == Atom::ATOMIC)
    error->all(FLERR, 2, "Cannot use fix ilves with non-molecular system");

  // parse required args: tol iter N constraint_specs

  tolerance = utils::numeric(FLERR, arg[3], false, lmp);
  max_iter = utils::inumeric(FLERR, arg[4], false, lmp);
  output_every = utils::inumeric(FLERR, arg[5], false, lmp);

  // check if any typelabels conflict with fix ilves arguments.

  bool allow_typelabels = (atom->labelmapflag != 0);
  if (allow_typelabels) {
    for (int i = Atom::ATOM; i < Atom::DIHEDRAL; ++i) {
      if ((atom->lmap->find_type("b", i) >= 0) || (atom->lmap->find_type("a", i) >= 0) ||
          (atom->lmap->find_type("t", i) >= 0) || (atom->lmap->find_type("m", i) >= 0))
        allow_typelabels = false;
    }
    if (!allow_typelabels && (comm->me == 0))
      error->warning(FLERR,
                     "At least one typelabel conflicts with a fix ilves option: support for "
                     "typelabels is disabled.");
  }

  // parse constraint specs: bond types (b keyword) and angle types (a keyword)

  int nb = atom->nbondtypes + 1;
  int na = atom->nangletypes + 1;
  bond_flag = new int[nb];
  bond_distance = new double[nb];
  angle_flag = new int[na];
  angle_distance = new double[na];

  for (int i = 0; i < nb; ++i) {
    bond_flag[i] = 0;
    bond_distance[i] = 0.0;
  }
  for (int i = 0; i < na; ++i) {
    angle_flag[i] = 0;
    angle_distance[i] = 0.0;
  }

  int next = 6;
  char mode = '\0';
  while (next < narg) {
    int i = -1;
    if (strcmp(arg[next], "b") == 0) {
      mode = 'b';
    } else if (strcmp(arg[next], "a") == 0) {
      mode = 'a';

      // break if known optional keyword

    } else if ((strcmp(arg[next], "kbond") == 0) || (strcmp(arg[next], "store") == 0)) {
      break;

      // get numeric types for b or a keywords.

    } else if (mode == 'b') {
      if (allow_typelabels)
        i = utils::expand_type_int(FLERR, arg[next], Atom::BOND, lmp);
      else
        i = utils::inumeric(FLERR, arg[next], false, lmp);
      if (i < 1 || i > atom->nbondtypes)
        error->all(FLERR, next, "Invalid bond type {} for fix ilves", arg[next]);
      bond_flag[i] = 1;
    } else if (mode == 'a') {
      if (allow_typelabels)
        i = utils::expand_type_int(FLERR, arg[next], Atom::ANGLE, lmp);
      else
        i = utils::inumeric(FLERR, arg[next], false, lmp);
      if (i < 1 || i > atom->nangletypes)
        error->all(FLERR, next, "Invalid angle type {} for fix ilves", arg[next]);
      angle_flag[i] = 1;
      has_angle = true;
    } else {
      error->all(FLERR, "Unknown fix {} option: {}", style, arg[next]);
    }
    next++;
  }

  // allocate bond distance array (indexed 1..nbondtypes)
  // statistics arrays (per bond type, indices 1..nbondtypes)

  if (output_every) {
    int nb = atom->nbondtypes + 1;
    b_count = new bigint[nb];
    b_count_all = new bigint[nb];
    b_ave = new double[nb];
    b_ave_all = new double[nb];
    b_max = new double[nb];
    b_max_all = new double[nb];
    b_min = new double[nb];
    b_min_all = new double[nb];

    int na = atom->nangletypes + 1;
    a_count = new bigint[na];
    a_count_all = new bigint[na];
    a_ave = new double[na];
    a_ave_all = new double[na];
    a_max = new double[na];
    a_max_all = new double[na];
    a_min = new double[na];
    a_min_all = new double[na];
  }

  // per-atom arrays, registered with Atom class for dynamic resizing

  n_constr = max_constr = 0;
  next_output = -1;

  FixIlves::grow_arrays(atom->nmax);
  atom->add_callback(Atom::GROW);

  comm_forward = 3;    // communicate 3 doubles (xshake x/y/z) per atom
}

/* ---------------------------------------------------------------------- */

FixIlves::~FixIlves()
{
  // unregister per-atom arrays with Atom class

  atom->delete_callback(id, Atom::GROW);

  // free per-atom arrays

  memory->destroy(ilves_flag);
  memory->destroy(xshake);

  // free constraint list

  memory->destroy(c_atom1);
  memory->destroy(c_atom2);
  memory->destroy(c_dist2);
  memory->destroy(c_lambda);

  // free input arrays

  delete[] bond_flag;
  delete[] bond_distance;
  delete[] angle_flag;
  delete[] angle_distance;

  // free statistics arrays

  if (output_every) {
    delete[] b_count;
    delete[] b_count_all;
    delete[] b_ave;
    delete[] b_ave_all;
    delete[] b_max;
    delete[] b_max_all;
    delete[] b_min;
    delete[] b_min_all;

    delete[] a_count;
    delete[] a_count_all;
    delete[] a_ave;
    delete[] a_ave_all;
    delete[] a_max;
    delete[] a_max_all;
    delete[] a_min;
    delete[] a_min_all;
  }
}

/* ---------------------------------------------------------------------- */

int FixIlves::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= POST_NEIGHBOR;
  mask |= POST_FORCE;
  return mask;
}

/* ----------------------------------------------------------------------
   Set bond equilibrium distances from force->bond.
   Called after force->bond is initialized.
------------------------------------------------------------------------- */

void FixIlves::init()
{
  // error if more than one ilves fix

  if (modify->get_fix_by_style(fmt::format("^{}", style)).size() > 1)
    error->all(FLERR, Error::NOLASTLINE, "More than one fix ilves instance");

  // cannot use during minimization (turns off bond forces)

  if (update->whichflag == 2)
    error->all(FLERR, Error::NOLASTLINE, "Fix ilves cannot be used with minimization");

  // error if box changes come before ilves fix

  bool boxflag = false;
  for (const auto &ifix : modify->get_fix_list()) {
    if (boxflag && utils::strmatch(ifix->style, fmt::format("^{}", style)))
      error->all(FLERR, Error::NOLASTLINE, "Fix ilves must come before any box-changing fix");
    if (ifix->box_change) boxflag = true;
  }

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Fix ilves is not compatible with run style respa");

  if (!force->bond)
    error->all(FLERR, Error::NOLASTLINE, "Bond style must be defined for fix ilves");

  // get equilibrium bond distances from the bond force style

  for (int i = 1; i <= atom->nbondtypes; i++)
    bond_distance[i] = force->bond->equilibrium_distance(i);

  // compute virtual bond distances for angle constraints
  // uses the law of cosines: d12 = sqrt(b1^2+b2^2-2*b1*b2*cos(theta0))

  if (has_angle) {
    if (!force->angle)
      error->all(FLERR, Error::NOLASTLINE,
                 "Angle style must be defined for fix ilves when using 'a' keyword");

    for (int i = 1; i <= atom->nangletypes; i++) {
      if (!angle_flag[i]) continue;

      // scan local atoms to find an example of this angle type
      // and extract the arm bond types from the central atom's bond list

      int b1type = 0, b2type = 0, flag = 0;

      for (int j = 0; j < atom->nlocal && !flag; j++) {
        for (int m = 0; m < atom->num_angle[j]; m++) {
          if (std::abs(atom->angle_type[j][m]) != i) continue;

          // find central atom (angle_atom2)
          int i2 = atom->map(atom->angle_atom2[j][m]);
          if (i2 == -1 || i2 >= atom->nlocal) continue;

          tagint ta1 = atom->angle_atom1[j][m];
          tagint ta3 = atom->angle_atom3[j][m];

          // search central atom i2's bonds for connections to ta1 and ta3
          int bt1 = 0, bt2 = 0;
          for (int n = 0; n < atom->num_bond[i2]; n++) {
            if (atom->bond_atom[i2][n] == ta1) bt1 = std::abs(atom->bond_type[i2][n]);
            if (atom->bond_atom[i2][n] == ta3) bt2 = std::abs(atom->bond_type[i2][n]);
          }
          if (bt1 == 0 || bt2 == 0) continue;

          int t1 = MIN(bt1, bt2);
          int t2 = MAX(bt1, bt2);
          if (b1type == 0) {
            b1type = t1;
            b2type = t2;
          } else if (t1 != b1type || t2 != b2type) {
            flag = 1;    // inconsistent bond types for this angle type
          }
          break;
        }
      }

      // error if inconsistent bond types found on any proc
      int flag_all = 0;
      MPI_Allreduce(&flag, &flag_all, 1, MPI_INT, MPI_MAX, world);
      if (flag_all)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix ilves angle type {} has inconsistent arm bond types", i);

      // ensure all procs have the bond types (from whichever proc found them)
      MPI_Allreduce(MPI_IN_PLACE, &b1type, 1, MPI_INT, MPI_MAX, world);
      MPI_Allreduce(MPI_IN_PLACE, &b2type, 1, MPI_INT, MPI_MAX, world);

      if (b1type == 0) {
        angle_distance[i] = 0.0;    // no angle of this type found locally
        continue;
      }

      const double theta0 = force->angle->equilibrium_angle(i);
      const double d1 = bond_distance[b1type];
      const double d2 = bond_distance[b2type];
      angle_distance[i] = sqrt(d1 * d1 + d2 * d2 - 2.0 * d1 * d2 * cos(theta0));
    }
  }

  // set/update timestep info

  dtv = update->dt;
  dtfsq = update->dt * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   Called at the start of each run to initialize.
   Mark constrained bond types negative and build initial constraint list.
------------------------------------------------------------------------- */

void FixIlves::setup(int vflag)
{
  // update local pointers

  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass = atom->mass;
  rmass = atom->rmass;
  type = atom->type;
  nlocal = atom->nlocal;

  // perform initial marking of constrained bonds as negative

  pre_neighbor();

  // build initial constraint list from the neighbor bondlist

  post_neighbor();

  // set up output schedule

  bigint ntimestep = update->ntimestep;
  if (output_every) {
    next_output = ntimestep + output_every;
    if (ntimestep % output_every != 0)
      next_output = (ntimestep / output_every) * output_every + output_every;
  } else {
    next_output = -1;
  }

  if (output_every) stats();

  // apply constraints at the current (time 0) positions

  post_force(vflag);
}

/* ----------------------------------------------------------------------
   Called before each neighbor list rebuild.
   Re-mark constrained bond types as negative in per-atom bond_type arrays
   so they are excluded from the regular bond force computation.
------------------------------------------------------------------------- */

void FixIlves::pre_neighbor()
{
  // update local copies of atom pointers (atom exchange may have occurred)

  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass = atom->mass;
  rmass = atom->rmass;
  type = atom->type;
  nlocal = atom->nlocal;

  // mark constrained bond types negative in per-atom arrays
  // (they may have been reset to positive during atom exchange)

  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;

  for (int i = 0; i < nlocal; i++) {
    for (int m = 0; m < num_bond[i]; m++) {
      int btype = std::abs(bond_type[i][m]);
      if (btype >= 1 && btype <= atom->nbondtypes && bond_flag[btype]) bond_type[i][m] = -btype;
    }
  }

  // mark constrained angle types negative and also mark their arm bonds negative
  // so those bonds are excluded from the regular force computation

  if (has_angle) {
    int *num_angle = atom->num_angle;
    int **angle_type = atom->angle_type;
    tagint **angle_atom1 = atom->angle_atom1;
    tagint **angle_atom2 = atom->angle_atom2;
    tagint **angle_atom3 = atom->angle_atom3;
    tagint *tag = atom->tag;

    for (int i = 0; i < nlocal; i++) {
      for (int m = 0; m < num_angle[i]; m++) {
        int atype = std::abs(angle_type[i][m]);
        if (atype < 1 || atype > atom->nangletypes || !angle_flag[atype]) continue;

        // mark the angle type negative (exclude from angle force computation)
        angle_type[i][m] = -atype;

        // mark arm bonds negative: the bonds connecting atom i to the
        // other two atoms in this angle via the arm-bond connections
        tagint ti = tag[i];
        tagint ta1 = angle_atom1[i][m];    // end atom 1 global tag
        tagint ta2 = angle_atom2[i][m];    // central atom global tag
        tagint ta3 = angle_atom3[i][m];    // end atom 3 global tag

        for (int n = 0; n < num_bond[i]; n++) {
          tagint bt = bond_atom[i][n];
          // atom i participates in this angle in one of three roles:
          //   (a) i is the central atom (ti == ta2): its arms connect to ta1 and ta3
          //   (b) i is end atom 1 (ti == ta1): its arm connects to central ta2
          //   (c) i is end atom 3 (ti == ta3): its arm connects to central ta2
          bool is_arm = (ti == ta2 && (bt == ta1 || bt == ta3)) || (ti == ta1 && bt == ta2) ||
              (ti == ta3 && bt == ta2);
          if (is_arm) bond_type[i][n] = -std::abs(bond_type[i][n]);
        }
      }
    }
  }

  // clear ilves_flag for all local atoms; will be re-set in post_neighbor

  for (int i = 0; i < nlocal; i++) ilves_flag[i] = 0;
}

/* ----------------------------------------------------------------------
   Called after each neighbor list rebuild.
   Rebuilds the flat constraint list from per-atom bond and angle data.

   Bonds marked negative (by pre_neighbor) are excluded from the neighbor
   bondlist, so we scan per-atom bond data directly.  For angles, we scan
   per-atom angle data using a minimum-global-tag ownership rule that
   ensures each constraint is added exactly once across all MPI ranks.
------------------------------------------------------------------------- */

void FixIlves::post_neighbor()
{
  // reset constraint count; keep allocated memory

  n_constr = 0;

  const int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  tagint *tag = atom->tag;

  // -----------------------------------------------------------------
  // Phase 1: bond constraints (b keyword).
  //
  // Bonds flagged by bond_flag were marked negative in pre_neighbor.
  // We scan per-atom bond data (not the neighbor bondlist, which
  // excludes negative-type bonds) to find them.
  //
  // Ownership rule: add constraint (i,j) when i < nlocal and either
  // j >= nlocal (j is a ghost) or i < j (i has smaller local index).
  // This ensures each bond is counted exactly once.
  // -----------------------------------------------------------------

  for (int i = 0; i < nlocal; i++) {
    for (int m = 0; m < num_bond[i]; m++) {
      const int btype = bond_type[i][m];
      if (btype >= 0) continue;    // not constrained (positive = normal)
      const int bt = -btype;
      if (!bond_flag[bt]) continue;

      int j = atom->map(bond_atom[i][m]);
      if (j == -1) continue;
      j = domain->closest_image(i, j);

      // Skip double-counting only when newton_bond=0 (bond stored at both
      // atoms): let the atom with the smaller local index own the constraint.
      // When newton_bond=1 the bond is stored at atom i only, so always add.
      if (!force->newton_bond && j < nlocal && j < i) continue;

      add_constraint(i, j, bond_distance[bt]);
    }
  }

  // -----------------------------------------------------------------
  // Phase 2: angle constraints (a keyword).
  //
  // Each constrained angle A(i1)-B(i2)-C(i3) is enforced via three
  // bond constraints:
  //   (a) B-A: equilibrium bond length of that bond type
  //   (b) B-C: equilibrium bond length of that bond type
  //   (c) virtual A-C: sqrt(b1^2+b2^2-2*b1*b2*cos(theta0))
  //
  // Ownership of bond (a) and (b) uses min-global-tag of the two
  // endpoint atoms, so each is added exactly once.  Bond (a) and (b)
  // are skipped when their bond type is already in bond_flag (the b
  // keyword constraint already covers it).
  //
  // Ownership of virtual (c): with newton_bond=1 (default) only the
  // central atom B stores this angle, so B always adds the virtual bond.
  // With newton_bond=0 all three atoms store the angle; the one with the
  // smallest global tag among all three adds it exactly once.  In both
  // cases the virtual bond connects the two end atoms A and C, not B.
  // -----------------------------------------------------------------

  if (has_angle) {
    const int *num_angle = atom->num_angle;
    int **angle_type = atom->angle_type;
    tagint **angle_atom1 = atom->angle_atom1;
    tagint **angle_atom2 = atom->angle_atom2;
    tagint **angle_atom3 = atom->angle_atom3;

    for (int i = 0; i < nlocal; i++) {
      for (int m = 0; m < num_angle[i]; m++) {
        const int atype = std::abs(angle_type[i][m]);
        if (atype < 1 || atype > atom->nangletypes) continue;
        if (!angle_flag[atype]) continue;
        if (angle_distance[atype] == 0.0) continue;

        tagint ta1 = angle_atom1[i][m];    // global tag of end atom 1
        tagint ta2 = angle_atom2[i][m];    // global tag of central atom
        tagint ta3 = angle_atom3[i][m];    // global tag of end atom 3

        // map to local indices (closest periodic images)
        int ia1 = atom->map(ta1);
        int ia2 = atom->map(ta2);
        int ia3 = atom->map(ta3);
        if (ia1 == -1 || ia2 == -1 || ia3 == -1) continue;
        ia1 = domain->closest_image(i, ia1);
        ia2 = domain->closest_image(i, ia2);
        ia3 = domain->closest_image(i, ia3);

        const tagint ti = tag[i];
        const tagint ti1 = tag[ia1];
        const tagint ti2 = tag[ia2];
        const tagint ti3 = tag[ia3];

        // --- real arm bond ia2-ia1 ---
        // owned by whichever of {ia2, ia1} has the smaller global tag
        if (ti == MIN(ti2, ti1)) {
          // search atom i's own bond list for the bond to the other atom
          tagint other_tag = (ti == ti2) ? ta1 : ta2;
          int btype_arm = 0;
          for (int n = 0; n < num_bond[i]; n++) {
            if (bond_atom[i][n] == other_tag) {
              btype_arm = std::abs(bond_type[i][n]);
              break;
            }
          }
          if (btype_arm > 0 && !bond_flag[btype_arm]) {
            int other = (ti == ti2) ? ia1 : ia2;
            add_constraint(i, other, bond_distance[btype_arm]);
          }
        }

        // --- real arm bond ia2-ia3 ---
        if (ti == MIN(ti2, ti3)) {
          tagint other_tag = (ti == ti2) ? ta3 : ta2;
          int btype_arm = 0;
          for (int n = 0; n < num_bond[i]; n++) {
            if (bond_atom[i][n] == other_tag) {
              btype_arm = std::abs(bond_type[i][n]);
              break;
            }
          }
          if (btype_arm > 0 && !bond_flag[btype_arm]) {
            int other = (ti == ti2) ? ia3 : ia2;
            add_constraint(i, other, bond_distance[btype_arm]);
          }
        }

        // --- virtual bond ia1-ia3 ---
        // Ownership rule differs by newton_bond setting:
        //   newton_bond=1: bond data stored at one atom only (central atom
        //     atom2 is the sole holder of this angle entry), so the central
        //     atom must always add the virtual bond (condition: ti == ti2).
        //   newton_bond=0: angle stored at all three atoms; add from the atom
        //     with the smallest global tag among all three (condition:
        //     ti == MIN(ti1, ti2, ti3)).
        // In both cases the virtual bond connects the two END atoms ia1 and
        // ia3, NOT atom i (which may be the central atom).
        const bool is_owner = force->newton_bond ? (ti == ti2) : (ti == MIN(MIN(ti1, ti2), ti3));
        if (is_owner) {
          // c_atom1 must be a local atom; prefer the end atom with the
          // smaller global tag as c_atom1.
          int c_va = (ti1 <= ti3) ? ia1 : ia3;
          int c_vb = (ti1 <= ti3) ? ia3 : ia1;
          if (c_va >= nlocal) {
            if (c_vb >= nlocal) continue;    // both end atoms are ghosts; skip
            std::swap(c_va, c_vb);
          }
          add_constraint(c_va, c_vb, angle_distance[atype]);
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixIlves::add_constraint(int i, int j, double dist)
{
  if (n_constr == max_constr) {
    max_constr += DELTA_CONSTR;
    memory->grow(c_atom1, max_constr, "ilves:c_atom1");
    memory->grow(c_atom2, max_constr, "ilves:c_atom2");
    memory->grow(c_dist2, max_constr, "ilves:c_dist2");
    memory->grow(c_lambda, max_constr, "ilves:c_lambda");
  }

  c_atom1[n_constr] = i;
  c_atom2[n_constr] = j;
  c_dist2[n_constr] = dist * dist;
  c_lambda[n_constr] = 0.0;
  n_constr++;

  ilves_flag[i] = 1;
  if (j < nlocal) ilves_flag[j] = 1;
}

/* ----------------------------------------------------------------------
   Main constraint enforcement using the ILVES Newton algorithm.

   The key difference from SHAKE:
   - All Newton corrections are computed simultaneously from the same
     positions, then applied all at once (Jacobi-style update).
   - This gives quadratic convergence for tightly coupled constraint networks
     and supports arbitrarily large constraint graphs.
------------------------------------------------------------------------- */

void FixIlves::post_force(int vflag)
{
  if (output_every && update->ntimestep == next_output) stats();

  // update local pointers in case they changed

  x = atom->x;
  v = atom->v;
  f = atom->f;
  mass = atom->mass;
  rmass = atom->rmass;
  type = atom->type;
  nlocal = atom->nlocal;

  // reset accumulated Lagrange multipliers for this timestep

  for (int k = 0; k < n_constr; k++) c_lambda[k] = 0.0;

  // -----------------------------------------------------------------
  // Step 1: Compute unconstrained positions (same as SHAKE).
  //         xshake[i] = x[i] + dt*v[i] + 0.5*dt^2*f[i]/m[i]
  // -----------------------------------------------------------------

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (ilves_flag[i]) {
        const double dtfmsq = dtfsq / rmass[i];
        xshake[i][0] = x[i][0] + dtv * v[i][0] + dtfmsq * f[i][0];
        xshake[i][1] = x[i][1] + dtv * v[i][1] + dtfmsq * f[i][1];
        xshake[i][2] = x[i][2] + dtv * v[i][2] + dtfmsq * f[i][2];
      } else {
        xshake[i][0] = xshake[i][1] = xshake[i][2] = 0.0;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (ilves_flag[i]) {
        const double dtfmsq = dtfsq / mass[type[i]];
        xshake[i][0] = x[i][0] + dtv * v[i][0] + dtfmsq * f[i][0];
        xshake[i][1] = x[i][1] + dtv * v[i][1] + dtfmsq * f[i][1];
        xshake[i][2] = x[i][2] + dtv * v[i][2] + dtfmsq * f[i][2];
      } else {
        xshake[i][0] = xshake[i][1] = xshake[i][2] = 0.0;
      }
    }
  }

  // communicate unconstrained positions to ghost atoms

  comm->forward_comm(this);

  // -----------------------------------------------------------------
  // Step 2: Newton iteration loop.
  //
  //   ILVES Newton step (diagonal Jacobian approximation, Eq. from paper):
  //
  //   g_k   = 0.5 * (|s_k|^2 - d_k^2)   [constraint residual]
  //   A_kk  = (1/m_a + 1/m_b) * (r_k . s_k)   [diagonal Jacobian]
  //   dl_k  = -g_k / A_kk               [Newton correction, note sign]
  //
  //   where r_k = x[a]-x[b]  (OLD positions -- fixed for entire loop)
  //         s_k = xshake[a]-xshake[b]  (current positions -- updated each iter)
  //
  //   Position update (simultaneous for all constraints):
  //   xshake[a] += dl_k * (1/m_a) * r_k   (dl_k < 0 for stretched bond)
  //   xshake[b] -= dl_k * (1/m_b) * r_k
  //
  //   Sign convention: when bond is stretched (g_k > 0) and r||s,
  //   A_kk > 0 => dl_k = -g/A < 0 => dl*invma*r points toward b (shortens) ✓
  // -----------------------------------------------------------------

  // dl array: Newton corrections for all constraints, computed simultaneously
  // allocated on stack for thread safety; will be dynamic for large n_constr

  auto dl = std::vector<double>(n_constr);

  for (int iter = 0; iter < max_iter; iter++) {

    // --- compute residuals and check convergence ---

    double max_err = 0.0;

    for (int k = 0; k < n_constr; k++) {
      const int a = c_atom1[k];
      const int b = c_atom2[k];

      const double sx = xshake[a][0] - xshake[b][0];
      const double sy = xshake[a][1] - xshake[b][1];
      const double sz = xshake[a][2] - xshake[b][2];
      const double s2 = sx * sx + sy * sy + sz * sz;

      const double g = 0.5 * (s2 - c_dist2[k]);
      const double err = std::abs(g) / c_dist2[k];
      if (err > max_err) max_err = err;
    }

    // global convergence check across all MPI ranks

    MPI_Allreduce(MPI_IN_PLACE, &max_err, 1, MPI_DOUBLE, MPI_MAX, world);
    if (max_err <= tolerance) break;

    // --- compute all Newton corrections simultaneously ---
    //     uses OLD positions (r_k) and CURRENT xshake (s_k)

    for (int k = 0; k < n_constr; k++) {
      const int a = c_atom1[k];
      const int b = c_atom2[k];

      // r_k = x[a] - x[b]  (OLD positions, unchanged throughout loop)
      const double rx = x[a][0] - x[b][0];
      const double ry = x[a][1] - x[b][1];
      const double rz = x[a][2] - x[b][2];

      // s_k = xshake[a] - xshake[b]  (current working positions)
      const double sx = xshake[a][0] - xshake[b][0];
      const double sy = xshake[a][1] - xshake[b][1];
      const double sz = xshake[a][2] - xshake[b][2];
      const double s2 = sx * sx + sy * sy + sz * sz;

      // constraint residual
      const double g = 0.5 * (s2 - c_dist2[k]);

      // inverse masses
      const double invma = rmass ? 1.0 / rmass[a] : 1.0 / mass[type[a]];
      const double invmb = rmass ? 1.0 / rmass[b] : 1.0 / mass[type[b]];

      // diagonal Jacobian entry: A_kk = (1/m_a + 1/m_b) * (r_k . s_k)
      const double rs = rx * sx + ry * sy + rz * sz;
      const double A_kk = (invma + invmb) * rs;

      // Newton step: dl = -g / A_kk  (negative sign is essential:
      // positive residual g means bond stretched, so dl < 0, which then
      // moves atoms toward each other via the += dl*invm*r update below).
      // Skip if denominator is too small.
      dl[k] = (std::abs(A_kk) > SMALL_DENOM) ? -g / A_kk : 0.0;
    }

    // --- apply ALL corrections simultaneously (Jacobi/Newton update) ---
    //     This is the key ILVES feature: corrections are computed from
    //     the SAME positions and applied all at once, unlike SHAKE's
    //     sequential (Gauss-Seidel) update which can slow convergence
    //     for coupled constraints.

    for (int k = 0; k < n_constr; k++) {
      if (dl[k] == 0.0) continue;

      const int a = c_atom1[k];
      const int b = c_atom2[k];

      const double rx = x[a][0] - x[b][0];
      const double ry = x[a][1] - x[b][1];
      const double rz = x[a][2] - x[b][2];

      const double invma = rmass ? 1.0 / rmass[a] : 1.0 / mass[type[a]];
      const double invmb = rmass ? 1.0 / rmass[b] : 1.0 / mass[type[b]];

      const double ca = dl[k] * invma;
      const double cb = dl[k] * invmb;

      xshake[a][0] += ca * rx;
      xshake[a][1] += ca * ry;
      xshake[a][2] += ca * rz;
      xshake[b][0] -= cb * rx;
      xshake[b][1] -= cb * ry;
      xshake[b][2] -= cb * rz;

      // accumulate total Lagrange multiplier for force application
      c_lambda[k] += dl[k];
    }

    // communicate updated ghost positions for the next iteration
    // (corrections to cross-boundary chains propagate one link per iteration)

    comm->forward_comm(this);
  }

  // -----------------------------------------------------------------
  // Step 3: Apply constraint forces to f[].
  //
  //   The constraint force for bond k is:
  //   f[a] += (c_lambda[k] / dtfsq) * r_k
  //   f[b] -= (c_lambda[k] / dtfsq) * r_k
  //   where r_k = x[a] - x[b] (OLD positions)
  //
  //   This is derived from: constraint force = position_correction * mass / dtfsq
  //   and position correction for a = c_lambda[k] * invma * r_k
  //   => f_constraint[a] = c_lambda[k] * invma * r_k * ma / dtfsq
  //                       = c_lambda[k] * r_k / dtfsq
  //
  //   Only apply to LOCAL atoms (ghosts don't store permanent forces).
  //   When both atoms of a constraint are local, both get their forces;
  //   when atom b is a ghost, only atom a (which is always local) is updated.
  // -----------------------------------------------------------------

  ev_init(0, vflag);

  for (int k = 0; k < n_constr; k++) {
    const int a = c_atom1[k];    // always local
    const int b = c_atom2[k];    // local or ghost

    const double rx = x[a][0] - x[b][0];
    const double ry = x[a][1] - x[b][1];
    const double rz = x[a][2] - x[b][2];

    const double fac = c_lambda[k] / dtfsq;

    // atom a is always local

    f[a][0] += fac * rx;
    f[a][1] += fac * ry;
    f[a][2] += fac * rz;

    // atom b: only apply if local

    if (b < nlocal) {
      f[b][0] -= fac * rx;
      f[b][1] -= fac * ry;
      f[b][2] -= fac * rz;
    }

    // accumulate virial if requested

    if (evflag) {
      double vir[6];
      vir[0] = fac * rx * rx;
      vir[1] = fac * ry * ry;
      vir[2] = fac * rz * rz;
      vir[3] = fac * rx * ry;
      vir[4] = fac * rx * rz;
      vir[5] = fac * ry * rz;

      int atomlist[2] = {a, b};
      int count = 0;
      if (a < nlocal) atomlist[count++] = a;
      if (b < nlocal) atomlist[count++] = b;
      double fpairlist[] = {fac};
      double dellist[][3] = {{rx, ry, rz}};
      int pairlist[][2] = {{a, b}};
      v_tally(count, atomlist, 2.0, vir, nlocal, 1, pairlist, fpairlist, dellist);
    }
  }
}

/* ----------------------------------------------------------------------
   Print ILVES statistics.
------------------------------------------------------------------------- */

void FixIlves::stats()
{
  if (!output_every) return;

  // reset statistics arrays

  int nb = atom->nbondtypes + 1;
  int na = atom->nangletypes + 1;
  for (int i = 1; i < nb; i++) {
    b_count[i] = 0;
    b_ave[i] = b_max[i] = 0.0;
    b_min[i] = BIG;
  }
  for (int i = 0; i < na; i++) {
    a_count[i] = 0;
    a_ave[i] = a_max[i] = 0.0;
    a_min[i] = BIG;
  }

  // collect per-bond-type bond length statistics

  for (int k = 0; k < n_constr; k++) {
    const int a = c_atom1[k];
    const int b = c_atom2[k];

    if (a >= nlocal && b >= nlocal) continue;    // skip if neither is local

    // determine bond type from per-atom array

    tagint btag = atom->tag[b];
    int btype = 0;
    for (int m = 0; m < atom->num_bond[a]; m++) {
      if (atom->bond_atom[a][m] == btag) {
        btype = std::abs(atom->bond_type[a][m]);
        break;
      }
    }
    if (btype == 0) {    // angle constraint

      utils::print(stderr, "Constraint {} has bond type {}\n", k, btype);
      continue;
    }

    const double dx = atom->x[a][0] - atom->x[b][0];
    const double dy = atom->x[a][1] - atom->x[b][1];
    const double dz = atom->x[a][2] - atom->x[b][2];
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz);

    b_count[btype]++;
    b_ave[btype] += r;
    b_max[btype] = MAX(b_max[btype], r);
    b_min[btype] = MIN(b_min[btype], r);
  }

  // MPI reduction

  MPI_Allreduce(b_count, b_count_all, nb, MPI_LMP_BIGINT, MPI_SUM, world);
  MPI_Allreduce(b_ave, b_ave_all, nb, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(b_max, b_max_all, nb, MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(b_min, b_min_all, nb, MPI_DOUBLE, MPI_MIN, world);

  // print output

  if (comm->me == 0) {
    utils::logmesg(lmp, "ILVES stats (type/ave/delta/count) at step {}\n", update->ntimestep);
    const int width = (int) log10(double(MAX(1, nb))) + 2;
    for (int i = 1; i < nb; i++) {
      if (b_count_all[i] == 0) continue;
      utils::logmesg(lmp, "  Bond:  {:>{}d}   {:<9.6} {:<11.6} {:>8d}\n", i, width,
                     b_ave_all[i] / b_count_all[i], b_max_all[i] - b_min_all[i], b_count_all[i]);
    }
  }

  next_output += output_every;
}

/* ----------------------------------------------------------------------
   Pack xshake positions for ghost communication.
------------------------------------------------------------------------- */

int FixIlves::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int m = 0;
  if (!pbc_flag) {
    for (int i = 0; i < n; i++) {
      const int j = list[i];
      buf[m++] = xshake[j][0];
      buf[m++] = xshake[j][1];
      buf[m++] = xshake[j][2];
    }
  } else {
    double dx, dy, dz;
    if (!domain->triclinic) {
      dx = pbc[0] * domain->xprd;
      dy = pbc[1] * domain->yprd;
      dz = pbc[2] * domain->zprd;
    } else {
      dx = pbc[0] * domain->xprd + pbc[5] * domain->xy + pbc[4] * domain->xz;
      dy = pbc[1] * domain->yprd + pbc[3] * domain->yz;
      dz = pbc[2] * domain->zprd;
    }
    for (int i = 0; i < n; i++) {
      const int j = list[i];
      buf[m++] = xshake[j][0] + dx;
      buf[m++] = xshake[j][1] + dy;
      buf[m++] = xshake[j][2] + dz;
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixIlves::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  const int last = first + n;
  for (int i = first; i < last; i++) {
    xshake[i][0] = buf[m++];
    xshake[i][1] = buf[m++];
    xshake[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   Resize per-atom arrays when atom count grows.
------------------------------------------------------------------------- */

void FixIlves::grow_arrays(int nmax)
{
  memory->grow(ilves_flag, nmax, "ilves:ilves_flag");
  memory->grow(xshake, nmax, 3, "ilves:xshake");
}

/* ----------------------------------------------------------------------
   Copy per-atom data from atom i to atom j (used during atom migration).
------------------------------------------------------------------------- */

void FixIlves::copy_arrays(int i, int j, int /*delflag*/)
{
  ilves_flag[j] = ilves_flag[i];
  xshake[j][0] = xshake[i][0];
  xshake[j][1] = xshake[i][1];
  xshake[j][2] = xshake[i][2];
}

/* ----------------------------------------------------------------------
   Initialize per-atom data for a newly created atom.
------------------------------------------------------------------------- */

void FixIlves::set_arrays(int i)
{
  ilves_flag[i] = 0;
  xshake[i][0] = 0.0;
  xshake[i][1] = 0.0;
  xshake[i][2] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixIlves::reset_dt()
{
  dtv = update->dt;
  dtfsq = update->dt * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   Report degrees of freedom removed by bond constraints.
   Each constrained bond removes 1 DOF from the system.
   Count each bond once (using the c_atom1 ownership invariant).
------------------------------------------------------------------------- */

bigint FixIlves::dof(int igroup)
{
  const int groupbit = group->bitmask[igroup];
  const int *mask = atom->mask;

  bigint n = 0;
  for (int k = 0; k < n_constr; k++) {
    const int a = c_atom1[k];    // always a local atom
    if (mask[a] & groupbit) n++;
  }

  bigint nall;
  MPI_Allreduce(&n, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  return nall;
}

/* ---------------------------------------------------------------------- */

double FixIlves::memory_usage()
{
  double bytes = 0.0;
  // per-atom arrays (ilves_flag + xshake)
  bytes += (double) atom->nmax * sizeof(int);
  bytes += (double) atom->nmax * 3 * sizeof(double);
  // constraint arrays
  bytes += (double) max_constr * (2 * sizeof(int) + 3 * sizeof(double));
  return bytes;
}
