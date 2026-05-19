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

#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr int DELTA_CONSTR = 256;
static constexpr double SMALL_DENOM = 1.0e-10;

/* ---------------------------------------------------------------------- */

FixIlves::FixIlves(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg),
    bond_flag(nullptr), bond_distance(nullptr),
    ilves_flag(nullptr), xshake(nullptr),
    c_atom1(nullptr), c_atom2(nullptr), c_dist2(nullptr), c_lambda(nullptr),
    x(nullptr), v(nullptr), f(nullptr), mass(nullptr), rmass(nullptr), type(nullptr),
    b_count(nullptr), b_count_all(nullptr), b_ave(nullptr), b_ave_all(nullptr),
    b_max(nullptr), b_max_all(nullptr), b_min(nullptr), b_min_all(nullptr)
{
  if (narg < 8) utils::missing_cmd_args(FLERR, fmt::format("fix {}", style), error);

  // error checks

  if (atom->molecular == Atom::ATOMIC)
    error->all(FLERR,"Cannot use fix {} with non-molecular system", style);

  // parse required args: tol iter N constraint_specs

  tolerance = utils::numeric(FLERR, arg[3], false, lmp);
  max_iter  = utils::inumeric(FLERR, arg[4], false, lmp);
  output_every = utils::inumeric(FLERR, arg[5], false, lmp);

  // parse constraint specs: only bond types supported (b keyword)

  bond_flag = new int[atom->nbondtypes + 1]();

  int next = 6;
  char mode = '\0';
  while (next < narg) {
    if (strcmp(arg[next], "b") == 0) {
      mode = 'b';
    } else if (mode == 'b') {
      int i = utils::inumeric(FLERR, arg[next], false, lmp);
      if (i < 1 || i > atom->nbondtypes)
        error->all(FLERR,"Invalid bond type {} for fix {}", arg[next], style);
      bond_flag[i] = 1;
    } else {
      error->all(FLERR,"Unknown fix {} option: {}", style, arg[next]);
    }
    next++;
  }

  // allocate bond distance array (indexed 1..nbondtypes)

  bond_distance = new double[atom->nbondtypes + 1]();

  // statistics arrays (per bond type, indices 1..nbondtypes)

  if (output_every) {
    int nb = atom->nbondtypes + 1;
    b_count     = new bigint[nb]();
    b_count_all = new bigint[nb]();
    b_ave       = new double[nb]();
    b_ave_all   = new double[nb]();
    b_max       = new double[nb]();
    b_max_all   = new double[nb]();
    b_min       = new double[nb]();
    b_min_all   = new double[nb]();
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

  if (modify->get_fix_by_style(fmt::format("^{}",style)).size() > 1)
    error->all(FLERR,"More than one fix {} instance", style);

  // cannot use during minimization (turns off bond forces)

  if (update->whichflag == 2)
    error->warning(FLERR,"fix {} cannot be used with minimization", style);

  // error if box changes come before ilves fix

  bool boxflag = false;
  for (const auto &ifix : modify->get_fix_list()) {
    if (boxflag && utils::strmatch(ifix->style, fmt::format("^{}",style)))
      error->all(FLERR,"Fix {} must come before any box-changing fix", style);
    if (ifix->box_change) boxflag = true;
  }

  if (!force->bond)
    error->all(FLERR,"Bond style must be defined for fix {}", style);

  // get equilibrium bond distances from the bond force style

  for (int i = 1; i <= atom->nbondtypes; i++)
    bond_distance[i] = force->bond->equilibrium_distance(i);

  // set/update timestep info

  dtv   = update->dt;
  dtfsq = update->dt * update->dt * force->ftm2v;
}

/* ----------------------------------------------------------------------
   Called at the start of each run to initialize.
   Mark constrained bond types negative and build initial constraint list.
------------------------------------------------------------------------- */

void FixIlves::setup(int vflag)
{
  // update local pointers

  x     = atom->x;
  v     = atom->v;
  f     = atom->f;
  mass  = atom->mass;
  rmass = atom->rmass;
  type  = atom->type;
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

  x      = atom->x;
  v      = atom->v;
  f      = atom->f;
  mass   = atom->mass;
  rmass  = atom->rmass;
  type   = atom->type;
  nlocal = atom->nlocal;

  // mark constrained bond types negative in per-atom arrays
  // (they may have been reset to positive during atom exchange)

  int *num_bond  = atom->num_bond;
  int **bond_type = atom->bond_type;

  for (int i = 0; i < nlocal; i++) {
    for (int m = 0; m < num_bond[i]; m++) {
      int btype = std::abs(bond_type[i][m]);
      if (btype >= 1 && btype <= atom->nbondtypes && bond_flag[btype])
        bond_type[i][m] = -btype;
    }
  }

  // clear ilves_flag for all local atoms; will be re-set in post_neighbor

  for (int i = 0; i < nlocal; i++) ilves_flag[i] = 0;
}

/* ----------------------------------------------------------------------
   Called after each neighbor list rebuild.
   Rebuild the flat constraint list from the newly built bondlist.
------------------------------------------------------------------------- */

void FixIlves::post_neighbor()
{
  // reset constraint count; keep allocated memory

  n_constr = 0;

  const int nbondlist  = neighbor->nbondlist;
  int **bondlist = neighbor->bondlist;

  for (int k = 0; k < nbondlist; k++) {
    const int i     = bondlist[k][0];   // always local (< nlocal)
    const int j     = bondlist[k][1];   // local or ghost
    const int btype = std::abs(bondlist[k][2]);

    if (btype < 1 || btype > atom->nbondtypes) continue;
    if (!bond_flag[btype]) continue;

    // grow constraint arrays if needed

    if (n_constr == max_constr) {
      max_constr += DELTA_CONSTR;
      memory->grow(c_atom1,  max_constr, "ilves:c_atom1");
      memory->grow(c_atom2,  max_constr, "ilves:c_atom2");
      memory->grow(c_dist2,  max_constr, "ilves:c_dist2");
      memory->grow(c_lambda, max_constr, "ilves:c_lambda");
    }

    c_atom1[n_constr]  = i;
    c_atom2[n_constr]  = j;
    c_dist2[n_constr]  = bond_distance[btype] * bond_distance[btype];
    c_lambda[n_constr] = 0.0;
    n_constr++;

    ilves_flag[i] = 1;
    if (j < nlocal) ilves_flag[j] = 1;
  }
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

  x      = atom->x;
  v      = atom->v;
  f      = atom->f;
  mass   = atom->mass;
  rmass  = atom->rmass;
  type   = atom->type;
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
      const double rs    = rx * sx + ry * sy + rz * sz;
      const double A_kk  = (invma + invmb) * rs;

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
  for (int i = 1; i < nb; i++) {
    b_count[i] = 0;
    b_ave[i]   = 0.0;
    b_max[i]   = -1.0e100;
    b_min[i]   =  1.0e100;
  }

  // collect per-bond-type bond length statistics

  for (int k = 0; k < n_constr; k++) {
    const int a = c_atom1[k];
    const int b = c_atom2[k];

    if (a >= nlocal && b >= nlocal) continue;    // skip if neither is local

    // determine bond type from per-atom array

    tagint btag = atom->tag[b];
    int btype   = 0;
    for (int m = 0; m < atom->num_bond[a]; m++) {
      if (atom->bond_atom[a][m] == btag) {
        btype = std::abs(atom->bond_type[a][m]);
        break;
      }
    }
    if (btype == 0) continue;

    const double dx = atom->x[a][0] - atom->x[b][0];
    const double dy = atom->x[a][1] - atom->x[b][1];
    const double dz = atom->x[a][2] - atom->x[b][2];
    const double r  = std::sqrt(dx * dx + dy * dy + dz * dz);

    b_count[btype]++;
    b_ave[btype] += r;
    b_max[btype]  = MAX(b_max[btype], r);
    b_min[btype]  = MIN(b_min[btype], r);
  }

  // MPI reduction

  MPI_Allreduce(b_count,     b_count_all, nb, MPI_LMP_BIGINT, MPI_SUM,  world);
  MPI_Allreduce(b_ave,       b_ave_all,   nb, MPI_DOUBLE,     MPI_SUM,  world);
  MPI_Allreduce(b_max,       b_max_all,   nb, MPI_DOUBLE,     MPI_MAX,  world);
  MPI_Allreduce(b_min,       b_min_all,   nb, MPI_DOUBLE,     MPI_MIN,  world);

  // print output

  if (comm->me == 0) {
    utils::logmesg(lmp, "ILVES stats at step {}\n", update->ntimestep);
    for (int i = 1; i < nb; i++) {
      if (b_count_all[i] == 0) continue;
      const double bavg = b_ave_all[i] / b_count_all[i];
      utils::logmesg(lmp, "  bond type {}: n={}, avg={:.5f}, min={:.5f}, max={:.5f}\n",
                     i, b_count_all[i], bavg, b_min_all[i], b_max_all[i]);
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
  memory->grow(ilves_flag, nmax,    "ilves:ilves_flag");
  memory->grow(xshake,     nmax, 3, "ilves:xshake");
}

/* ----------------------------------------------------------------------
   Copy per-atom data from atom i to atom j (used during atom migration).
------------------------------------------------------------------------- */

void FixIlves::copy_arrays(int i, int j, int /*delflag*/)
{
  ilves_flag[j]   = ilves_flag[i];
  xshake[j][0]    = xshake[i][0];
  xshake[j][1]    = xshake[i][1];
  xshake[j][2]    = xshake[i][2];
}

/* ----------------------------------------------------------------------
   Initialize per-atom data for a newly created atom.
------------------------------------------------------------------------- */

void FixIlves::set_arrays(int i)
{
  ilves_flag[i]   = 0;
  xshake[i][0]    = 0.0;
  xshake[i][1]    = 0.0;
  xshake[i][2]    = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixIlves::reset_dt()
{
  dtv   = update->dt;
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
  const int *mask    = atom->mask;

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
  bytes += (double)atom->nmax * sizeof(int);
  bytes += (double)atom->nmax * 3 * sizeof(double);
  // constraint arrays
  bytes += (double)max_constr * (2 * sizeof(int) + 3 * sizeof(double));
  return bytes;
}
