.. index:: fix ilves

fix ilves command
=================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID ilves tol iter N b bond-type ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ilves = style name of this fix command
* tol = accuracy tolerance of ILVES solution
* iter = max # of iterations in each Newton solve
* N = print ILVES statistics every this many timesteps (0 = never)
* one or more *b* bond-type pairs must be appended

  .. parsed-literal::

       *b* values = one or more bond types

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ilves 0.0001 200 10 b 4
   fix 1 polymer ilves 1.0e-5 500 0 b 1 2 3
   fix 1 all ilves 0.0001 200 100 b 2

Description
"""""""""""

Apply bond constraints to specified bond types in the simulation using
the ILVES (Iterative Linear-Velocity-Equation Solver) algorithm
:ref:`(Garcia-Risueno) <Garcia-Risueno>`.  This typically enables a longer
timestep than unconstrained dynamics.

The key advantage of ILVES over :doc:`fix shake <fix_shake>` is that it
handles **arbitrarily large and connected constraint networks**.  In LAMMPS,
fix shake can only be applied to small clusters of atoms (one central atom
bonded to 1, 2, or 3 others, with an optional angle constraint).  This
means long polymer chains -- for example all the C-C backbone bonds in a
protein -- cannot be constrained with fix shake.  Fix ilves removes this
restriction entirely: all bonds of the specified types are constrained,
regardless of how the constraint network is connected or how many atoms it
involves.

.. versionadded:: TBD

**Algorithm overview:**

1. At each timestep the unconstrained positions after the Verlet half-step
   are computed (identical to fix shake).
2. A Newton iteration loop corrects these positions:

   a. Compute the residual for every constrained bond:
      :math:`g_k = \tfrac{1}{2}(|\mathbf{s}_k|^2 - d_k^2)` where
      :math:`\mathbf{s}_k = \mathbf{x}'_a - \mathbf{x}'_b` and :math:`d_k`
      is the equilibrium distance.
   b. Compute a Newton correction using the diagonal entry of the Jacobian
      :math:`A_{kk} = (m_a^{-1}+m_b^{-1})\,\mathbf{r}_k\cdot\mathbf{s}_k`
      where :math:`\mathbf{r}_k = \mathbf{x}_a - \mathbf{x}_b` uses the
      *old* (pre-step) positions.
   c. **Apply all corrections simultaneously** (Jacobi-style update).
      Unlike SHAKE which applies corrections one bond at a time (a
      Gauss-Seidel iteration), ILVES computes all corrections from the
      same positions and then applies them together.  This simultaneous
      update is the fundamental algorithmic difference and gives
      quadratic convergence for tightly coupled constraint networks.
   d. Communicate ghost-atom positions via MPI forward communication.
   e. Repeat until the maximum relative constraint violation
      :math:`\max_k |g_k| / d_k^2 \leq` *tol*, or until *iter* iterations
      have been performed.

3. Constraint forces are computed from the accumulated Lagrange multipliers
   and added to the force array.

**MPI parallelism note:**

When a constraint network spans multiple MPI ranks, the corrections to
atoms on neighboring ranks propagate through ghost-atom communication.
Each Newton iteration propagates information by one rank boundary, so
convergence for chains that cross many ranks may require more iterations
than for purely local chains.  The *iter* parameter should be increased
accordingly for large parallel constraint networks.

**Comparison with fix shake:**

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Feature
     - fix shake
     - fix ilves
   * - Max cluster size
     - 4 atoms
     - unlimited
   * - Connected chains / rings
     - no
     - yes
   * - Update style
     - sequential (Gauss-Seidel)
     - simultaneous (Jacobi/Newton)
   * - Convergence rate (coupled)
     - linear
     - quadratic
   * - Angle constraints
     - yes
     - no (bonds only)

The *b* keyword followed by bond type integers specifies which bond types
to constrain.  At least one *b* specification is required.  The constraint
distance for each bond type is read from the bond potential (same mechanism
as fix shake).

Statistics output
"""""""""""""""""

If *N* is positive, ILVES prints statistics every *N* timesteps.  For each
active bond type the number of bonds, mean length, minimum, and maximum are
reported.

Restrictions
""""""""""""

This fix can only be used with molecular systems (not ``atom_style atomic``).

A bond style must be defined before using this fix (to obtain equilibrium
bond distances).

Only bond constraints are supported.  Angle constraints are not implemented;
use :doc:`fix shake <fix_shake>` for angle constraints.

Only the *b* keyword is recognized for specifying constraints.

This fix may not be used with minimization.

**Related commands:**

:doc:`fix shake <fix_shake>`, :doc:`fix rattle <fix_shake>`,
:doc:`fix restrain <fix_restrain>`

Default
"""""""

None.

----------

.. _Garcia-Risueno:

**(Garcia-Risueno)** P. Garcia-Risueno, L. Perez-Segui, J. Buitrago,
J. M. Cela, "ILVES: An Efficient Parallel Constraint Solver for Molecular
Dynamics", J. Chem. Theory Comput. (2025),
https://doi.org/10.1021/acs.jctc.5c01376
