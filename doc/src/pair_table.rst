.. index:: pair_style table
.. index:: pair_style table/gpu
.. index:: pair_style table/kk
.. index:: pair_style table/omp
.. index:: pair_style table/mod

pair_style table command
========================

Accelerator Variants: *table/gpu*, *table/kk*, *table/omp*

pair_style table/mod command
============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style table style N keyword ...

* style = *lookup* or *linear* or *spline* or *bitmap* = method of table evaluation
* N = use N values in *lookup*, *linear*, *spline* tables
* N = use 2\^N values in *bitmap* tables
* zero or more keywords may be appended
* keyword = *ewald* or *pppm* or *msm* or *dispersion* or *tip4p* or *rlinear* or *rsquared*
   *ewald* = flag compatibilty with kspace style ewald (and compatible)
   *pppm* = flag compatibility with kspace style pppm (and compatible)
   *msm* = flag compatibility with kspace style msm (and compatible)
   *dispersion* = flag compatibility with a dispersion kspace style
   *tip4p* = flag compatibility with a tip4p kspace style
   *rlinear* = use internal tables with linear spacing in r
   *rsquared* = use internal tables with squared spacing in r

Examples
""""""""

.. code-block:: LAMMPS

   pair_style table lookup 1000 rsquared
   pair_style table linear 2000 pppm
   pair_style table spline 500 rlinear
   pair_style table bitmap 12 rsquared
   pair_coeff * 3 morse.table ENTRY1
   pair_coeff * 3 morse.table ENTRY1 7.0

Description
"""""""""""

Pair style *table* uses tabulated data to compute pair-wise forces and
energies instead of analytical functions.  When performing dynamics or
minimization internal tables are used to evaluate energy and forces for
pairwise interactions between particles.  The internal tabulated data is
imported from table files where the energy and force are listed as a
function of distance into internal tables.  The pair style supports
multiple "evaluation styles" (*lookup*, *linear*, *spline*, and
*bitmap*) that differ in how the values from the internal tables are
looked up and evaluated, resulting in different computational effort and
accuracy.  The internal tables use equidistant spacing in either *r*
(*rlinear*) or :math:`r^2` (*rsquared*) that also result in different
computational cost and accuracy across the range of the tabulated data.
A detailed discussion of these options follows below.

.. versionchanged:: TBD

The internal tables can be created with linear spacing for improved
accuracy and consistency at small distances at a slightly increased
compatational cost.  Previously only squared spacing was available.

The internal tables are in most cases created as a pre-computation by
fitting cubic splines to the values in the provided table files and
interpolating energy and force values at each of requested *N* distances
according to the *rlinear* or *rsquared* spacing style.  Only if the
data in the table file exactly matches the requirements of the requested
evaluation and spacing style, it is copied without interpolation.

During a simulation, the internal table data is used used to evaluate
energy and force values as needed for each pair of particles separated
by a distance :math:`r_{ij}`.  This evaluation is done in one of 4
styles: *lookup*, *linear*, *spline*, or *bitmap*\ .

For the *lookup* style, the distance :math:`r_{ij}` is used to find the
closest table entry which is then used directly.  This is fast but incurs
a significant error depending on the granularity of the data.  The *lookup*
evaluation style is not available for *rlinear* spacing.

For the *linear* style, the distance :math:`r_{ij}` is used to find the
2 surrounding table values from which an energy or force is computed by
linear interpolation.  This is significantly more accurate than *lookup*
at only a small additional computational cost.

For the *spline* style, cubic spline coefficients are pre-computed and
stored for each of the *N* values in the internal table, one set of
splines for energy, another for force.  Note that these splines are
fitted to the **internal** tabulated data which may have already be
determined from a spline interpolation of the orignal tabulated data
read from a file (see above).  The distance :math:`r_{ij}` is used to
find the appropriate set of spline coefficients which are used to
evaluate a cubic polynomial which computes the energy or force.
This is the most accurate evaluation style available but also the one
with the largest computational cost.

The *bitmap* style is similar to *linear*, only that the lookup is done
using a fast bit-mapping technique due to :ref:`(Wolff) <Wolff2>`, which
requires the number of tabulated data points in the internal table to be
a power of 2.  Thus a value of *N* will result in :math:`2^N` data
points in the internal table.  The *bitmap* style is not available for
*rlinear* spacing.

The performance and accuracy of the tabulation depends on multiple
factors like the number of values, the evaluation style, the spacing
style and the precision of the imported data.  A table with more values
is usually more accurate, but less efficient since it requires more
memory.  A tabulation with *rlinear* (internal) spacing is more accurate
at small distances (where the curvature of typical atomic potential
functions tends to be larger) but incurs additional computational cost
compared to the *rsquared* style.  Spline interpolation is usually more
accurate than a linear interpolation, and that more accurate than
lookup.  How large the difference is depends on the curvature of the
potential.  The added cost of a spline interpolation may be offset by
needing fewer tabulation values to represent the potential with the
desired accuracy.  The optimal choice typically requires some
experimentation and testing.

The following coefficients must be defined for each pair of atoms
types via the :doc:`pair_coeff <pair_coeff>` command as in the examples
above.

* filename
* keyword
* cutoff (distance units)

The filename specifies a file containing tabulated energy and force
values.  The keyword specifies a section of the file.  The cutoff is an
optional coefficient.  If not specified, the outer cutoff in the table
itself (see below) will be used to build the internal table for forces
and energy.  If the custoff is specified, only data from the file to the
cutoff is used to create the internal table.  The format of this file is
described below.

If your tabulated potential(s) are designed to be used as the
short-range part of one of the long-range solvers specified by the
:doc:`kspace_style <kspace_style>` command, then you must use one or
more of the optional keywords listed above for the pair_style command.
These are *ewald* or *pppm* or *msm* or *dispersion* or *tip4p*\ .  This
is so LAMMPS can insure the short-range potential and long-range solver
are compatible with each other, as it does for other short-range pair
styles, such as :doc:`pair_style lj/cut/coul/long <pair_lj_cut_coul>`.
Note that it is up to you to insure the tabulated values for each pair
of atom types has the correct functional form to be compatible with the
matching long-range solver.

----------

Here are some guidelines for using the pair_style table command to
best effect:

* Vary the number of table points; you may need to use more than you think
  to get good resolution.  It is strongly recommended to first experiment
  with tables for a Lennard-Jones or Morse potential, for which LAMMPS has
  analytical versions, to get a sense for what kind of accuracy is possible
  for different tabulation settings and tabulated data.  Another option
  would be to use :doc:`pair style python <pair_python>` to define an
  analytic version of the desired potential, if the analytic version is
  otherwise not available in LAMMPS.
* Use the :doc:`pair_write <pair_write>` command to produce a file with
  data of the final evaluated potential so it can be compared to the
  original table data.  During the import of data from a file there may
  be artifacts resulting from the spline interpolation used to create
  the internal tabulated data.  This typically happens at small
  distances and is more likely with *rsquared* spacing than with
  *rlinear* spacing.
* Start with the *linear* style; it is the evaluation style least likely to cause problems.
* The additional interpolation and its artifacts can be avoided, if the
  data in the imported table matches the choice of cutoff, number of
  tabulation values, evaluation style, and spacing style.
* Make sure that your tabulated forces and tabulated energies are
  consistent (dE/dr = -F) over the entire range of r values.  LAMMPS
  will warn if this is not the case.  This test is approximate though
  and may produce false positives if the imported data has low precision,
  uses a different spacing or has much fewer data points.
* Use as large an inner cutoff as possible.  This avoids fitting splines
  to very steep parts of the potential.  Use *rlinear* spacing to reduce
  artifacts, if necessary.

----------

Suitable tables in the correct format for use with these pair styles can
be created by LAMMPS itself using the :doc:`pair_write <pair_write>`
command.  In combination with the :doc:`pair style python <pair_python>`
this can be a powerful mechanism to implement and test tables for use
with LAMMPS.  Another option to generate tables is the Python code in
the ``tools/tabulate`` folder of the LAMMPS source code distribution.

The format of a tabulated file has an (optional) header followed by a
series of one or more sections, defined as follows (without the
parenthesized comments). The header must start with a `#` character
and the DATE: and UNITS: tags will be parsed and used.  Using the
UNITS: tag is strongly encouraged:

.. parsed-literal::

   # DATE: 2020-06-10  UNITS: real  CONTRIBUTOR: ... (header line)
   # Morse potential for Fe   (one or more comment or blank lines)

   MORSE_FE                   (keyword is first text on line)
   N 500 R 1.0 10.0           (N, R, RSQ, BITMAP, FPRIME parameters)
                              (blank)
   1 1.0 25.5 102.34          (index, r, energy, force)
   2 1.02 23.4 98.5
   ...
   500 10.0 0.001 0.003

A section begins with a non-blank line whose first character is not a
"#"; blank lines or lines starting with "#" can be used as comments
between sections.  The first line begins with a keyword which
identifies the section.  The line can contain additional text, but the
initial text must match the argument specified in the pair_coeff
command.  The next line lists (in any order) one or more parameters
for the table.  Each parameter is a keyword followed by one or more
numeric values.

The parameter "N" is required and its value is the number of table
entries that follow.  Note that this may be different than the *N*
specified in the :doc:`pair_style table <pair_style>` command.  Let
Ntable = *N* in the pair_style command, and Nfile = "N" in the tabulated
file.  What LAMMPS does is a preliminary interpolation by creating
splines using the Nfile tabulated values as nodal points.  It uses these
to interpolate energy and force values at Ntable different points to
create the internal tables.  The resulting internal tables of length
Ntable are then used as described above, when computing energy and force
for individual pair distances.  This means that if you want the
interpolation tables of length Ntable to match exactly what is in the
tabulated file (with effectively no preliminary interpolation), you
should set Ntable = Nfile, and use the "RSQ" or "BITMAP" parameter with
the *rsquared* spacing setting.  This is because the internal table
abscissa is RSQ (separation distance squared) in this case, for
efficient lookup.  If the *rlinear* spacing is used, then the "R"
parameter is needed to avoid interpolation during import of the table.

All other parameters are optional.  If "R" or "RSQ" or "BITMAP" does not
appear, then the distances in each line of the table are used as-is to
perform spline interpolation.  In this case, the table values can be
spaced in *r* uniformly or however you wish to position table values in
regions of large gradients.

If used, the parameters "R" or "RSQ" are followed by 2 values *rlo* and
*rhi*\ .  If specified, the distance associated with each energy and
force value is computed from these 2 values (at full floating point
precision), rather than using the (potentially lower precision) value
listed in each line of the table.  The distance values in the table file
are ignored in this case.  For "R", distances uniformly spaced between
*rlo* and *rhi* are computed; for "RSQ", squared distances uniformly
spaced between *rlo\*rlo* and *rhi\*rhi* are computed.

.. note::

   If you use "R" or "RSQ", the tabulated distance values in the file
   are effectively ignored, and replaced by new values as described in
   the previous paragraph.  If the distance value in the table is not
   very close to the new value (i.e. round-off difference), then you
   will be assigning energy/force values to a different distance, which
   is probably not what you want.  LAMMPS will warn if this is
   occurring.

If used, the parameter "BITMAP" is also followed by 2 values *rlo* and
*rhi*\ .  These values, along with the "N" value determine the ordering
of the N lines that follow and what distance is associated with each.
This ordering is complex, so it is not documented here, since this
file is typically produced by the :doc:`pair_write <pair_write>` command
with its *bitmap* option.  When the table is in BITMAP format, the "N"
parameter in the file must be equal to 2\^M where M is the value
specified in the pair_style command.  Also, a cutoff parameter cannot
be used as an optional third argument in the pair_coeff command; the
entire table extent as specified in the file must be used.

If used, the parameter "FPRIME" is followed by 2 values *fplo* and
*fphi* which are the derivative of the force at the innermost and
outermost distances listed in the table.  These values are needed by the
spline construction routines.  If not specified by the "FPRIME"
parameter, they are estimated (less accurately) by the first 2 and last
2 force values in the table.  This parameter is not used by BITMAP
tables.

Following a blank line, the next N lines list the tabulated values.  On
each line, the first value is the index from 1 to N, the second value is
r (in distance units), the third value is the energy (in energy units),
and the fourth is the force (in force units).  The r values must
increase from one line to the next (unless the BITMAP parameter is
specified).

Note that one file can contain many sections, each with a tabulated
potential.  LAMMPS reads the file section by section until it finds one
that matches the specified keyword.

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support mixing.  Thus, coefficients for all
I,J pairs must be specified explicitly.

The :doc:`pair_modify <pair_modify>` shift, table, and tail options are
not relevant for this pair style.

This pair style writes the settings for the "pair_style table" command
to :doc:`binary restart files <restart>`, so a pair_style command does
not need to be specified in an input script that reads a restart file.
The coefficient information, however, is not stored in the restart file,
since it is tabulated in the potential files.  Thus, pair_coeff commands
do need to be specified in the restart input script.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""
none

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`pair_write <pair_write>`

Default
"""""""

*rsquared*

----------

.. _Wolff2:

**(Wolff)** Wolff and Rudd, Comp Phys Comm, 120, 200-32 (1999).
