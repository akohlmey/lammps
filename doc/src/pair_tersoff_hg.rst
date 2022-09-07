.. index:: pair_style tersoff/hg

pair_style tersoff/hg command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style keywords values

* style = *tersoff/hg

Examples
""""""""

.. code-block:: LAMMPS

   pair_style tersoff/hg
   pair_coeff * * Si.tersoff.mod Si Si

Description
"""""""""""

The *tersoff/hg* style computes a bond-order type interatomic potential

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair style is part of the MANYBODY package.  It is only enabled
if LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`

Default
"""""""

none
