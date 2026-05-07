.. index:: fix ave/specden

fix ave/specden command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID ave/specden Nevery Nrepeat Nfreq value1 value2 ... keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* ave/specden = style name of this fix command
* Nevery = use input values every this many timesteps
* Nrepeat = number of time samples in one spectral density window
* Nfreq = calculate spectral density every this many timesteps
* one or more input values can be listed
* value = c_ID, c_ID[N], f_ID, f_ID[N], v_name

  .. parsed-literal::

       c_ID = global scalar calculated by a compute with ID
       c_ID[I] = Ith component of global vector calculated by a compute with ID, I can include wildcard (see below)
       f_ID = global scalar calculated by a fix with ID
       f_ID[I] = Ith component of global vector calculated by a fix with ID, I can include wildcard (see below)
       v_name = global value calculated by an equal-style variable with name
       v_name[I] = Ith component of a vector-style variable with name, I can include wildcard (see below)

* zero or more keyword/arg pairs may be appended
* keyword = *ave* or *start* or *prefactor* or *file* or *append* or *overwrite* or *title1* or *title2* or *title3*

  .. parsed-literal::

       *ave* args = *one* or *running*
         one = compute fresh spectral density each Nfreq steps
         running = accumulate and average spectral density continuously
       *start* args = Nstart
         Nstart = start accumulating spectral density on this timestep
       *prefactor* args = value
         value = prefactor to scale all spectral density data by
       *file* arg = filename
         filename = name of file to output spectral density data to
       *append* arg = filename
         filename = name of file to append spectral density data to
       *overwrite* arg = none = overwrite output file with only latest output
       *title1* arg = string
         string = text to print as 1st line of output file
       *title2* arg = string
         string = text to print as 2nd line of output file
       *title3* arg = string
         string = text to print as 3rd line of output file

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all ave/specden 5 100 1000 c_myTemp file temp.specden
   fix 1 all ave/specden 1 512 5120 &
             c_thermo_press[1] c_thermo_press[2] c_thermo_press[3] &
             ave running file press.specden
   fix 1 all ave/specden 2 256 5120 c_myVacf ave running overwrite &
             file vacf.specden title1 "# VDOS from VACF"

Description
"""""""""""

.. versionadded:: TBD

Use one or more global scalar values as inputs every few timesteps,
collect a rolling window of :math:`N_\text{repeat}` samples, and
compute the power spectral density (PSD) of each input signal via the
Discrete Fourier Transform (DFT).  Results are written to a file and
are accessible via :doc:`variables <variable>` or other
:doc:`output commands <Howto_output>`.

The group specified with this command is ignored.  However, note that
specified values may represent calculations performed by computes and
fixes which store their own "group" definitions.

Each listed value must be a global scalar produced by a
:doc:`compute <compute>`, :doc:`fix <fix>`, or an equal-style or
vector-style :doc:`variable <variable>`.  Per-atom or local quantities
cannot be used directly; use :doc:`compute reduce <compute_reduce>` to
convert them to a global scalar first.

----------

The :math:`N_\text{every}`, :math:`N_\text{repeat}`, and
:math:`N_\text{freq}` arguments control the sampling and output
schedule.  Input values are sampled every :math:`N_\text{every}`
timesteps and stored in a rolling ring buffer of length
:math:`N_\text{repeat}`.  On every timestep that is a multiple of
:math:`N_\text{freq}`, the DFT of the buffered samples is computed and
output.  :math:`N_\text{freq}` must be a multiple of
:math:`N_\text{every}`, and :math:`N_\text{repeat} \cdot
N_\text{every} \le N_\text{freq}` must hold.

----------

The power spectral density for input signal :math:`x` is computed as:

.. math::

   S(k) = \frac{1}{N} \left| \sum_{n=0}^{N-1} x_n \,
     e^{-2\pi i k n / N} \right|^2

where :math:`N = N_\text{repeat}` is the DFT window length, and
:math:`x_n` are the :math:`N` most recently sampled values.  The
nominal DFT length :math:`N_\text{repeat}` is used for the frequency
axis even if the ring buffer has not yet accumulated
:math:`N_\text{repeat}` samples (missing earlier samples are treated
as zero).

The output contains :math:`N_\text{repeat}/2 + 1` frequency bins
(positive frequencies including DC and Nyquist).  The frequency
corresponding to bin :math:`k` is:

.. math::

   f_k = \frac{k}{N_\text{repeat} \cdot N_\text{every} \cdot \Delta t}

where :math:`\Delta t` is the simulation timestep in LAMMPS time units.

----------

For input values from a compute, fix, or variable, the bracketed
index :math:`I` can be specified using a wildcard asterisk to
effectively specify multiple values.  This takes the form "\*" or
"\*n" or "m\*" or "m\*n".  If :math:`N` is the size of the vector,
then an asterisk with no numeric values means all indices from 1 to
:math:`N`.  A leading asterisk means all indices from 1 to n
(inclusive).  A trailing asterisk means all indices from m to
:math:`N` (inclusive).  A middle asterisk means all indices from m to
n (inclusive).

----------

The *ave* keyword controls how spectral density estimates from
successive :math:`N_\text{freq}` windows are combined.

If *ave* is set to *one* (default), the spectral density is computed
from the :math:`N_\text{repeat}` samples in the current window only,
and the accumulator is reset after each output.  This gives an
independent PSD estimate at each :math:`N_\text{freq}` timestep.

If *ave* is set to *running*, the PSD estimates from all windows since
the start (or since *start* Nstart) are averaged together.  The
output at timestep :math:`t` is the average of all PSD estimates
computed so far.  This is analogous to Welch's method and reduces
statistical noise in the spectral estimate.

----------

The *start* keyword specifies the timestep at which accumulation begins.
Before this timestep the fix does nothing.  Default is timestep 0.

The *prefactor* keyword multiplies all output PSD values by a constant
factor.  Default is 1.0.

The *file* keyword causes output to be written to the named file.  If
the file already exists it is truncated.  If no file keyword is given,
no file is written.

The *append* keyword is the same as *file* except the file is opened
in append mode and new data are appended to any existing content.

The *overwrite* keyword causes the output file to be rewound to the
position after the header before each write, so only the most recent
data are retained in the file.  This keyword requires *ave running*
and is most useful in combination with that setting.

The *title1*, *title2*, and *title3* keywords allow custom header
lines to replace the defaults.  The default headers are:

.. parsed-literal::

   # Spectral density data for fix ID
   # Timestep Number-of-freq-bins
   # Index Frequency Ncount specden(val1) specden(val2) ...

----------

Output info
"""""""""""

This fix produces a global array with :math:`N_\text{repeat}/2 + 1`
rows (frequency bins) and :math:`N_\text{values} + 2` columns.  The
columns are:

* column 1: frequency in cycles per LAMMPS time unit
* column 2: number of DFT windows accumulated (*noutput*)
* columns 3, 4, ...: PSD for each listed input value

The array values can be accessed by any command that uses global
values from a fix, e.g., :doc:`thermo_style custom <thermo_style>`.
See the :doc:`Howto output <Howto_output>` doc page for an overview.

The file (if requested) has the same format as the array.  Each output
block starts with a line giving the current timestep and the number of
frequency bins, followed by one line per frequency bin in the format:

.. parsed-literal::

   index frequency ncount specden1 specden2 ...

Restrictions
""""""""""""

This fix is part of the EXTRA-FIX package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` doc page for more info.

This fix only accepts global scalar inputs.  Per-atom or local
quantities cannot be used.

Related commands
""""""""""""""""

:doc:`fix ave/correlate <fix_ave_correlate>`,
:doc:`fix ave/correlate/long <fix_ave_correlate_long>`,
:doc:`fix ave/histo <fix_ave_histo>`,
:doc:`fix ave/time <fix_ave_time>`

Default
"""""""

The option defaults are ave = one, start = 0, prefactor = 1.0,
no file output, and no overwrite.
