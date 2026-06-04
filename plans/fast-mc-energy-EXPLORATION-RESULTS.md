# Fast Monte Carlo trial energies in `fix gcmc`: exploration results

*Prototype + benchmarks, single-node CPU, for discussion. Not a finished feature.*

## Motivation

With a long-range solver, a tail correction, or a non-`single()` pair style,
`fix gcmc` runs in `full_energy` mode: **every trial move** re-decomposes the
domain, rebuilds ghosts and the neighbor list, and recomputes the whole-system
energy. For inserting/deleting one small molecule this is dominated by the
neighbor-list rebuild (originally measured ~93% of a trial). The exploration
asks how far a *local* trial-energy evaluation can go.

## What was prototyped (opt-in `fast` keyword, default off)

1. **`ENERGY_ONLY` in the CPU kspace styles** (ewald, ewald/disp, pppm ik/ad,
   pppm/cg, pppm/stagger, pppm/disp; pppm/tip4p inherits). When only the global
   energy is requested, the force-producing work (inverse FFTs / E-field /
   Ewald force loop) is skipped; energy + virial come from the structure
   factors / pre-FFT tally. Bit-for-bit vs the full path (unit test).
2. **Local trial energy change** for molecule **insertion and deletion**: the
   moved molecule's real-space pair energy via `Pair::single()` (no neighbor
   rebuild), a reciprocal-energy recompute via (1), the analytic tail delta, and
   a self-calibrated internal constant. No `comm->exchange()` inside a trial.
3. **Private cell list** so the local pair energy is O(N_local) not O(N_all).
4. **`validate` mode** (run both paths, report max relative disagreement) and a
   per-segment wall-time breakdown printed at end of run.

Correct and committed; **molecule translation/rotation are *not* included**
(an unresolved multi-atom/periodic-boundary discrepancy in the move delta), and
the local path requires **pair cutoff <= 1/2 the shortest box length**
(minimum image; otherwise it falls back to full energy with a warning).

## Benchmarks (1 MPI rank, `run 2000`, molecule exchange only)

| System | atoms | kspace | `fast off` | `fast yes` | overall | /trial | validate |
|---|---|---|---|---|---|---|---|
| N2 in CHO framework (3x2x2) | 100,536 | pppm | 376.6 s | 55.0 s | **6.85x** | ~55x | 3.8e-12 |
| CO2 fluid (5x5x5) | 3,000 | ewald | 75.5 s | 34.2 s | **2.2x** | ~13x | 1.0e-12 |
| H2O fluid + shake (5x5x5) | 3,000 | ewald | 84.8 s | 36.1 s | **2.35x** | ~13x | 3.7e-9 |

Where the GCMC time goes (`fast yes`, fraction of the fix's own time):

| System | full-energy baseline | **kspace recompute** | local pair (cell list) | ghost |
|---|---|---|---|---|
| N2 100k (pppm) | 56.6% | 38.7% | **0.9%** | 3.9% |
| CO2 3k (ewald) | 15.6% | **80.0%** | 1.9% | 2.5% |
| H2O 3k (ewald) | 14.1% | **81.3%** | 2.0% | 2.5% |

## Findings

- The local path is **correct** on all real systems (1e-9..1e-12 relative) and
  reproduces the full-path loading bit-for-bit; it removes the per-trial
  neighbor rebuild, giving **~13x-55x per trial** and **2.2x-6.85x overall**.
- **Speedup scales with the eliminated rebuild cost** (system size): large
  framework -> 6.85x; small fluids -> ~2.2x (there the run is MD-bound).
- The **cell list works**: local pair energy is <=2% everywhere.
- The **dominant remaining MC cost is the per-trial reciprocal recompute** -
  for Ewald it is ~80%. A full O(N.k) recompute is redone every trial even
  though one molecule changes the reciprocal energy by an O(1) amount.
- The **per-sweep full-energy baseline** is the other redundant cost (largest,
  56.6%, for the big PPPM system): the MD just computed that energy the step
  before; it is recomputed because the MD moved the gas in between.
- The shipped `examples/mc` GCMC inputs do **not** exercise the path as-is
  (`in.gcmc.lj` is LJ-only; `in.gcmc.co2`/`in.gcmc.h2o` use a 10 A box with a
  14 A cutoff -> cutoff > 1/2 box -> falls back). They only engage once enlarged
  by `replicate` (as benchmarked above).

## Architectural directions suggested by the data

1. **Incremental kspace per trial** is the highest-leverage change (the ~80%
   Ewald cost). PPPM cannot be made incremental cheaply (global FFT), but
   **Ewald's structure factors `S(k)` can**: cache `S(k)`, update one molecule
   in `O(k_count)` (`delta S`, `delta E`), commit on accept. This is the
   standard dedicated-GCMC approach.

2. **Hybrid PPPM-MD + Ewald-increments.** Keep PPPM for the MD and the periodic
   full energy (good O(N log N) scaling); use Ewald *only* for the per-trial
   increments. Cost model: per-trial increment `O(M.k_count)` replaces an
   `O(N log N)` FFT; the gating cost is one `S(k)` build per sweep
   `O(N.k_count)`, amortized over the sweep's trials. `k_count ~ g_ewald^3 V`,
   so this is cheapest exactly when the **real-space cutoff is large** (coarse
   reciprocal mesh) - i.e. the framework/adsorption regime where the per-trial
   PPPM recompute is most wasteful. Caveat: PPPM and Ewald are different
   approximations agreeing only to ~the chosen accuracy, so the MC effectively
   samples the Ewald Hamiltonian while the MD uses PPPM - a deliberate, bounded,
   *verifiable* (via `validate`) inconsistency of order the kspace accuracy.
   Requires a shared `g_ewald` and consistent self/`qsum`/boundary terms.

3. **Reuse the MD energy** instead of the per-sweep full re-evaluation.

4. **Amortize ghosts + cell list to accept-only** (rebuild on accept, reuse
   across rejected trials) to remove the per-trial O(N) overhead that remains.

5. The **minimum-image restriction** (cutoff <= 1/2 box) excludes single-unit-
   cell studies with large cutoffs; a self-image-aware local energy would lift
   it. The **molecule-move delta** needs a component-level (evdwl/ecoul/ebond/
   elong) breakdown to localize its bug before translation/rotation can ship.

## Reproduction

Branch `kspace-energy-only`. Inputs: `build-test/gcmc_bench.inp` (framework;
data files alongside) and `replicate`d copies of `examples/mc/in.gcmc.{co2,h2o}`
with `fast ${fm}` appended. Run with `-var fm off|yes|validate`; the breakdown
and validation summary print at end of run.
