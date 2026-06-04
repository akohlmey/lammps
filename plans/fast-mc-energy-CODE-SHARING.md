# Reducing duplication in the MC package: base class vs mixin vs composition

*Design discussion companion to fast-mc-energy-EXPLORATION-RESULTS.md.*

## Problem

`fix gcmc`, `fix widom`, `fix atom_swap`, `fix mol_swap`,
`fix charge_regulation` (and partly `fix sgcmc`) already duplicate a lot:
`energy_full()`, the `full_energy` decision, the `neigh_modify exclude` setup,
the `thermo_pe` hookup, and the `ENERGY_GLOBAL|ENERGY_ONLY` plumbing. The
proposed fast path (cell list, local Delta-energy, caches, `validate`, timing)
multiplies that duplication. Goal: share it so bugs are fixed once.

## The LAMMPS-specific constraint: the `Pointers` diamond

`Fix`/`Compute`/`Pair`/... all derive from `Pointers` (holds `lmp` and the
`atom`/`force`/`comm`/... references). A naive `class MCEnergy : public Pointers`
mixed into `class FixGCMC : public Fix, public MCEnergy` yields **two
`Pointers` subobjects**. Changing `Fix` to inherit `Pointers` virtually is too
invasive. LAMMPS already navigates this two ways:
- **stateless interface mixin** - `ElectrodeKSpace` (pure virtuals, no data, no
  `Pointers`): no diamond.
- **stateful mixin that owns its own `lmp`** - `ThrOMP` (does *not* derive from
  `Pointers`; stores `LAMMPS *lmp`, reaches subsystems via `lmp->atom`): no
  diamond. Used as `class PairLJCutOMP : public PairLJCut, public ThrOMP`.

So every option below is feasible; the choice is coupling vs testability vs
migration risk.

## What is shared vs divergent (drives the answer)

- Shared: an **energy engine + setup** (`energy_full`, full_energy decision,
  exclusion group, kspace ENERGY_ONLY plumbing, the fast-path machinery).
- Divergent: the **trial loop and acceptance rule** (grand-canonical
  insert/delete/move; Widom no-accept; type/identity swaps with fixed positions;
  protonation reaction ensemble). `sgcmc` is an EAM special case.

Key point: the shared part is *computation*, not a common *control flow*. That
is the opposite of the `FixNH` family (shared integrate loop, hooks for
thermostat specifics) - which is exactly where a base class earns its keep. Here
the skeleton differs per fix and the engine is shared, which favors a composed
engine over a template-method base.

## Options

**A. `FixMC` base class** (single inheritance; `FixNH` precedent
`FixNVT/FixNPT : public FixNH : public Fix`).
+ shared state/lifecycle in one place; natural is-a; virtual hooks; shared
  option parsing / setmask / restart conventions.
- single-inheritance lock-in; god-base risk; converting working fixes changes
  their base (restart/behavior compat); `sgcmc` must be carved around; tight
  coupling (base change ripples to all). Best when a trial-loop skeleton is
  shared - which these mostly do not.

**B. MI mixin** (`ThrOMP` form: `class FixGCMC : public Fix, public MCEnergy`,
mixin owns its `lmp`).
+ orthogonal add-on without touching the primary base; idiomatic; members
  directly visible in the derived (no `helper->`); diamond-free.
- two sources of truth for subsystem access (inherited `atom` vs `lmp->atom`);
  MI lookup/ctor-order friction; harder to unit-test in isolation; ambiguity if
  two energy-like mixins are ever combined. Best for exactly one capability when
  avoiding delegation boilerplate matters.

**C. Composition** (each fix owns an `MCEnergy *`, a `Pointers` subclass built
with `lmp`; the plan's original proposal).
+ cleanest separation; no inheritance/diamond; **independently unit-testable**
  (build an `MCEnergy`, check `delta_insert` vs brute force) - decisive given
  this engine already hit 3+ subtle bugs (tail/volume, stale cache, move-Delta);
  decoupled; **additive, low-risk migration** (add a member + delegate, no base
  change, no restart impact); a fix can own one or several helpers.
- delegation boilerplate; trial data (molecule template, params) passed in.
  Best when the shared thing is a computational engine and correctness matters.

**D. Hybrid** - thin `FixMC` base for the Fix-lifecycle commonalities
(full_energy detection, exclusion group, `c_pe`, shared keywords, `post_run`
timing) + composed `MCEnergy` helper for the engine. Inheritance for the family,
composition for the reusable core. Usually the right real-world answer.

**E. CRTP / static mixin** - zero-overhead dispatch is irrelevant (energy
compute dominates), non-idiomatic in LAMMPS, bloats. Dismiss.

## Comparison

| | fix-once | coupling | unit-testable engine | migration risk | MI/diamond |
|---|---|---|---|---|---|
| A base class | yes | high | no | high (changes base, restart) | n/a |
| B MI mixin | yes | medium | hard | medium (changes each fix's bases) | ok if ThrOMP-style |
| C composition | yes | low | **yes** | **low (additive)** | none |
| D hybrid | yes | low-med | yes | low-med | none |

## Recommendation

Lead with **C (composed `MCEnergy` helper)**: the duplication *and* the bugs
concentrate in the energy engine, and composition makes it a unit under test
with a permanent `validate` harness. Prove it on `fix gcmc`, then migrate
`widom`/`atom_swap`/`mol_swap`/`charge_regulation` one at a time with no base
change. Add a thin **`FixMC` base (-> D)** later only if Fix-level setup
duplication proves annoying. Do **not** lead with a base-class refactor, and do
**not** use the naive `MCEnergy : public Pointers` mixin (the one shape that
hits the diamond) - if MI is chosen it must follow the ThrOMP own-`lmp` form.
