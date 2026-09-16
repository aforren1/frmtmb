# Handing a round to a new session

Written 2026-09-15, at the frmtmb.sample generics release. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`. Read `dev/machine-library.md` BEFORE you run
anything: this machine has lost its R library five times, and the cause
is OPEN rather than settled.

## Where the tree stands

Main carries the Phase 2 release. Four lanes merged: `wt-frailty`,
`wt-eamhier`, `wt-learnhier` and `wt-coh`. The release commit and the
docs rebuild are separate commits, in that order. **The user pushes; no
session pushes for them.**

Versions: frmtmb **0.57.0**, frmtmb.sample **0.5.0**, frmtmb.eam 0.8.1,
frmtmb.learn 0.4.1, frmtmb.spline 0.5.1, frmtmb.coupling 0.3.1, frmtmb.ode
0.4.0, frmtmb.latent 0.3.0. frmtmb.sample floors on core 0.57.0 as a HARD
requirement: it cannot load against 0.56.0, because its `.onLoad()` calls
`frm_install_generics()`, which 0.56.0 does not export.

**Core is not a lane by default.** The user does brms compatibility
there separately. Phase 2.5 is the exception and they asked for it
directly; see below.

### Verified at the release commit

- Suite: 221 files, 13159 assertions, no failures and no errors.
- Gated, NO SKIPS: 23 of 23 files, 2491 assertions. Every row identical
  to the previous release except `test-stan-identity.R` at 71 from 55,
  which is the correlated-block Stan programs item 2.2 added, exactly
  the +16 the total moved by.
- Scale tier and `R CMD check --as-cran` and pkgdown: see the release
  commit message for the numbers.
- One baseline count fell, `frmtmb.spline/test-surface.R` 48 to 47, and
  it was reconciled in review before release: three assertions removed,
  two added.

**Read the comment at the top of each `dev/release/` script before you
change it.** Two flags are forbidden there for reasons this round paid
for again.

## What is next, in order

**Phase 2.5 is done except items 2.5e and 2.5f.** Neither frmtmb nor
frmtmb.sample owns a generic it did not define any more. Both adopt from
the real owners through ONE implementation, core's
`frm_install_generics(pkgname, owners)`, exported on the extension API.
Across 14 load orders, 67 of 67 brms methods lost at base became 0 of 67,
and 60 calls on a brms fit that differed on 12 now differ on 0.

**Next, and the user has decided both:**

- **2.5e**, core's argument spellings. `fitted(fit, re_formula = NA)` is
  swallowed by `...` and silently returns the conditional fit.
- **2.5f**, frmtmb.sample matching brms. `rhat()` and `neff_ratio()` answer
  a DIFFERENT QUESTION from brms, and 10 of 28 methods take arguments in a
  different order, 2 of which answer silently. **Match brms on all of it.**
- Then Phase 2.6, pass brms's own test suite; tmbstan 1.2.1's verification;
  and `frmtmb_control(vectorize = FALSE)`.

**Settled, do not reopen:** the non-generic name collisions with brms stay,
since `::` is sufficient in both load orders. And a gratia older than 0.9.0
loaded before frmtmb.sample now stops it loading; the user accepted that and
the declared floor as sufficient.

## How a round runs here

Lanes in manual git worktrees off main, one item or one coherent group
per lane. A worker writes. A reviewer whose job is to FALSIFY rather
than confirm reads the same worktree against a shared reference build
of the base commit. Punch rounds go back to the worker, and the SAME
reviewer re-checks.

    git worktree add ../frmtmb-wt-<name> -b wt-<name> <sha>

Two punch rounds is the cap, a third only for a blocker. **Build ONE
reference library for the round.** This round did not build one at all:
the previous release's library was already the base commit at the right
versions, so it served as the reference and four lanes read it
read-only. Check whether that is true again before paying for a build.

Consolidation: commit each lane on its branch, merge, set versions,
roxygenise, THEN install (that order), verify each installed NAMESPACE
against the tree, run the tiers, commit the release and the docs
separately, remove the worktrees, prune branches, regenerate
`dev/suite-baseline.tsv`.

**Run `R CMD check` on a QUIET machine.** Its examples-timing NOTE
measures load here, established with a control: a fixed arithmetic
control swings a factor of 4.7 on identical work, and the base build
itself crossed the five second threshold at 6.75 s under load.

## What the user has settled

- **The user pushes.** No session pushes, and no git operation without
  their say so.
- **Nobody uses this package yet.** Break backward compatibility
  freely: ship the refusal, bump, and say plainly in NEWS what stops
  working. No disclosure is owed beyond that, which is why item 2.2's
  correction to what `?rlddm` says needs a NEWS bullet and nothing more.
- **StanHeaders is pinned OUTSIDE `%LOCALAPPDATA%`**, in `pinlib` at
  2.32.10. It survived all five library losses. The user library keeps
  2.39.1 and is not to be downgraded.
- **The R user library stays under `%LOCALAPPDATA%`** by the user's
  decision, knowing it will recur.
- Core is the user's lane except where they ask otherwise.

## What this round is evidence for

Four lanes, four reviewers, six punch rounds. Every row's claim held,
and three of the six rows found something the row was not sent for.

**The thing worth carrying is how often an early number dissolved.**
Two signals that looked like findings did not survive their own
replicates: a coverage of 3 of 6 that became 0.85 at 20, and an `se/sd`
of 0.694 at 14 that became 0.984 at 60. The second is the instructive
one, because the lane first explained it as "a ratio of two spreads is
noisy" and the review showed the real reason: 0.694 sat at the 0.3rd
PERCENTILE of 20,000 random 14-subsets, so it was an unlucky prefix,
and the coverage count had agreed with it at 14. Both instruments were
too few. A dissolved signal is a result about the count and belongs on
the page, not in the bin.

**And a target can be the biased quantity.** Item 2.2 spent most of a
long mechanism hunt asking why an estimate sat 22 percent below a
target, refuting three candidate explanations including two of its own,
before finding that `qlogis(ndt/floor) = log(ndt) - log(m)` identically
and the target was 94 percent nuisance. Nobody had asked whether the
target was right. When an estimator looks biased against one quantity
and unbiased against eight others, ask what is different about the
quantity before asking what is wrong with the estimator.

**Guards improved.** Two failed CLOSED on their first spelling, after
three rounds in which every guard failed open. Both were written with
their inverse case at the same time as the assertion rather than after
it, and one was built deliberately to FAIL when a filed defect is
fixed, with a comment saying to flip it rather than delete it.

**The provenance failures were all one failure.** A count taken from
launches, a guard whose condition had never been observed true, and a
count taken from the plan rather than from the files. In each the thing
being measured was not the thing being reported. The third is the one
that generalizes: a lane's verifier checks files, and that was a
SENTENCE. It was closed structurally, by generating counts into the
document, and the fix found a real error on its first run.

## Worktrees

`wt-frailty`, `wt-eamhier`, `wt-learnhier` and `wt-coh` are merged and
removed, with their evidence committed on main under `dev/`. Create
fresh worktrees off the current main.
