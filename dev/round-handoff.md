# Handing a round to a new session

Written 2026-09-09, at frmtmb 0.55.1. Read this, then
`dev/extension-gaps-plan.md`, then `dev/organizer-rules.md` and
`dev/lane-rules.md`.

## Where the tree stands

Main is at the 0.55.1 release: `9033a26` plus its docs rebuild
`fad7a9f`. Twelve commits sit ahead of `origin/main`. **The user pushes;
no session pushes for them.**

Versions: frmtmb 0.55.1, frmtmb.eam 0.6.0, frmtmb.learn 0.3.0,
frmtmb.sample 0.4.0, frmtmb.spline 0.5.0, frmtmb.coupling 0.3.0,
frmtmb.ode 0.2.0, frmtmb.latent 0.2.2. Every extension floors on frmtmb
0.55.1; frmtmb.learn also floors on frmtmb.eam 0.6.0.

Verified at that commit, one test file per R process: core 134 files
8482 passing, eam 1408, sample 1118, coupling 428, spline 407, learn
333, latent 228, ode 223, no failures and no crashes. Gated tiers with
no skips: brms likelihood 404, agreement 187, methods 950, priors 103,
port 18, BCM 363 over thirteen chapters, reinforcement-learning identity
102, fuzz green, sample 327, learn 388. The scale tier green in all
three configurations. `R CMD check --as-cran` with four packages at OK
and four carrying only the environmental V8 NOTE.

## What is next, in order

`dev/extension-gaps-plan.md` is the plan and it is edited in place
rather than appended to. Phase 0 is done and Phase 1 is closed except
for the three items Phase 0 promoted into it:

- **1.0a**, the per-subject non-decision-time bound in frmtmb.eam. This
  is the most consequential item in the backlog, because it fixes a
  silent wrong answer that SHIPPED in 0.6.0 with a disclosure and no
  fix. Its acceptance criterion is measured: reach at least -7027.4,
  against -7148.8 today, with a positive definite Hessian and `sd(ndt)`
  recovered. `extensions/frmtmb.eam/tests/testthat/test-scale.R`
  currently asserts the DEFECT, with the recovery assertion commented
  beneath it, so it flips when this lands.
- **1.0b**, the same bound in `rlddm()`. It WAITS ON 1.0a, because
  frmtmb.learn takes its diffusion parameterization from frmtmb.eam and
  inherits both the bug and the fix.
- **1.0c**, `frm_lincmt()`, the analytic one- to three-compartment
  solution. Independent of the other two and can run beside 1.0a. It was
  Phase 5 until Phase 0 measured one `frm(se = TRUE)` on the plan's ode
  design at 3948 seconds. Its acceptance criterion is an identity with
  `frm_ode()` to 1e-8 across the schedule space, and that half is what
  makes the speed half worth anything.

After those, Phase 2 opens: recovery tables and third-party identities
at the scale Phase 0 measured. It was gated on Phase 0's table, which
now exists in `dev/scale-findings.md`.

## How a round runs here

Lanes in manual git worktrees off main, one item or one coherent group
per lane. A worker writes; a reviewer whose job is to FALSIFY rather
than confirm reads the same worktree against a reference build of the
base commit; punch rounds go back to the worker; the same reviewer
re-checks, because it holds the harness it already built. Two to four
rounds is normal. Nobody but the consolidating session commits.

    git worktree add ../frmtmb-wt-<name> -b wt-<name> <sha>

Consolidation: commit each lane on its branch, merge in an order that
puts renames first, resolve, set versions, roxygenise, install all eight
into one private library, run the suites and the gated tiers and the
checks and the docs, commit the release and the docs separately, remove
the worktrees, prune merged branches, update memory.

Two artifacts help and are worth keeping current:

- `dev/suite-baseline.tsv` and its note. A file-count audit cannot see a
  file that RAN and asserted less, which happens when a `test_that()`
  block throws. Regenerate it at each release from the release run.
- `tests/testthat/test-ci-siblings.R`. An extension depending on a
  sibling needs its workflow to install that sibling and to list it in
  `paths:`. frmtmb.learn's CI failed exactly that way.

## What the last two rounds cost, and what to carry

`dev/lane-rules.md` and `dev/organizer-rules.md` carry the operational
rules. The two that generalize furthest:

**Every guard built in the 0.55.1 round failed open on its first try.**
A loop that hit `next` on every name and asserted nothing. A hash check
that passed on a deleted file. An argument whose affirmative value
silently disabled the guard it gated. A count computed and never
asserted. When reviewing a guard, construct the case where the thing it
guards is ABSENT, not merely wrong.

**Instruments need checking as much as code.** Two timing claims
evaporated under replication, three test counts turned out to be capped
or stale, two contention artifacts were nearly filed as regressions, and
one reviewer's own patch had three holes. A number that will not
reproduce usually means the construction has not been found yet, not
that the record is wrong: twice a lane called a shipped figure false and
twice it was the lane that was wrong.

## Open decisions that belong to the user, not to a session

- Whether to push. Twelve commits are waiting.
- Any change that alters what a shipped parameter MEANS, or that refuses
  a model which currently fits. The standing policy is to break
  backward compatibility freely, since nothing external depends on this
  yet, but the user has wanted to hear about each one, and each needs a
  NEWS bullet saying plainly what stops working.
- Where a disclosure goes. The 0.6.0 `ndt` disclosure went in
  frmtmb.eam's NEWS only, at the user's direction, rather than anywhere
  more visible.

## Worktrees

Two were created for 1.0a and 1.0c and then removed unused when this
handoff was written. Create fresh ones off the current main.
