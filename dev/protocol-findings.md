# PROTOCOL-SEAMS lane, 2026-09-05

Worktree `C:/Users/adf44/source/r/frmtmb-wt-protocol`, branch `wt-protocol`,
from 68d6782. Private library
`.../scratchpad/ps-lib`, core plus `frmtmb.latent`, `frmtmb.eam` and
`frmtmb.sample` installed from this worktree.

## Reading, before any edit

`dev/structured-family-protocol.md` (404 lines, status COMPLETE through
step 10), `dev/rl-findings.md:236-330` (the rewritten seam section),
`dev/reviews/2026-09-05-reinforcement-learning.md:289-420` and `:586-612`,
`dev/reviews/2026-09-05-rdm-gng.md:510-585` (the `vint()` finding),
`dev/reviews/2026-09-05-core-seams.md:556-570` and `:780-806` (the
withdrawn row 5 and the residual cache item), `inst/rl/rw-delta.R`,
`extensions/frmtmb.latent/R/hmm.R`, `.../lca.R`,
`extensions/frmtmb.sample/R/loo.R`, `extensions/frmtmb.eam/R/ddm-shared.R`.

## What the code says, measured against the briefs

SEAM 1. `R/importance.R:660-661` is

    ll  <- wts * row_lpdf(fam, yraw, yraw, dpv, av, extra)
    agg <- smat %*% RTMB::matrix(ll, n, nd)

so per-row densities are collapsed into per-group sums on the next line.
`imp_group_map()` (`R/importance.R:333-343`) already refuses a row that
reaches two grouping levels. The refusal that blocks a structured family
is `R/importance.R:124-133`, and it fires on the PRESENCE of
`structure$loglik`, not on any property of it.

`imp_verify()` (`R/importance.R:733`) checks the stacked per-group
values against the plain objective, per group and in total. That is the
safety net that makes admitting a family-supplied per-group slot safe:
a family that mishandles the `n * nd` stacking is caught loudly at the
first freeze rather than silently.

SEAM 2. `core_aterms` (`R/parse.R:21-22`) is `weights trials cens trunc
se vint vreal mi`; everything else an `frm()` formula accepts comes from
`frmtmb_aterm_registry`. `subset`, `rate`, `index`, `thres` and `cat`
are brms aterms that frmtmb does not parse at all, so no declaration can
name them; `parse_response()` already refuses them by name.
`required_aterms` is a conjunction only (`R/frame.R:1186-1211`); nothing
refuses a term the density never reads.

SEAM 3. `frmtmb_compat_features_tbl()` (`R/compat.R:505`) rebuilds a
107-row data frame per call, ~16 ms. `frmtmb_register_aterm()` reads it
once at `R/compat.R:339` and then mutates
`frmtmb_compat_contrib$features` at `:357`. Plain invalidation would
still cost one build per registration, because each registration's own
read comes before its own write. Four registrations cost one build only
if the cache is EXTENDED with the rows the registration adds. That is
what this lane implements.

## Design decisions taken before writing code

1. Two new optional structure slots, `loglik_row` and `loglik_group`,
   both with `loglik`'s signature `(y, dpars, aterms, weights, block,
   extra)`. `unit` is untouched and keeps its meaning: the leave-one-out
   declaration.
2. A reserved block name, `group`: the family's own independent-unit
   label per row. It is what lets core derive `n`, the replicate count
   and the grouping-alignment check without any family convention.
   Required when `loglik_group` is declared.
3. Stacking is part of the contract, not an accident: the core may call
   either slot on the design stacked `nrep = length(y) / n` times, row
   `j` being original row `((j - 1) %% n) + 1`. `loglik_row` returns
   `n * nrep`; `loglik_group` returns `ng * nrep` in draw-major order,
   which is the `RTMB::matrix(ll, n, nd)` layout the correction already
   uses.
4. Deviance residuals need a saturated per-row comparison that no
   conditional log-density supplies on its own. Rather than a third
   slot, `loglik_row` may attach `attr(x, "saturated")`. Without it,
   deviance stays refused, by a sentence that names the omission.

## Landed so far

SEAM 3, done and measured. `R/compat.R`: `frmtmb_compat_features_tbl()`
is now a cache lookup over `compat_features_build()`, validated against
`frmtmb_compat_contrib$features` by `identical()` rather than by a dirty
flag, and EXTENDED (not dropped) by both registration paths. Plain
invalidation would not have helped: a registration reads the vocabulary
before it writes to it, so each one would still pay a build. Measured on
this machine, registry reset between cells:

| | before | after |
|---|---|---|
| bare build | 16.0 ms | 16.0 ms |
| cached read | - | 0.010 ms |
| four `frmtmb_register_aterm()` | 4 builds, 110 ms | 1 build, 70 ms |

Five tests in `tests/testthat/test-compat-register.R`, each comparing
the cached table against a fresh `compat_features_build(NULL)` rather
than against a stored expectation.

SEAM 2, done. `frmtmb_family(accepts_aterms =)`, `NULL` by default so a
custom family written before it keeps every term. Read together with
`required_aterms` (`accepted_aterm_names()`), spelled in formula names
without parentheses (`aterm_base()` maps `trunc_ub` to `trunc`, `vint2`
to `vint`, `reward1` to `reward`). The refusal runs at `R/frame.R` after
`valid_y`, LAST of the addition-term guards, so every existing specific
refusal still fires first and its message is unchanged; what it catches
is only what nothing else did.

Measured, on this worktree:

    gaussian: the addition term `vint()` is not one this family reads,
    so writing it would change nothing about the fit. This family takes
    `cens()`, `mi()`, `se()`, `trunc()`, `weights()`.

`frmtmb.eam`: the five families declare theirs and
`ddm_refuse_dec()` is deleted from `R/ddm-shared.R` with both call
sites. `rdm x vreal()` moves from "works, carried without effect" to
"refused", which is the compat row agreeing with the declaration rather
than describing the silence.

SEAM 1, core landed. Slots `loglik_row` and `loglik_group` on
`frmtmb_structure()`, reserved block name `group`, block validation at
frame assembly, the importance admission plus the grouping-alignment
check, and deviance residuals off `loglik_row` with its `saturated`
attribute. `inst/rl/rw-delta.R` declares both slots off ONE recursion
(`rw_terms()`), so the three quantities cannot drift.

Measured at the fitted values, 12 subjects x 40 trials:

    loglik total    -206.984938563
    sum(rows)       -206.984938563   (480 values)
    sum(groups)     -206.984938563   (12 values)
    sum dbinom(y, 1, fitted(fit))    -206.984938563
    deviance sum of squares 413.96988 = -2 * loglik exactly

and `frm(importance = 64)` on that model FITS, logLik -222.26 against
the Laplace -222.43, mcse 0.148. `imp_verify()` passed, which is the
check that matters: it compares the stacked per-group values against
the plain objective per group and in total, so the family's handling of
the `n * nd` stacking is verified rather than trusted.

## The aterm declaration table, measured against the compat rows

`accepted_aterm_names()` for every constructible family in
`family_registry` plus the four `frmtmb.eam` families, crossed with the
nine `aterm`-kind features in the compat vocabulary, against
`frm_compat()`'s verdict for the same pair:

| declaration | compat "works" | "refused" | "untested" | "conditional" |
|---|---|---|---|---|
| accepts | 62 | 0 | 1 | 2 |
| refuses | 0 | 143 | 143 | 0 |

The audit found ONE contradiction and it was a defect in my first
declaration, not in the table: `poisson` was given `cens()` because it
carries an `lcdf`, and `cens()` is refused for a DISCRETE response one
guard earlier. Fixed; `poisson` takes `trunc()` and `weights()`. The
two "conditional" cells (`cox` and `poisson` with `trunc()`) are not
contradictions: conditional means the pair works with a condition,
which is what an allow-list entry says too.

RESIDUAL, measured and not acted on: 143 pairs are declared refused and
read "untested" in the table. That is not a disagreement in the table's
own vocabulary - untested means nobody exercised the pair - but the
declaration now makes those answers computable, and the table could
derive its family-by-aterm rows from `accepts_aterms` instead of
leaving them to a kind-level default. That is a change to the shared
rule machinery in `R/compat.R` and belongs to whoever owns that next.

## The Laplace caveat, measured

40 subjects, the vignette's own model and truth, 6 seeds per design,
`importance = 100`. Laplace fit 0.4-7.7 s, corrected fit 14-121 s.

| | 40 x 100 trials | 40 x 20 trials |
|---|---|---|
| replicates the correction completed | 4 of 6 | 5 of 6 |
| mean shift, `alpha_(Intercept)` | -0.004 | +0.034 |
| mean shift, `alpha_conditiontrt` | +0.010 | +0.014 |
| mean shift, `beta_(Intercept)` | -0.002 | +0.048 |
| mean shift, `log sd(alpha)` | +0.357 | +0.480 |
| largest single shift in `log sd(alpha)` | 0.545 | 1.386 |
| mcse | 0.13-0.18 | 0.02-0.70 |
| smallest ESS per draw | 0.53 | 0.03 |

The fixed effects barely move, which is the vignette's existing
conclusion confirmed by a second route. The variance component moves a
lot, and the diagnostics say how far to trust that.

WHY THREE REPLICATES DID NOT COMPLETE, run individually. None is a
failure of the seam: two are the correction's own convergence guard
("the corrected negative log-likelihood ROSE ... check whether
`alpha: 1 | id + beta: 1 | id [ID]` has enough rows per group to
identify every variance component") and one is an NA gradient on a
Laplace fit that itself reported false convergence with a maximum
gradient of 1.8e9.

DRAW-COUNT SENSITIVITY, one 20-trial dataset (seed 4206). Laplace
`log sd(alpha)` -1.2024; at 100 draws the correction returns +0.1831,
a shift of 1.39 that takes sd(alpha) from 0.30 to 1.20 against a truth
of 0.5; at 400 draws it REFUSES, with the same unidentified-covariance
guard. So the corrected answer at 20 trials is not an answer: with
enough draws to see its own weights, the correction declines.

That is the honest finding, and it is not the one the first table
suggests: the vignette must not claim the variance-component bias is
Laplace error. What it can claim is that the question is now ASKABLE,
that at 100 trials the correction moves the variance component up
substantially where it converges, and that at 20 trials it refuses,
which is a reason to treat a 20-trial subject-level SD as a lower bound
rather than evidence about which of the two causes is at work.

## Is the instability the family's pieces or the model's covariance?

The vignette's model puts a correlated pair of subject effects on the
learning rate and the temperature (`(1 | p | id)` twice). Re-run with
one scalar random intercept instead, same family, `importance = 200`,
40 subjects x 100 trials:

| seed | Laplace `log sd(alpha)` | corrected | shift | mcse | min ESS/draw |
|---|---|---|---|---|---|
| 5301 | -1.816 | (refused) | | | |
| 5302 | -0.527 | -0.406 | +0.121 | 0.053 | 0.94 |
| 5303 | -0.833 | -0.651 | +0.182 | 0.037 | 0.89 |
| 5304 | -0.663 | -0.527 | +0.136 | 0.070 | 0.66 |

With a well-identified single variance component and 200 draws the
diagnostics are good and the shift is +0.12 to +0.18, against +0.36 for
the correlated pair at 100 draws. So the Laplace error in this family's
variance component is REAL and modest, and the larger shifts measured
on the vignette's model are partly the correction's own noise on a
weakly identified 2x2 covariance. The one refusal (seed 5301) is a
Laplace fit whose sd(alpha) had already collapsed to 0.16.

That is why the vignette says "about a third of a log unit, the same
order as the whole bias" and does not report a decomposition. The
scalar-model table is the cleaner measurement and it is here rather
than in the vignette, which is about the correlated model.

## On the vignette's own dataset, 30 subjects x 80 trials, seed 2024

| | sd(alpha) | mcse | min ESS/draw | seconds |
|---|---|---|---|---|
| Laplace | 0.4606 | | | 8 |
| importance = 100 | 0.6955 | 0.189 | 0.40 | 55 |
| importance = 200 | 0.5743 | 0.090 | 0.81 | 88 |

Truth is 0.5. The vignette's live chunk uses 200 draws, and keeps
`error = TRUE`, because the correction can decline on this model and a
vignette that fails to build over that would be worse than one that
shows the refusal.

## After the 2026-09-05 restart

The machine restarted at about 22:35 and killed every process. The
worktree (35 changed paths) and the private library survived; all four
packages load from it. Everything below was re-run one process at a
time.

APPLIED FOR wt-car-jacobian, at the coordinator's instruction.
`R/compat.R:1307` (the `car()` rule, not :1231 as numbered in their
findings) said `con_sd` "does scale that type's prediction standard
errors and ranef() conditional SDs as con_sd^2 in the variance". That
lane's exact Jacobian makes it false, and their replacement sentence is
in now, verbatim. NOTE for the merge: this worktree does not contain
their code, so on this branch the sentence is ahead of the behavior;
it becomes true when wt-car-jacobian lands. No test pins the rule text.

FOUND WHILE RE-CHECKING THE LIBRARY, and not mine to fix.
`frmtmb.sample` cannot be loaded in a session that has already loaded
`frmtmb.latent`:

    .onLoad failed in loadNamespace() for 'frmtmb.sample': ...
    frmtmb_register_compat(expects =) names 'hmm', and the registry
    already has that feature, so the declaration exempts nothing.

`frmtmb.latent` supplies `hmm` as a feature; `frmtmb.sample` declares
`expects = "hmm"`, and the registry refuses an expectation that is
already met. I checked it is NOT the vocabulary cache: with the cache
defeated so that every read is a fresh build, the failure is identical,
so this is the pre-existing `expects =` semantics meeting a load order
nothing had tried. The two fixes available are for `expects =` to
tolerate an already-present name (it is a forward reference, and one
that arrived early is not an error) or for `frmtmb.sample` to stop
expecting `hmm`. Neither is this lane's file. It also means an
extension suite must be run one package per process, which is what the
verification below does.

## Verification after the restart

Every suite re-run one process per file, `NOT_CRAN=true`.

| suite | files | pass | fail | error | skip |
|---|---|---|---|---|---|
| core | 117 | 6621 | 0 | 0 | 91 |
| frmtmb.latent | 4 | 227 | 0 | 0 | 0 |

The core row is the suite re-run AFTER the two late fixes below, on
the final source: 6621 passing, seven more than the pre-fix run,
which is the two new tests.
| frmtmb.eam | 17 | 1289 | 0 | 0 | 0 |
| frmtmb.sample | 10 | 888 | 0 | 0 | 2 |

The core file list was audited by name against `ls tests/testthat`:
117 files on disk, 117 result lines, no name missing. The brief says
116 at 68d6782; the tree has 117 and I added no file, so the count in
the brief is one low rather than a file having appeared.

Two runs needed repeating and neither was a test failure:
`test-rdm-gng.R` and `test-variability.R` produced no result line while
three suites ran concurrently (killed processes, not assertions), and
both pass when run alone (169 and 140).

TWO LATE FIXES to my own work, made after the first check run was
started and before it was re-run:

* `structure_unit_deviance()` clamped a negative unit deviance to zero.
  A saturated log-density below the fitted one is a family that has its
  two quantities the wrong way round, and clamping would report a
  residual of ZERO at exactly those rows, which reads as a perfect fit.
  It now refuses, naming the row and both values.
* The allow-list refusal named only the first offending term. It names
  all of them now, so a formula with two unread terms is one round trip
  rather than two.

Both have tests in `test-structure.R`, which is at 113 passing.

## test-rl-example.R, gated and not

NOT gated: 96 passing, 2 skipped (the Stan tier), 0 failures.

GATED (`FRMTMB_BRMS_FIT_TESTS=true`, `NOT_CRAN=true`, the Stan cache
copied from the main checkout to `dev/stan-cache`): 97 passing and 2
errors, both in the two Stan-identity tests, and both are Stan MODEL
COMPILATION failures rather than assertion failures.

That is environmental and I checked it rather than assuming it. The
copied cache holds 54 programs and NOT this one: the RL program's key
under rstan 2.32.7 is `3087656b89a076a3e34b80556035865f`, which is
absent, so the tier needs a fresh compile. A fresh compile of a
two-line Stan program

    parameters { real x; } model { x ~ normal(0, 1); }

fails on this machine too, after 99 seconds, with the same
`make ... Error 1` out of StanHeaders and RcppEigen. So the toolchain
cannot build any Stan model here, the tier cannot run, and nothing in
this lane touches `rl_stan_code()`, `helper-brms.R` or `helper-rl.R`.
The identity itself is unaffected by anything I changed: the family's
total log-likelihood is now the sum of its per-subject pieces, and
`test-rl-example.R` asserts that sum against the independent scalar
reference (`the taped recursion equals an independent scalar
reference`), which passes.

## Did folding `loglik` into the pieces cost the hot path anything?

`rw_loglik()` is now `sum(rw_loglik_group(...))` rather than a scalar
accumulated one trial at a time, so the vignette's shape-and-cost table
is a claim about code I changed. Both shapes measured in ONE process,
interleaved, five rounds, 40 subjects x 100 trials, gradients over
batches of 20 (the machine was running the as-cran check at the time,
so read the columns against each other and not against the vignette's
absolute numbers):

| `loglik` shape | tape build | one gradient | objective |
|---|---|---|---|
| new: sum of the pieces | 0.140 s | 122.0 ms | 1982.6691425636 |
| old: scalar accumulator | 0.160 s | 128.0 ms | 1982.6691425636 |

Identical to 1e-12, and the new shape is not slower on either column.
So the vignette's table stands, and the reason it stands is worth
stating: the pieces were always there. The recursion produced one
vector per trial either way; the old shape summed each vector as it
went and the new one sums them at the end, which is the same number of
additions on the tape.

## R CMD check --as-cran, core

`_R_CHECK_CRAN_INCOMING_=false`, pandoc from the RStudio quarto tools
directory on PATH, checked against the private library. The first run
stopped at `checking package dependencies ... ERROR: Package suggested
but not available: 'frmtmb.spline'`, which is a sibling extension the
library did not hold; installing it from `extensions/` and re-running
gave:

    Status: 1 ERROR, 1 WARNING, 2 NOTEs

Every one of the four, named:

* ERROR, `checking PDF version of manual without index`, and
* WARNING, `checking PDF version of manual`. One cause, and it is not
  an Rd problem despite the message saying it typically is:
  `frmtmb.Rcheck/Rdlatex.log` says `pdflatex is not available`, and
  there is no pdflatex on this machine. Every Rd stage the check can
  run passed: Rd files, Rd metadata, Rd line widths, Rd
  cross-references, Rd \usage sections, Rd contents, all OK.
* NOTE, `checking HTML version of manual`: "Skipping checking math
  rendering: package 'V8' unavailable".
* NOTE, `checking for non-standard things in the check directory`:
  `frmtmb-manual.tex`, the file the failed PDF build left behind.

What matters passed: R code for possible problems (52 s) OK, examples
(42 s) OK, examples with --run-donttest (74 s) OK, tests (294 s) OK,
re-building of vignette outputs (358 s) OK. The vignette rebuild is the
one that exercises the new live `importance` chunk.

## The cache, measured again on a shared machine

Six sibling R processes were running, so absolute times are inflated
and the robust numbers are the ratio and the count:

| | measured |
|---|---|
| one build | 16.0 ms quiet, 18.3 ms and 23.5 ms under load |
| one cached read | 0.0010 ms (200000 reads in 0.20 s) |
| a build costs | about 23500 cached reads |
| four `frmtmb_register_aterm()` calls | 1 build, was 4 |

The build count is the claim that cannot move with load, and
`test-compat-register.R` pins it.

## What was left out, and why

* A per-row slot for `hmm()`. A row's emission density is not its
  contribution to the likelihood, because the state that emitted it was
  reached through every earlier row; there is no per-row factor to
  declare. `deviance` stays refused there in the words it already used.
* A factorization slot for `lca()`. Its likelihood IS rowwise, one row
  per subject, so a slot would be a second definition of numbers the
  core already has. The constructor refuses one, and lca.R says why.
* An exclusivity rule for `required_aterms` alternative groups. That is
  what would close `wiener()` accepting `dec()` and `vint()` together
  with the second unread, which the allow-list cannot: both spellings
  are legitimately accepted, and the fault is supplying both. `gddm()`
  is why the rule is not a one-liner - it reads `dec()` and `vint()`
  together, with the condition index moving between `vint1` and
  `vint2` depending on whether `dec()` is present.
  DONE 2026-09-07 as `frmtmb_family(exclusive_aterms =)`, opt-in for
  that reason; see `dev/aterms-findings.md`.
* Deriving the compat table's 143 family-by-aterm "untested" cells from
  the new declarations. Now computable, but it is a change to the
  shared rule machinery rather than to a declaration.
  DONE 2026-09-07. 148 cells, not 143, because `dec()` joined the
  vocabulary; audited against every hand-written row first, and no row
  disagreed with any declaration.
* `structure_unit()` is still dead code in core (`R/structure.R`), the
  residual `dev/rl-findings.md` records. Nothing this lane added reads
  it; `frmtmb.sample` still inlines the expression deliberately.
  DONE 2026-09-07: `structure_generic()` reads it, which is the caller
  the slot was documented for.
* Fixing the `frmtmb.sample` load-order defect found above. Not this
  lane's file, and the fix is a semantics decision about `expects =`.

# Punch round, 2026-09-06

Review: `dev/reviews/2026-09-06-protocol.md`, verdict PUNCH, five items.
All five addressed; none disputed. The reviewer also showed my gated
blocker was my own toolchain, and it was: `R_MAKEVARS_USER` pointed at
`dev/stan-cache/makevars-cxx17.mk` compiles the two-line program in
79 s here, and the gated RL file then passes 102 / 0 / 0 / 0, the
reviewer's count exactly.

**Item 1, the weights convention.** `R/predict.R:2417` no longer
multiplies the returned unit deviance by the row weights; the family
applies them inside the slot, which is what `R/importance.R:756`
already did with the same slot. The convention is now written down at
`R/structure.R:168-178`, beside the saturated-value paragraph, and the
saturated paragraph says the attribute carries the same weight so the
deviance comes out weighted once.

Measured, and the pin verified by reinstating the old behavior in the
namespace: with weights 1, 2, 4, 3 the ratio of the toy family's
deviance residuals to `gaussian()`'s is constant at 1.04680 (spread
7.97e-14) after the fix and spreads from 1.04680 to 2.09360 before it,
the two differing by exactly 1.0000, 1.4142, 2.0000, 1.7321 = sqrt(w).
Test: "the family applies the row weights, and the core does not",
`tests/testthat/test-structure.R`.

While writing it I hit a related edge the review did not name: the
structured deviance path lives inside the `fitted_mean` branch of
`residuals()`, so a structure with `loglik_row` and `deviance = TRUE`
but NO `fitted_mean` falls through to the rowwise path and is refused
with "family 'x' has no standard unit deviance", which is a misleading
message for a family that supplied both halves. The toy declares a
`fitted_mean`, which is what a real family wanting deviance residuals
would do; the misleading fall-through is recorded, not fixed.

**Item 2, the loo unit.** Both halves. `loo_matrix()`
(`extensions/frmtmb.sample/R/loo.R:405-419`), the one funnel `loo()`,
`waic()` and `psis()` pass through, now emits a message naming the unit
and the column count when the matrix carries one; and the claim at
`loo.R:220` is corrected to say the attribute never reaches
`loo::loo.matrix()` and that the printed elpd carries no mark of being
leave-one-unit-out. Test: "a group-unit matrix says so, because loo()
cannot", `test-loo.R`, 74 passing.

**Item 3, the motivating example.**
`dev/structured-family-protocol.md:516-541` now carries the reviewer's
measurement: `wiener()` CAN read `vint()` as its boundary, so
`rt | vint(upper)` and `rt | dec(upper)` agree bit for bit at
-130.566406836; the real defect is the two supplied together, and the
allow-list cannot close it because both are legitimately on the list.
The doc now agrees with this lane's own "left out" entry instead of
contradicting it.

**Item 4, the 400-draw refusal.** Named and bounded in
`vignettes/reinforcement-learning.Rmd`. Re-measured before writing:
seed 4206 at 40 x 20 gives Laplace `log sd(alpha)` -1.2024, +0.1831 at
100 draws (mcse 0.336, min ESS/draw 0.329) and a refusal at 400 draws,
reproducing my earlier numbers exactly. The vignette now names that
dataset, and says in the next paragraph that other 20-trial datasets
are fine, quoting the reviewer's two seeds. "The instability is plain"
is gone; the claim is now that the correction is unreliable at this
design rather than uniformly broken.

**Item 5, the unused level.** `droplevels()` in
`structure_group_codes()` (`R/structure.R:757-763`) so the core's own
codes are contiguous whatever arrives, AND a named refusal in
`check_structure_block()` (`R/structure.R:816-831`). droplevels alone
would have been silent, and silence is wrong here: a family that
returns one value per DECLARED level would then be one value long and
aligned to the wrong groups. The refusal names the family, the response
and the unused levels. Test: "an unused level in the block's grouping
is refused by name".

Findings 5 and 8 need no action per the review and got none. Finding 3
(`frmtmb.sample` after `frmtmb.latent`) is pre-existing and not this
lane's; the reviewer confirmed it at 68d6782 with a base-built library
and named the one-line fix.

## Punch-round verification

Fresh process per file, `NOT_CRAN=true`, counts by name:

| file | pass | fail | error | skip |
|---|---|---|---|---|
| core test-structure.R | 122 | 0 | 0 | 0 |
| core test-importance.R | 201 | 0 | 0 | 0 |
| core test-rl-example.R (ungated) | 96 | 0 | 0 | 2 |
| core test-rl-example.R (GATED) | 102 | 0 | 0 | 0 |
| core test-method-residue.R | 64 | 0 | 0 | 0 |
| core test-message-uniqueness.R | 6 | 0 | 0 | 0 |
| sample test-loo.R | 74 | 0 | 0 | 1 |
| sample test-message-uniqueness.R | 6 | 0 | 0 | 0 |
| latent, 4 files | 227 | 0 | 0 | 0 |

test-structure.R rose 113 to 122 and test-loo.R 70 to 74, which is the
three new tests.

GATED TIER, run by me this time. `R_MAKEVARS_USER` pointing at a copy
of `dev/stan-cache/makevars-cxx17.mk` under my own prefix,
`FRMTMB_STAN_CACHE` at a private copy of the cache,
`FRMTMB_BRMS_FIT_TESTS=true`: the two-line Stan program compiles in
79 s and `test-rl-example.R` passes 102 / 0 / 0 / 0. My earlier report
that the toolchain could not build any Stan model was true of my
environment as I had it and false of the machine; the reviewer was
right and the fix was one environment variable.

ROXYGEN across all four packages, run twice: the second pass changes
nothing. `git status` over `man/` and `NAMESPACE` is identical between
passes and the two regenerated files hash the same.

AS-CRAN, with `RSTUDIO_PANDOC` and TinyTeX's `bin/windows` prepended to
PATH:

    Status: 1 NOTE

The one NOTE is `checking HTML version of manual ... Skipping checking
math rendering: package 'V8' unavailable`, an absent optional package.
Everything else is OK, including the three that were not last time:
`checking PDF version of manual ... OK` (TinyTeX supplies the pdflatex
whose absence produced the previous ERROR and WARNING), and `checking
for non-standard things in the check directory ... OK` (no leftover
`frmtmb-manual.tex`, because the PDF build now succeeds). Examples 26 s
OK, examples with --run-donttest 26 s OK, tests 128 s OK, re-building
of vignette outputs 150 s OK.
