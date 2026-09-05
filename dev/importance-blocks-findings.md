# wt-imp-blocks findings: several blocks over one grouping factor

Extends `frm(importance = N)` (0.50.0, R/importance.R) from exactly one
random-effect block to any number of blocks over ONE grouping factor.
Blocks over different factors stay refused.

The scratchpad is wiped without warning, so every number below is
transcribed here as it was measured. Re-measure before quoting.

## ENVIRONMENT
Private library `<scratchpad>/ib-lib`, a mirror of the user library
(356 packages: RTMB, testthat, pkgload, lme4, GLMMadaptive, numDeriv,
MASS all present and loadable). The shared user library
`C:/Users/adf44/AppData/Local/R/win-library/4.6` was HEALTHY at the
start of this lane and is used only as a read fallback in `R_LIBS`.
Every command sets `R_LIBS` = ib-lib; user library, and `R_LIBS_USER` =
ib-lib.

## THE FRAME FACT THE WHOLE LANE RESTS ON (measured, not assumed)
`(1 | g)` in mu and `(1 | g)` in sigma are TWO blocks, never one. Only
an `|ID|` key merges terms into a single block (R/frame.R:2040). The
probe:

    n blocks: 2
      label=1 | g         dpar=mu    cs=us dim=1 nlev=12 group=g
        b_idx: 1-12   c_idx: 1-12   theta_idx: 1
      label=sigma: 1 | g  dpar=sigma cs=us dim=1 nlev=12 group=g
        b_idx: 13-24  c_idx: 13-24  theta_idx: 2
    lp y.mu    Z: 96x24, nonzero cols 1-12
    lp y.sigma Z: 96x24, nonzero cols 13-24

So each linear predictor's `Z` spans the FULL coefficient vector and
fills only its own block's columns, which is why
`build_importance_objective()`'s `Z %*% u` needed no change at all.
`b_idx` and `c_idx` are identical for every structure in
`imp_covstructs`; they diverge only for `rr`, which the whitelist
refuses, and `imp_plan()`'s count check would catch it anyway because
the two spaces have different lengths.

`(1 + x || g)` also produces two blocks (diag, dim 1 each) over one
factor, so the double-bar spelling is covered by this lane too.

## THE SCATTERED INDEX (imp_layout)
`b` is level-major WITHIN a block, but blocks are CONCATENATED, so with
more than one block a group owns one run per block and the runs sit
`n_levels * dim` apart. Group k's positions are

    idx[, k] = for each block m: c_idx_m[(k - 1) * q_m + seq_len(q_m)]

built as one reshape per block, `matrix(c_idx, q_m, ng)`, whose column
k IS that run. `imp_layout()` returns `idx` (q_total x n_group) and
`pos_group` (one group label per position, the inverse lookup).

Measured layouts:
- `(1|g) + (0+x|g)`, 20 groups: `idx[, 1] = (1, 21)`, `idx[, 20] =
  (20, 40)`.
- `(1|g)` mu + `(1|g)` sigma, 12 groups: `idx[, 1] = (1, 13)`.
- `(1+x|g)` mu + `(1|g)` sigma, 12 groups: `idx[, 1] = (1, 2, 25)`,
  `idx[, 2] = (3, 4, 26)`. The mu block's two coefficients are
  adjacent; the sigma block's sits 2 * ng away.

### THE SIX SITES IT REPLACED
The previous lane's note named five; the reviewer named four
(R/importance.R:308, 344, 364, 477) and did not name the verification
swap. The union is six, and all six now read the same object, so they
cannot drift apart the way six spellings could:

1. `imp_plan()` block-diagonality gate: `lev <- rep(seq_len(ng), each =
   q)` became `lev <- lay[["pos_group"]]`. With two blocks the Hessian
   has genuine cross-block entries WITHIN a group, and the gate must
   not read them as cross-group.
2. `imp_plan()` scalar branch: it assumed the layout implicitly
   (`u <- uhat + sdv * z` over the whole vector). Now addressed through
   `p1 <- idx[1L, ]`. It is still the q = 1 fast path (one scalar block
   and nothing else), and it is bit-identical there because `idx[1, ]`
   is then `seq_len(ng)`.
3. `imp_plan()` general branch: `ii <- (k - 1L) * q + seq_len(q)`
   became `ii <- idx[, k]`. This drives BOTH the Hessian slice
   `hess[ii, ii]` and the draw placement `u[ii, ]`, which is what makes
   the proposal the group's joint conditional over every block.
4. `imp_plan()` zsq: `z[(seq_len(ng) - 1L) * q + j, ]` became
   `z[idx[j, ], ]`. The reviewer flagged this one specifically; missed,
   the half-norms would be read off the wrong rows with NO error
   raised, because the shapes still match.
5. `build_importance_objective()` u_slices: same stride, same fix,
   `plan[["u"]][plan[["idx"]][j, ], ]`.
6. `imp_verify()` group swap: `idx <- (g - 1L) * q + seq_len(q)` became
   `idx <- plan[["idx"]][, g]`. The per-group check swaps the WHOLE
   group and nothing else, or the difference stops isolating one
   group's term.

A sweep for any other contiguity assumption (`(k-1)*q`, `each = q`,
`%/% q`, `%% q`) across all of R/ turns up nothing else. The two
remaining `each =` uses are correct: `pos_group`'s construction, and
`col_level` inside one block, where level-major contiguity DOES hold.

## THE OTHER CHANGES
- `imp_prior_terms()` is now a SUM over blocks. Each block's per-level
  Gaussian is extracted from its own registry density at its own
  `theta_idx`, on its own slice of the group's coefficients. The blocks
  are independent given theta, so the group's precision is block
  diagonal and no cross terms appear. The extraction body moved to
  `imp_block_prior()`; it needs no `ADoverload()` of its own, because
  the arithmetic there is `*` and `+` on advectors and
  `covstruct_registry[[k]]$nll` installs its own overloads
  (R/covstruct.R:25 and elsewhere). Stated in a comment rather than
  assumed, since the previous lane lost overloads exactly this way.
- `imp_group_map()` gathers (row, level) pairs across EVERY block and
  every linear predictor before checking, so the same pass that catches
  a multi-membership term also catches two blocks that disagree about
  which level a row belongs to. Its duplicate key is now formed in
  DOUBLE arithmetic: `n * (ng + 1)` overflows an integer at a few
  million rows, and an NA key would silently drop the clash it looks
  for. (A pre-existing latent bug in a line this lane rewrote anyway.)
- `imp_ess_warning()` names the GROUPING FACTOR rather than a term
  label: one proposal covers a level's coefficients from every block,
  so there is one effective sample size per level, not one per block.
  Its signature is unchanged because R/fit.R:942 calls it and R/fit.R
  is outside this lane.
- `imp_frozen_proposal()` returns `layout` where it returned `block`.
  R/confint.R passes that object around opaquely, so nothing outside
  R/importance.R reads the renamed field (checked).

## THE REFUSALS
Zero blocks, previously folded into the "exactly one" message, now has
its own:

> `importance` needs a random-effect block to correct, and this model
> has none. The correction reweights an integral over random effects,
> and a model without any has no such integral and no approximation to
> improve on. Use importance = 0

Different grouping factors, naming the blocks AND the factors:

> `importance` takes several random-effect blocks only when they share
> ONE grouping factor, and this model spreads 2 blocks over 2 factors
> (`1 | g` over g, `1 | h` over h). Crossed or nested factors make the
> marginal likelihood one integral over every factor at once, which
> does not split into the per-group integrals the correction resamples.
> Fit with importance = 0

One factor, different level sets:

> `importance` needs every block over `g` to carry the same grouping
> levels, and `1 | g` has 12 where `sigma: 1 | g` has 11 (first
> difference: '3'). The proposal draws one level's coefficients from
> every block at once, so a level carried by only some of them has no
> joint Gaussian to be drawn from. Fit with importance = 0

Nesting is refused by the SAME message as crossing, because a nested
factor is still a second factor: `(1 | g) + (1 | g:h)` spreads two
blocks over two factors.

The level-set refusal has NO user-facing route: two blocks over the
same column take their levels from that column and always agree. It is
tested on a doctored frame through `check_importance_scope()` directly,
because a level set that silently disagreed would put a group's draws
on another group's rows rather than raise anything.

## THE BRUTE-FORCE REFERENCE (validation b)
Written in tests/testthat/test-importance.R, sharing NO code with the
package: Gauss-Hermite nodes by Golub-Welsch (eigenvalues of the Jacobi
matrix of the Hermite recurrence; the squared first components of the
eigenvectors are the weights, already normalized to sum to 1), then a
tensor 2-D quadrature over `(b_mu, b_sigma)` per group against
`dnorm()` directly.

Two mistakes worth recording, both caught by measurement:

1. WEIGHT NORMALIZATION. Dividing the Golub-Welsch weights by
   `sqrt(pi)` double-counts the normalizer. On a 50-group design the
   reference came out 28.62 too high, which is exactly
   `50 * log(sqrt(pi))`. The correct normal-expectation weight is the
   squared eigenvector component itself.
2. NODE COUNT. Non-adaptive Gauss-Hermite converges slowly when the
   conditional is tight relative to the prior, and at nq = 60 the
   reference was wrong by an amount that LOOKED like estimator bias
   (10 to 20 times the seed-to-seed spread). Convergence on the design
   finally used:

   | nq | reference | change |
   |---|---|---|
   | 60 | 193.16131561 | - |
   | 100 | 193.16184589 | +5.30e-04 |
   | 150 | 193.16182912 | -1.68e-05 |
   | 200 | 193.16182817 | -9.46e-07 |
   | 300 | 193.16182816 | -1.32e-08 |

   The test uses nq = 150, converged to 1.1e-06, which is 37000 times
   below the Monte Carlo standard error it is compared at. The
   reference is vectorized (one `dnorm` call per b_mu node against an
   m x nq sigma matrix), so nq = 150 costs 2.7 s rather than minutes.

PIN ON THE REFERENCE ITSELF: on a single-block gaussian model, where
the Laplace approximation is exact, the 1-D form of the same reference
reproduces the Laplace objective to 1.8e-12 (470.1021040398 both). That
is what fixes the parameterization: `theta` is the log standard
deviation of a dim-1 `us` block, `betad` the log residual SD.

### THE AGREEMENT
Design: 12 groups of 10 rows, `y ~ x + (1 | g)` with
`sigma ~ 1 + (1 | g)`, gaussian, data seed 19, random-effect SDs 0.8
(mu) and 0.4 (sigma). All values at the Laplace optimum, so the
proposal is anchored where the fit would report from.

| quantity | value | error vs reference | in MCSE |
|---|---|---|---|
| GHQ reference (nq = 150) | 193.161829 | - | - |
| importance, 2000 draws, seed 1 | 193.137101 | -0.024728 | 0.61 |
| Laplace | 193.396458 | +0.234629 | 5.81 |

Reported MCSE 0.04039, min ESS 0.395, median ESS 0.882.

The Laplace error is real here and not a rounding difference, which is
the point of choosing a sigma block: sigma depends on its random effect
through a log link, so the joint density is NOT quadratic in the random
effects and there is something to correct. (On the gaussian mu-only
designs the Laplace approximation is exact and the correction has
nothing to do.)

Over ten proposal seeds at 2000 draws the worst miss was 1.30 MCSE, so
the test's `3 * mcse` tolerance has a 2.3x margin over what was
actually observed:

| seed | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| miss / MCSE | 0.61 | 0.09 | 0.30 | 0.73 | 0.97 | 0.40 | 0.07 | 0.81 | 0.65 | 1.30 |
| min ESS | 0.395 | 0.772 | 0.655 | 0.780 | 0.431 | 0.540 | 0.817 | 0.578 | 0.611 | 0.313 |

The MCSE itself ranged over 0.0296 to 0.0469, and min ESS over 0.313
to 0.817, all above the 0.25 floor, so this design never trips the
coverage warning.

DESIGNS REJECTED, and why (all measured over 10 seeds against the
converged reference):
- 15 groups of 8, SDs (0.9, 0.7): min ESS falls to 0.004 on some
  seeds. The estimator stays honest (worst miss 1.57 MCSE) but the
  proposal degenerates, which makes a poor validation design.
- 15 groups of 10, SDs (0.9, 0.4): worst miss 2.24 MCSE, min ESS
  0.011. Too close to a 3 MCSE tolerance to be a stable test.

## VALIDATION (c): A DIM-2 BLOCK PLUS A SCALAR ONE
`y ~ x + (1 + x | g)` with `sigma ~ 1 + (1 | g)`, so a group's
coefficients are scattered AND its Hessian slice is a genuine 3 x 3
with cross-block entries.

`imp_verify()` passes. Gradient of the corrected objective against
numDeriv at 500 draws:

| shift | mean relative difference |
|---|---|
| 0.0 | 1.58e-06 |
| 0.3 | 6.87e-06 |

The TOLERANCE IS SET BY numDeriv, not by the AD gradient: two numDeriv
settings (its default, and `r = 6, d = 1e-4, eps = 1e-6`) disagree with
EACH OTHER by 3.0e-04 at the optimum and 1.0e-02 at the shifted point,
which is larger than the gap being measured. The test therefore uses
1e-4 rather than the single-block file's 1e-6, and says why.

Coverage at the optimum, 500 draws: min ESS 0.172, median 0.881,
MCSE 0.117.

## TIMINGS versus the single-block case at the same draw count
Same data (12 groups of 10 rows, 120 rows) so only the block structure
differs. `one` is `(1|g)`, `two` is `(1|g)` mu + `(1|g)` sigma, `three`
is `(1+x|g)` mu + `(1|g)` sigma. Each figure is the mean of as many
repetitions as fit in 2 seconds, after a warm-up call. (An earlier
attempt at these numbers was meaningless: single calls sat below the
15 ms Windows clock resolution and the first design measured absorbed
the JIT warm-up, which made the two-block plan look 10x FASTER than the
one-block one.)

| N | case | q | plan+build | tape | fn | gr |
|---|---|---|---|---|---|---|
| 500 | one | 1 | 3.3 ms | 0.48 s | 0.04 ms | 4.61 ms |
| 500 | two | 2 | 13.2 ms | 0.50 s | 0.04 ms | 4.82 ms |
| 500 | three | 3 | 17.5 ms | 0.67 s | 0.04 ms | 7.12 ms |
| 2000 | one | 1 | 19.6 ms | 2.05 s | 0.05 ms | 21.05 ms |
| 2000 | two | 2 | 37.0 ms | 7.39 s | 0.15 ms | 42.08 ms |
| 2000 | three | 3 | 137.3 ms | 5.60 s | 0.08 ms | 48.78 ms |

**CORRECTED, punch round 2026-09-05. No ratio should be read off the
gr column above either.** The claim it carried, "about 1.05x at 500
draws and 2.00x at 2000", does not reproduce, and the table says so on
its own terms: its one-block gr row scales 4.6x from N = 500 to 2000
while its two-block row scales 8.7x, and the tape is linear in N for
both arms, so the 42.08 ms is an outlier and the 2.00x is that outlier
divided by a good number. The figures came from unpaired means, one
process per arm, which is not a measurement of a 10 percent effect on
a loaded machine.

Re-measured twice, every arm timed inside ONE batch loop so machine
drift lands on all of them alike, both arms warmed before either is
timed, and every gradient call at a different parameter vector so that
nothing can be served from a cache. Median of 15 batches.

| comparison | N | one | two | ratio |
|---|---|---|---|---|
| `(1\|g)` mu -> `(1\|g)` mu + `(1\|g)` sigma | 500 | 7.03 / 3.85 ms | 8.96 / 4.83 ms | **1.27x / 1.25x** |
| `(1\|g)` mu -> `(1\|g)` mu + `(1\|g)` sigma | 2000 | 29.00 / 15.88 ms | 34.78 / 19.12 ms | **1.20x / 1.20x** |
| `(0+x+z\|g)` -> `(0+x\|g) + (0+z\|g)` | 500 | 5.17 / 5.63 ms | 5.06 / 5.09 ms | **0.98x / 0.90x** |
| `(0+x+z\|g)` -> `(0+x\|g) + (0+z\|g)` | 2000 | 24.02 / 25.76 ms | 19.91 / 20.99 ms | **0.83x / 0.82x** |

The two figures in each cell are the two runs. The absolute values move
between runs by up to 1.8x, which is the machine, and the ratios do
not, which is why they have to be measured paired. Every arm now scales
3.9x to 4.7x from N = 500 to N = 2000, against the 4x the draw count
implies: the internal inconsistency is gone.

So the second block costs **roughly 1.2x per gradient here, and the
reviewer's independent paired measurement of the same comparison gave
1.14x and 1.09x. Call it 1.1 to 1.3x**, not 2x.

**And that ratio is not a cost of this lane.** The first comparison
changes the MODEL: it adds a sigma linear predictor. The last two rows
are the controlled one, which the lane should have run and did not:
same data, same predictors, same `q_total`, same number of random
effects, varying ONLY the block count, one dim-2 block against two
dim-1 blocks. It is 0.98x / 0.90x and 0.83x / 0.82x, so two blocks are
never slower than the one block that holds the same coefficients, and
at N = 2000 they are faster, because a dim-2 `us` block carries a
correlation parameter that two dim-1 blocks do not. **The
several-blocks machinery costs nothing per gradient.**

The plan build figure and its stated mechanism survive, measured as
plan + build together. `(1|g)` -> `(1|g)` + sigma `(1|g)` is 2.5x /
2.2x at N = 500 and 1.9x / 2.0x at N = 2000, the low end of the "2 to
7x" claimed above. The controlled comparison, where BOTH arms have
`qt = 2` and so both take the per-group Cholesky loop, is 1.08x / 1.06x
and 0.96x / 0.83x. That is the mechanism confirmed from both sides: the
2x is the cost of leaving the `q == 1` vectorized branch, not the cost
of a second block. It is milliseconds against a tape measured in
seconds, so it still does not matter.

The TAPE column is allocator-noisy (three at 2000 draws taping faster
than two is not reproducible) and no ratio should be read off it.

End to end through `frm()`, 12 groups of 10, two blocks:
Laplace 2.1 s, importance at 500 draws 13.9 s (5 rounds, not capped,
grad 4.0e-04, min ESS 0.716, median 0.860, MCSE 0.068). Laplace nll
193.39646, corrected 193.11008.

## WHAT DID NOT CHANGE
The estimator, the stacking, the reporting, the diagnostics, the round
loop, the divergence guard and the frozen-proposal profile warning are
untouched. The single-block path is bit-identical: the whole existing
test file (123 assertions) passes unchanged apart from the call sites
that now pass a layout instead of a block, and the two refusal messages
that were reworded.

## DELIBERATE OMISSIONS
- BLOCKS OVER DIFFERENT FACTORS stay refused, by name. This is not a
  gap to be filled later by a draw count: the marginal likelihood of a
  crossed or nested model is one integral over every factor at once and
  has no per-group factorization to resample.
- `R/fit.R`'s `@param importance` roxygen still says "exactly one
  random-effect block over a grouping factor". It is now wrong, and
  R/fit.R is outside this lane's declared surface (two sibling lanes
  own regions of it), so it was NOT edited. The integrator should
  replace that clause with: "any number of random-effect blocks over
  ONE grouping factor, of any dimension and any covariance structure
  that is Gaussian within a level and independent between levels". The
  vignette, NEWS and `frm_compat("importance")` all state the new scope
  correctly.
- No test of three or more blocks. The index is built by a loop over
  blocks with no special case at two, and the dim-2-plus-scalar design
  already exercises unequal block dimensions, so a third block would
  add cost without adding a code path.
- The per-group proposal is the JOINT conditional over all blocks,
  including its cross-block entries. A cheaper block-diagonal proposal
  (one Cholesky per block per group) was not tried: it would be a
  different, worse proposal, and the cross terms are exactly the
  coupling the group's rows create.

## VERIFICATION

### The five named files, one file per process, load_all, NOT_CRAN=true
| file | failed | errors | warnings | skipped | passed |
|---|---|---|---|---|---|
| test-importance.R | 0 | 0 | 0 | 0 | 173 |
| test-compat.R | 0 | 0 | 0 | 0 | 260 |
| test-message-uniqueness.R | 0 | 0 | 0 | 0 | 6 |
| test-bracket-access.R | 0 | 0 | 0 | 0 | 7 |
| test-autoscale.R | 0 | 0 | 0 | 0 | 46 |

test-importance.R was 123 assertions before this lane and is 173 after.
The 123 that existed pass unchanged; the only edits to them were the
call sites that now pass a layout instead of a block, and two refusal
messages that were reworded (the zero-block case, and two blocks over
different factors).

### R CMD check --as-cran --no-manual
`_R_CHECK_CRAN_INCOMING_=false`, RSTUDIO_PANDOC on PATH, private
library. **Status: OK. 0 NOTEs, 0 WARNINGs, 0 ERRORs.**
  examples [79s] OK
  examples with --run-donttest [85s] OK
  tests [292s] OK
  re-building of vignette outputs [309s] OK

A FIRST run of this check reported 2 WARNINGs, and they were MINE, not
the package's: building the tarball with `--no-build-vignettes` leaves
no `inst/doc`, so check reports "Files in the 'vignettes' directory but
no files in 'inst/doc'" and lists all SEVEN vignettes, including six
this lane never touched. Rebuilt with vignettes and the check is clean.
Worth recording so the next lane does not chase it: that warning pair
is a build-flag artifact, not a finding.

Note the standing caveat from the previous lane: R CMD check runs the
suite WITHOUT NOT_CRAN, and test-importance.R opens with
`skip_on_cran()`, so the check never exercises the correction itself.
The per-file suite under NOT_CRAN=true is what tests it.

### FULL CORE SUITE, one file per process
**104 files: failed = 0, errors = 0, warnings = 0, skipped = 21,
passed = 5473.** No file failed to run.

The 21 skips are all opt-in tiers, not capability gaps, and none is
this lane's: test-brms-likelihood.R (18) and test-brms-agreement.R (2)
sit behind `FRMTMB_BRMS_FIT_TESTS=true` (tests/testthat/helper-brms.R),
and test-fuzz.R (1) behind `FRMTMB_FUZZ=true`. The previous lane's
figure was 102 files and 4 skips at a233c3c; HEAD has since added the
brms log-density tier (564e185), which is where the extra files and
the extra skips come from.

### SCOPE
`git diff --name-only 564e185`:
  NEWS.md, R/compat.R, R/importance.R,
  tests/testthat/test-importance.R, vignettes/diagnostics.Rmd
Untracked: dev/importance-blocks-findings.md
Nothing outside the lane's declared surface. R/frame.R, R/fit.R,
R/objective.R and R/priors.R are untouched. No man/ regeneration was
needed: every function this lane changed is `@noRd`, and no .Rd file
embeds the `frm_compat()` note text (checked).

The lane made NO commits; the worktree is still at 564e185. The main
checkout is clean and was never written to. (It has moved on its own
to df9d6bc since this lane started, by another lane's integrator, not
by this one.)

## Punch round 2026-09-05 (F1-F9 from dev/review-importance-blocks.md)

Worked in `C:/Users/adf44/source/r/frmtmb-wt-imp-blocks` on branch
`wt-imp-blocks`, still uncommitted and still based at 564e185. This
round made no commits. The main checkout was read twice and never
written to: it was at 5dfdd84 and clean when this round started and at
84d328b and clean when it ended, moved by another lane's integrator.
Every
number below was measured in this round in a private library holding
this worktree's own install; where the reviewer measured the same
quantity, both figures are given.

### F1. `R/fit.R:158-162` scope sentence - FIXED
The `@param importance` roxygen said "exactly one random-effect block
over a grouping factor", which the lane made false. Replaced with the
draft from the review's section 9: "any number of random-effect blocks,
of any dimension, provided they share ONE grouping factor and one set
of levels, in any covariance structure that is Gaussian within a level
and independent between levels". `R/fit.R:149` took the plural in the
same pass ("the blocks may have any dimension"). The paragraph the
review recommended is added after the scope list at `R/fit.R:170-174`,
naming distributional regression as the case users will actually meet
and repeating that crossed or nested factors are refused.

`man/frm.Rd` regenerated. The edit is confined to the two importance
paragraphs; the `mo()` region of the same file's roxygen was not
touched.

### F2. `R/fit.R:1313-1325` `@param importance_ess` - FIXED
It promised "the hardest design in the test suite holds `0.43` at its
own optimum". That figure was retired from `R/importance.R` by the
previous lane and this lane retired it a second way. Rewritten to
describe the two regimes the threshold separates, with every figure
carrying its design and its draw count, all re-measured here:

| design | draws | worst group |
|---|---|---|
| probe, 60 groups of 8 Bernoulli rows, correlated slope | 1000 | 0.947 |
| the same probe, both log SDs displaced by 0.5 | 1000 | 0.033 |
| gaussian response, where the correction has nothing to correct | any | 1.000 |
| `(1 + x \| g)` mu + `(1 \| g)` sigma, 12 groups of 10 | 500 | 0.172 |
| the same design | 1000 | 0.268 |
| the same design | 2000 | 0.400 |

The 0.172 is the review's number, reproduced to the digit. No single
retired figure is promised now: the text says what moves the number
(draws) and which regime each design sits in. `man/frmtmb_control.Rd`
regenerated.

### F3. `dev/importance-blocks-findings.md` TIMINGS - CORRECTED IN PLACE
The claim "about 1.05x at 500 draws and 2.00x at 2000" does not
reproduce, and the old table contradicts itself: its one-block gr row
scales 4.6x from N = 500 to 2000 where its two-block row scales 8.7x,
on a tape that is linear in N for both arms.

Re-measured twice here, all four arms timed inside ONE batch loop so
drift lands on every arm alike, both arms warmed before either is
timed, median of 15 batches, and **every gradient call at a different
parameter vector**. That last point matters: a first run of this probe
called `gr()` at the same parameter every time and produced a 4x draw
count that cost only 2x, which is not possible for this tape. With the
parameter moving, every arm scales 3.9x to 4.7x from N = 500 to 2000,
against the 4x the draw count implies, and the internal inconsistency
is gone.

| comparison | N | ratio, run 1 / run 2 | reviewer |
|---|---|---|---|
| `(1\|g)` mu -> `(1\|g)` mu + `(1\|g)` sigma | 500 | 1.27x / 1.25x | 1.14x |
| the same | 2000 | 1.20x / 1.20x | 1.09x |
| `(0+x+z\|g)` -> `(0+x\|g) + (0+z\|g)` | 500 | 0.98x / 0.90x | 0.97x |
| the same | 2000 | 0.83x / 0.82x | 0.98x |

So roughly 1.1x to 1.3x across both of us for the second block, not 2x.
The controlled rows are the ones that answer the lane's own question:
same data, same predictors, same `q_total`, varying only the block
count. Two dim-1 blocks are never slower than the one dim-2 block that
holds the same coefficients, and at N = 2000 they are faster, since the
dim-2 `us` block carries a correlation parameter they do not. The
several-blocks machinery costs nothing per gradient.

Plan build, measured as plan + build together: 2.5x / 2.2x at N = 500
and 1.9x / 2.0x at N = 2000 for one block -> two, and 1.08x / 1.06x and
0.96x / 0.83x for the controlled pair where both arms already have
`qt = 2`. The lane's stated MECHANISM is therefore confirmed from both
sides: the cost is leaving the `q == 1` vectorized branch for the
per-group Cholesky loop, not the second block. Absolute times moved by
up to 1.8x between the two runs while the ratios did not, which is the
whole argument for pairing.

### F4. `R/importance.R:1022` and `R/fit.R:942` one-block spellings - FIXED
`imp_record()` and the `imp_ess_warning()` call site now read the
layout, which already carries `levels` and `group_name` for exactly
this purpose. `importance_fit()` hands `lay` out with its result
(`R/importance.R:957-962`), `imp_record(imp)` lost its `frame` argument
(`R/importance.R:1014-1024`), and `imp_ess_warning()`'s second argument
is named `lay` rather than `bk` (`R/importance.R:1036-1042`), so the
name says where the fields come from rather than which block happened
to be first. `frame[["re_blocks"]][[1L]]` no longer appears anywhere in
`R/importance.R` or in `R/fit.R`; the only package-wide survivor is the
doc example at `R/sampling-api.R:167`, which is unrelated.

### F5. `R/importance.R:1112` unread `layout` field - DROPPED
`imp_frozen_proposal()` returned `layout = lay` and nothing read it.
Dropped rather than used: its two consumers, `imp_ess_at()` and the
profile-ESS warning, need `io`, `plan` and `template` only, and naming
a level in that warning would be a behavior change the punch list did
not ask for.

### F6. `R/importance.R:464` mismatched delimiters - FIXED
The Cholesky-failure message opened the term label with a backtick and
closed it with an apostrophe. It now closes with a backtick.

### F7. `R/compat.R:460-461` the `double_bar` row - DRAW COUNT STATED
Measured on the row's own design, 20 groups of 6 gaussian rows,
`(1 + x || g)`, which really does build two `diag` blocks over one
factor with group 1 at positions 1 and 21:

| draws | rounds | capped | min ESS | MCSE | corrected - Laplace |
|---|---|---|---|---|---|
| 200 | 5 | yes, warns | 1.000000 | 2.6e-09 | 0.018441 |
| 500 | 2 | no | 1.000000 | 0 | 0.002083 |
| 1000 | 2 | no | 1.000000 | 0 | 0.000060 |
| 2000 | 2 | no | 1.000000 | 0 | 0.000114 |
| 4000 | 1 | no | 1.000000 | 0 | 0.000192 |

The Laplace value is 195.118267. The note now says the claim was
verified at 1000 draws, gives the 6e-05 agreement and the two rounds,
and adds that the weights are equal at every draw count while the ROUND
LOOP is not: at 200 draws the design exhausts its five rounds and
warns. That reproduces the reviewer's finding exactly.

### F8. the `imp_ess_warning()` rename is now pinned - ADDED
`tests/testthat/test-importance.R:240-252`, inside "the effective
sample size diagnostic separates the two regimes", which runs on a
SINGLE-block design: that is the case the rename also changed, so the
pin belongs there. The message is captured and matched on
``groups of `g` poorly`` (the grouping factor), asserted NOT to contain
the term label `x | g`, and matched on the level-and-value spelling
`'<level>' (`. The two call sites in that test now pass `p$lay`, which
is what the fit passes since F4, and the helper's `bk` field is gone
with them, so the test file no longer spells `re_blocks[[1L]]` either.
That block went from 5 assertions to 8.

### F9. the non-gaussian multi-block value check - ADDED
`tests/testthat/test-importance.R:664-767`, section (b2), with its own
data helper and its own reference integrator, placed after the gaussian
brute-force section so it can reuse that section's `imp_gh()` nodes and
`imp_lse()`.

30 groups of 8 Bernoulli rows, `y ~ x + (1 | g) + (0 + x | g)`, two
blocks over one factor with group 1 owning positions 1 and 31,
integrated by tensor Gauss-Hermite quadrature over
`(u_intercept, u_slope)` against `dbinom()` directly. Seed 5, so it is
deterministic; the file-level `skip_on_cran()` already gates it.

| quantity | value | error vs reference | in MCSE |
|---|---|---|---|
| reference, nq = 70 (nq = 50 agrees to 8e-09, nq = 100 to 1e-11) | 140.180284248 | - | - |
| Laplace | 141.040485564 | **+0.860201** | 56 at 8000 draws |
| importance, 2000 draws | 140.175337990 | -0.004946 | 0.17 |
| importance, 8000 draws | 140.184317279 | +0.004033 | 0.26 |

Min ESS 0.866 at 2000 draws and 0.880 at 8000, so the joint proposal
covers this design and the agreement is not a degenerate cancellation.
The correction removes 99.5 percent of a Laplace error that no
algebraic identity fixes. The test asserts the correction within
`3 * mcse` of the reference and the Laplace value beyond `10 * mcse` of
it, at BOTH draw counts, plus the scattered index and the coverage.

The test also pins the reference before trusting it: `nq = 50` must
give the same number as `nq = 70`.

Two designs were rejected before this one, and are worth recording. At
seed 37 the fitted standard deviations are 1.62 and 2.35 and min ESS
falls to 0.054 at 8000 draws, so the proposal does not cover it; at
seed 101 the Laplace error is only 0.249, which is a weaker
demonstration. Seed 5 has a large Laplace error AND a well-covered
proposal, which is the combination the check needs.

Reviewer's own version of this design landed within 0.03 MCSE at 8000
draws against a Laplace error of 0.45. Mine lands within 0.26 MCSE
against a Laplace error of 0.86. Both say the same thing; the MCSE
ratio differs because the designs are not the same draw.

### Two things the punch round did NOT do
`NEWS.md` gained one sentence, not a new bullet: the existing bullet
claims validation "against a brute-force reference that shares no code
with the package", and after F9 there are two such references, so the
sentence names the Bernoulli one and its numbers. No NEWS entry was
written for F1 or F2: F1 documents a feature that has not shipped yet,
and F2 is a correction to prose about a released argument whose
behavior did not change.

`R/importance.R`'s `imp_ess_floor` comment was left alone. It carries
the previous lane's own measurements, F2 was about `R/fit.R` only, and
two of its figures do not reproduce for me at the settings it states
(probe 0.947 where it says 0.93, displaced 0.033 where it says 0.009).
Recording that here rather than editing it: the numbers are close
enough to be a different draw count or a different displacement, the
comment is out of this round's scope, and rewriting it on a guess would
put a third unverified figure in the file. **For a later round.**

### VERIFICATION

Every run below is against this worktree, one test file per process,
`NOT_CRAN=true`, in the private library `scratchpad/ibp-lib` holding
this worktree's own install, with counts audited by name.

**A harness defect of mine, recorded so the next lane does not chase
it.** The first run of this suite used
`pkgload::load_all(export_all = FALSE)` and reported 7 failures and 9
errors in `test-compat.R` and 9 errors in `test-autocor.R`. Those files
reference package internals by BARE NAME, and `test_file()` on a path
runs them under an environment whose parent is the search path rather
than the namespace, so the internals were not findable. Nothing was
wrong with the package. With the default `export_all = TRUE` both files
are clean. That run was discarded, and only the numbers below stand.

**The four named files, one process each:**

| file | blocks | passed | failed | errors | warnings | skipped | secs |
|---|---|---|---|---|---|---|---|
| test-importance.R | 31 | 185 | 0 | 0 | 0 | 0 | 311.6 |
| test-compat.R | 29 | 260 | 0 | 0 | 0 | 0 | 15.2 |
| test-message-uniqueness.R | 1 | 6 | 0 | 0 | 0 | 0 | 7.2 |
| test-bracket-access.R | 3 | 7 | 0 | 0 | 0 | 0 | 1.2 |

`test-importance.R` was 30 blocks and 173 assertions when this round
started and is 31 and 185 now: F9 adds one block with 9 assertions and
F8 adds 3 to "the effective sample size diagnostic separates the two
regimes", which went from 5 to 8. 173 + 9 + 3 = 185. Every other block
holds the count it had. The F9 block itself runs in **4.8 s**.

`test-compat.R` is 260 assertions, the same number the lane recorded,
so the reworded `double_bar` note broke none of the registry
invariants; `test-message-uniqueness.R` is 6 and still passes with the
reworded Cholesky-failure message from F6.

**The vignette:** `vignettes/diagnostics.Rmd` renders clean through
`rmarkdown::render()` with pandoc from the RStudio quarto tools
directory. (The first attempt failed with "pandoc version 2.8 or higher
is required" because that directory was not on PATH in the driver
script, not because of anything in the vignette.)

**The full core suite, one process per file:**

    files=104  blocks=942  passed=5485  failed=0  errors=0
    warnings=0  skipped=21

No file failed to run. The lane recorded 104 files, 941 blocks and 5473
passed; this round adds exactly the one block and the twelve assertions
F8 and F9 introduced (941 + 1 = 942, 5473 + 12 = 5485). The 21 skips
are the same opt-in tiers the lane recorded and none is in this lane's
scope: `test-brms-likelihood.R` (18) and `test-brms-agreement.R` (2)
behind `FRMTMB_BRMS_FIT_TESTS`, `test-fuzz.R` (1) behind `FRMTMB_FUZZ`.

**`R CMD check --as-cran`** with `_R_CHECK_CRAN_INCOMING_=false`, the
tarball built WITH vignettes (which is what avoids the reviewer's two
`inst/doc` warnings), pandoc 3.8.3 from the RStudio quarto tools
directory, private library:

    Status: OK        0 ERRORs, 0 WARNINGs, 0 NOTEs

    checking whether package 'frmtmb' can be installed ... [12s] OK
    checking R code for possible problems ............... [32s] OK
    checking for code/documentation mismatches .............. OK
    checking Rd \usage sections ............................. OK
    checking Rd line widths ................................. OK
    checking examples .................................. [34s] OK
    checking examples with --run-donttest .............. [24s] OK
    checking tests .................................... [103s] OK
    checking re-building of vignette outputs .......... [100s] OK

No V8 NOTE arose. The doc-mismatch and Rd stages are the ones F1 and F2
could have broken, and both pass.

The standing caveat from the previous lanes still holds: `R CMD check`
runs the suite WITHOUT `NOT_CRAN`, and `test-importance.R` opens with
`skip_on_cran()`, so the check never exercises the correction. The
per-file suite above is what tests it, and that is also why F9's 4.8 s
does not land on CRAN.

### SCOPE
`git diff --name-only 564e185`:

    NEWS.md  R/compat.R  R/fit.R  R/importance.R
    man/frm.Rd  man/frmtmb_control.Rd
    tests/testthat/test-importance.R  vignettes/diagnostics.Rmd

Untracked: `dev/importance-blocks-findings.md` (this file),
`dev/review-importance-blocks.md` (the reviewer's).

The three files beyond the lane's own five are exactly what F1 and F2
require: `R/fit.R`, which the lane deliberately did not own, and the
two `.Rd` files roxygen regenerates from it. `roxygen2::roxygenise()`
was run twice and the second run wrote nothing, so the generated
documentation is idempotent. The `R/fit.R` edits sit in the
`@param importance` paragraphs, the `@param importance_ess` paragraph
and the two `imp_record()` / `imp_ess_warning()` call sites; the `mo()`
region another lane is editing was not touched.

### WHAT THE PUNCH ROUND CHANGED ABOUT THE VERDICT
Nothing. The review's GO WITH FIXES stands, and none of F1 to F9 was a
correctness defect in the shipped path. Two of them turned out to be
worth more than "nit": F3, because an unpaired timing measurement had
put a wrong number in the permanent record and the same mistake was
easy to repeat (this round made it once, with a cached gradient, and
caught it only by checking that a 4x draw count costs 4x); and F9,
because it is now the only multi-block check whose answer is fixed by
neither an identity nor a gaussian.
