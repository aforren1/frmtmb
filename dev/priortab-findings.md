# Lane priortab: prior tables and prior messages against brms

Round of 2026-11-11, base e8abed0a (frmtmb 0.81.0), branch
`wt-r82-priortab`. Brief: `round-briefs/r082-priortab.md`. Nothing is
committed.

## How every number here is made

| what | script | output (gitignored) |
| --- | --- | --- |
| brms 2.23.0 on every case, no compile | `dev/priortab-brms.R` | `dev/priortab-out/brms.txt` |
| frmtmb's side of the same cases | `dev/priortab-frm.R` (`PRIORTAB_LIB` unset = rellib-r20) | `frm-before.txt`, `frm-after2.txt` |
| messages, fits and round trips | `dev/priortab-frm2.R` | `frm2-before.txt`, `frm2-after.txt` |
| install into the lane library | `dev/priortab-install.sh` | `inst-*.log` |
| suites, one file per process, 3 at once | `dev/priortab-run-par.sh` | `priortab-<tier>.log` |

All data: seed 1, `n = 120`, `g` with 10 levels (the scripts' header).

## Per item: brms, frmtmb before, frmtmb after

STATUS: see the section "Status" at the end; this table is the record.

| item (backlog line) | brms 2.23.0 | frmtmb 0.81.0 | lane |
| --- | --- | --- | --- |
| 1. `lb` on `sd`/`sigma` (3232, 4250, 4530, 5081) | `sd` class row `lb = 0`, group rows blank; `sigma` 0; `nu` 1; `shape`, `phi`, `kappa`, `beta` 0; `zi`, `hu`, `zoi`, `coi`, `quantile` 0 and 1; `ndt` 0 and `min_Y`; `xi`, `alpha` blank | every `lb`/`ub` NA | as brms on both routes; `ndt`'s `ub` only on `route = "sample"` (frm()'s log link has no top) |
| 2. mixture `theta` row (1115, 3229) | `dirichlet(1) theta` row | none | NOT DONE: needs a decision (below) |
| 3. table `+` set_prior (4293, 4648) | brmsprior of both, either order | "Incompatible methods", then "non-numeric argument" | a `frmtmb_priorlist` of both, either order, route kept |
| 4. `rescor` with `coef` (4643) | refused (`Lrescor_y__y2`); `group` refused too | accepted, applied to the whole matrix | refused, `coef` and `group` |
| 5a. `cor` with `dpar` (1109) | refused (`L_mub`, `L_sigma`); `nlpar` refused (`L_a`) | accepted | refused, `dpar` and `nlpar` |
| 5b. `cor` rows of a dpar block (3201) | one `cor g` row | `cor g` and `cor g sigma` | one `cor g` row |
| 5c. mv `cor` rows carry `resp` (4638) | one `cor g` row, no resp; `resp = "y2"` refused | `cor g y`, `cor g y2`; resp accepted | one `cor g` row; `resp` still accepted (user decision of 2026-09-24) |
| 6. `Y_a` printed (3161) | `Ya`, `Yb` | `Y_a`, `Y_b` | `Ya`, `Yb`; both spellings still accepted |
| 7. one `prior_sd_g` per id (4366) | two statements `prior_sd_1`, `prior_sd_2`, both renamed `prior_sd_g` | one statement per group | NOT DONE (below) |
| 8. `car()` rho flat on its logit (4360) | `car` row `lb = 0`, `ub = 1`, (flat) | flat on logit(rho) | NOT DONE: M, a seam (below) |
| 9. 4283 one-column table | "undefined columns selected" | "Duplicated ... 'b' is given 5 times" | refused for the dropped columns |
| 9. 4312 `list(sp1)` | "Argument 'prior' must be a 'brmsprior' object" | "a list of 1 set_prior() specifications ... sp1 + sp2" | "Pass the specification itself"; `list()` and `list("flat")` get their own |
| 9. 4316 `frm(prior = "flat")` | (no such string) | generic message | points at frm_sample() and at leaving `prior =` out |
| 9. 4266 coerced sampling table refused | n/a | refused | per row since the fix round: `""` rows read, printed `"(flat)"` rows on a default slot refused and named (the first version accepted the whole frame; review MAJOR 3) |
| 9. 4246 `b (flat)` line | n/a | printed whatever the prior | "(flat) where your prior sets none" when class b has one |
| 9. 4307 `ar` order 2/3 under `"flat"` | no such model | `uniform(-2, 2)` | NO CHANGE (below) |
| 9. 4925 `uniform(0, 0.2)` on a slope | warns, samples unbounded | warns, then "NA/NaN gradient evaluation" | warns, fits at the bound 0.2, as with `lb = 0, ub = 0.2`; since the fix round only a density positive at its end is boxed (review MAJOR 2) |

## What changed, and where

The review of the first version (`dev/reviews/2026-11-11-priortab.md`,
FIX) changed the support box, the drop of bounds-only rows and the
coerced-table check; where this section and "Fix round" below
disagree, "Fix round" is the record.

Core `R/priors.R` (owned by this lane):

- `prior_table()`: the class-wide `sd` row of each prefix carries
  `lb = 0`; a distributional parameter's own class row carries
  `dpar_row_bounds()`, the ends of its link's image at the ends brms
  declares (`prior_class_bounds()`, brms's `dpar_bounds()`), plus on
  `route = "sample"` a bound the registered default carries. `cor`
  rows are keyed by group alone. `resp` is written as
  `brms_stan_name()`. New helpers `dpar_row_bounds()`, `link_image()`.
- `+.frmtmb_priorlist` reads a table as `c()` reads it; the same
  function object is registered as `+.frmtmb_prior_rows`, so R's Ops
  group dispatch finds two identical methods and no longer warns
  "Incompatible methods".
- The resolver refuses `cor` with `dpar`/`nlpar` and `rescor` with
  `coef`/`group`, and boxes a class `b` coefficient to its density's
  support where no bound is given (`support_box()`,
  `prior_dist_support()`).
- `drop_unbounding_slots()`: a bounds-only specification that writes
  no finite bound on the model is dropped before the fit. Without it,
  every table passed back (now carrying `lb = 0` on `sd`/`sigma`)
  would have made `fit$prior` non-NULL, which `vcov_cluster()` refuses
  and gp-start and the fit label read. It also removes a pre-existing
  case: a table of any model with a smooth (`sds` row, `lb = 0`).
- `check_table_slot_columns()`; `prior_input_refusal()` messages for
  `list(sp)`, `list()`, `"flat"` and a density string;
  `prior_spec_slot()` keys `resp` by its Stan name, so `Y_a` and `Ya`
  are one slot.

Outside `R/priors.R`, each a small hunk:

- `R/fit.R` (shared with quadseam, boundfam, nlstop, gpstop, ordstart):
  three lines after `prior_drop_unmatched()` in `frm()`, calling
  `drop_unbounding_slots()`.
- `R/sampling-api.R`: `registered_prior_defaults()` attaches the
  default's `lb`/`ub` as attributes of its density string (stripped in
  the table by `as.vector()`).
- frmtmb.sample `R/sample.R`: `announce_default_priors(user_pl =)` and
  its call; `table_lists_defaults()` before the coerced-table refusal.
- `NAMESPACE` (`S3method("+", frmtmb_prior_rows)`), `man/set_prior.Rd`,
  `man/default_prior.Rd` (roxygen).
- Tests: `tests/testthat/test-prior-table-brms.R` (new),
  `extensions/frmtmb.sample/tests/testthat/test-prior-table-sample.R`
  (new).
- `NEWS.md`, `extensions/frmtmb.sample/NEWS.md`, `dev/test-backlog.md`
  (CLOSED lines only).

## Does any fit move? (`dev/priortab-frm2.R`, seed 1)

The brief asks, since `lb = 0` on `sd` is a change of what a fit does
if the row is used.

| case | before | after |
| --- | --- | --- |
| `frm(y ~ z + (1 \| g), prior = default_prior(...))` objective | 35.0222221 | 35.0222221, `identical(opt$par)` TRUE |
| the same, student | 59.9320008 | 59.9320008, identical par |
| the same, zero_inflated_poisson | 182.0708116 | 182.0708116, identical par |
| sample route `(1 \| g)`: -log prior of the table vs the defaults, at the mode | 8.303773143 vs 8.303773143 | 8.303773143 vs 8.303773143 |
| sample route student: same | 55560143.88 both | 55560143.88 both |

The restated bounds map to infinite internal bounds (`internal_bound()`),
so no objective and no box changes. On the sampling route a table now
carries a default's own bound (`nu` above 1 under a log link, `ndt`
below `min_Y`), which a table passed back used to lose: that is a fix
of the round trip, not a change of the default.

`uniform(0, 0.2)` on `z` with no bounds (`y2 = y + 0.5 z`): before,
"NA/NaN gradient evaluation"; after, `z = 0.2`, logLik
-126.7129011913, code 0, the same fit as with `lb = 0, ub = 0.2`
(both arms of that control: 0.9906694047, 0.2). `uniform(0, 2)`
(interior) gives 0.7535962762 on both arms: inside its support the
uniform is a constant, and the box changes nothing there
(`test-prior-table-brms.R` checks it against the plain ML fit to a
relative 1e-5). `frm_sample()` of the bounded case
(`dev/priortab-probe2.R`, seed 3, 200 draws): before, "the sampler
returned no draws"; after, every draw of `b_z` in [0.1544289,
0.1999861].

## Decisions this lane did not take, with a recommendation

### Q1. The mixture weights' row (items 1115, 3229)

brms lists `dirichlet(1)` on class `theta` for a mixture's weights.
frmtmb's class `"theta"` is the internal covariance vector, and
`set_prior("dirichlet(...)", class = "theta")` is refused with a message
that says the word names two sets of parameters. A row that cannot be
passed back would break the table's round trip, so none was added.

Recommendation: rename frmtmb's internal class (no existing users; for
example `thetacov`), then list `theta` as brms does, and route
`dirichlet(alpha)` on class `theta` to the weights on both routes
(`frm()` a penalty on the softmax chart, `frm_sample()` in place of
`mix_theta_nlp()`'s `dirichlet(1)`). That is a lane of its own: every
reader of class `"theta"` in core, the extensions and the tests.

### Q2. `cor` with `resp` (item 4638, second half)

brms refuses it; frmtmb keeps it by the user's decision of 2026-09-24.
This lane refused `cor` with `dpar` and `nlpar` (item 1109, brms
refuses both), so `resp` is now the one narrowing a correlation prior
takes. Recommendation: refuse it too, for one rule (brms's: `cor` is
keyed by group alone). Not done, because the decision is the user's.

## Not done, and why

- Item 7 (4366), one `prior_sd_g` per id: brms makes two Stan
  statements (`prior_sd_1`, `prior_sd_2`, `dev/priortab-out/brms.txt`
  case 7b) and its renaming gives both the name `prior_sd_g`
  (`dev/priordraws-findings.md`, section 8). Two draws columns of one
  name cannot be told apart by name, so brms's output is not one to
  copy (the project tiebreaker); both draw from the same prior. The code
  is frmtmb.sample `R/prior-draws.R`, lane sampfront's.
- Item 8 (4360), `car()`'s rho: M, not S. The fix needs no new
  density: `logistic(0, 1)` on logit(rho) is exactly brms's (flat) on
  (0, 1). What it needs is a covariance-structure seam slot that names
  a logit-charted position (core `frmtmb_register_covstruct()`), its
  use in frmtmb.spatial's car registration (lane arspat's package) and
  in frmtmb.sample's defaults and `prior = "flat"`. Three packages and a
  seam; left filed.
- Item 4307, `ar` of order 2 or 3 under `prior = "flat"`: no change.
  The `uniform(-2, 2)` printed is the density sampled, and the only
  spelling that round-trips: a `(flat)` row with `lb`/`ub` on class
  `ar` of order above 1 is refused by the resolver ("A bound needs order
  one"). brms has no such model.
- Item 4266 is closed for the coerced sampling table (three of the four
  reported cases); a data frame written by hand with frmtmb's printed
  `(flat)` on a slot brms gives a default is still refused, by design:
  brms writes `""`, never `"(flat)"`.

## Tests

New files, each seen failing on rellib-r20 first:

- `tests/testthat/test-prior-table-brms.R`: rellib-r20 pass=20 fail=24
  err=4 warn=1 (`dev/priortab-out/t-ptb-before.txt`; the `lb` columns
  NA, `+` "non-numeric argument", the four spellings accepted, `Y_a`
  listed, the old messages, "NA/NaN gradient evaluation", `Y_a` and
  `Ya` taken as two slots); lane pass=71 fail=0.
- `extensions/frmtmb.sample/tests/testthat/test-prior-table-sample.R`:
  rellib-r20 pass=2 fail=6 err=1 (`t-pts-before.txt`: `lb`/`ub` NA, the
  coerced table refused, the old slopes line, no support box); lane
  pass=12 fail=0.

Suites on the lane build, one file per process, 3 at a time
(`dev/priortab-run-par.sh`, `NOT_CRAN=true`, every log's `lib:` line
names `wt-priortab-lib` for frmtmb):

| suite | files | pass | fail | err | skip | warn |
| --- | --- | --- | --- | --- | --- | --- |
| frmtmb (core, `priortab-core1.log`) | 318 | 23558 | 0 | 0 | 239 | 0 |
| frmtmb.sample (`priortab-sample1.log`) | 87 | 4050 | 0 | 0 | 18 | 0 |
| frmtmb.coupling (`priortab-ext1.log`) | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1848 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 687 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.spatial | 8 | 478 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 15 | 556 | 0 | 0 | 0 | 0 |
| frmtmb.sample, new file on the final build | 1 | 12 | 0 | 0 | 0 | 0 |

The extensions other than frmtmb.sample are rellib-r20's builds, on the
lane's core (their `lib:` lines). The frmtmb.sample suite above ran
before two last edits (`prior_spec_slot()`'s Stan-name key in core,
the same key in `table_lists_defaults()`); the sessions below run the
final build.

Sessions on the final build (`dev/priortab-run-sessions.sh`, a copy of
`dev/release/run-sessions.sh` writing under `dev/priortab-out/`, log
`priortab-sess1.log`); every LEAK line is tagged `[with load of ...]`:

| session | files | pass | fail | err | skip | warn |
| --- | --- | --- | --- | --- | --- | --- |
| frmtmb.sample forward | 83/88 | 4062 | 0 | 0 | 18 | 0 |
| frmtmb.sample reverse | 83/88 | 4062 | 0 | 0 | 18 | 0 |
| frmtmb.sample preload | 83/88 | 4062 | 0 | 0 | 18 | 0 |
| frmtmb forward | 305/318 | 23558 | 0 | 0 | 239 | 0 |
| frmtmb reverse | 305/318 | 23723 | 0 | 0 | 239 | 0 |
| frmtmb preload | 305/318 | 23826 | 0 | 0 | 241 | 0 |

Gated tier (`FRMTMB_BRMS_FIT_TESTS=true`) on every file that names a
prior table, `cor`, `rescor` or a uniform (core) or a prior (sample),
log `priortab-gated1.log`, no skip:

| suite | files | pass | fail | err | skip | warn |
| --- | --- | --- | --- | --- | --- | --- |
| frmtmb gated | 57 | 7015 | 0 | 0 | 0 | 0 |
| frmtmb.sample gated | 38 | 2092 | 0 | 0 | 0 | 0 |

The ported brms suite: `dev/priortab-brmsport-run.R` (a copy of
`dev/brmsport-run.R` that takes `lib;rellib-r20` and loads every
extension from both) records `test-brms-suite-priors.R` (pass=36,
67 recorded rows) and `test-brms-suite-brmsformula.R` (pass=16, 19
rows) on both arms, all eight extensions loaded. The record files are
identical before and after (`rec-before-*.tsv` against
`rec-after-*.tsv`): no ledger row flips, so no manual verdict file is
edited. The only bin-2 priors row (`priors:130`, the intercept default
under an offset) reads a sampling default this lane does not change.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
`OPENBLAS_NUM_THREADS=4`; Rblas 74 symbols and Rlapack 616 forwarded;
`dev/priortab-out/ob-*.txt`; `dev/ciharden-out/` deleted after):
test-prior-table-brms.R 71/0/0, test-prior-table-sample.R 12/0/0,
test-prior-mv-resp.R 36/0/0, test-prior-location-dpar.R 55/0/0,
test-lkj.R 213/0/0, test-prior-one-slot.R 90/0/0 (pass/fail/err; no
skip, no warning). 0.3.32 not run: no new assertion sits near a
rounding-decided branch (below).

The new assertions that a rounding could move, over 20 responses
perturbed by `y * (1 + 2e-15 k)`, k = 1 to 20 (`dev/priortab-ulps.R`,
`ulps.txt`, seed 1): the support box and the explicit box give
identical `opt$par` 20 of 20; the slope sits at exactly 0.2 20 of 20;
the interior uniform's slope against plain ML differs by at most
1.236e-13 relative (the test allows 1e-5); the round-trip objective
differs by 0 relative 20 of 20, and `fit$prior` is NULL 20 of 20.

## Status

Done: items 1, 3, 4, 5 (`cor` with `dpar`/`nlpar` refused, `cor` rows
by group), 6, and from item 9: 4283, 4312, 4316, 4266 (coerced
sampling table), 4246, 4925. Found and fixed on the way: `Y_a` and
`Ya` were two slots (both accepted, one dropped in silence; on both
arms, `dev/priortab-probe2.R`), and a table passed back made an ML fit
a penalized one by label.

Not done: item 2 and the `resp` half of 4638 (decisions Q1 and Q2),
item 7 (brms's duplicate names; sampfront's file), item 8 (M, a seam
across three packages), 4307 (no change, with the reason).

Bump: frmtmb minor (BREAKING table rows, refusals, and the drop of
no-op bounds). frmtmb.sample minor (it now accepts tables it refused).
Floor: frmtmb.sample's code calls only `spec_spelling()` and
`brms_stan_name()`, both in 0.81.0, but its new test reads the lane's
`lb` columns and support box, so its `frmtmb (>= )` floor is the core
version this round releases.

`R CMD check --as-cran` (`dev/priortab-check.sh`, R_LIBS = lane
library, rellib-r20, user library; pandoc and TinyTeX on PATH):
frmtmb `Status: 1 NOTE` (the HTML manual's V8 math-rendering NOTE,
expected).
frmtmb.sample `Status: 1 NOTE` (the same V8 NOTE). Tests inside the
checks: frmtmb 46 min OK, frmtmb.sample 12 min OK. The check trees
under `dev/priortab-check/` were deleted after reading.

## Fix round (review `dev/reviews/2026-11-11-priortab.md`, FIX)

Every finding was reproduced first: the reviewer's probes 2, 3 and 4,
rerun on the first lane build, gave output identical to the review's
logs (`dev/priortab-out/rev-probe*-lane.txt` against
`dev/reviews/priortab-review-probe*-lane.txt`, no difference outside the
library line). The fixed build's runs of the same probes are
`dev/priortab-out/fix-probe*-lane.txt`; the lane's own cases are
`dev/priortab-fixprobe.R` (`fixprobe-lane.txt`, seed 1, n = 120).

| finding | fix | evidence on the fixed build |
| --- | --- | --- |
| MAJOR 1, a user's bounds-only prior dropped | `drop_unbounding_slots()` drops only rows a table restated: `priorlist_from_rows(mark = TRUE)` and the brms-table reader set `from_table` on a density-free row whose `source` is `"default"`; a row by hand (no `source`) or marked `"user"` is kept | `set_prior("", class = "sd", lb = 0)` and the same for `sigma`: `fit$prior` kept; `frm_sample(fit)` and `frm_sample(prior = )` resolve the same rows, `"" sd lb 0` (and `"" sigma lb 0`); a passed-back table still leaves `fit$prior` NULL |
| MAJOR 2, closed box at a zero-density end | see "MAJOR 2: what was built" below | beta(2, 2) on `yneg`: no bound, "optimizer failed" as on base; `lb = 1e-8` 129.8616523 at z = 0.006822082; `lb = 1e-10` 129.8616523 at 0.006822074; `lb = 0` code 1 with the new warning, objective Inf |
| MAJOR 3, 4266 decided per frame | per row: `printed_flat` on every printed `"(flat)"` row of a table read as a prior; frmtmb.sample's `table_flat_rows()` refuses those whose slot carries a proper brms default and names each (`table_lists_defaults()` and `table_flat_lost()` removed) | H1 refused naming `(class=Intercept)`, `(class=sd)`; H4 refused naming `(class=Intercept)`; H3 and the 4266 cases (`""` on `sd`, `+ set_prior("", class = "sigma")`) read; a coerced fit table relabelled `source = "user"` still refused (test-flat-prior.R) |
| MINOR 4, `Ya` rows did not supersede | `spec_target()` keys `resp` by `brms_stan_name()` | probe4 T1 and T2: 8 effective rows (were 12) |
| MINOR 5, softplus, squareplus | `link_image()`: two far doublings that each shrink the value by a third or more go to 0 | probe6: `sigma lb` 0 under both, as brms; identity stays NA |
| MINOR 6, box classes | the box reads each entry's own scale through `internal_bound()`: slopes, an uncentered intercept, a dpar's own class; a centered intercept is not covered (documented, filed) | `uniform(0.6, 1)` on `sigma`: `opt$par` identical to `lb = 0.6, ub = 1`, sigma 0.6; `uniform(0, 0.5)` on `0 + Intercept`: identical to the explicit box, 0.5 |
| MINOR 7, slopes line | `flat_slope_count()` counts the slopes no prior reaches (probe rows from `default_prior()`); no line when 0 | class `b` on every slope: no line; `coef = "z"` only: "(flat) where your prior sets none"; no prior: the old line |
| MINOR 8, 1109 | the CLOSED line names the one spelling fixed; `coef = "Intercept"` on `y ~ 1` is refused already (checked, `dev/priortab-out/p8.R`); a bound with `coef` stays open as frmtmb's documented box | `dev/test-backlog.md` |
| NEWS | the drop names `frm_sample(fit)`'s effect, the box is stated as built, the guard has its own bullet; frmtmb.sample's NEWS gets a "Needs" line | |

### MAJOR 2: what was built, and where it departs from the instruction

The transform route (lb + exp(u) on the fit route) was not built: it
changes how every bounded parameter is optimized and reported (a box
at 0.2 for `uniform(0, 0.2)` is reached exactly today, and a transform
can only approach it), which is more than a few functions; no
measurement of interior-MAP fits under it was needed.

The fallback the coordinator named was "box only uniform(); refuse the
other bounded densities without lb/ub". The refusal was built, run,
and taken out, because it refuses a fit that base gets right:
`beta(2, 2)` on `b z` with no bound on `y` (slope near 0.26) fits on
rellib-r20 with code 0, objective 59.80059518, z = 0.2567826375
(`dev/reviews/priortab-review-probe2-base.txt`, B4), and the refusal
stopped it (`fix-probe2-lane.txt` of the refusing build). The lane rule
"a check that fires on a correct model is worse than none" decides it.
What is built instead:

- A parameter with no bound is boxed to its density's support only
  where the density is finite at that end (`uniform()`, and
  `exponential()` at 0; `support_box()` evaluates the density there).
  `beta(2, 2)`, `lognormal()`, `gamma(0.5, 1)` get no box: an inside
  mode is found as on base, and a slope the data push out fails loudly
  as on base ("NA/NaN gradient evaluation").
- A fit whose final stop sits on a finite bound with a non-finite
  objective is reported as not converged (`nonfinite_not_converged()`
  on the final `opt` in `fit_assembled()`, code 1, the message names
  the objective). First placed in `run_optimizer()`, it marked a
  nested-quadrature stage that ends at -Inf on its trust box and is
  recovered by later rounds (test-families-extreme.R, the GEV fit:
  final objective 214.7792, code 0 on base), so it moved to the final
  stop; without the bound condition it would also have warned on
  test-brms-families.R's NaN-everywhere family.
- The explicit `lb = 0` defect (a closed box at a zero-density end) is
  filed in `dev/test-backlog.md`, "Filed in the round of 2026-11-11",
  with the reproducer; the transform is named there as the cure.

If the coordinator still wants the refusal, it is the four lines
removed from `support_box()`; the case above is what it costs.

### Files touched in the fix round

`R/priors.R`, `R/fit.R` (the final-stop guard, one line before
`fit <- structure(` in `fit_assembled()` and the helper after
`run_optimizer()`; `run_optimizer()` itself is unchanged), frmtmb.sample
`R/sample.R`, `NEWS.md`, frmtmb.sample `NEWS.md`, the two new test
files, `dev/test-backlog.md` (CLOSED lines reworded; two items filed),
and one line of frmtmb.sample's `tests/testthat/test-flat-prior.R`
(lane sampfront's file): its expected text of this lane's refusal,
which now names rows, `"(class=Intercept), (class=sigma)"`.

### Tests of the fix round

Changed test files, each on base (rellib-r20) and the lane:

| file | rellib-r20 | lane |
| --- | --- | --- |
| test-prior-table-brms.R | pass=30 fail=30 err=5 warn=2 | pass=96, 0 fail |
| test-prior-table-sample.R | pass=7 fail=8 err=1 | pass=25, 0 fail |
| test-flat-prior.R (one expectation) | (message text of this lane) | pass=132, 0 fail |

Of the new assertions, two pass on base because base was right there
and the first lane build was not: a user's bounds-only prior kept on
the fit (MAJOR 1; the first build failed it, review probe2 D1 to D3),
and the H1/H4 refusals (MAJOR 3; review probe3 H1, H4).

Suites, one file per process, final build (`priortab-fix2.log`; every
`lib:` line names the lane library for frmtmb and frmtmb.sample):

| suite | files | pass | fail | err | skip | warn |
| --- | --- | --- | --- | --- | --- | --- |
| frmtmb | 318 | 23583 | 0 | 0 | 239 | 0 |
| frmtmb.sample | 88 | 4075 | 0 | 0 | 18 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1848 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 12 | 687 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.spatial | 8 | 478 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 15 | 556 | 0 | 0 | 0 | 0 |

Gated tier on the same 95 files (`priortab-gated2.log`, no skip):
frmtmb 57 files pass=7040, frmtmb.sample 38 files pass=2105, 0 fail,
0 err. Sessions (`priortab-sess2.log`; every LEAK line tagged
`[with load of ...]`): frmtmb.sample forward, reverse and preload
83/88 files pass=4075 each; frmtmb forward 305/318 pass=23583, reverse
pass=23748, preload pass=23851 (skip 241); all 0 fail, 0 err, 0 warn.

OpenBLAS 0.3.26 with LAPACK (`OPENBLAS_NUM_THREADS=4`, Rblas 74 and
Rlapack 616 symbols forwarded; `ob2-*.txt`; `dev/ciharden-out/`
deleted): test-prior-table-brms.R 96, test-prior-table-sample.R 25,
test-families-extreme.R 124, test-brms-families.R 1257,
test-flat-prior.R 132, test-prior-one-slot.R 90; 0 fail, 0 err, 0 warn
each.

The fix round's rounding-sensitive assertions over 20 ulp-perturbed
responses (`dev/priortab-ulps2.R`, `ulps2.txt`, seed 1): the box and
the explicit bounds give identical `opt$par` on `sigma` 20 of 20, with
sigma exactly 0.6; the uncentered intercept exactly 0.5 and the
exponential's slope exactly 0, 20 of 20; beta(2, 2) at `lb = 1e-8`
against `1e-10` differs by at most 3.06e-14 relative (the test allows
1e-8), with z at least 0.006822082 (the test asks for 1e-5); the
`lb = 0` stop is code 1, 20 of 20.

`R CMD check --as-cran` on the fix build: frmtmb `Status: 1 NOTE`
(tests 37 min OK), frmtmb.sample `Status: 1 NOTE` (tests 587 s OK); the
NOTE is the HTML manual's V8 math rendering on both.

Merge check, line endings normalized (`git merge-file` of each shared
file against every lane that changed it, inputs in the scratchpad):
`R/fit.R` with gpstop, names, nlstop, ordstart, quadseam; `R/sampling-api.R`
with names; frmtmb.sample `R/sample.R` with names, sampfront; `NAMESPACE`
with names, sampfront, silent: 0 conflicts each. `NEWS.md`,
frmtmb.sample `NEWS.md` and `dev/test-backlog.md` conflict with every
lane that adds a "development version" heading or a "Filed in the
round of 2026-11-11" section at the same line, as expected; they are
unions for the consolidation.
