# Lane bin2port: bin 2 of brms's suite reaches 0.71.0's families

Date 2026-10-14. Worktree `frmtmb-wt-bin2port`, branch `wt-bin2port`,
base 47cf659d (frmtmb 0.71.0). Build judged: the read-only base
library `C:/Users/adf44/source/r/rellib-r9` (frmtmb 0.71.0,
frmtmb.sample 0.19.0, frmtmb.eam 0.12.0), because this lane changed no
package code: the private library `wt-bin2port-lib` stays empty, and
"before" and "after" differ only in the files under `dev/`. brms
2.23.0, R 4.6.1 on Windows.

This is a measurement lane. Nothing under `R/`, `extensions/*/R`,
`tests/`, `man/` or `vignettes/` changed. Every table in a generated
block below is pasted by `dev/bin2port-fill.R` from the logs.

STATUS: after fix round 2 of 2 of the review
(`dev/reviews/2026-10-14-bin2port.md`, MERGEABLE AFTER FIXES; the
re-check found every round-1 finding fixed). The sections "Fix round
1" and "Fix round 2" say what changed.

## Result

The 22 rows of scopes 1 and 2 each have a true verdict. Of them, 17
statements (19 assertions) pass, 2 are divergences (domain), and 3 are
decisions (no quantile function). The multinomial rows `log_lik:504`
and `posterior_predict:349` stay gaps, blocked on lane catmulti. Of
the 25 rows of scope 3, 2 pass through a new wiener port, 15 had a
false class and
are now gaps (13 SAR and FCOR rows, 2 ARMA leave-one-out rows), and 8
stay "not ported (cost)" with the exact blocker in the row. The value
checks put 40 more shape passes beside a value comparison with brms's
own function, and find one frmtmb defect. The counts are in the
generated blocks "Rows moved" and "The bin-2 summary after".

## Fix round 2

1. **The newdata defect is a pointer, not a to-file item.** The
   re-check found it wider than this record said, with a silent wrong
   answer for binomial; lane ndaterm is fixing it in package code this
   round. "To file" item 3 is now a short pointer with the reviewer's
   summary.
2. **`posterior_predict:420` passes (shape).** The adapter answers
   brms's helper signature, so it now returns that helper's container:
   for `fun == "wiener"`, `brms_adapt_posterior_predict()` returns
   `matrix(out, ncol = 1)`, with a comment that brms's
   `posterior_predict_wiener()` alone returns ndraws x 1. No frmtmb
   value changes; the row's `nrow()` counts frmtmb.eam's draws, as a
   length check does for every other family. The boundary
   conditioning stays a note (in `dev/brmsport-bin2-verdicts-2p.R` and
   below). The row has no manual verdict, and the `container` artifact
   class is gone (nothing else used it). This supersedes fix round 1,
   item 7.

Generated after this round: bin 2 has 179 passes (37 value, 142
shape; 40 of the shape passes carry a value check that holds, 102
check only a shape), and bins 1 and 2 pass 570 of 823; the rate is 570
of 706 (80.7%), 570 of 709 (80.4%) over 0.71.0's fixed denominator
(blocks "The bin-2 summary after" and "rate709"). One row moved
against round 1: `posterior_predict:420`, cannot transfer (container)
to pass (shape). Two driver runs give identical outputs (block "md5"),
and the tier passes in every file (block "Tier after").

## Fix round 1

1. **`log_lik:130`, `:133` are gaps** (finding 1, accepted). The survey
   class of the ARMA log_lik block is `absent`, with the reason the
   review gives, and the `no-pointwise` class is gone from `ct_kind`.
   The rate's denominator is now 706 (generated: "Over every
   assertion but the 117 ... 569 of 706 pass"). The backlog text is
   under "To file".
2. **Per-entry bound in the value check** (finding 2, accepted).
   `brms_value_compare()` now also requires
   `d <= sqrt(eps) * pmax(abs(ref), 1)` for every entry. The 40 checks
   that held still hold. `dev/bin2port-tol.R` (the reviewer's cases,
   seed 534, block "tol") shows the old rule and the new rule on the
   logistic-normal block: one entry of |ref| 0.263 moved by 1e-6, 1e-4
   and 0.01 holds under `all.equal()` alone and fails under the new
   rule, as does every entry moved by 1e-6 (the `catlpdf` mutation the
   old rule let through at `log_lik:534`). The unmoved values hold
   under both.
3. **A value verdict names its message** (finding 3, accepted).
   `dev/brmsport-bin2-valuechecks-manual.tsv` has a `msg_has` column,
   the text the verdict was written for. The ledger stops when the
   recorded message does not contain it, and the generator passes it
   to `brms_value_check(msg_has = )`, so the tier fails too. Both are
   seen to fire (blocks "guard msg" and "guard ledger"): `log_lik:462`
   given the text it would have after the cumulative fix gives 83
   pass, 1 fail in the tier, and the ledger stops with "VALUE CHECK
   'log_lik:462' no longer matches its verdict's message" and writes
   nothing (the ledger's md5 is unchanged). Consolidation note under
   "For the consolidation".
4. **hurdle_cumulative `y = 0`** (finding 4, accepted). The defect
   record, the manual reason of `log_lik:462` and the reproduction
   now name it. `dev/bin2port-defect-cumulative.R` (block "defect"): at
   thresholds (0, 0, 1), mu 0, hu 0.2, y = 0, 1, 4 frmtmb gives NaN,
   NaN, -1.536405 where the exact values are -1.609438, -0.916291,
   -1.536405.
5. **Shape only when the shape call is the first argument** (finding
   5, accepted). `bin2_shape_re` is anchored: an `expect_equal()`,
   `expect_true()`, `expect_identical()` or `expect_equivalent()`
   whose first argument is a `dim()`, `length()`, `names()`,
   `dimnames()`, `nrow()`, `ncol()`, `is.numeric()` or
   `is.character()` call, or `expect_is()`, `expect_output()`,
   `expect_length()`, `expect_named()`, optionally inside `SW()`. The
   three rowSums rows are now value passes: 37 value and 141 shape
   (was 34 and 144), and 101 check only a shape. On the before ledger
   the new rule changes no row (30 value, 128 shape).
6. **K = 5 in probe1** (finding 6, accepted). A logistic-normal case
   with four log ratios and six distinct correlations agrees with brms
   to 2.9e-15 relative; the same correlations in column-major order,
   the control, differ by up to 0.86 relative (block "probe1"). The
   correlation path, `brms_adapt_cat_lncor()` and its use of frmtmb's
   `us_theta_cor()`, is checked by `dev/bin2port-probe1.R` only:
   brms's blocks set `prep$lncor`, which brms does not read, so no
   tier row reaches it. The categorical-like `posterior_predict` rows
   are shape and value-property checks only, so a category swap in
   the simulator passes the tier (the review's mutation: 55 pass, 0
   fail).
7. **`posterior_predict:420`** (finding 7; superseded by fix round 2,
   item 2, which makes it a pass). Measured with
   `dev/bin2port-probe6.R` and `-6b.R` (blocks "probe6", "probe6b"):
   frmtmb drops no dimension. The draws-level `posterior_predict()`
   gives 150 x 200 on a wiener fit and 150 x 1 on one gaussian row,
   and `sim_draw()`, which the adapter calls, returns a plain vector
   of one value per draw for wiener as for gaussian, which is its
   contract. The adapter's `as.vector()` therefore drops nothing
   either. The row fails because brms's INTERNAL helper
   `posterior_predict_wiener(i, prep)` returns an ndraws x 1 matrix
   (`rcontinuous(types = "q")`), unlike brms's other helpers. Neither
   branch of the coordinator's instruction applies (no frmtmb defect,
   no adapter drop), and reshaping in the adapter is not allowed, so
   the row is now "cannot transfer" with a new artifact class,
   `container`: the assertion tests the container of a brms-internal
   helper that no frmtmb function returns. It stays in the
   denominator. The reason states the boundary conditioning as a
   separate difference the row does not test. Found on the way:
   `predict(newdata =)` and `posterior_predict(newdata =)` on a wiener
   fit fail even when newdata carries the boundary column, while
   `fitted(newdata =)` reads it (now "To file" item 3, a pointer).
8. **`distributions:92`, `:102` are divergences (domain)** (finding 8).
9. **Wording** (finding 9). The NEWS bullet no longer claims the
   multinomial mapping works; "labeled"; "Not done" names
   `log_lik:580` (betareg) with the extraDistr rows.

I disagree with nothing in the review. On finding 7 the evidence
supported neither outcome first offered; fix round 2 adopts the
reviewer's recommendation.

## What changed

All in `dev/`, all files lane-env.md names for this lane, and one
NEWS bullet:

- `dev/brmsport-bin2-adapters.R`: the categorical-like adapter
  (`brms_adapt_cat_family()`, `_dpars()`, `_lncor()`, `_run()`,
  `_log_lik()`, `_posterior_predict()`); the G1 shims `dfrechet`,
  `pfrechet`, `rfrechet`, `dgen_extreme_value`, `pgen_extreme_value`,
  `rgen_extreme_value` and their helpers `brms_g1_p()` and
  `brms_g1_frechet_dp()`; the G2 shims `dcategorical` and
  `inv_link_categorical` (`brms_g2_cat_shims`); the wiener route (the
  `dec` term, frmtmb.eam's namespace, `brms_adapt_brms_skip()`); the
  value checks (`brms_value_mark()`, `brms_value_env()`,
  `brms_value_compare()`, `brms_value_check()`).
- `dev/brmsport-bin2-adapters-spec.R`: the new blocks (2 categorical
  maps, 2 wiener), the new shim names, the substitutions for RWiener,
  the shared shape pattern `bin2_shape_re` and the value-check
  candidate rule `bin2_value_check_ok()`.
- `dev/brmsport-bin2-gen.R`: emits `brms_value_mark()` and
  `brms_value_check()`, with the expected outcome and message text
  read from `dev/brmsport-bin2-valuechecks.tsv`.
- `dev/brmsport-bin2-ledger.R`: reads the value records, stops on a
  value check of a shape pass that does not hold without a manual
  verdict, on a stale one, and on one whose message no longer matches;
  writes `dev/brmsport-bin2-valuechecks.tsv`; adds the columns
  `value_check`, `value_verdict` and `value_class`, the summary
  section "Value checks beside the shape passes".
- `dev/brmsport-bin2-verdicts-2p.R`: the verdict moves of this lane and
  the manual verdicts of the value checks that do not hold (written to
  `dev/brmsport-bin2-valuechecks-manual.tsv`).
- `dev/brmsport-bin2-survey.R`: the G2 and G3 counterpart texts, and a
  class and reason of their own for the 11 G3 and G4 blocks left to
  the survey (`bin2port_g34`).
- Regenerated by `sh dev/brmsport-bin2-all.sh`: the ledger, summary,
  verdicts, survey and groups tables, the three generated 2P test files.
- `NEWS.md`: one bullet under `# frmtmb (development version)`.
- New scripts: `dev/bin2port-probe0.R` to `-probe6.R`, `-probe4b.R`,
  `-probe6b.R`, `dev/bin2port-tol.R`, `dev/bin2port-defect-cumulative.R`,
  `dev/bin2port-diff.R`, `dev/bin2port-fill.R`. Logs in
  `dev/bin2port-out/` and `dev/bin2port-before-log/` (gitignored). The
  reviewer's `dev/bin2port-rev-*.R` are not mine and were not edited.

## Scope 1 and 2: the rows ported

Each port evaluates frmtmb's own family. The adapter builds the family
with the bare constructor brms's block names, resolved through the
family's own `defer` slot on `prep$data$ncat` categories, which is the
route `frm()` takes. It maps brms's `mu1`, `mu2`, ... (and
logistic_normal's `sigma*`) by position, which is how brms's `get_Mu()`
reads them, onto frmtmb's per-category dpars, and evaluates
`row_lpdf()` or `sim_draw()`. Only frmtmb's own chart inverse,
`us_theta_cor()`, maps a logistic-normal correlation into frmtmb's
coordinates. The values agree with brms's own `log_lik_<family>()` on
the blocks' preps (`dev/bin2port-probe1.R`, seed 495, block
"probe1"): categorical 4.4e-16, dirichlet 5.3e-15 (refcat 1) and
3.6e-15 (refcat 2), dirichlet2 2.7e-15, logistic_normal 4.9e-12
relative, K = 4 with three correlations per draw 8.1e-15 relative,
K = 5 with six 2.9e-15 relative. 20000 logistic-normal draws have the
category means of brms's simulator to 2e-3 and the log-ratio
correlation 0.603 against 0.608.

| row | verdict | value or shape |
|---|---|---|
| `log_lik:495` categorical | pass | shape (dim); value check holds |
| `log_lik:516` dirichlet | pass | shape (dim); value check holds |
| `log_lik:523` dirichlet2 | pass | shape (dim); value check holds |
| `log_lik:534` logistic_normal | pass | shape (dim); value check holds |
| `posterior_predict:344` categorical | pass | shape (dim) |
| `posterior_predict:360`, `:368`, `:380` | pass | shape (dim) |
| `posterior_predict:361`, `:369`, `:381` | pass | value (row sums equal 1) |
| `distributions:527` dcategorical (loop, 2) | pass | value (closed form) and dim |
| `distributions:549` inv_link_categorical (loop, 2) | pass | value (closed form) and dim |
| `distributions:95` pfrechet | pass | shape (length); value check holds, 1.1e-16 |
| `distributions:134` dgen_extreme_value | pass | shape; value check holds |
| `distributions:137` pgen_extreme_value | pass | shape; value check holds |
| `distributions:141` rgen_extreme_value | pass | shape (draws, not checked) |
| `distributions:92` dfrechet(shape = 1) | divergence (domain) | |
| `distributions:102` rfrechet(shape = 1:10) | divergence (domain) | |
| `distributions:97`, `:100` qfrechet, `:139` qgen_extreme_value | cannot transfer: no-quantile (decision) | |
| `log_lik:504`, `posterior_predict:349` multinomial | cannot transfer: absent (gap), blocked on lane catmulti | |

`distributions:92` and `:102`: frmtmb's frechet(), like brms's frechet
family, takes the mean `mu = scale * gamma(1 - 1/nu)`, which is
infinite at shape 1. The shim refuses the call rather than hand the
density `mu = Inf`, where brms's length-only assertion would hold on
ten NaN values. The q rows take the class of the other q rows
(`:8`, `:68`, `:71`): SPEC.md:153-154 puts quantile functions outside
frmtmb, and the review accepts this (finding 8).

The multinomial rows fail at the block's own `multinomial()`, which
frmtmb 0.71.0 refuses ("needs the number of categories"), and the
adapter builds the family from the same bare constructor. With a
stand-in for lane catmulti's change (a multinomial whose `defer` reads
K from the response), the same adapter answers both rows and agrees
with brms's `log_lik_multinomial()` to 2.8e-14 (`dev/bin2port-probe5.R`,
block "probe5"). The NEWS bullet does not claim this mapping yet.

A brms test bug, found on the way: brms's logistic-normal blocks set
`prep$lncor`, which brms 2.23.0 never reads (`get_Sigma()` reads
`prep$dpars$lncor`), so brms evaluates them with uncorrelated log
ratios. brms's value is identical with and without `prep$lncor`, and
moves by up to 3.3e7 log units when the same draws are put in
`prep$dpars$lncor` (block "probe1"). The adapter reads what brms reads.

## Scope 3: the G3 and G4 rows

A run on 0.71.0 contradicts the survey's group reason, "frmtmb has
each of these models, so the row is transferable", for 15 of the 25
rows (`dev/bin2port-probe4.R`, block "probe4"):

- frmtmb has no `sar()` and no `fcor()`: `frm()` refuses both ("is a
  brms autocorrelation term that frmtmb does not support"). The 13 SAR
  and FCOR rows (`log_lik:150`, `:152`, `:157`, `:159`, `:174`, `:176`,
  `:230`, `posterior_predict:111`, `:113`, `:118`, `:120`, `:133`,
  `:135`) are gaps.
- frmtmb fits `arma(cov = TRUE)`, and `log_lik()` of its draws refuses
  it ("needs a likelihood that factors into one term per
  observation"), pointing to the `cov = FALSE` form, which it
  supports. brms's `log_lik_gaussian_time()` and
  `log_lik_student_time()` are the leave-one-out conditional densities
  of one observation given the rest of its group. `log_lik:130` and
  `:133` are gaps (fix round 1, finding 1); see "To file".

Ported, through frmtmb.eam's `wiener()` with brms's `dec` as its
`dec()` term:

- `log_lik:565` passes (shape); its value check holds against brms's
  `log_lik_wiener()` to 8.9e-16 (`dev/bin2port-probe3.R`).
- `posterior_predict:420` passes (shape, fix round 2): the adapter
  returns brms's `posterior_predict_wiener()` container, ndraws x 1,
  around frmtmb.eam's draws, and `nrow()` counts them. Note, a
  difference the row does not test: frmtmb.eam simulates the response
  time given the boundary the row ended at, and brms draws the
  boundary too; frmtmb.eam documents this (?frmtmb.eam::wiener,
  "Post-fit quantities against brms"). brms 2.23.0 skips its own block
  for an RWiener fault, which frmtmb.eam does not use, so the port
  runs it in place of the skip.

Still not ported (cost), 8 rows, each with its blocker in the ledger:

- `log_lik:107`, `:109`, `posterior_predict:69`, `:72` (multivariate
  with residual correlation) and `posterior_predict:93`, `:97`
  (`arma(cov = TRUE)` draws): frmtmb has the models, and `log_lik()`
  of draws answers per row for the multivariate one, but the density
  (`rescor_row_loglik()`) and the structured simulators
  (`sim_autocor_rows()`, the rescor draw) read a fit: the correlation
  parameters from `fit$estimates`, the structure from `fit$frame`. No
  exported function takes brms's per-draw matrices, so a port must
  build a fit per draw or call internals.
- `log_lik:603`, `posterior_predict:440` (custom family): brms
  dispatches to a user's prep-level `log_lik_<name>(i, prep)` found by
  name; frmtmb's `custom_family()` is `frmtmb_family()`, which takes
  `lpdf(y, dpars, aterms)` at construction. A port must restate the
  user's function in that signature, a 2T-style restatement.

## Scope 5: value checks

Where brms's assertion checks only a shape and the line before it (or,
for wiener, the call inside it) computed a log density or a CDF
through the adapter, the generator adds a check after the assertion.
It evaluates the same line, as brms wrote it, against brms's own
function (`brms:::log_lik_<family>()`, or brms's d or p function under
the name the header binds to a shim) on the same arguments, from the
same random number state. It holds when `all.equal()` holds
(testthat's tolerance), every entry is within
`sqrt(eps) * max(|ref|, 1)` of brms's, and the NA and infinite
positions are identical. It is recorded as kind "value" under the
brms row's id, never in that row's verdict, and the tier asserts its
recorded outcome and, for one that does not hold, the message text
its verdict names. Simulators are not checked: two simulators agree in
distribution, not draw by draw, so a draw-by-draw comparison would
test the order of the random calls.

Of the 142 shape passes, 47 carry a value check: 40 hold, 1 finds a
defect, 2 find divergences, and 4 cannot run their reference (brms
needs extraDistr or betareg, which are not installed). Of the 179
bin-2 passes, 37 compare a value as brms wrote them, 40 more carry a
value check that holds, and 102 check only a shape (before: 30 of 158
compared a value). The table "Every value check" lists all 51 checks,
4 of them on rows that do not pass. Of the 40 that hold, the largest
relative difference is 1.26e-9 (`log_lik:475`, acat under the probit
link); every other is 3.9e-14 or smaller.

The three that differ:

- `log_lik:462`, **defect** (low). frmtmb's cumulative log probability
  of an end category is NaN on doubles when the threshold pair the row
  is clamped to is tied, where the category has positive probability,
  and so is hurdle_cumulative's `y = 0` row at a tie of the lowest
  pair. brms's prep has all thresholds 0, and brms's values for
  categories 1 and 4 are finite. Reproduction
  `dev/bin2port-defect-cumulative.R` (block "defect"): at thresholds
  (0, 0, 0), (-1, -1, 1) and (-1, 1, 1) every end category next to the
  tie is NaN where the exact value is finite, and hurdle_cumulative at
  (0, 0, 1) gives NaN at `y = 0` and `y = 1`. With distinct thresholds
  frmtmb is exact. On the RTMB tape the same sum is exact (-1.650688),
  because RTMB folds the masked `0 * -Inf` term (rtmb-pitfalls item
  22), so fits are not affected; numeric evaluation is, which is what
  `log_lik()` of draws runs. Cause: `ord_cumulative_logpmf()`
  (R/families.R) evaluates the interior pair `ord_log_interior(a, b)`
  on every row and masks it with `(1 - i1 - iK) *`, and
  `fam_hurdle_cumulative()` calls it on `pmax(y, 1)` and masks with
  `(1 - i0) *`; at a tie the pair is `log(0)`. The fix and the test
  the fixing lane needs are under "To file". After the fix the check
  stays "differs" with a new message (frmtmb 0 NA, brms 200), so the
  tier and the ledger stop on it and it must be re-judged.
- `log_lik:271`, divergence (precision): discrete_weibull at y = 12.
  brms's pmf underflows to -Inf in 17 of 25 draws (brms's block says
  so); frmtmb's log-scale density is finite there and matches the
  same closed form taken on the log scale (`dev/bin2port-probe2.R`).
- `log_lik:478`, divergence (domain): hurdle_cumulative at a
  zero-width category, frmtmb -Inf (the exact log) and brms NaN.

## To file

For `dev/test-backlog.md` at consolidation, items 1 and 2 (item 3 is
a pointer to a lane that is fixing it):

1. **log_lik() of draws for residual covariance blocks.** frmtmb fits
   `arma(cov = TRUE)`, `cosy()` and `unstr()`, and `log_lik()` of
   their draws refuses the pointwise matrix ("needs a likelihood that
   factors into one term per observation"). brms returns the
   leave-one-out conditional density of each observation given the
   rest of its group (the non-factorized LOO of Burkner, Gabry and
   Vehtari), the conditional normal or Student t of
   `log_lik_gaussian_time()` and `log_lik_student_time()`; brms's own
   suite checks the FCOR Student t case in closed form
   (`tests.log_lik.R:179-229`). The refusal's reason ("leaving one out
   would drop part of a joint density") does not hold for a
   conditional density. Bin-2 rows `log_lik:130`, `:133` are gaps
   until it exists. Neither SPEC.md, the compat table nor
   `vignettes/brms-migration.Rmd` lists the refusal.
2. **The cumulative tied-threshold NaN.** In `ord_cumulative_logpmf()`
   give a masked row a legal pair from data alone,
   `a <- q(disc * (tau[ym] - eta)) + (i1 + iK)`, which adds 0 on the
   interior rows; the review measured it exact at every tie for logit,
   probit and cloglog, bitwise unchanged at distinct thresholds, and a
   cumulative fit unchanged (`dev/bin2port-rev-cumfix.R`). The fixing
   lane needs a test, seen to fail on rellib-r9, that pins cumulative
   `y = 1` and `y = K` and hurdle_cumulative `y = 0` at tied
   thresholds, and the OpenBLAS run lane-env asks of numerical code.
3. **Pointer, being fixed in lane ndaterm: `predict(newdata =)` and
   `posterior_predict(newdata =)` drop every addition term but
   `trunc()`, `rate()` and `thres(gr = )`.** Not a to-file item. The
   reviewer's summary (`dev/reviews/2026-10-14-bin2port.md`, "The
   newdata defect"; `dev/bin2port-rev-newdata.R`,
   `dev/bin2port-rev-newdata2.R`, seed 7):
   - The serious case is silent: on
     `frm(bf(y | trials(N) ~ x), binomial())` with newdata `N = 1000`
     and `N = 2`, `predict(newdata =, summary = FALSE)` gives draw
     means 0.59 and 0.56 with maximum 1, so every draw is Bernoulli.
     frmtmb.sample's `posterior_predict(newdata =)` gives 0.587 and
     0.547. `simulate(newdata =)` (575.6 and 1.17) and
     `posterior_epred()` (562.2 and 1.12) are right.
   - The wiener case this lane found raises an error instead: "needs
     the boundary each row ended at, which is missing", though newdata
     carries the `dec()` column, which `fitted(newdata =)` reads
     (blocks "probe6", "probe6b", seed 420).
   - The three sites: core's `predict_simulate()`
     (`R/predict-brms.R:599-610`), frmtmb.sample's
     `posterior_predict.frmtmb_draws()`
     (`extensions/frmtmb.sample/R/methods-draws.R:1081-1092`), and
     core's `arma_cond_fill_dpars()` (`R/autocor.R:1362`), which reads
     `trunc()` only.
   - The fix: take `aterms_for_newdata(rspec, newdata)` whenever the
     response has addition terms, as `simulate(newdata =)` and
     `fitted()` do, pinned by a binomial `trials()` and a wiener
     `dec()` newdata test in core and frmtmb.sample.

## Evidence the instruments work

- Before: the driver on the unchanged tree reproduces the tracked
  ledger exactly (`git status` clean after the run), and its logs are
  kept in `dev/bin2port-before-log/`.
- The tier of the generated files passes on the new verdicts (block
  "Tier after"). It fails when a recorded value-check outcome is wrong
  in either direction (block "guard", round 0: `log_lik:85` set to
  "differs" and `log_lik:462` set to "holds" give exactly two
  failures, in those two blocks), and when a verdict's message text no
  longer matches (block "guard msg", fix round 1). The ledger stops on
  the same mismatch without writing (block "guard ledger").
  Construction: `sed` of the expectation into a scratch copy of the
  generated file, run with `dev/brmsport-bin2-run.R`; for the ledger,
  the manual tsv edited, the ledger run, the tsv restored.
- Two full driver runs after fix round 1 give identical ledgers,
  summaries, tables and generated files (md5, block "md5"). The
  value-check messages carry differences near 1e-16, which could
  change by an ulp between the two core types of this machine
  (rtmb-pitfalls item 21); the outcomes cannot.

## Not done, and why

- No package code changed, so no package suite, OpenBLAS pass or
  `R CMD check` was run: no path under `R/` or `extensions/` changed,
  and the generated bin-2 files live under `dev/`, outside every
  package's tests. Bin 1 was not regenerated (lane-env: once, at
  consolidation).
- The cumulative defect is filed, not fixed. The newdata defect is
  being fixed in lane ndaterm.
- The 8 cost rows: each needs a fit per draw or internals, or a
  restatement of a user's function.
- No value checks on simulators (above), or on the 2T methods rows,
  whose objects are frm_sample() draws against brms's fixtures.
- The 4 value checks whose reference cannot run: `distributions:123`,
  `:177` and `log_lik:436` need extraDistr, and `log_lik:580` needs
  betareg. Installing them would let the checks run (and the brms row
  `log_lik:267`, class reference, also needs extraDistr).

## For the consolidation

- Overlaps: `dev/brmsport-bin2-*` (this lane is the named exception)
  and `NEWS.md`.
- Lane catmulti's bare `multinomial()` should flip `log_lik:504` and
  `posterior_predict:349` to pass: remove their manual verdicts in
  `dev/brmsport-bin2-verdicts-2p.R`. Its value check `log_lik:504`
  then moves from "error" to "holds", which the regeneration records.
  Its `categorical(refcat = )` is not needed by brms's rows (refcat
  1); the adapter refuses another refcat for categorical by name.
- After lane catmulti's cumulative fix, `log_lik:462` must be
  re-judged: its message no longer contains "NA in 200 values from the
  adapter", so the ledger stops. Expected: a divergence (domain), as
  `log_lik:478` (frmtmb -Inf, brms NaN at the zero-width categories),
  or "holds" if the fix also changes those. Write the new verdict and
  its `msg_has` in `dev/brmsport-bin2-verdicts-2p.R`.
- The ledger has three new columns, so a script that reads
  `dev/brmsport-bin2-ledger.tsv` by position would break; the ones in
  `dev/` read by name.
- The rate's denominator moves from 709 to 706: the 3 q rows are
  decisions.
- Lane ndaterm's newdata fix reaches no bin-2 row: the adapter calls
  `sim_draw()` with brms's addition terms, not the `newdata` route.

## Generated blocks

### The bin-2 summary before (rellib-r9, the unchanged tree)

<!-- begin before -->

### Bin 2 per file and bin

| file (bin) | assertions | pass (value / shape) | defect | divergence | cannot transfer: decision / gap / cost / artifact |
|---|---|---|---|---|---|
| `tests.brmsfit-helpers.R` (2P) | 8 | 0 (0 / 0) | 0 | 0 | 8 / 0 / 0 / 0 |
| `tests.brmsfit-methods.R` (2T) | 58 | 31 (1 / 30) | 0 | 8 | 0 / 17 / 0 / 2 |
| `tests.brmsterms.R` (2P) | 18 | 0 (0 / 0) | 0 | 0 | 18 / 0 / 0 / 0 |
| `tests.distributions.R` (2P) | 113 | 43 (12 / 31) | 0 | 2 | 40 / 15 / 13 / 0 |
| `tests.log_lik.R` (2P) | 65 | 40 (17 / 23) | 0 | 2 | 0 / 4 / 17 / 2 |
| `tests.posterior_predict.R` (2P) | 66 | 44 (0 / 44) | 0 | 0 | 0 / 3 / 19 / 0 |
| `tests.priors.R` (2T) | 1 | 0 (0 / 0) | 0 | 1 | 0 / 0 / 0 / 0 |
| **total** | 329 | 158 (30 / 128) | 0 | 13 | 66 / 39 / 49 / 4 |

### Bin 2 outcome by class

| bin | outcome | class | kind | assertions |
|---|---|---|---|---|
| 2T | pass | - | - | 31 |
| 2T | divergence | policy | - | 9 |
| 2T | cannot transfer | absent | gap | 17 |
| 2T | cannot transfer | fixture | artifact | 1 |
| 2T | cannot transfer | upstream | artifact | 1 |
| 2P | pass | - | - | 113 |
| 2P | pass | per-line | - | 14 |
| 2P | divergence | domain | - | 3 |
| 2P | divergence | parameterization | - | 1 |
| 2P | cannot transfer | absent | gap | 7 |
| 2P | cannot transfer | adapter | artifact | 1 |
| 2P | cannot transfer | brms-backends | decision | 3 |
| 2P | cannot transfer | brms-internal | decision | 48 |
| 2P | cannot transfer | no-cdf | gap | 15 |
| 2P | cannot transfer | no-density | decision | 12 |
| 2P | cannot transfer | no-quantile | decision | 3 |
| 2P | cannot transfer | not-ported | cost | 49 |
| 2P | cannot transfer | reference | artifact | 1 |

Bin 2T passes: 31 of 59 (52.5%), judged on the target object.
Bin 2P: 178 of 270 assertions run through the adapter, 127 of them pass (the rest it refused or they failed); the other 92 are classified from the survey.
Of the 158 bin-2 passes, 30 compare a value and 128 check only a shape, length, names or type.
The 58 methods rows of 2T on the maximum-likelihood fit (the restatement, not counted): divergence 52, pass 6.

### Bins 1 and 2 together

| outcome | bin 1 | bin 2T | bin 2P | bins 1 and 2 |
|---|---|---|---|---|
| pass | 391 | 31 | 127 | 549 |
| defect | 1 | 0 | 0 | 1 |
| divergence | 51 | 9 | 4 | 64 |
| cannot transfer: no counterpart (decision) | 48 | 0 | 66 | 114 |
| cannot transfer: absent (gap) | 2 | 17 | 22 | 41 |
| cannot transfer: not ported (cost) | 0 | 0 | 49 | 49 |
| cannot transfer: port artifact | 1 | 2 | 2 | 5 |
| **total** | **494** | **59** | **270** | **823** |

Bins 1 and 2 pass: 549 of 823 (66.7%); bin 1 alone 391 of 494, as dev/brmsport-ledger.tsv records it (regenerated once at consolidation).
Over every assertion but the 114 that frmtmb has no counterpart for by decision, 549 of 709 pass (77.4%). The 41 gaps, 49 rows not ported (cost) and 5 port artifacts stay in that denominator.
The rate's denominator excludes exactly the "cannot transfer" rows of these classes: brms-backends 3, brms-internal 63, mcmc 4, no-density 12, no-draws 3, no-quantile 3, stan 26.
<!-- end before -->

<!-- begin rate709before -->
Over the 709 assertions of 0.71.0's denominator (every assertion but its 114 decisions): 549 of 709 pass (77.4%).
<!-- end rate709before -->

### The bin-2 summary after (rellib-r9, this lane's dev files)

<!-- begin after -->

### Bin 2 per file and bin

| file (bin) | assertions | pass (value / shape) | defect | divergence | cannot transfer: decision / gap / cost / artifact |
|---|---|---|---|---|---|
| `tests.brmsfit-helpers.R` (2P) | 8 | 0 (0 / 0) | 0 | 0 | 8 / 0 / 0 / 0 |
| `tests.brmsfit-methods.R` (2T) | 58 | 31 (1 / 30) | 0 | 8 | 0 / 17 / 0 / 2 |
| `tests.brmsterms.R` (2P) | 18 | 0 (0 / 0) | 0 | 0 | 18 / 0 / 0 / 0 |
| `tests.distributions.R` (2P) | 113 | 51 (16 / 35) | 0 | 4 | 43 / 15 / 0 / 0 |
| `tests.log_lik.R` (2P) | 65 | 45 (17 / 28) | 0 | 2 | 0 / 13 / 3 / 2 |
| `tests.posterior_predict.R` (2P) | 66 | 52 (3 / 49) | 0 | 0 | 0 / 9 / 5 / 0 |
| `tests.priors.R` (2T) | 1 | 0 (0 / 0) | 0 | 1 | 0 / 0 / 0 / 0 |
| **total** | 329 | 179 (37 / 142) | 0 | 15 | 69 / 54 / 8 / 4 |

### Bin 2 outcome by class

| bin | outcome | class | kind | assertions |
|---|---|---|---|---|
| 2T | pass | - | - | 31 |
| 2T | divergence | policy | - | 9 |
| 2T | cannot transfer | absent | gap | 17 |
| 2T | cannot transfer | fixture | artifact | 1 |
| 2T | cannot transfer | upstream | artifact | 1 |
| 2P | pass | - | - | 134 |
| 2P | pass | per-line | - | 14 |
| 2P | divergence | domain | - | 5 |
| 2P | divergence | parameterization | - | 1 |
| 2P | cannot transfer | absent | gap | 22 |
| 2P | cannot transfer | adapter | artifact | 1 |
| 2P | cannot transfer | brms-backends | decision | 3 |
| 2P | cannot transfer | brms-internal | decision | 48 |
| 2P | cannot transfer | no-cdf | gap | 15 |
| 2P | cannot transfer | no-density | decision | 12 |
| 2P | cannot transfer | no-quantile | decision | 6 |
| 2P | cannot transfer | not-ported | cost | 8 |
| 2P | cannot transfer | reference | artifact | 1 |

Bin 2T passes: 31 of 59 (52.5%), judged on the target object.
Bin 2P: 184 of 270 assertions run through the adapter, 148 of them pass (the rest it refused or they failed); the other 86 are classified from the survey.
Of the 179 bin-2 passes, 37 compare a value and 142 check only a shape, length, names or type.
The 58 methods rows of 2T on the maximum-likelihood fit (the restatement, not counted): divergence 52, pass 6.

### Value checks beside the shape passes

A value check compares the value the line before a shape assertion computed through the adapter with brms's own function on the same arguments (brms:::log_lik_<family>(), or brms's d or p function). It is recorded beside the brms row and does not change that row's verdict. Simulators are not checked: two of them agree in distribution, not draw by draw.

| file | shape passes | with a value check | holds | defect | divergence | reference cannot run |
|---|---|---|---|---|---|---|
| `tests.brmsfit-methods.R` | 30 | 0 | 0 | 0 | 0 | 0 |
| `tests.distributions.R` | 35 | 19 | 17 | 0 | 0 | 2 |
| `tests.log_lik.R` | 28 | 28 | 23 | 1 | 2 | 2 |
| `tests.posterior_predict.R` | 49 | 0 | 0 | 0 | 0 | 0 |
| **total** | 142 | 47 | 40 | 1 | 2 | 4 |

Of the 179 bin-2 passes, 37 compare a value as brms wrote them, 40 more check a shape and carry a value check that holds against brms's own function, and 102 check only a shape.
Value checks that do not hold: `distributions:123` (cannot transfer, reference); `distributions:177` (cannot transfer, reference); `log_lik:271` (divergence, precision); `log_lik:436` (cannot transfer, reference); `log_lik:462` (defect, tied-thresholds); `log_lik:478` (divergence, domain); `log_lik:580` (cannot transfer, reference).

### Bins 1 and 2 together

| outcome | bin 1 | bin 2T | bin 2P | bins 1 and 2 |
|---|---|---|---|---|
| pass | 391 | 31 | 148 | 570 |
| defect | 1 | 0 | 0 | 1 |
| divergence | 51 | 9 | 6 | 66 |
| cannot transfer: no counterpart (decision) | 48 | 0 | 69 | 117 |
| cannot transfer: absent (gap) | 2 | 17 | 37 | 56 |
| cannot transfer: not ported (cost) | 0 | 0 | 8 | 8 |
| cannot transfer: port artifact | 1 | 2 | 2 | 5 |
| **total** | **494** | **59** | **270** | **823** |

Bins 1 and 2 pass: 570 of 823 (69.3%); bin 1 alone 391 of 494, as dev/brmsport-ledger.tsv records it (regenerated once at consolidation).
Over every assertion but the 117 that frmtmb has no counterpart for by decision, 570 of 706 pass (80.7%). The 56 gaps, 8 rows not ported (cost) and 5 port artifacts stay in that denominator.
The rate's denominator excludes exactly the "cannot transfer" rows of these classes: brms-backends 3, brms-internal 63, mcmc 4, no-density 12, no-draws 3, no-quantile 6, stan 26.
<!-- end after -->

<!-- begin rate709 -->
Over the 709 assertions of 0.71.0's denominator (every assertion but its 114 decisions): 570 of 709 pass (80.4%).
<!-- end rate709 -->

### Rows moved (dev/bin2port-diff.R)

<!-- begin diff -->
Rows moved: 39 statements, 41 assertions

| row | n | before | after |
|---|---|---|---|
| `distributions:92` | 1 | cannot transfer (not-ported) | divergence (domain) |
| `distributions:95` | 1 | cannot transfer (not-ported) | pass (shape) |
| `distributions:97` | 1 | cannot transfer (not-ported) | cannot transfer (no-quantile) |
| `distributions:100` | 1 | cannot transfer (not-ported) | cannot transfer (no-quantile) |
| `distributions:102` | 1 | cannot transfer (not-ported) | divergence (domain) |
| `distributions:134` | 1 | cannot transfer (not-ported) | pass (shape) |
| `distributions:137` | 1 | cannot transfer (not-ported) | pass (shape) |
| `distributions:139` | 1 | cannot transfer (not-ported) | cannot transfer (no-quantile) |
| `distributions:141` | 1 | cannot transfer (not-ported) | pass (shape) |
| `distributions:527` | 2 | cannot transfer (not-ported) | pass (value) |
| `distributions:549` | 2 | cannot transfer (not-ported) | pass (value) |
| `log_lik:130` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `log_lik:133` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `log_lik:150` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `log_lik:152` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `log_lik:157` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `log_lik:159` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `log_lik:174` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `log_lik:176` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `log_lik:230` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `log_lik:495` | 1 | cannot transfer (not-ported) | pass (shape) |
| `log_lik:516` | 1 | cannot transfer (not-ported) | pass (shape) |
| `log_lik:523` | 1 | cannot transfer (not-ported) | pass (shape) |
| `log_lik:534` | 1 | cannot transfer (not-ported) | pass (shape) |
| `log_lik:565` | 1 | cannot transfer (not-ported) | pass (shape) |
| `posterior_predict:111` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `posterior_predict:113` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `posterior_predict:118` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `posterior_predict:120` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `posterior_predict:133` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `posterior_predict:135` | 1 | cannot transfer (not-ported) | cannot transfer (absent) |
| `posterior_predict:344` | 1 | cannot transfer (not-ported) | pass (shape) |
| `posterior_predict:360` | 1 | cannot transfer (not-ported) | pass (shape) |
| `posterior_predict:361` | 1 | cannot transfer (not-ported) | pass (value) |
| `posterior_predict:368` | 1 | cannot transfer (not-ported) | pass (shape) |
| `posterior_predict:369` | 1 | cannot transfer (not-ported) | pass (value) |
| `posterior_predict:380` | 1 | cannot transfer (not-ported) | pass (shape) |
| `posterior_predict:381` | 1 | cannot transfer (not-ported) | pass (value) |
| `posterior_predict:420` | 1 | cannot transfer (not-ported) | pass (shape) |

Moves by kind (assertions):

| before | after | assertions |
|---|---|---|
| cannot transfer (not-ported) | cannot transfer (absent) | 15 |
| cannot transfer (not-ported) | cannot transfer (no-quantile) | 3 |
| cannot transfer (not-ported) | divergence (domain) | 2 |
| cannot transfer (not-ported) | pass (shape) | 14 |
| cannot transfer (not-ported) | pass (value) | 7 |

Value checks on shape passes: 47 statements; cannot transfer 4, defect 1, divergence 2, pass 40 
<!-- end diff -->

### Every value check (dev/brmsport-bin2-valuechecks.tsv)

<!-- begin values -->
| row | brms row | value check | verdict | message |
|---|---|---|---|---|
| `distributions:50` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 3.33e-16, max rel 1.77e-15 |
| `distributions:61` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 3.04e-18, max rel 3.69e-15 |
| `distributions:80` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 0, max rel 0 |
| `distributions:92` | divergence | error | (row does not pass) | NOT CHECKED: the adapter's value does not exist |
| `distributions:95` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 1.11e-16, max rel 1.22e-16 |
| `distributions:109` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 0, max rel 0 |
| `distributions:111` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 0, max rel 0 |
| `distributions:123` | pass | error | cannot transfer (reference) | REFERENCE ERROR: Please install the 'extraDistr' package. |
| `distributions:134` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 4.16e-17, max rel 4.72e-16 |
| `distributions:137` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 1.11e-16, max rel 2.17e-16 |
| `distributions:148` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 6.94e-18, max rel 2.97e-16 |
| `distributions:162` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 5.55e-17, max rel 1.29e-16 |
| `distributions:167` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 4.16e-16, max rel 2.03e-15 |
| `distributions:172` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 5.55e-17, max rel 5.26e-16 |
| `distributions:177` | pass | error | cannot transfer (reference) | REFERENCE ERROR: Please install the 'extraDistr' package. |
| `distributions:183` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 2.22e-16, max rel 1.48e-16 |
| `distributions:193` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 0, max rel 0 |
| `distributions:198` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 5.55e-17, max rel 1.89e-16 |
| `distributions:203` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 0, max rel 0 |
| `distributions:208` | pass | holds | pass | HOLDS: 10 values (0 NA, 0 infinite), max abs diff 0, max rel 0 |
| `log_lik:85` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 1.98e-16 |
| `log_lik:88` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 1.78e-15, max rel 3.08e-16 |
| `log_lik:271` | pass | differs | divergence (precision) | DIFFERS: 25 values (0 NA, 17 infinite), max abs diff 1.73e-06, max rel 2.37e-09 |
| `log_lik:389` | pass | holds | pass | HOLDS: 15 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 3.91e-16 |
| `log_lik:392` | cannot transfer | error | (row does not pass) | NOT CHECKED: the adapter's value does not exist |
| `log_lik:414` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 3.55e-15, max rel 2.08e-16 |
| `log_lik:417` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 5.77e-15, max rel 2.89e-15 |
| `log_lik:420` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 0, max rel 0 |
| `log_lik:423` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 4.44e-16, max rel 1.27e-16 |
| `log_lik:426` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 0, max rel 0 |
| `log_lik:429` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 1.42e-14, max rel 5.43e-16 |
| `log_lik:433` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 4.44e-16, max rel 1.23e-16 |
| `log_lik:436` | pass | error | cannot transfer (reference) | REFERENCE ERROR: Please install the 'extraDistr' package. |
| `log_lik:440` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 2.22e-16, max rel 1.68e-16 |
| `log_lik:443` | pass | holds | pass | HOLDS: 50 values (0 NA, 0 infinite), max abs diff 0, max rel 0 |
| `log_lik:462` | pass | differs | defect (tied-thresholds) | DIFFERS: NA in 200 values from the adapter, 200 from brms |
| `log_lik:465` | pass | holds | pass | HOLDS: 400 values (0 NA, 0 infinite), max abs diff 9.95e-14, max rel 1.61e-14 |
| `log_lik:468` | pass | holds | pass | HOLDS: 400 values (0 NA, 0 infinite), max abs diff 3.73e-14, max rel 3.89e-14 |
| `log_lik:471` | pass | holds | pass | HOLDS: 400 values (0 NA, 0 infinite), max abs diff 3.55e-15, max rel 2.26e-14 |
| `log_lik:475` | pass | holds | pass | HOLDS: 400 values (0 NA, 0 infinite), max abs diff 4.46e-08, max rel 1.26e-09 |
| `log_lik:478` | pass | differs | divergence (domain) | DIFFERS: NA in 0 values from the adapter, 50 from brms |
| `log_lik:495` | pass | holds | pass | HOLDS: 400 values (100 NA, 0 infinite), max abs diff 4.44e-16, max rel 6.14e-15 |
| `log_lik:504` | cannot transfer | error | (row does not pass) | NOT CHECKED: the adapter's value does not exist |
| `log_lik:510` | cannot transfer | error | (row does not pass) | NOT CHECKED: the adapter's value does not exist |
| `log_lik:516` | pass | holds | pass | HOLDS: 400 values (0 NA, 0 infinite), max abs diff 3.55e-15, max rel 7.68e-16 |
| `log_lik:523` | pass | holds | pass | HOLDS: 400 values (0 NA, 0 infinite), max abs diff 4.44e-15, max rel 1.26e-14 |
| `log_lik:534` | pass | holds | pass | HOLDS: 400 values (0 NA, 0 infinite), max abs diff 7.45e-09, max rel 1.21e-14 |
| `log_lik:547` | pass | holds | pass | HOLDS: 90 values (0 NA, 0 infinite), max abs diff 2.22e-16, max rel 1.41e-16 |
| `log_lik:550` | pass | holds | pass | HOLDS: 90 values (60 NA, 0 infinite), max abs diff 4.44e-16, max rel 2.03e-16 |
| `log_lik:565` | pass | holds | pass | HOLDS: 5 values (0 NA, 0 infinite), max abs diff 8.88e-16, max rel 3.26e-16 |
| `log_lik:580` | pass | error | cannot transfer (reference) | REFERENCE ERROR: Please install the 'betareg' package. |

51 value checks: differs 3, error 8, holds 40.
<!-- end values -->

### Tier before

<!-- begin tierbefore -->
```
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-distributions.R pass=56 fail=0 err=0 skip=0 blocks=15
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=52 fail=0 err=0 skip=0 blocks=11
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-methods.R pass=58 fail=0 err=0 skip=0 blocks=13
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=54 fail=0 err=0 skip=0 blocks=11
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-priors.R pass=1 fail=0 err=0 skip=0 blocks=1
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb.sample 0.19.0 
RESULT frmtmb.sample test-brms-suite-bin2-methods.R pass=58 fail=0 err=0 skip=0 blocks=13
```
<!-- end tierbefore -->

### Tier after

<!-- begin tierafter -->
```
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-distributions.R pass=78 fail=0 err=0 skip=0 blocks=17
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=84 fail=0 err=0 skip=0 blocks=12
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-methods.R pass=58 fail=0 err=0 skip=0 blocks=13
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=55 fail=0 err=0 skip=0 blocks=12
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-priors.R pass=1 fail=0 err=0 skip=0 blocks=1
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb.sample 0.19.0 
RESULT frmtmb.sample test-brms-suite-bin2-methods.R pass=58 fail=0 err=0 skip=0 blocks=13
```
<!-- end tierafter -->

### guard (round 0, flipped outcomes)

<!-- begin guard -->
```
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=82 fail=2 err=0 skip=0 blocks=12
FAILING BLOCK:
  log_lik of aysm_laplace models runs without errors
  log_lik for ordinal models runs without erros
```
<!-- end guard -->

### guard msg (fix round 1, a stale message text in the tier)

<!-- begin guardmsg -->
```
lib: C:/Users/adf44/source/r/rellib-r9 frmtmb 0.71.0 frmtmb 0.71.0 
RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=83 fail=1 err=0 skip=0 blocks=12
FAILING BLOCK:
  log_lik for ordinal models runs without erros
```
<!-- end guardmsg -->

### guard ledger (fix round 1, a stale message text in the ledger)

<!-- begin guardledger -->
```
VALUE CHECK 'log_lik:462' no longer matches its verdict's message 'NA in 0 values from the adapter, 200 from brms': DIFFERS: NA in 200 values from the adapter, 200 from brms
Error: 1 problem(s); the ledger was not written
Execution halted
```
<!-- end guardledger -->

### md5

<!-- begin md5 -->
```
a6815f1d25f76422e944d385ef9093e7 *dev/brmsport-bin2-ledger.tsv
bae8caecf7ac268030ef5542b3211b96 *dev/brmsport-bin2-summary.md
ef85fbc71811fd0e24d81343929a17bd *dev/brmsport-bin2-verdicts.tsv
779e1278543fb1ef7f1590c350d28d84 *dev/brmsport-bin2-valuechecks.tsv
eef46a0cef462e793adc85db50a4ab15 *dev/brmsport-bin2-valuechecks-manual.tsv
3841f11196cdec316c167c9d46337e03 *dev/brmsport-bin2-verdicts-manual-2p.tsv
7930f070a6c9322a468dc61437e50248 *dev/brmsport-bin2-survey.tsv
55f71870a4b3335980d689d0018fd372 *dev/brmsport-bin2-groups.tsv
5be80200dee759687b5727d52992c2d9 *dev/brmsport-bin2-split.tsv
b5aefc5642b469c6ff8cc30db8c479b5 *dev/brmsport-bin2/frmtmb.sample/test-brms-suite-bin2-methods.R
477a5db996108215b45629d3ea270826 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-distributions.R
45cd99301ba12d2be76f26f280f11b28 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R
ac508cf19ee138d8a1d1568f7e7478a8 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-methods.R
b0a81a86bfc1c484028d65235788adf0 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-posterior_predict.R
8f7817700cd09f7ee5047dd1ae204a51 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-priors.R
```

Second run:

```
a6815f1d25f76422e944d385ef9093e7 *dev/brmsport-bin2-ledger.tsv
bae8caecf7ac268030ef5542b3211b96 *dev/brmsport-bin2-summary.md
ef85fbc71811fd0e24d81343929a17bd *dev/brmsport-bin2-verdicts.tsv
779e1278543fb1ef7f1590c350d28d84 *dev/brmsport-bin2-valuechecks.tsv
eef46a0cef462e793adc85db50a4ab15 *dev/brmsport-bin2-valuechecks-manual.tsv
3841f11196cdec316c167c9d46337e03 *dev/brmsport-bin2-verdicts-manual-2p.tsv
7930f070a6c9322a468dc61437e50248 *dev/brmsport-bin2-survey.tsv
55f71870a4b3335980d689d0018fd372 *dev/brmsport-bin2-groups.tsv
5be80200dee759687b5727d52992c2d9 *dev/brmsport-bin2-split.tsv
b5aefc5642b469c6ff8cc30db8c479b5 *dev/brmsport-bin2/frmtmb.sample/test-brms-suite-bin2-methods.R
477a5db996108215b45629d3ea270826 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-distributions.R
45cd99301ba12d2be76f26f280f11b28 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R
ac508cf19ee138d8a1d1568f7e7478a8 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-methods.R
b0a81a86bfc1c484028d65235788adf0 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-posterior_predict.R
8f7817700cd09f7ee5047dd1ae204a51 *dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-priors.R
```
<!-- end md5 -->

### tol (dev/bin2port-tol.R, seed 534)

<!-- begin tol -->
```
|ref| range: 2.63e-01 5.84e+06 
control: the adapter's values                all.equal alone: TRUE   new rule: TRUE   HOLDS: 400 values (0 NA, 0 infinite), max abs diff 9.31e-10, max rel 2.02e-14
one entry (|ref| 0.263) + 1e-06              all.equal alone: TRUE   new rule: FALSE  DIFFERS: 400 values (0 NA, 0 infinite), max abs diff 1e-06, max rel 3.8e-06
one entry (|ref| 0.263) + 0.0001             all.equal alone: TRUE   new rule: FALSE  DIFFERS: 400 values (0 NA, 0 infinite), max abs diff 0.0001, max rel 0.00038
one entry (|ref| 0.263) + 0.01               all.equal alone: TRUE   new rule: FALSE  DIFFERS: 400 values (0 NA, 0 infinite), max abs diff 0.01, max rel 0.038
one entry (|ref| 0.263) + 1                  all.equal alone: FALSE  new rule: FALSE  DIFFERS: 400 values (0 NA, 0 infinite), max abs diff 1, max rel 3.8
every entry + 1e-6 (the catlpdf mutation)    all.equal alone: TRUE   new rule: FALSE  DIFFERS: 400 values (0 NA, 0 infinite), max abs diff 1e-06, max rel 3.8e-06
```
<!-- end tol -->

### probe0 (dev/bin2port-probe0.R)

<!-- begin probe0 -->
```

Attaching package: 'frmtmb'

The following object is masked from 'package:stats':

    ar

 [1] "family"           "dpars"            "links"            "lpdf"            
 [5] "valid_y"          "init_dpars"       "type"             "post"            
 [9] "sim"              "sim_ctx"          "sim_refusal"      "primary_dpars"   
[13] "lcdf"             "lccdf"            "required_aterms"  "accepts_aterms"  
[17] "se_dpar"          "exclusive_aterms" "family_finalize"  "extra_pars"      
[21] "drop_intercept"   "structure"        "defer"           
NULL
[1] "logit"
[1] "mu2" "mu3"
[1] "mu2" "mu3" "phi"
List of 5
 $ kind: chr "dirichlet"
 $ cats: chr [1:3] "1" "2" "3"
 $ ref : int 1
 $ pred: int [1:2] 2 3
 $ mu  : chr [1:2] "mu2" "mu3"
[1] "mu1" "mu2" "mu3"
[1] "mu2"    "mu3"    "sigma2" "sigma3"
List of 6
 $ kind : chr "logistic_normal"
 $ cats : chr [1:3] "1" "2" "3"
 $ ref  : int 1
 $ pred : int [1:2] 2 3
 $ mu   : chr [1:2] "mu2" "mu3"
 $ sigma: chr [1:2] "sigma2" "sigma3"
[1] "y"      "dpars"  "aterms" "extra" 
[1] FALSE
[1] "mu1" "mu3" "phi"
[1] "multinomial() needs the number of categories, e.g. multinomial(K = 3)"
[1] "mu2" "mu3"
[1] "mu" "nu"
 [1] "family"           "dpars"            "links"            "lpdf"            
 [5] "valid_y"          "init_dpars"       "type"             "post"            
 [9] "sim"              "sim_ctx"          "sim_refusal"      "primary_dpars"   
[13] "lcdf"             "lccdf"            "required_aterms"  "accepts_aterms"  
[17] "se_dpar"          "exclusive_aterms" "family_finalize"  "extra_pars"      
[21] "drop_intercept"   "structure"       
[1] "mu"    "sigma" "xi"   
[1] TRUE
[1] 0.3144855
     [,1] [,2]
[1,]  1.0  0.3
[2,]  0.3  1.0
```
<!-- end probe0 -->

### probe1 (dev/bin2port-probe1.R, seed 495 and 1234)

<!-- begin probe1 -->
```
categorical log_lik                dims 50x8 / 50x8  max abs diff 4.44e-16  max rel 7.48e-16  NA 100/100
multinomial log_lik (bare constructor): multinomial() needs the number of categories, e.g. multinomial(K = 3) 
multinomial(K = 3) row_lpdf        dims 50x8 / 50x8  max abs diff 4.26e-14  max rel 6.66e-15  NA 0/0
dirichlet log_lik                  dims 50x8 / 50x8  max abs diff 5.33e-15  max rel 7.65e-16  NA 0/0
dirichlet log_lik, refcat 2        dims 50x8 / 50x8  max abs diff 3.55e-15  max rel 7.66e-16  NA 0/0
dirichlet2 log_lik                 dims 50x8 / 50x8  max abs diff 2.66e-15  max rel 7.41e-16  NA 0/0
logistic_normal log_lik            dims 50x8 / 50x8  max abs diff 2.24e-08  max rel 4.86e-12  NA 0/0
brms with and without prep$lncor identical: TRUE 
brms with prep$dpars$lncor differs: max abs change 33051061 
logistic_normal, prep$dpars$lncor  dims 50x8 / 50x8  max abs diff 1.19e-05  max rel 2.76e-12  NA 0/0
logistic_normal K = 4, 3 cors      dims 50x2 / 50x2  max abs diff 2.33e-10  max rel 8.12e-15  NA 0/0
row-major  0.30 -0.20  0.25  0.10 -0.15  0.05 
col-major  0.30 -0.20  0.10  0.25 -0.15  0.05 
logistic_normal K = 5, 6 cors      dims 50x2 / 50x2  max abs diff 7.45e-09  max rel 2.86e-15  NA 0/0
control: column-major fill (must differ) dims 50x2 / 50x2  max abs diff 7.49e+05  max rel 0.859  NA 0/0
pp categorical dim 50 8 range 1 3 
pp dirichlet dim 50 3 max |rowSums - 1| 1.11e-16 
pp dirichlet2 dim 50 3 max |rowSums - 1| 1.11e-16 
pp logistic_normal dim 50 3 max |rowSums - 1| 2.22e-16 
LN 20000 draws: frmtmb colMeans 0.3154 0.4590 0.2256  brms 0.3166 0.4575 0.2259 
LN cor of log ratios: frmtmb 0.603  brms 0.608 
```
<!-- end probe1 -->

### probe2 (dev/bin2port-probe2.R, the blocks' seeds 233 and 446)

<!-- begin probe2 -->
```
block line of log_lik:271: 233 
  withr::local_seed(233L) 
y = 12 
        brms
frmtmb   -Inf finite
  finite   17      8
         mu   shape       frmtmb      brms
9  0.140263 3.07648 -4.10458e+03      -Inf
10 0.223063 2.33133 -4.92165e+02 -492.1647
24 0.254391 7.93428 -4.99907e+08      -Inf
22 0.287050 4.01884 -2.71210e+04      -Inf
13 0.316813 2.19742 -2.70332e+02 -270.3322
11 0.329131 4.10036 -2.95710e+04      -Inf
20 0.351132 4.24050 -3.94502e+04      -Inf
6  0.359817 2.44790 -4.47969e+02 -447.9691
17 0.397232 8.56491 -1.61584e+09      -Inf
23 0.411833 6.29807 -5.55585e+06      -Inf
25 0.421937 3.54522 -5.77955e+03      -Inf
14 0.451127 2.22165 -1.98829e+02 -198.8292
2  0.461094 4.49630 -5.50994e+04      -Inf
3  0.465745 1.29849 -1.93812e+01  -19.3812
1  0.486797 5.68580 -9.84658e+05      -Inf
8  0.518993 2.82402 -7.31894e+02 -731.8941
5  0.559828 4.65239 -6.08545e+04      -Inf
15 0.565098 4.40172 -3.21149e+04      -Inf
21 0.620510 4.77986 -6.87145e+04      -Inf
7  0.637821 3.25042 -1.44783e+03      -Inf
18 0.695010 2.02208 -5.53469e+01  -55.3469
12 0.696879 3.85979 -5.28563e+03      -Inf
16 0.776636 1.99841 -3.62591e+01  -36.2591
4  0.812921 3.81311 -2.69937e+03      -Inf
19 0.827293 9.51994 -3.56096e+09      -Inf
closed form on the log scale (no underflow) vs frmtmb: max abs diff 5.72e-06 
           mu      sh    naive     logscale       frmtmb     brms
[1,] 0.140263 3.07648     -Inf -4.10458e+03 -4.10458e+03     -Inf
[2,] 0.223063 2.33133 -492.165 -4.92165e+02 -4.92165e+02 -492.165
[3,] 0.254391 7.93428     -Inf -4.99907e+08 -4.99907e+08     -Inf
[4,] 0.287050 4.01884     -Inf -2.71210e+04 -2.71210e+04     -Inf
[5,] 0.316813 2.19742 -270.332 -2.70332e+02 -2.70332e+02 -270.332
[6,] 0.329131 4.10036     -Inf -2.95710e+04 -2.95710e+04     -Inf
[7,] 0.351132 4.24050     -Inf -3.94502e+04 -3.94502e+04     -Inf
[8,] 0.359817 2.44790 -447.969 -4.47969e+02 -4.47969e+02 -447.969
cumulative, all thresholds 0, Y = 1 2 3 4 1 2 3 4 
category 1 : frmtmb NaN NaN NaN  brms -0.720118 -0.008825 -0.966682 
category 2 : frmtmb -Inf -Inf -Inf  brms NaN NaN NaN 
category 3 : frmtmb -Inf -Inf -Inf  brms NaN NaN NaN 
category 4 : frmtmb NaN NaN NaN  brms -0.67654 -0.01798 -0.78723 
      brms
frmtmb finite
  -Inf      0
hurdle_cumulative, y = 3 : frmtmb -Inf -Inf -Inf  brms NaN NaN NaN 
< table of extent 1 x 0 >
cumulative, thresholds -1, 0, 1: max abs diff 6.22e-14 
hurdle_cumulative, thresholds -1, 0, 1: max abs diff 5.95e-14 
```
<!-- end probe2 -->

### probe3 (dev/bin2port-probe3.R, seed 555)

<!-- begin probe3 -->
```
frmtmb.eam 0.12.0 
brmsfamily('wiener') dpars: mu bs ndt bias  accepts: dec vint weights cens trunc ndt_group  finalize: TRUE 
wiener(link_ndt = 'log') dpars: mu bs ndt bias  accepts: dec vint weights cens trunc ndt_group 
[1] "y"      "dpars"  "aterms"
[1]  -6.1405819  -0.8861875 -36.9778983   0.1477658  -0.7221328
[1]  -6.1405819  -0.8861875 -36.9778983   0.1477658  -0.7221328
[1]  -6.1405819  -0.8861875 -36.9778983   0.1477658  -0.7221328
```
<!-- end probe3 -->

### probe4 (dev/bin2port-probe4.R, seed 4)

<!-- begin probe4 -->
```
frmtmb 0.71.0 frmtmb.sample 0.19.0 
sar(W) lag : ERROR: sar() is a brms autocorrelation term that frmtmb does not support: sar(W, type = "lag"). The supported residual correlation terms are ar(), ma(), arma(), cosy(), unstr() and car() 
sar(W) error : ERROR: sar() is a brms autocorrelation term that frmtmb does not support: sar(W, type = "error"). The supported residual correlation terms are ar(), ma(), arma(), cosy(), unstr() and car() 
fcor(M) : ERROR: fcor() is a brms autocorrelation term that frmtmb does not support: fcor(M). The supported residual correlation terms are ar(), ma(), arma(), cosy(), unstr() and car() 
Warning message:
Optimizer did not report convergence: singular convergence (7). A Newton step from the stop gains 3.14e-09 in log-likelihood, within frmtmb_control(grad_tol = 0.001), so the estimates are at the maximum and their standard errors are judged as for a converged fit 
arma(cov = TRUE) fit converged: FALSE 
log_lik() of draws on arma(cov = TRUE): log_lik() needs a likelihood that factors into one term per observation, and the response 'y' does not: it carries arma(time = t, gr = g, p = 1, q = 1, cov = TRUE), an R-side residual correlation MATRIX, whose smallest independent unit is a whole group, so a column of the matrix would be a GROUP and leaving one out would drop part of a joint density. brms's default cov = FALSE form of arma() does factor, one term per row given that group's observed earlier residuals, and log_lik() supports it: refit without cov = TRUE for the pointwise density. Compare these models with AIC() on the ML fits, or with frm_bootstrap() 
Warning message:
Large maximum absolute gradient at the optimum (363 at (Intercept)): the curvature there is unusable, so how much log-likelihood is left could not be measured. The fit may not have converged. diagnose() reports all three numbers; see the 'Convergence problems' section of vignette('diagnostics') for the remedies 
log_lik() of draws on arma(cov = FALSE): no refusal 
mvbind rescor fit converged: TRUE 
log_lik() of draws on a rescor model: no refusal 
Error: wiener: simulating a response time needs the boundary each row ended at, which is missing. See ?wiener.
Execution halted
```
<!-- end probe4 -->

### probe4b (dev/bin2port-probe4b.R, seed 409)

<!-- begin probe4b -->
```
wiener sim_draw with dec: class numeric  dim 
[1] 0.7722316 0.4612562 0.4581513 0.3785003 0.5905511
brms posterior_predict_wiener: class matrix array  dim 5 1 
         [,1]
.x1 2.4445742
.x2 1.2034155
.x3 6.4946227
.x4 0.5870558
.x5 0.5930424
```
<!-- end probe4b -->

### probe5 (dev/bin2port-probe5.R, seed 504)

<!-- begin probe5 -->
```
log_lik multinomial dim 50 8  max abs diff vs brms 2.84e-14 
posterior_predict multinomial dim 50 3  row sums == trials[3]: TRUE 
```
<!-- end probe5 -->

### probe6 (dev/bin2port-probe6.R, seed 420)

<!-- begin probe6 -->
```
frmtmb 0.71.0 frmtmb.eam 0.12.0 frmtmb.sample 0.19.0 
frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
  Intercept          student_t(3, 0.6, 2.5)
  b                  (flat), as brms leaves slopes
posterior_predict(draws), all rows: dim 150 200 
newdata columns: rt upper 
posterior_predict(draws, newdata = 1 row) : ERROR wiener: simulating a response time needs the boundary each row ended at, which is missing. See ?wiener. 
posterior_predict(draws, newdata = 2 rows) : ERROR wiener: simulating a response time needs the boundary each row ended at, which is missing. See ?wiener. 
posterior_predict(draws, newdata = all rows) : ERROR wiener: simulating a response time needs the boundary each row ended at, which is missing. See ?wiener. 
predict(fit, newdata = 1 row, summary = FALSE) : ERROR wiener: simulating a response time needs the boundary each row ended at, which is missing. See ?wiener. 
predict(fit, summary = FALSE), all rows : class matrix/array  dim 1000 200 
simulate(fit, nsim = 3) : class data.frame  dim 200 3 
frm_sample(): default priors (brms 2.23 defaults; prior = "flat" opts out)
  Intercept          student_t(3, 0, 2.5)
  sigma              student_t(3, 0, 2.5)  [natural scale]
  b                  (flat), as brms leaves slopes
gaussian posterior_predict(draws, newdata = 1 row) : class matrix/array  dim 150 1 
sim_draw(): class numeric  dim NULL  length 5 
sim_draw(gaussian): class numeric  dim NULL 
```
<!-- end probe6 -->

### probe6b (dev/bin2port-probe6b.R, seed 420)

<!-- begin probe6b -->
```
fitted(fit, newdata = 2 rows):
      Estimate  Est.Error      Q2.5     Q97.5
[1,] 0.7483912 0.02485278 0.6996807 0.7971018
[2,] 0.7483912 0.02485278 0.6996807 0.7971018
predict(fit, newdata = 2 rows): wiener: simulating a response time needs the boundary each row ended at, which is missing. See ?wiener. 
```
<!-- end probe6b -->

### defect (dev/bin2port-defect-cumulative.R)

<!-- begin defect -->
```
frmtmb 0.71.0 from C:/Users/adf44/source/r/rellib-r9 
thresholds 0 0 0 
  y = 1  frmtmb NaN NaN NaN   exact -0.403186 -0.693147 -1.541008
  y = 2  frmtmb -Inf -Inf -Inf   exact -Inf -Inf -Inf
  y = 3  frmtmb -Inf -Inf -Inf   exact -Inf -Inf -Inf
  y = 4  frmtmb NaN NaN NaN   exact -1.103186 -0.693147 -0.241008
thresholds -1 -1  1 
  y = 1  frmtmb NaN NaN NaN   exact -0.854355 -1.313262 -2.395545
  y = 2  frmtmb -Inf -Inf -Inf   exact -Inf -Inf -Inf
  y = 3  frmtmb -0.867555 -0.771937 -1.095314   exact -0.867555 -0.771937 -1.095314
  y = 4  frmtmb -1.867786 -1.313262 -0.554355   exact -1.867786 -1.313262 -0.554355
thresholds -1  1  1 
  y = 1  frmtmb -0.854355 -1.313262 -2.395545   exact -0.854355 -1.313262 -2.395545
  y = 2  frmtmb -0.867555 -0.771937 -1.095314   exact -0.867555 -0.771937 -1.095314
  y = 3  frmtmb -Inf -Inf -Inf   exact -Inf -Inf -Inf
  y = 4  frmtmb NaN NaN NaN   exact -1.867786 -1.313262 -0.554355
thresholds -1  0  1 
  y = 1  frmtmb -0.854355 -1.313262 -2.395545   exact -0.854355 -1.313262 -2.395545
  y = 2  frmtmb -1.41622 -1.46508 -2.09523   exact -1.41622 -1.46508 -2.09523
  y = 3  frmtmb -1.72965 -1.46508 -1.55404   exact -1.72965 -1.46508 -1.55404
  y = 4  frmtmb -1.867786 -1.313262 -0.554355   exact -1.867786 -1.313262 -0.554355
taped sum over rows y = 1, 1, 4 at tau = (0, 0, 1): -1.650688 
exact: -1.650688 
on doubles: NaN 
hurdle_cumulative, hu 0.2, y 0 1 4 , thresholds 0 0 1 
  frmtmb       NaN       NaN -1.536405 
  exact  -1.6094379 -0.9162907 -1.5364052 
hurdle_cumulative, hu 0.2, y 0 1 4 , thresholds -1  0  1 
  frmtmb -1.609438 -1.536405 -1.536405 
  exact  -1.609438 -1.536405 -1.536405 
```
<!-- end defect -->
