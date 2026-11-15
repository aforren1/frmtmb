# Lane priors2, round of 2026-11-13 (base ec8178d3, frmtmb 0.82.0)

Branch `wt-r83-priors2`; brief `round-briefs/r083-priors2.md`. Nothing
is committed. This part is the record of this round's lane. The record
of the earlier lane of the same name (2026-09, the dpar-class flip)
follows under "Earlier lane".

STATUS: fix round 3 done (section "Fix round 3"); nothing pending but
the user's answer on the mixture-weights Jacobian.

## Fix round 3 (re-check after fix round 2, short review)

This section is the record where it and the earlier ones disagree.

1. Silent local mode. `open_end_probe()` now also checks each
   infinite-density end the fit ended away from: the objective at the
   nearest double past the end, the other parameters at the fit
   (`obj$fn`, not the optimizer's counted wrapper). Where it is `-Inf`
   or below the fit's by max(1, sqrt(eps) |objective|), the fit is kept
   and `check_map_unbounded()` warns, class
   `frmtmb_prior_unbounded_local`, naming the parameter and the end. No
   refusal. Where no end is near, the objective's state is restored as
   saved (not settled), so a fit whose end is not lower keeps its state
   and evaluation count as before. `open_end_rescue()`'s early return
   goes through the probe too. Before/after (`fix2-probe-after.txt`,
   `fix3-probe-after.txt`):
   - J, `y ~ x + z`, `gamma(0.5, 1)` lb 0 on x: code 0 at 0.3929088467,
     163.395461346473, silent; now the same fit with the warning.
   - D `y` (60.0663831131304) and C (29.5619832612003): the same fits,
     now with the warning (the same posterior).
   - A to I: every other line identical (`diff --strip-trailing-cr`),
     including H with no new warning (its posterior vanishes at 0).
   - J with `lognormal(0, 1)` and `gamma(2, 2)`: the fit as before, no
     warning.
   Tests (`test-prior-bounds-ends.R`): D `y` and the short review's P3
   warn with the class; case H at 1e7 does not; `gamma(2, 2)` on P3 is
   silent. NEWS: folded into the infinite-density bullet. The
   hard-bounds section of `?set_prior` (`R/priors.R`) still said "where
   it runs to the end ... frm() refuses", stale since fix round 2; it
   now says what the probe does, so `man/set_prior.Rd` changed and core
   was checked again.
2. Filed LOW in `dev/test-backlog.md`: the cauchit `linkfun` (tan())
   maps bounds outside (0, 1) to finite values, a box at 0.5, with the
   short review's numbers.

Runs: `fix3.log`, lane library: `test-prior-bounds-ends.R` pass 85,
`test-prior-mixture-theta.R` 55, `test-prior-table-brms.R` 101,
`test-conditions-census.R` 11; fail, error, skip, warn 0. 33-fit survey
(`survey-hex-fix3.txt`): 29 of 33 bitwise identical, the same 4 others.
R CMD check --as-cran, core (`check-frmtmb-fix3.txt`): Status: 1 NOTE
(V8); tests 28 min OK, examples OK. The check tree was deleted.

## Fix round 2 (short review `dev/reviews/2026-11-13-priors2-short.md`, FIX)

This section is the record where it and "Fix round" disagree. F2, F4,
F5, L1 to L3 and the famend merge are untouched. Outputs:
`fix2-probe-after.txt` (`dev/priors2-fix-probe.R`, now with case H, the
end label and `zi lb`), `survey-fix2.tsv`, `survey-hex-fix2.txt`,
`review-bitwise-fix2.tsv`, `short-bitwise-fix2.tsv`, `fix2.log`,
`ob0.3.32-fix2/`, `check-frmtmb-fix2.txt`.

### 1. The refusal is a proof (`open_end_probe()`)

`check_map_unbounded()` no longer reads position. Where the fit ends
within the slack of an infinite-density end, `open_end_probe()`
evaluates the objective along that coordinate at the end plus (slack
times 100^-k), down to the nearest double past the end (2^-1074 at 0).
It refuses only where the objective falls at every point (a `-Inf` only
as the tail) and its last value is below the best objective found (the
box run, the refit, the first probe point) by max(1, sqrt(eps) |best|).
Otherwise the best point is kept; where that is a probe point, it is
returned at code 1 with a warning ("not shown to rise without bound").
The unit test drives all three branches with stand-in objectives.

| case (`fix-probe-after.txt` / `fix2-probe-after.txt`) | before (fix round) | after |
| --- | --- | --- |
| A gamma(0.5,1), poisson identity | code 0, 120.129961301298, b = 1.383301921, NA/NaN warning | the same fit, no warning |
| B beta(0.5,2) | code 0, 88.9031549557894, 0.6478533563, NA/NaN warning | the same, no warning |
| C gaussian (review's) | code 0, 29.5619832612003 | the same |
| D `y` | code 0, 60.0663831131304 | the same |
| D `yneg` | refused | refused, by the probe |
| D' `yneg`, beta(0.5,2), lb 0, ub 1 | refused, named "ub = 1" | refused, named "lb = 0" |
| E lognormal, beta(2,2), gamma(2,2) (pushed) | 133.228272690687, 129.86165226729, 130.273893123621, NA/NaN warning | the same, no warning |
| F famend start | code 0, 129.86165226729 | the same |
| G, G' sd | code 0, 36.5172433865701, 35.868587845369 | the same |
| H gamma, 1e5, S 1 / 3 | code 0 at b* | the same, no NA/NaN warning |
| H gamma, 1e7, S 1 | refused (b* = 4.854e-10) | code 0, b = 4.853823773e-10 (b* 4.853823771e-10) |
| H gamma, 1e7, S 3 | refused | code 0, 2.426911886e-09 (b* the same) |
| H gamma, 1e9, S 1 / 3 | refused | code 0, 4.853823776e-12 / 2.426911934e-11 (b* 4.853823776e-12 / 2.426911888e-11) |
| H beta(0.5,2), 1e5 | code 0 | the same |
| H beta(0.5,2), 1e7, S 1 / 3 | refused, named "ub = 1" | code 0, 4.853823772e-10 / 2.426911886e-09 |
| H beta(0.5,2), 1e9, S 1 / 3 | refused | code 0, 4.853823728e-12 / 2.426911911e-11 |

H at 1e7 and 1e9 still warns "Standard errors are not available ... a
bound holds it": the SE check's positional slack, as with `lb = b*/1e4`
in the short review's own table. Filed LOW in the backlog; not a
refusal, and 0.82.0 failed on these fits. Case H at 1e7 is in
`test-prior-bounds-ends.R`.

The end label: `box_end_density()` keeps one label per end
(`end_labels$lower`, `$upper`), and `open_end_flags()` one per position
and side, so the refusal names the end the probe ran to.

### 2. The box run's "NA/NaN function evaluation"

`open_end_hold()` holds that warning from the box run (its runs,
restarts and `natural_finish()`) where a box end is open, and only
there. `open_end_rescue()` marks a kept refit (`opt$open_refit`); the
held warnings are dropped then, and raised again where the box run's
own result is returned or the rescue raises. A, B, E and H give no
NA/NaN warning now (`fix2-probe-after.txt`). The `allow_warnings()` for
it in `test-prior-bounds-ends.R` (two sites) are gone, so the test fails
if it returns.

### 3. A NaN lower bound above the support

`internal_bound()` reads a NaN lower bound as "constrains nothing" only
where it is at or below the infimum of the inverse link's range
(`min(linkinv(-xmax), linkinv(xmax))`); above it, the bound is refused
as an empty box. Through `frm()`: `set_prior("", class = "zi", lb = 2)`
was code 0 at the unbounded fit; now "lb = 2 is outside the support of
the parameter this prior is about, so it describes an empty box".
`lb = -1` still gives the unbounded fit, bitwise. Both are tested
through `frm()`. NEWS has a BREAKING bullet for it.

### Merge note correction

For `test-prior-table-brms.R`, do a three-way merge (base ec8178d3,
priors2, famend), not famend's whole file: priors2 also changed the
`cor`/`resp` test (base line 180), which famend's copy keeps as
accepted. The at0 block takes the resolution under "Merge note for
consolidation" below; the short review's scratch merge did this and
passed 105 of 105.

### Runs

- One file per process, lane library (`fix2.log`):
  `test-prior-bounds-ends.R` pass 80, `test-prior-mixture-theta.R` 55,
  `test-prior-table-brms.R` 101, `test-conditions-census.R` 11; fail,
  error, skip and warn 0.
- 33-fit survey, hex against rellib-r21 (`survey-hex-fix2.txt`): 29 of
  33 bitwise identical, the same 4 others. The review's 15 fits: `cmp`
  identical to `priors2-review-bitwise-base1.tsv`. The short review's 5
  unbounded fits: `cmp` identical to `priors2-short-bitwise-base.tsv`.
- OpenBLAS 0.3.32 with LAPACK (`ob0.3.32-fix2/`):
  `test-prior-bounds-ends.R` 80, `test-prior-mixture-theta.R` 55, fail,
  error and warn 0. `dev/ciharden-out/` was then deleted.
- R CMD check --as-cran, core (`check-frmtmb-fix2.txt`): Status: 1 NOTE
  (V8 math rendering); tests 37 min OK, examples 49 s and 50 s
  (donttest) OK. The build predates the last NEWS rewording only. The
  check tree was deleted.

## Fix round (review `dev/reviews/2026-11-13-priors2.md`, FIX)

Where this section and the item sections below disagree, this section
is the record. New scripts: `dev/priors2-fix-probe.R` (output
`fix-probe-{before,after}.txt`), `dev/priors2-survey-hex.R`
(`survey-hex.txt`); the review's `priors2-review-bitwise.R` rerun on
the lane (`review-bitwise-lane2.tsv`).

### F2: the closed box first (decided)

`optimize_obj()` is again 0.82.0's: the box run, its restarts and
`natural_finish()`. Then `open_end_rescue()` (`R/fit-open-ends.R`)
refits in the open chart only where that run stopped on a box end
where the prior's density is 0 or infinite (within the
`grad_bound_active()` slack), or ended at a non-finite objective beside
one, or failed beside one (`open_end_failed()` holds the error and
raises it if the refit fails too). A box run at an interior maximum is
returned untouched.

- The 33-fit survey (`survey-cmp.txt` of the fix round,
  `survey-hex.txt`): 29 of 33 bitwise identical to rellib-r21 in hex,
  including the four interior fits that moved 4.7e-7 to 7.6e-6 in the
  first version. The other 4: `zi_unif_on` errored on 0.82.0 (5027),
  and the three pushed fits, code 1 at objective Inf on 0.82.0, now
  converge (code 0) at 207.2657, 206.0200, 205.5936.
- The review's 15 fits (`priors2-review-bitwise.R lane`): `cmp` of the
  lane's tsv against `priors2-review-bitwise-base1.tsv` is identical.
- Pushed fits still reach their mode (`fix-probe-after.txt`, E):
  `lognormal(0, 1)` 133.228272690687 at z = 0.02012771547 (lb = 1e-8:
  133.228272690687, 0.02012771984); `beta(2, 2)` 129.86165226729 at
  0.006822071213 (lb = 1e-8: 129.861652267295, 0.006822082412);
  `gamma(2, 2)` 130.273893123621 at 0.006775969295 (130.273893123626,
  0.006775981223). 20 of 20 ulp perturbations converge on each
  (`ulps2.txt`).

### F1: an infinite-density end is open, refused only after the fit

`box_end_density()` lists infinite-density ends as open, and also as
`inf_*`; there is no pre-fit refusal. `check_map_unbounded()` now
refuses after the fit, and only where the open refit ran to an
infinite-density end (`opt$open_unbounded`).

| case (`fix-probe-*.txt`) | 0.82.0 | first lane build | fix round |
| --- | --- | --- | --- |
| A: poisson("identity"), `gamma(0.5, 1)`, `lb = 0` | error, NA/NaN gradient | refused | code 0, 120.129961301298, b = 1.383301921 (lb = 1e-10: 1.383301902) |
| B: the same, `beta(0.5, 2)`, `lb = 0`, `ub = 1` | error, NA/NaN gradient | refused | code 0, 88.9031549557894, b = 0.6478533563 (lb = 1e-10: 0.6478533564) |
| C: the review's gaussian slope | code 0, local mode 0.3022615629 | refused | the same as 0.82.0 to the digits printed (29.5619832612003): the box run ends inside |
| D: gaussian `y ~ z`, slope inside | code 1, -Inf | refused | code 0, 60.0663831131304, z = 0.2498450905 |
| D: gaussian `yneg ~ z`, slope below 0 | code 1, -Inf | refused | refused: the refit ran to the end |
| G: `class = "sd"`, `gamma(0.5, 1)`, `lb = 0` | code 0, 36.5172433865701 | code 0 (the same) | code 0 (the same) |

G was never refused, on any build, and must not be: a bound at 0 on
class `sd` is no box (its log is held, `log(0) = -Inf`), and the
density on the log, `gamma(0.5, 1)` on the sd times the Jacobian sd, is
`sd^0.5 exp(-sd)`, which vanishes at sd = 0. So the posterior has its
maximum inside, as G' (no group effect in the data) shows at log sd
-3.459.

The gaussian-slope description is fixed in NEWS: 0.82.0 returned the
local mode where the box run ended inside (C) and stopped at -Inf where
it ran to the end (D).

### F3: `bound_inside()`

The step doubles until it lands inside and a bisection then takes the
inside value nearest the original; a bound no value maps inside is
refused by name. The review's 301 logit ends `1 - 10^-seq(2, 5, 0.01)`
and `ub = 0.99981802991413904` are now in `test-prior-bounds-ends.R`
(both sides inside on all), with the refusal (`lb = 2` through
`plogis`).

### F4: NEWS and the floor

"Through 0.82.0" and "Through 0.29.0" are gone; the infinite-density
refusal is BREAKING and described as what remains after F1. The
frmtmb.sample floor must rise at release: its DESCRIPTION still says
`frmtmb (>= 0.82.0)`, and with core 0.82.0 its class `theta` tests fail
(core refuses the row).

### F5, L1 to L3

- F5: `R/compat.R:2167` names class `"thetacov"`.
- L1: `class = "theta"` with a coef `^theta(ac|r)?_[0-9]+$` is refused
  with "Write class = "thetacov", coef = ...", with a density, a bound
  or neither; an empty density on a narrowed weights row prints as
  `theta_a ~ (flat)`. Tested in `test-prior-mixture-theta.R`.
- L2: the refusal says that `frm_sample()` takes the bound on the
  class row, not on one coefficient.
- L3: filed (LOW) with the three sites, and said in `open_end_flags()`.

### Merge note for consolidation (famend)

`tests/testthat/test-prior-table-brms.R`, base lines 327 to 346: famend
keeps "the caller's start at the box's end: a stop at an infinite
objective is reported as not converged", `at0 <- allow_warnings(frm(yneg
~ z, data = d, prior = sp(lb = 0, ub = 1), start = list(beta = c(1,
0))), c("the objective is Inf", "gradient at the reported"), require =
"the objective is Inf")` with `expect_identical(at0$opt$convergence,
1L)` and `expect_false(is.finite(at0$opt$objective))`. After the merge
priors2's rescue converges there (`fix-probe-after.txt`, F:
code 0, 129.86165226729, z = 0.006822067572). Resolution (corrected
in fix round 2): merge the file three ways, not by taking famend's
copy, and resolve the one conflict, the at0 block, with

    at0 <- frm(yneg ~ z, data = d, prior = sp(lb = 0, ub = 1),
               start = list(beta = c(1, 0)))
    expect_identical(at0$opt$convergence, 0L)
    expect_lte(abs(at0$opt$objective - m8$opt$objective),
               1e-8 * abs(m8$opt$objective))

and drop priors2's own version of that test's last block (identical in
effect). The guard `nonfinite_not_converged()` keeps its unit test in
`test-prior-bounds-ends.R`. No double move: famend's
`prior_start_inside()` moves a cold start whose prior density is not
finite before the optimizer runs; the box run then ends inside, and
`open_end_rescue()` does nothing. A caller's start, which famend
leaves, is the only case the rescue moves, and it moves it after the
box run, once.

### Fix-round runs

- One file per process, the whole core suite and every changed
  extension file (`priors2-fix-files/`, `priors2-fix.log`, 338 of 338,
  every `lib:` line the lane library): frmtmb 332 files, pass 24565,
  fail 1, error 0, skip 240, warn 0; frmtmb.sample 5 files (the
  changed ones), pass 433, fail 0; frmtmb.spatial `test-sar.R`, pass
  157, fail 0, skip 4. The one failure, `test-conditions-census.R`:
  bare `stop()` in `R/fit-open-ends.R`; now `frm_stop()`, reinstalled
  and rerun: pass 11, fail 0. `test-prior-bounds-ends.R` (pass 66) and
  `test-prior-mixture-theta.R` (pass 55) were rerun on that build.
- OpenBLAS with LAPACK, 0.3.26 (rebuilt here) and 0.3.32 (the
  reviewer's build), `dev/priors2-openblas-run.sh`, `ob0.3.26/`,
  `ob0.3.32/`: on each, `test-prior-bounds-ends.R` 66,
  `test-prior-mixture-theta.R` 55, `test-predfix.R` 92,
  `test-priors-bounds-grcov.R` 54, `test-mixture-starts.R` 135,
  `test-prior-table-brms.R` 101, `test-priors-autocor-classes.R` 186;
  fail, error and warn 0. (The build before the `frm_stop()` change;
  that change touches only a failing optimizer's error class.)
- New assertions over 20 ulp perturbations: 20 of 20 on each
  (`ulps2.txt`).
- One session, forward, on the final build (`fix-sessions/`): frmtmb
  files 319/332, pass 24566, fail 0, error 0, skip 240, warn 0;
  frmtmb.sample files 86/91, pass 4248, fail 0, error 0, skip 18, warn
  0. No untagged LEAK line.
- R CMD check --as-cran on the final build (`dev/priors2-check.sh`,
  `check-frmtmb-fix.txt`, `check-sample-fix.txt`): frmtmb Status: 1
  NOTE (V8 math rendering), tests 52 min OK, examples 94 s and 106 s
  (donttest) OK; frmtmb.sample Status: 1 NOTE (V8), tests 13 min OK.
  The check trees were deleted after these lines were copied.

### Mixture-weights Jacobian

Unchanged, pending the user's answer.

## How every number here is made

| what | script | output (gitignored, `dev/priors2-out/`) |
| --- | --- | --- |
| brms 2.23.0 on every case, no compile | `dev/priors2-brms.R` | `brms.txt` |
| frmtmb's side of the same cases | `dev/priors2-frm.R` (`PRIORS2_LIB=none` = rellib-r21) | `frm-before.txt`, `frm-after.txt` |
| bounds: 5005, 5019, 5027, 5646 | `dev/priors2-bounds.R` | `bounds-before.txt`, `bounds-after.txt` |
| 33 bounded fits on both builds | `dev/priors2-survey.R`, `dev/priors2-survey-cmp.R` | `survey-*.tsv`, `survey-cmp.txt` |
| the four moved fits against a tight optimum | `dev/priors2-survey-tight.R` | `survey-tight.txt` |
| link round trips | `dev/priors2-ulp.R` | `ulp.txt` |
| new assertions over 20 ulp perturbations | `dev/priors2-ulps.R` | `ulps.txt` |
| mixture weights' SE in summary() | `dev/priors2-theta-se.R` | (stdout) |
| install; one test file | `dev/priors2-install.sh`; `dev/priors2-test1.sh` | `inst-*.log`; `t1-<arm>-<file>.txt` |
| suites, one file per process | `dev/priors2-run-par.sh` | `priors2-all.log`, `priors2-all-files/` |
| suites, one session per package | `dev/priors2-run-sessions.sh` | `priors2-sess.log` |

Data and seeds are in each script's header lines.

## Item 1: the `theta` rename and the mixture-weights row (1132, 3279)

### brms 2.23.0 (`brms.txt`)

- A mixture with no weight formula and no fixed weight has one row
  `dirichlet(1) theta` (M1, M2, M5 ordinal, M7 poisson, M40 with
  `sigma1 ~ x`), one per response in a multivariate model (M6,
  `resp = y`, `y2`). With `theta1 ~ x` (M3) or `theta1 = 0.4` (M4)
  there is no `theta` row.
- Accepted: `set_prior("dirichlet(2)", class = "theta")` (M8,
  `con_theta` 2 2), a vector per component (M9, 1 2 3), brms's
  `prior(dirichlet(2), class = theta)` (M24), `""` (M29, keeps
  `dirichlet(con_theta)` at 1), concentrations below 1 (M28), `resp`
  in a multivariate model (M21), an ordinal mixture (M27).
- Refused: any other density, "Currently 'dirichlet' is the only valid
  prior for simplex parameters" (M11); `theta1` (M12, M13); `dpar`,
  `coef`, `group`, `nlpar` (M14 to M17, M30: "theta_theta1",
  "theta_a", "theta_g"); a weight with a formula (M19) or a fixed value
  (M25); a non-mixture (M20); no `resp` in a multivariate model (M22);
  a concentration vector of the wrong length (M41, "Invalid Dirichlet
  prior. Expected input of length 3."). `lb = 0` is accepted and
  silently dropped (M18); frmtmb refuses it instead, as it refuses a
  bound on every simplex class (a silent drop is the worse answer).

### What changed

- Core class `"theta"` (the internal covariance vector) is now
  `"thetacov"` everywhere a class is named: `frmtmb_prior_classes`, the
  resolver branch, `default_prior()` rows and their ordering
  (`prior_rows_sd_order()`, `prior_table()`), every refusal message,
  `?set_prior`, `?default_prior`, `?frm`, `?par_template`,
  `?frm_simulate`, the inputs and brms-migration vignettes, and
  frmtmb.sample's `refuse_coef_bounds()` and docs. The internal
  component and parameter names (`theta`, `theta_12`, `thetaac_1`) do
  not change. Historical `dev/` scripts of earlier lanes keep their
  spelling: they record runs on earlier builds. No live `dev/` reader
  (release drivers, `equiv-*`, `brmsport-*`) names the class.
- Class `"theta"` is the weights' simplex: `mix_weight_simplexes()`
  (which responses have one, their `K`, reference and `betad`
  positions; brms's rule, plus frmtmb.latent's `lca()` excluded because
  its weights are location parameters), `mix_weight_entry()` (the
  resolved Dirichlet), `set_prior()`'s simplex rules, `resp` required in
  a multivariate model (`resp_keyed_prior_classes`),
  `brms_direct_prior_classes` (a brms table row carries over).
- frm(): the entry is `dirichlet_dist(alpha, K, "softmax")` with the
  reference position; `dirichlet_logdens()` puts the 0 log ratio where
  the family's reference is. That is the density of the log ratios,
  the Dirichlet times the softmax Jacobian, the same density
  frm_sample() samples, so the MAP fit is that posterior's mode. A
  decision of this lane: the brief says "a penalty on the softmax
  chart"; with the chart's Jacobian the MAP is the sampled posterior's
  mode (internal consistency, and frmtmb's MAP convention for `simo`,
  `sd`, dpar classes), and every concentration above 0 has a maximum,
  so brms's `dirichlet(0.5)` (M28) is not refused. Without it,
  concentrations below 1 would have had to be refused on frm().
- default_prior(): one `theta` row per simplex, `(flat)` on
  `route = "fit"` (frm() applies no density without a row, as for
  `simo`), brms's `dirichlet(1)` on `route = "sample"`.
- frmtmb.sample: `mix_theta_nlp(fit, entries)` steps aside for a
  simplex a user entry covers (`mix_theta_own()`), the announcement
  leaves those weights out, and `prior_plan()` takes the row's density
  for the prior draws.
- `frm_simulate()` refuses a `dirichlet()` entry with a reason (a
  simplex of two weights reached `rprior_one()` and stopped with "No
  sampler for prior kind 'dirichlet'"; filed).

### Before and after (`frm-before.txt`, `frm-after.txt`)

- M1: 0.82.0 lists no `theta` row; the lane lists `(flat) theta`
  after `sigma2`, brms's place; M1s (sample route) lists
  `dirichlet(1) theta`. M6 lists one per response.
- M8b to M9, M21, M28, M31: 0.82.0 refuses every `dirichlet()` on
  class `theta` ("dirichlet() is the density of a simplex; ... not
  class = "theta""); the lane fits. Weights, two gaussian components,
  160 rows, seed 3: ML 0.305628 / 0.694372; `dirichlet(2)` 0.311037 /
  0.688963 (objective 311.9146428 to 313.2132785); `dirichlet(c(20,
  1))` 0.661547 / 0.338453 at another mode of the mixture (objective
  329.7702481); three components under `dirichlet(c(1, 2, 3))`
  0.310659, 0.199513, 0.489828.
- M30: `frm(prior = default_prior(fit))` is bitwise the ML fit on both
  builds.
- R1: the `toep()` rows are class `thetacov` (were `theta`).

### Tests

- New `tests/testthat/test-prior-mixture-theta.R`: rows, refusals in
  brms's words, the MAP against an independent optimizer of the same
  posterior (two components at three concentrations, three components
  at `c(9, 1, 3)`: within 1e-3 of the move the prior made, objectives
  equal to 1e-8 relative), the round trip, brms's spelling,
  `frm_simulate()`, `thetacov`. rellib-r21: pass 5, fail 15, error 5
  (first version of the file); lane: pass 52, fail 0.
- frmtmb.sample `test-default-priors-dpar.R`: brms's table and the
  sampling table now agree on the mixture (was "theta||||: brms
  dirichlet(1), frmtmb NA"); a `theta` row replaces `dirichlet(1)` on
  the tape (density checked to 1e-12) and in `prior_plan()`.
  rellib-r21: pass 56, fail 1, error 1; lane: pass 64.
- frmtmb.sample `test-prior-draws.R`: `sample_prior = "only"` under
  `dirichlet(c(30, 1))` draws theta1 with mean 30/31 and sd 0.0306,
  each within 4 Monte Carlo errors. rellib-r21: error ("dirichlet() is
  the density of a simplex"); lane: pass 95.
- Renamed in 12 test files (`class = "thetacov"`); the mixture
  expectations of `test-get-prior-route.R`, `test-simplex-pars.R` and
  `test-prior-compat.R` now expect brms's row and spellings.

## Item 2: the step-start finding is a message

`mixture_starts_signal()` calls `frm_message(class =
"frmtmb_mixture_step_start")`. `?mixture`, `?frmtmb_control`
(`starts`) and the source docs say so. `test-mixture-starts.R` captures
the message class, asserts no warning carries it, and gains a direct
test (on rellib-r21 the signal is a warning, so the new test fails
there). Four files that allowed it as a warning no longer do
(`test-ordinal-mixture-brms.R`, `test-ordinal-mixture.R`,
`test-ordinal-thres-names.R`, which now muffles the message class, and
`test-predict-draw-edge.R`, a file of lane predraw2's area: one hunk).
The autoscale pre-fit runs no starts (`ctl$starts <- 0L`), so no
message leaks from it, and `autoscale_choose()` replays messages too.

## Item 3: `cor` with `resp` is refused

brms 2.23.0 refuses it with and without `group`, and on a univariate
model (`brms.txt`, C2 to C5: "L_g__y2", "L_y2", "L_y"). The resolver's
`cor` branch refuses `resp` beside `dpar` and `nlpar`, opening with
brms's sentence and the parameter brms names. Both routes reach the
resolver. `cov_class_refusal()` no longer suggests a `cor` call with
`resp`. 0.82.0 fitted C2 and C3 (MAP, log prior 0.205304); the lane
refuses both. `test-prior-table-brms.R`: rellib-r21 accepted the row
(the old test asserted it); now refused through `validate_prior()` and
`frm()`, with the group spelling reaching both responses' blocks.
Decision recorded in the backlog CLOSED line (2026-11-13, reversing
2026-09-24).

## Item 4: prior bounds on the fit route (5005, 5019, 5027)

### Measured first (`survey-cmp.txt`, `survey-tight.txt`)

33 bounded fits (bounds on `b`, `sd`, `sigma`, `zi`, `phi`, `shape`,
`nu`, an uncentered intercept, `nlpar`, a dpar's `b`, `thetacov`; with
the maximum inside, on the bound, and pushed past a zero-density end):

- 25 bitwise unchanged.
- 4 moved, all interior maxima behind a zero-density end that the open
  chart now takes (`lognormal`, `gamma(2, 2)`, `beta(2, 2)`,
  `inv_gamma(3, 1)` with `lb = 0`): by 4.67e-7, 1.34e-6, 2.1e-6 and
  7.59e-6 relative, the objective by -1.3e-10, -1.4e-10, +4.8e-10 and
  -6.2e-9. Against a tight optimum (nlminb, rel.tol 1e-15) the box of
  0.82.0 was 9.58e-7, 9.16e-7, 0 and 6.88e-6 off and the lane is
  6.26e-7, 4.28e-7, 2.1e-6 and 7.09e-7 off: the lane is nearer in three
  of four. Three of the four exceed the brief's 1e-6; each is within
  the optimizer's own stopping tolerance, and none loses objective
  beyond 4.8e-10.
- 1 base error now fits: `uniform(0, 0.1)` on `zi` (5027, "NA/NaN
  gradient evaluation" on 0.82.0).
- 3 pushed past a zero-density end: 0.82.0 stopped on the bound, code
  1, objective Inf; the lane converges (code 0) at the mode.

### What changed

- 5027: `internal_bound()` moves a bound carried through a link inward
  an ulp at a time until `linkinv()` lands inside (`bound_inside()`).
  frmtmb's logit misses 37 of 99 ends (22 above, 15 below), the log 11
  of 100 (`ulp.txt`). `uniform(0, 0.3)` on `zi`: code 1 at 186.1733365
  with "false convergence (8)" and a 0.22 Newton gain on 0.82.0; code 0
  at the bound at 185.9492756 on the lane (`bounds-*.txt`, C2, C3).
  `uniform(0.35, 5)` on `sigma` with the data's sigma under 0.35:
  "NA/NaN gradient evaluation" on 0.82.0, code 0 at 0.35 on the lane;
  0.36, which round-trips, is the same on both (C6).
- 5005: the transform the brief prefers, where it is needed. A user's
  bound at an end of the density's support where the density is zero
  (`box_end_density()`, `prior_end_density()` read off the density's
  form, since `lognormal()` at 0 is NaN on doubles) is marked open;
  `optimize_obj()` then steps in `u` with `lb + exp(u)`, `ub - exp(u)`
  or a logit between two (`open_end_coords()`, `R/fit-open-ends.R`),
  no Jacobian, and finishes with one run in the box's own coordinates
  from there (`open_end_finish()`). Every other bound keeps its closed
  box, so the maximum on a bound, the SE check's "a bound holds it",
  `grad_bound_active()` and the restarts are untouched. Not a transform
  of every bounded coefficient: that would move every boundary maximum
  off its bound by the optimizer's tolerance and reach every reader of
  the box. On `yneg ~ z` (seed 1, 120 rows): `lognormal(0, 1)`,
  `beta(2, 2)`, `gamma(2, 2)` with `lb = 0` stopped at z = 0 with the
  objective Inf, code 1, on 0.82.0; the lane converges at z =
  0.0068221 (beta, objective 129.8616523, the mode the record gives at
  `lb = 1e-8`).
- An end where the density is infinite (`gamma(0.5, 1)`, `beta(a < 1)`)
  has no maximum: frm() refuses it (`check_map_unbounded()`). 0.82.0
  stopped there with the objective -Inf, also on data whose maximum
  likelihood slope is inside (case A gamma(0.5, 1) on `y`). The
  refusal is a proof: the likelihood is finite at the end, and the
  prior's log density rises without bound. A bound inside the support
  (`lb = 1e-8`) and a density finite at its end (`exponential(1)`,
  `beta(1, 2)`) keep the closed box (tested). frm_sample() is not
  refused.
- 5019 (a uniform on a centered intercept): not fixed; see item 5.

### Tests

New `tests/testthat/test-prior-bounds-ends.R`. rellib-r21: pass 19,
fail 26; lane: pass 54, fail 0. `dev/priors2-ulps.R`: every new fit
assertion holds on 20 of 20 ulp perturbations (`ulps.txt`).

## Item 5: a bound on class `Intercept` (5646)

brms 2.23.0 bounds the temporary intercept for centered predictors:
`real<lower=3> Intercept;` (I1, I2), and `Intercept_sigma` under
`sigma ~ x` (I4); a bound with `coef`, the only spelling for the
intercept at zero, is refused (I3, I7). frmtmb's `Intercept` density
reads the centered intercept too. 0.82.0 boxed the intercept at zero:
`x = z + 10`, `lb = 3`, intercept 3.0000 and slope 0.2078 where lm()
gives 1.0155 and 0.4055 and brms's bound binds nothing (I0, I1).

Decision: the centered bound is a linear constraint on frm()'s
parameters, which its box cannot hold without teaching every reader of
the box (filed, MEDIUM). So `lb`/`ub` on class `Intercept` is refused
where the centering offset is material (`center_offset_material()`: a
column mean beyond `sqrt(eps)` of its sd, so `scale()`d predictors and
no predictors keep the bound, which is then brms's own), on both
routes, with the two working spellings named: center the predictors,
or `bf(center = FALSE)` with `class = "b", coef = "Intercept"`. Tested
(`test-prior-bounds-ends.R`, its guard-absent cases), and
`test-priors-bounds-grcov.R` and `test-predfix.R` moved to the
uncentered spellings.

## Decisions taken that the user may want to revisit

1. frm()'s class `theta` density carries the softmax Jacobian (above).
2. `lb`/`ub` on class `theta` is refused, where brms drops it silently.
3. The open chart only where a closed box was wrong (zero-density
   ends), not for every bounded coefficient.
4. A centered-intercept bound is refused rather than held.

## Not done, and filed (backlog, "From lane priors2")

- MEDIUM: a box on the centered intercept (5646's full fix, and 5019).
- LOW: `summary()` of a mixture prints its weights with NA standard
  errors (pre-existing on rellib-r21; `vcov()` has 0.184 on the log
  ratio).
- LOW: frm_simulate() draws from no `dirichlet()` prior.

## Bump and floor

Core: minor (BREAKING class rename, refusals, message). frmtmb.sample:
minor, floor the new core (it needs class `theta` resolved as weights
and the reference position in `dirichlet_logdens()`). Other extensions:
no change of code; frmtmb.spatial's tests name the renamed class.


## Files touched

- Core R: `R/priors.R` (lane-owned), `R/fit.R` (three small hunks in
  `fit_assembled()` and `optimize_obj()`, two `grad` signatures in
  `opt_coords()`, docs), new `R/fit-open-ends.R`, `R/simulate-new.R`,
  `R/mixture-starts.R`, `R/families.R` (`?mixture` doc lines),
  `R/par-template.R` (doc line); man pages regenerated
  (`default_prior`, `frm`, `frm_simulate`, `frmtmb_control`, `mixture`,
  `par_template`, `set_prior`); `NEWS.md`; vignettes `inputs.Rmd`,
  `brms-migration.Rmd`.
- frmtmb.sample: `R/sample.R`, `R/prior-draws.R`, `man/frm_sample.Rd`,
  `NEWS.md`.
- Tests, core: new `test-prior-mixture-theta.R`,
  `test-prior-bounds-ends.R`; edited `test-autocor-seam.R`,
  `test-brms-formula-priors.R`, `test-get-prior-route.R`,
  `test-latres.R`, `test-lkj.R`, `test-mixture-starts.R`,
  `test-nonfinite-gradient.R`, `test-ordinal-mixture-brms.R`,
  `test-ordinal-mixture.R`, `test-ordinal-thres-names.R`,
  `test-predfix.R`, `test-predict-draw-edge.R` (predraw2's area, one
  hunk), `test-prior-compat.R`, `test-prior-location-dpar.R`,
  `test-prior-mv-resp.R`, `test-prior-one-slot.R`,
  `test-prior-table-brms.R`, `test-priors-autocor-classes.R`,
  `test-priors-bounds-grcov.R`, `test-simplex-pars.R`, `test-v19.R`.
- Tests, extensions: frmtmb.sample `test-default-priors-autocor.R`,
  `test-default-priors-dpar.R`, `test-flat-prior.R`,
  `test-me-prior-draws.R`, `test-prior-draws.R`; frmtmb.spatial
  `test-sar.R`.
- dev: `test-backlog.md` (six CLOSED lines, section "Filed in the round
  of 2026-11-13"), this file, `priors2-*.R`, `priors2-*.sh`.
- Fix round: `R/compat.R` (one row's text, F5), `R/fit-open-ends.R`
  (rewritten: box first), `R/fit.R` (`optimize_obj()` back to 0.82.0's
  run plus one rescue call and the failed-run hand-off;
  `fit_assembled()` refusal after the fit), `R/priors.R`
  (`box_end_density()`, `bound_inside()`, the L1 refusal, docs),
  `NEWS.md`, frmtmb.sample `NEWS.md`; tests `test-prior-bounds-ends.R`
  and `test-prior-mixture-theta.R`; dev `priors2-fix-probe.R`,
  `priors2-survey-hex.R`, `priors2-openblas-run.sh` (rewritten), one
  more LOW backlog item and the 5005 CLOSED line reworded.
- Fix round 2: `R/fit-open-ends.R` (`open_end_probe()`,
  `open_end_steps()`, `open_end_hold()`, per-side labels, the refusal
  and warning text), `R/fit.R` (`optimize_obj()`: the held warning),
  `R/priors.R` (`box_end_density()` labels, `internal_bound()`),
  `NEWS.md`; `test-prior-bounds-ends.R`; dev `priors2-fix-probe.R`
  (cases H, D', I), `priors2-openblas-run.sh` (`OB_FILES`, `OB_TAG`),
  one more LOW backlog item.

## Runs

All on the lane library (`wt-priors2-lib`: core and frmtmb.sample from
this worktree), rellib-r21 behind it; every log's `lib:` line names the
lane library (531 of 531).

### One file per process, every suite (`priors2-all-files/`, 2 at a time)

| package | files | pass | fail | error | skip | warn |
| --- | --- | --- | --- | --- | --- | --- |
| frmtmb | 332 | 24539 | 8 | 1 | 240 | 0 |
| frmtmb.sample | 91 | 4248 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 8 | 529 | 0 | 0 | 6 | 0 |
| frmtmb.latent | 13 | 697 | 0 | 0 | 2 | 0 |
| frmtmb.survival | 16 | 575 | 0 | 0 | 0 | 0 |
| frmtmb.eam | 33 | 1851 | 0 | 0 | 5 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |

The 9 core failures are 5 files, each an old expectation of a behavior
this lane changed, each edited and rerun alone (`t1-after-*.txt`):

- `test-priors-autocor-classes.R` (4): `identical(ub, log(3))` on
  sigma; the bound is now one ulp inside, `exp(log(3))` being above 3.
  Rerun: pass 186, fail 0. (frmtmb.sample `test-me-prior-draws.R` had
  the same on `sdme`, found by its session; rerun: pass 41.)
- `test-prior-one-slot.R` (2), `test-prior-location-dpar.R` (1 error):
  `cor` with `resp` accepted and suggested. Rerun: pass 91; pass 56.
- `test-ordinal-mixture.R` (1), `test-nonfinite-gradient.R` (1): a
  required warning that, on rellib-r21, was the step-start warning
  (`dev/priors2-stepmsg.R`: on both builds the only condition naming
  the family or "degenerate boundary" is the step-start finding, a
  warning there and a message here, logLik -388.457855925 and
  -543.287982475 on both). The tests now accept the warning or the
  message. Rerun: pass 99; pass 5.
- also found on the OpenBLAS arm: `test-prior-table-brms.R`'s 5005 case
  (an Inf stop at a closed `lb = 0`) now reaches the mode; the guard
  `nonfinite_not_converged()` keeps a unit test of its own in
  `test-prior-bounds-ends.R`. Rerun: pass 101.

### OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`)

`dev/priors2-openblas-run.sh`, output `ob0.3.26/`: dgemm 2000 in 0.54 s
(OpenBLAS). `test-prior-bounds-ends.R` pass 58, `test-prior-mixture-theta.R`
52, `test-priors-bounds-grcov.R` 54, `test-predfix.R` 92,
`test-mixture-starts.R` 135, `test-prior-table-brms.R` 101; fail 0 in
each. 0.3.32 was not run: no change of this lane sits near a
rounding-decided branch of BLAS or LAPACK (`bound_inside()` is double
arithmetic, decided the same on every BLAS).

### One session per package

frmtmb.sample and frmtmb.spatial, `dev/priors2-run-sessions.sh`
(`priors2-sess-ext.log`), three orders each:

- frmtmb.sample forward, reverse, preload: files 86/91, pass 4246,
  fail 1 (the `sdme` log(3) expectation above, since edited and rerun),
  error 0, skip 18, warn 0. Every LEAK line carries "[with load of
  ...]": an extension a test loads, not a test's own registration.
- frmtmb.spatial forward, reverse, preload: files 8/8, pass 529, fail 0,
  skip 6; one tagged LEAK (GDAL_DATA, with the load of sf).


Core, `priors2-sess-core.log` (after the five test edits above):

- forward: files 319/332, pass 24555, fail 0, error 0, skip 240, warn 0
- reverse: files 319/332, pass 24720, fail 0, error 0, skip 240, warn 0
- preload: files 319/332, pass 24823, fail 0, error 0, skip 242, warn 0

Every LEAK line carries "[with load of ...]" (frmtmb.sample,
frmtmb.survival, frmtmb.eam loaded by a test); none is untagged.

The other six suites were not run in one session: this lane's test
edits register nothing and those suites changed no file; they ran one
file per process above.

### Gated tier (`FRMTMB_BRMS_FIT_TESTS=true`, `gated/`)

`test-ordinal-mixture-brms.R` pass 24, `test-brms-priors.R` 102,
`test-prior-compat.R` 205, `test-brms-formula-priors.R` 211,
`test-prior-table-brms.R` 101, `test-get-prior-route.R` 83; fail 0,
skip 0 in each.

### R CMD check --as-cran (`dev/priors2-check.sh`, in `dev/priors2-check/`)

- frmtmb: Status: 2 NOTEs. Examples over 5 s (`frm_hazard_reads`
  11.69 s elapsed, `dharma_residuals` 7.52 s: load, per
  `dev/lane-rules.md`); the V8 math-rendering NOTE. Tests 95 min, OK.
- frmtmb.sample: Status: 1 NOTE (V8). Tests 14 min, OK.
- frmtmb.spatial (tests changed only): Status: 1 NOTE (V8).

The check trees were deleted after these lines were copied; the
one-file logs, session logs and probe outputs stay in
`dev/priors2-out/`. `dev/ciharden-out/` was deleted in the first round;
the reviewer rebuilt it, the fix round left it alone, and fix round 2
deleted it after its OpenBLAS run.

---

# Earlier lane: PRIORS-CONSISTENCY (wt-priors2, 2026-09)


The decision this lane implements, verbatim: "re: sigma, I'd like to
maintain consistency with brms". That is option 3 of
`dev/priors-findings.md` ("flip the meaning of the existing
spelling"), which the review at `dev/reviews/2026-09-05-priors.md`
recommended against on cost grounds. The maintainer overrode it.

## brms's own acceptance rules, re-measured here

`brms::make_stancode()`, brms 2.23.0, no compilation. Each cell is the
`lprior` line brms writes, or its refusal.

| model | `class = <dpar>` | `class = "Intercept", dpar = <dpar>` |
| --- | --- | --- |
| `y ~ x`, gaussian | `normal_lpdf(sigma \| 0, 3)` | ERROR "do not correspond to any model parameter: Intercept_sigma" |
| `bf(y ~ x, sigma ~ 1)` | ERROR "... : sigma" | `normal_lpdf(Intercept_sigma \| 0, 3)` |
| `bf(y ~ x, sigma ~ x)` | ERROR "... : sigma" | `normal_lpdf(Intercept_sigma \| 0, 3)` |
| `cnt ~ x`, negbinomial | `normal_lpdf(shape \| 0, 3)` | ERROR "... : Intercept_shape" |
| `cnt ~ x`, zero_inflated_poisson | `normal_lpdf(zi \| 0, 3)` | (not measured) |
| `bf(cnt ~ x, zi ~ 1)`, zip | ERROR "... : zi" | `normal_lpdf(Intercept_zi \| 0, 3)` |
| `p ~ x`, Beta | `normal_lpdf(phi \| 0, 3)` | (not measured) |
| `bf(p ~ x, phi ~ 1)`, Beta | (not measured) | `normal_lpdf(Intercept_phi \| 0, 3)` |

Two things this pins that the review left implicit:

1. the rule generalizes past `sigma`: `shape`, `zi` and `phi` behave
   the same way, so it is a rule about dpars and not about one name;
2. **`sigma ~ 1` counts as a predictor.** brms refuses `class =
   "sigma"` on `bf(y ~ x, sigma ~ 1)` exactly as it does on
   `sigma ~ x`. The distinction is "did the user write a formula",
   not "does the design have more than an intercept".

## Can frmtmb make that distinction?

Yes, from the spec alone (`resolve_priorlist()` gets `spec` and
`frame`, never `bform`: `R/fit.R:722` and `R/par-template.R:157` call
it with a two-element shim). Measured field sets of
`fit$spec$responses[[1]]$dpars$sigma`:

| model | fields |
| --- | --- |
| `bf(y ~ x)` | `name link fixed re rhs smooth constant` |
| `bf(y ~ x, sigma ~ 1)` | `name link constant fixed re smooth mo miterms csterms gpterms carterms spdeterms acterms rhs` |
| `bf(y ~ x, sigma ~ x)` | same 14 as above |

`plain_dpar()` (`R/parse.R:1244`) builds the first; a written formula
goes through `parse_linpred()`, which adds the special-term slots. The
design matrix is `(Intercept)` only in the first two, so `ncol(X)` does
NOT separate them. The written-formula marker does.

## What landed

| # | change | file:line |
| --- | --- | --- |
| 1 | `set_prior()` takes a distributional parameter's own name as a class and marks it `natural` | `R/priors.R:404-405`, `:452-470` |
| 2 | the class vocabulary is open: frmtmb's ten own names, the brms classes refused BY NAME (now from the native path too), then any dpar name | `R/priors.R:479` (`frmtmb_prior_classes`), `:493` (`check_dpar_prior_class`) |
| 3 | the brms route is a pass-through: one rule, written once | `R/priors.R:838` (`brms_prior_route`), `:706-709` |
| 4 | the shape gate, brms's own acceptance rule, by name in both directions | `R/priors.R:1921` (`dpar_shape_refusal`), called at `:1646` |
| 5 | `dpar_has_predictor()` reads the spec's record of who wrote the formula | `R/priors.R:1896` |
| 6 | `set_prior()` refuses an unhonored `coef` on `sd`/`cor` (review residual R1) | `R/priors.R:444-445` |
| 7 | a natural spec prints and reports as the class it was WRITTEN with | `R/priors.R:1338` (`spec_spelling`), print at `:963`, `spec_target()` at `:1313` |
| 8 | `get_prior()` lists a predictor-free dpar under its own class, as brms does | `R/priors.R:1144-1148` |
| 9 | the sample route keys its defaults by the spelling the table lists | `R/sampling-api.R:335-338`, exported at `:231` |

No `scale =` or `natural =` argument was added to `set_prior()`: after
the flip the placement follows from the class, which is brms's own
rule, so there is nothing left for such an argument to override. The
`natural` field survives as an internal marker only, because
`frmtmb.sample:::natural_dpar_prior()` writes it directly.

## Measured, before and after

Same data, same seeds, same starts; `p2-before2.R` / `p2-after.R`.

| shape | spelling | before | after |
| --- | --- | --- | --- |
| `y ~ x`, gaussian, n = 200 (ML sigma 2.010268914) | `set_prior("normal(0, 0.5)", class = "sigma")` | error (class refused) | **sigma 1.942946001** |
| the same | `set_prior(..., class = "Intercept", dpar = "sigma")` | sigma 1.996511359 | **refused by name** |
| the same | `brms::prior(normal(0, 0.5), class = "sigma")` | sigma 1.942946001 | 1.942946001 (unmoved) |
| `bf(y ~ x, sigma ~ 1)` | `class = "Intercept", dpar = "sigma"` | sigma 1.996511359 | 1.996511359 (unmoved) |
| `bf(y ~ x, sigma ~ 1)` | `class = "sigma"` | (unreachable) | **refused by name** |

The two routes now agree to every digit on the same model, which is
the point: `set_prior(class = "sigma")` and
`brms::prior(..., class = "sigma")` are one code path.

`get_prior()` tables, compared against `brms::get_prior()` row for row:

| model | brms | frmtmb before | frmtmb after |
| --- | --- | --- | --- |
| `y ~ x` | b, b/x, Intercept, **sigma** | b, b/x, Intercept, **Intercept+dpar=sigma** | b, b/x, Intercept, **sigma** |
| `bf(y ~ x, sigma ~ x)` | ... Intercept+dpar, b+dpar, b/x+dpar | same | same (unmoved) |

## The eam sign flip the earlier review priced does NOT happen

The review's case against the flip rested on two costs, and one of
them is not a cost of the flip as the maintainer specified it.

`extensions/frmtmb.eam/tests/testthat/test-surface.R:133` writes
`set_prior("normal(0, 0.1)", class = "b", dpar = "mu")`. That is class
`"b"`, which is link-scale on the SLOPES in brms and here, so this
change does not touch it. Measured on the same fixture (`ddm_fit()`,
seed 404, n = 400):

| quantity | before | after |
| --- | --- | --- |
| `mu.cond` unpenalized | 0.9165381216 | 0.9165381216 |
| `mu.cond` with the class "b" prior | **+0.2654012893** | **+0.2654012893** |
| `shrink` ratio | 0.2895692858 | 0.2895692858 |

The review's `-0.3301564534` is what a flip that ALSO retargeted class
`"b"` on a dpar would produce. brms does not do that, so neither does
this. The `(0, 1)` ratio assertion is left exactly as it was, and it
still passes; the block gains a positive test of the new spelling
(`class = "bs"`, `normal(1, 0.05)`, pulling bs from 1.382781771 to
1.27555188) and of the refusal of the link spelling on the same model.

## The other test that had to move, and it is not the shrink one

`extensions/frmtmb.eam/tests/testthat/test-surface.R:126` asserted
`all(c("bs", "ndt") %in% gp$dpar)`. Those two dpars have no predictor
in that model, so `get_prior()` now lists them under their own CLASS
with an empty `dpar`, which is brms's own row. The assertion moves to
`gp$class`. This is in the eam lane's file and outside the one block
the brief named, because the brief could not have known which block
`get_prior()`'s row change would reach.

`extensions/frmtmb.sample/tests/testthat/test-prior-route.R:26` moved
for the same reason and in the same way: it read the sample route's
sigma default off `gp$class == "Intercept" & gp$dpar == "sigma"`, and
the row is now `gp$class == "sigma"`. The DEFAULT still attaches
(`student_t(3, 0, 2.5)`), because `registered_prior_defaults()` keys by
the spelling the table lists.

## The substantive-differences table, refreshed

Only the rows this lane moves are restated; every other row of
`dev/priors-findings.md`'s table stands as written there.

| brms row | frmtmb after this lane | substantive difference remaining | measured consequence |
| --- | --- | --- | --- |
| `class = "sigma"` and every other dpar class, from `brms::prior()` | a density on the parameter itself, through its inverse link | none | unmoved: sigma 1.942946001 on the gaussian shape before and after |
| the same class from frmtmb's OWN `set_prior()` | **the same density, the same code path** | **none. The lane's row 7 divergence is gone** | the two routes give one objective at 1e-12 on six shapes covering Intercept, sd, cor, sigma, shape, zi and phi (`test-brms-priors.R`) |
| `class = "Intercept", dpar = <name>` | the link-scale intercept of that dpar's predictor | none | brms writes `Intercept_sigma` for exactly the models frmtmb accepts it on |
| either spelling on the WRONG model shape | refused by name, naming the other | none | brms errors "do not correspond to any model parameter"; frmtmb errors naming the spelling that applies |
| `class = "b", dpar = <name>` | the dpar's slopes, link scale | none | unmoved: the eam `mu.cond` stays +0.2654012893 |
| `get_prior()`'s dpar row | its own class where the dpar has no predictor | none | matches `brms::get_prior()` row for row on `y ~ x`, `sigma ~ 1` and `sigma ~ x` |

### Every spelling brms honors that frmtmb still cannot

Swept 17 families through `brms::get_prior(y ~ x, family = )` and
compared the non-`b`, non-`Intercept` classes brms writes against the
dpars frmtmb's own family gives the same model:

`sigma`, `nu`, `alpha`, `beta`, `shape`, `phi`, `zi`, `hu`, `kappa`,
`ndt`, `quantile` - **every one is a dpar frmtmb has**, so every one is
now a class `set_prior()` takes. gaussian, student, lognormal,
skew_normal, exgaussian, Gamma, weibull, Beta, negbinomial, poisson,
bernoulli, zero_inflated_poisson, hurdle_gamma, cumulative, von_mises,
shifted_lognormal and asym_laplace all resolve with no class left over.
The families that carry a class frmtmb cannot name are the ones whose
CLASS was already refused for a structural reason (`simo`, `sds`,
`sdgp`, `lscale`, `sdcar`, `car`, mixture `theta`), and that list is
unchanged by this lane.

So the answer to the brief's question is: **none**. There is no dpar
spelling left whose brms meaning frmtmb cannot honor. What remains
refused is refused because frmtmb holds the QUANTITY somewhere else,
not because it spells the prior differently.

## What this cost, and what it did not

| item | cost |
| --- | --- |
| `print(prior_summary())` on a natural spec | the line changed: `class=Intercept dpar=sigma scale=natural` becomes `class=sigma scale=natural`. Deliberate: the old line printed a spelling the same model now refuses, and `print.frmtmb_priorlist()` promises what it prints can be pasted back into `set_prior()` |
| `frm_sample()`'s defaults | unchanged in substance and in sampling. `default_priors_for()` still builds `class = "Intercept", dpar = "sigma"` + `natural` for a sigma with no predictor and the plain link-scale spec for one with a predictor, which is brms's own distinction, and `registered_prior_defaults()` keys it onto the table's `sigma` row |
| `frm_sample()`'s disclosure MESSAGE | still labels that default `Intercept (sigma)  [natural scale]`. `announce_default_priors()` (`extensions/frmtmb.sample/R/sample.R:728`) reads `s$class` directly. Now inconsistent with the print and with `get_prior()`. Left to the wt-sample-ce lane, whose file it is; `spec_spelling()` is exported for it (`R/sampling-api.R:224`) and the fix is one line |

## Files touched

Core, mine by the brief:

- `R/priors.R` - the class vocabulary, the normalization, the shape
  gate, `get_prior()`'s dpar row, the printed spelling, the roxygen.
- `R/sampling-api.R` - the default-key spelling, plus `spec_spelling`
  added to the exported seam so the sampling package can label a
  natural default the way the table does.
- `tests/testthat/helper-brms-priors.R`, `test-brms-priors.R`,
  `test-prior-compat.R`, `test-get-prior-route.R`,
  `test-simulate-ergonomics.R`.
- `NEWS.md`, `vignettes/brms-migration.Rmd`, `SPEC.md`,
  `vignettes/inputs.Rmd`, `man/` via roxygen, `NAMESPACE`.

Outside the brief's list, each forced by the change and each minimal:

- `R/simulate-new.R` - two lines of ONE roxygen example, which wrote
  `class = "Intercept", dpar = "sigma"` on a model with no sigma
  formula and would now stop under `R CMD check`.
- `vignettes/frmtmb.Rmd` - one prior-predictive chunk, the same
  spelling, same reason. It now writes `exponential(1)` on
  `class = "sigma"`, a positive density on a positive parameter.
- `extensions/frmtmb.eam/tests/testthat/test-surface.R` - the
  `get_prior()` column assertion (see above), plus two assertions on
  the new spelling in the same block. The shrink assertion the brief
  named is untouched, because it did not move.
- `extensions/frmtmb.sample/tests/testthat/test-prior-route.R` - one
  line, the same `get_prior()` column.

## Verification

One `test_file()` per process, private library
`.../scratchpad/p2-lib` built from this worktree, my own copy of
`dev/stan-cache` (54 programs), `R_MAKEVARS_USER` at its
`makevars-cxx17.mk`. Every Stan program was a cache hit; nothing
compiled.

| run | result |
| --- | --- |
| `test-brms-priors.R` (gated) | **103 pass**, 0 fail (was 83 before the lane; +20 from the route-agreement block and the retargeted row 5) |
| `test-brms-likelihood.R` (gated) | **373 pass**, 0 fail, unmoved |
| `test-prior-compat.R` | **185 pass**, 0 fail (was 161) |
| `test-setprior.R` | 27, unmoved |
| `test-get-prior-route.R` | **37** (was 36; one row assertion became two) |
| `test-priors-autocor-classes.R` | 63, unmoved |
| `test-priors-bounds-grcov.R` | 49, unmoved |
| `test-simulate-ergonomics.R` | **50 pass** (was 45; two new assertions and the block that used to stop) |
| `test-sample-direct.R` (frmtmb.sample) | 136, unmoved |
| `test-message-uniqueness.R` | 6, unmoved |
| `test-bracket-access.R` | 8, unmoved |
| **full core suite**, one file per process | **117 of 117 files ran**, list `diff`ed against `ls tests/testthat/test-*.R` (identical, 0 duplicates), **6562 pass, 0 fail**, 93 skip, 2 errors |

The 2 errors are `test-influence-plot.R:68` and `:130`, both
`local_mocked_bindings()` without a `.package` argument, which needs
pkgload and so fails under a bare `test_file()` harness and not under
`test_check()`. `R CMD check`'s own `checking tests` runs the same file
green, which is the confirmation. The same artifact hits
`frmtmb.sample`'s `test-parallel-chains.R:54`. Neither touches priors.

### as-cran

`R CMD build` (vignettes built) then
`R CMD check --as-cran --no-multiarch`, `_R_CHECK_CRAN_INCOMING_=false`,
pandoc 3.8.3 on `PATH`.

| check | result |
| --- | --- |
| `checking tests` | **[31m] OK**, `FAIL 0 | WARN 1 | SKIP 98 | PASS 6554` |
| `checking examples` / `--run-donttest` | [76s] OK / [117s] OK |
| `checking re-building of vignette outputs` | **[427s] OK** |
| `checking files in 'vignettes'` / `package vignettes` | OK / OK |
| all seven Rd checks, `code/documentation mismatches`, `Rd contents` | OK |
| Status | 1 ERROR, 1 WARNING, 2 NOTEs |

The four remaining items are the environmental ones this machine always
has, each with a verified cause: `Rdlatex.log` says **"pdflatex is not
available"** (the manual WARNING, the ERROR, and the leftover
`frmtmb-manual.tex` NOTE) and **"package 'V8' unavailable"** (the HTML
NOTE). `Sys.which("pdflatex")` and `Sys.which("texi2dvi")` are both
empty here. Nothing that reads the package's code, documentation,
tests, examples or vignettes is anything but OK.

`_R_CHECK_FORCE_SUGGESTS_=false` was needed: `brokenstick` and
`frmtmb.spline` are Suggests and neither is installed on this machine.
Both uses are guarded (`skip_if_not_installed()`,
`requireNamespace()`), and main has since dropped `frmtmb.spline` from
core's Suggests (`612cde1`), so only `brokenstick` will remain.

The one WARNING inside the test run is
`test-prior-compat.R:505`, "Optimizer did not report convergence:
singular convergence (7)", on the `(x|g) + (z|h)` shape the previous
lane added. It is PRE-EXISTING: the same warning came out of that file
before any edit of mine.

### Merge

Main advanced to `b131fe1` during the run (`612cde1`, `de9d639`,
`b131fe1`). Those three commits touch `.github/`, `DESCRIPTION`,
`dev/build-docs.R`, `docs/**` and `vignettes/case-studies.Rmd`, and
this lane touches **none of them**, not even `NEWS.md`. There is no
overlap with main at all; the usual `NEWS.md` collision will only
appear against whichever sibling lane merges first.

## Where frmtmb's meaning still differs substantively from brms

The complete list after this lane, each row with what it costs. Only
the last two rows are about a prior SPELLING at all; the rest are
about the quantity underneath it.

| # | what | difference | measured consequence |
| --- | --- | --- | --- |
| 1 | a prior is a PENALTY | `frm()` is MAP, not a posterior, so the reported log likelihood, AIC and `anova()` are penalized quantities | the estimates land close where priors are weak (kidney `sd(patient)` 0.38 against brms's posterior mean 0.40) but they are not the same quantity |
| 2 | `class = "sd"` | brms writes the NORMALIZED half-t (density minus one `lccdf` per element); frmtmb evaluates the unfolded density plus the log-Jacobian | exactly `log(2)` per lower-bounded PRIORED element. A constant: it moves no mode. S1 `dT = 2 log 2` |
| 3 | `class = "cor"`, `"cortime"`, `"rescor"` | the same LKJ density on the same correlation, carried onto frmtmb's coordinate with that map's exact Jacobian, which is not Stan's coordinate | `(eta + (d-1)/2) log(1 - rho^2)`, measured -0.0021842 on S3: a coordinate change, not a different density |
| 4 | `class = "ar"`, `"ma"`, `"cosy"` | the same, for the stationarity map | ANALYTIC, not measured against Stan: brms declares `vector<lower=-1,upper=1>[Kar] ar` under `cov = TRUE` and frmtmb uses its own map |
| 5 | `simo` (`dirichlet`) | REFUSED. frmtmb holds a `mo()` simplex as free softmax coordinates and puts no density on it | a `mo()` table stops the call, naming the row |
| 6 | `sds`, `sdgp`, `sdcar` | REFUSED, presentational. The parameters exist as random-effect blocks and are reachable by hand (`class = "sd"`, `group = "s(x)"`), but `as_priorlist()` has no model in hand to read the label from | one named edit round per table |
| 7 | `lscale`, `car` | REFUSED. frmtmb keeps them in the raw internal covariance vector, so they are reachable only on the INTERNAL scale through `class = "theta"` | the density a user writes is not the density they get unless they rewrite it for the internal scale |
| 8 | mixture `theta`/`theta1`/`theta2` | REFUSED as a bare class. A modeled proportion with its own predictor is `class = "Intercept", dpar = "theta1"`, which is honored | brms writes no Stan statement for its reference component either |
| 9 | `horseshoe()`, `R2D2()`, `lasso()` | REFUSED. frmtmb has no hierarchical-shrinkage prior at all | the nearest honest substitute is a tight `normal()` on class `"b"`, which is a different estimator |
| 10 | `constant()`, `dirichlet()`, `uniform()` | REFUSED. `constant()` is brms's way of FIXING a parameter, whose frmtmb counterpart is a constant dpar in `bf()`; `dirichlet` needs a simplex target frmtmb does not have | those tables stop |
| 11 | a distributional parameter's class | **none any more.** Every class brms writes for a dpar is a class `set_prior()` takes, with brms's meaning and brms's model-shape rule | 17 families swept: `sigma`, `nu`, `alpha`, `beta`, `shape`, `phi`, `zi`, `hu`, `kappa`, `ndt`, `quantile`, all honored |
| 12 | `class = "theta"` | frmtmb-only. brms has no spelling for a raw internal covariance parameter | an addition, not a divergence |

Rows 2 to 4 are coordinate or normalization differences that no mode
can see. Rows 5 to 10 are refusals: frmtmb holds the quantity
somewhere else or does not hold it, and every one of them is refused
BY NAME from both routes now, which is the property this lane added to
the ones that were only refused on the translated route.

### Extension suites, against this lane's core

| package | files | result |
| --- | --- | --- |
| frmtmb.eam | 17 of 17 | **1292 pass, 0 fail, 0 error**, `test-surface.R` at **50** (was 47) |
| frmtmb.sample | 10 of 10 | **899 pass, 0 fail**, 1 error, `test-sample-direct.R` at **136** and `test-prior-route.R` at **9** |

The one frmtmb.sample error is `test-parallel-chains.R:54`, the
`local_mocked_bindings()`/pkgload harness artifact described above, not
a prior.

# Punch round, 2026-09-06

Against `dev/reviews/2026-09-05-priors2.md` (verdict PUNCH). The
review's own verification of the placement (4e-15 across five links)
and of the acceptance rule (60 of 60 cells on the two spellings this
round is about) is not repeated here; what follows is the nine items.

## 1. F1, blocker. `get_prior()` no longer names a class `set_prior()` refuses

`R/priors.R:1144-1157`. The own-class branch now skips any dpar whose
class `brms_prior_class_refusal()` rejects, which is the same gate
`set_prior()` uses, so the two cannot disagree by construction.

| model | class column before | after |
| --- | --- | --- |
| `bf(y ~ x) + mixture(gaussian(), gaussian())` | Intercept, b, b, sigma1, sigma2, **theta1** | Intercept, b, b, sigma1, sigma2 |

The row is DROPPED rather than sent back to `Intercept` + `dpar`,
because that spelling is refused on that model too (theta1 has no
predictor): falling back would have restored a quieter version of the
same lie. A mixture proportion frmtmb cannot prior is now absent from
the table, and `?set_prior` says so.

The regression test is `tests/testthat/test-get-prior-route.R:242-317`,
"every row get_prior() lists is one set_prior() accepts". It is
stronger than the review asked for: for **every row** of **seven
shapes** (ordinary, distributional, intercept-only dpar, negbinomial,
zero-inflated, Beta, mixture) it re-spells the row through
`set_prior()` AND resolves it against the fit, so a row that parses
but addresses nothing fails too. That is what makes it a test of the
promise rather than of one column.

## 2. F2, blocker. `class = "b", dpar = X` on an intercept-only predictor is refused

`R/priors.R:1976-1993`, the third case in `dpar_shape_refusal()`.
`dpar ~ 1` gives the parameter a predictor (so it takes the Intercept
spelling) but no population-level slopes (so class `"b"` addresses an
empty set). Measured before: the penalty was bit-identical to no prior
at all, with no message. Measured after: refused, naming
`class = "Intercept", dpar = X`.

Tests: `tests/testthat/test-prior-compat.R:703-716` for the refusal and
for `sigma ~ x` still being accepted, and `:718-737`, a new block that
asserts **brms refuses the same row**
(`brms::validate_prior()` answers "do not correspond to any model
parameter") and accepts it once sigma has slopes. The rule is brms's,
asserted against brms, not restated from the review.

## 3. F3. `frm_sample()` announces the spelling the table lists

`extensions/frmtmb.sample/R/sample.R:724-735`: `lab` is built from
`spec_spelling(s)`. Taken rather than deferred, as the review asked;
the file is this lane's now that wt-sample-ce has merged.

```
before:   Intercept (sigma)  student_t(3, 0, 2.5)  [natural scale]
after:    sigma              student_t(3, 0, 2.5)  [natural scale]
```

which is the word `get_prior(route = "sample")` prints and the word
`set_prior()` accepts.

## 4. F7. `coef` and `group` on a natural class are refused

`R/priors.R:726-763`, `unhonored_coef_refusal()` gains the natural
case, and both call sites pass `group` (`:444`, `:693-694`). The call
already saw the class the user WROTE, so no move was needed. `resp` is
still accepted, because it is the one narrowing such a class takes.
Test: `tests/testthat/test-prior-compat.R:772-786`.

## 5. F4. Which mode, in the migration vignette

`vignettes/brms-migration.Rmd:80-95`, next to the `log(2)` sentence.
Re-measured here rather than copied, because the review's data is not
in the tree; gaussian `y ~ x`, n = 150, `normal(2.5, 0.4)` on sigma:

| quantity | mine | the review's |
| --- | --- | --- |
| frmtmb penalized sigma | **2.133344** | 2.170587 |
| mode of `L * p` in sigma coordinates | **2.126533** | 2.163707 |
| mode in log-sigma coordinates | **2.133358** | 2.170589 |
| unpenalized ML sigma | 2.091048 | - |

Same phenomenon, same gap to three digits (0.0068 against 0.0069), on
different draws. The vignette carries MY pair, because a number in the
repository should be one the repository can reproduce; the review's
pair is recorded here as the independent confirmation.

## 6. F8 and F6, documented in the same pass

`vignettes/brms-migration.Rmd:97-105`. A resp-less natural class
broadcasts across a multivariate model where brms refuses the row:
**kept, not refused**, because it is frmtmb's own convention for every
class-wide prior (`class = "sd"` with no `group` covers every block,
`class = "b"` covers every slope), and making one class family the
exception would be a worse inconsistency than the difference with
brms. Documented instead, with `resp` named as the way to narrow.
`cumulative()` has no `disc` dpar at all, so the acceptance rule is
never reached there; also documented.

## 7. F5. `?frm_simulate` says which scale a dpar column is on

`R/simulate-new.R:497-506`.

## 8. NEWS

`NEWS.md:30-56`: the `get_prior()` bullet gains the dropped mixture
row, and three bullets are added for the `b`/dpar refusal, the
`coef`/`group` refusal, and the changed `print()` / `spec_target()`
text.

## 9. Re-measured: as-cran on a tarball built WITH vignettes

`R CMD build` (no `--no-build-vignettes`), then
`R CMD check --as-cran --no-multiarch`, `_R_CHECK_CRAN_INCOMING_=false`,
`_R_CHECK_FORCE_SUGGESTS_=false`, pandoc on `PATH`, and **no
`--library=`** this time. Status: **1 ERROR, 1 WARNING, 2 NOTEs**, and
every one of the four is named with its cause below.

The review's two vignette WARNINGs are gone, which confirms they were
its own `--no-build-vignettes` artifact: `checking files in
'vignettes'`, `checking package vignettes` and `checking installed
files from 'inst/doc'` are all OK, and `frmtmb.Rcheck/frmtmb/doc/`
holds the eight built `.html` files.

| item | text | measured cause |
| --- | --- | --- |
| `checking PDF version of manual` WARNING | "LaTeX errors when creating PDF version" | `frmtmb.Rcheck/Rdlatex.log` says **"pdflatex is not available"**, twice. `Sys.which("pdflatex")` and `Sys.which("texi2dvi")` are both `""`, and `R_LATEXCMD` is unset |
| `checking PDF version of manual without index` ERROR | "Re-running with no redirection of stdout/stderr" | the same retry, the same missing pdflatex |
| `checking HTML version of manual` NOTE | "Skipping checking math rendering: package 'V8' unavailable" | `requireNamespace("V8")` is FALSE on this machine |
| `checking for non-standard things in the check directory` NOTE | `'frmtmb-manual.tex'` | the leftover of the two items above: `Rd2pdf` writes the `.tex`, `texi2dvi` then fails, so nothing cleans it up. There is **no `frmtmb-manual.pdf`** beside it |

One INFO, not counted in the status: `checking package dependencies`
reports "Package suggested but not available for checking:
'frmtmb.spline'". It is a Suggests of core at this base and is not
installed here; main dropped it from core's Suggests at `612cde1`, so
it disappears on merge. `brokenstick`, the other Suggests, IS
installed and is not reported.

### The reviewer's doubt about pdflatex, settled

The review inferred that pdflatex must be available because
`frmtmb-Ex.pdf` was produced. That file is not a LaTeX product:

```
$ grep -a -o "Producer[^)]*)" frmtmb.Rcheck/frmtmb-Ex.pdf
Producer (R 4.6.1)
$ grep -a -o "Creator[^)]*)"  frmtmb.Rcheck/frmtmb-Ex.pdf
Creator (R)
```

`frmtmb-Ex.pdf` is R's own `pdf()` graphics device output from running
the examples, which needs no TeX at all. The LaTeX product would be
`frmtmb-manual.pdf`, and it does not exist in either check directory,
while `frmtmb-manual.tex` does. Together with `Rdlatex.log`'s explicit
"pdflatex is not available" and the two empty `Sys.which()` results,
that is three independent measurements agreeing. The attribution
stands; the inference from `frmtmb-Ex.pdf` does not.

### One correction to my own earlier report

My first as-cran attempt died on "Packages suggested but not
available: 'brokenstick', 'frmtmb.spline'", and I recorded brokenstick
as uninstalled. That was wrong, and the cause was my own
`--library=<private lib>` on `R CMD check`: it also narrows where
check LOOKS for suggested packages, hiding the user library where
brokenstick lives. The same flag was what killed 11 files of the first
full-suite run mid-flight, because check installs the package into
that library. Dropped here, as instructed.
