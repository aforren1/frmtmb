# Lane ordinal: disc, threshold structures and acat links

Base: 57c25589 (frmtmb 0.66.0), branch `wt-ordinal`, uncommitted.
Reference: brms 2.23.0, read from its source (`dev/ordinal-brms-read.R`,
`-read2.R`, `-read3.R`, `-grep.R`) and from its generated Stan code and
data (`dev/ordinal-brms-code.R`); outputs are the matching
`dev/ordinal-log-brms-*.txt`.

## Punch round 1 (review `dev/reviews/2026-09-30-ordinal.md`)

Verdict of round 0: NOT MERGEABLE, three blockers and seven minors.
What changed, with the evidence. The earlier sections are round 0's
record; where this section withdraws one of them, that section says so.

**B1. Loop closures in `draws_ordinal_cols()` (frmtmb.sample).** The
maps read `K` and `dd` from the loop frame, so every block took the
last block's. Fixed by building each block in `draws_ordinal_block()`,
whose closures own their values. New test: two ordinal responses of
different families and threshold counts (5 and 3 categories), under
flexible, equidistant, and equidistant beside sum_to_zero;
`posterior_epred()` at the first and last draw against brms's R-side
density at the stored columns, and each stored `delta_<resp>` against
the stored spacing. On the round-0 lane build it fails
(`dev/ordinal-p1-log-ts-draws-before.txt`: pass=32 fail=3 err=1; the
flexible epred is off by 0.121, `delta_y` by 0.791, the mixed model
stops in sampling). After the fix: pass=44 fail=0
(`dev/ordinal-p1-log-ts-draws-after.txt`). The review's own scripts on
the fixed build: `dev/ordinal-p1-log-mvsilent.txt` (max |delta_y -
spacing| = 0, was 1.48) and `dev/ordinal-p1-log-mvdraws-lane.txt` (all
four shapes give `posterior_epred()`).

**B2. The family-default `disc` is hidden.** `ord_family_tail()` sets
`fam$hidden_fixed_dpars <- "disc"` for all five families.
`lp_hidden_fixed()` is TRUE for a constant linear predictor whose dpar
the family hides and whose value is the family's own default; such a
predictor is left out of the `Links:` string (`family_link_str(fam,
shown =)`, with `lp_shown_dpars()` naming a modeled disc), the
`Fixed dpar:` block and `summary()$fixed_dpars`, the summary's
coefficient blocks, `fixef(flatten = TRUE)` (and so `frm_bootstrap()`'s
default statistic and `allFit`'s table) and `coef()`. The mapped
coefficient stays in the objective. A disc the user fixes at another
value (`bf(disc = 2)`) still shows. `test-thres-refit.R` is back to the
base file; `test-ordinal.R` and `test-xbeta-zibb-hurdle-cum.R` pin the
hidden forms.

Against rellib-r4, `dev/ordinal-p1-bitwise.R` (the review's 26-model
script, copied; data seed 20261002) and
`dev/ordinal-p1-bitwise-compare.R`, log
`dev/ordinal-p1-log-bitwise-compare.txt`: 570 of 581 outputs
`identical()`. `summary()`, `print()`, `fixef()`, `fixef(flatten =
TRUE)`, `coef()`, `vcov()` and `confint()` are identical on every
flexible ordinal model (cumulative logit, probit with `(1 | g)`,
cloglog with `thres(gr = )`, ordered factor, sratio with `cs()`,
sratio cauchit, cratio cloglog with `(1 | g)`, cratio probit_approx,
acat with `cs()`, acat with `thres(gr = )`, the ordinal-gaussian
multivariate model) and on the hurdle with disc modeled. The 11 that
differ, each intended:

- `default_prior()` of the two `thres(gr = )` models: the per-level
  `Intercept` rows of m7.
- the hurdle with disc at 1 (9 outputs): `summary()` and `print()` lose
  `disc = log` and `Fixed dpar: disc = 1`, `fixef(flatten = TRUE)`,
  `coef()` and the bootstrap lose the constant `disc_(Intercept)`
  (B2 for the hurdle, which on the base showed them); the draws carry
  `b_Intercept[k]` where the base had `tau_raw_k` (round 0's hurdle
  map), so `variables()`, the draws matrix and the draws summary
  differ by name, and `posterior_epred()` of the draws by at most
  1.1e-16, from the round trip through the stored thresholds; the
  sampled columns they share are identical.

So the coordinator's "identical on the hurdle" holds for every output
but the ones B2 asks to change there.

**B3. sum_to_zero with one threshold per vector.** The refusal is gone:
`thres_structure_check()` now refuses equidistant alone. A response
whose every threshold vector has one threshold holds them all at 0 and
has an empty `tau_raw`. That needed empty-block safety where names were
built with `paste0(cp, "_", seq_along(v))`, which recycles an empty
index to one name (`R/par-template.R`, `R/confint.R`,
`R/brms-names.R`, `R/interop.R`, `R/priors.R`, frmtmb.sample's
`draws-brms.R` and `methods-draws.R`, now `recycle0 = TRUE`), and the
empty block now reports its thresholds (`b_Intercept[a,1] = 0`) as brms
does. `dev/ordinal-p1-stz1.R` (seed 20261001,
`dev/ordinal-p1-log-stz1.txt`): with every `thres(gr = g)` level at two
categories and `disc ~ 0 + z`, cumulative, sratio and acat fit, the
logLik is brms's density at thresholds 0 to at most 7.7e-16 relative,
`fitted()` to 1.7e-16, `summary()`, `simulate()` and `predict(newdata
=)` run; the ungrouped two-category model equals
`glm(y == 1 ~ 0 + I(-x))` to 4.6e-16. New tests pin both.

**m1. Disc at a fixed point, and the inverse map.** New ungated tests:
each density (4 families, 2 links each, ungrouped and under `thres(gr
= g)`) is evaluated with `obj$fn()` at the fit's parameters with the
disc coefficient set to 0.7, against brms's R-side density there; the
threshold maps are inverted for every structure, ordered and not, on a
three-group layout; and frmtmb.sample's sum-to-zero `posterior_epred()`
is checked against brms at the stored columns. The review's four
escaping mutants, rebuilt from the current source by
`dev/ordinal-p1-mutants.R` into `dev/ordinal-p1-mutlib/` (removed after
the runs; the script rebuilds them), now fail ungated
(`dev/ordinal-p1-log-mut-*.txt`):

| mutant | test-ordinal-disc-thres.R | frmtmb.sample draws test |
|---|---|---|
| M2 acat logit data path, no disc | fail=1 | |
| M2b grouped sequential, no disc | fail=4 | |
| M2c cratio, no disc | fail=2 | |
| M3b sum-to-zero inverse, unordered | fail=1 | fail=4 |

**m2. Grouped flexible thresholds are bitwise again.** Their reported
thresholds go through the per-slice `ord_tau_from_raw()` map of the base
again; the running sum of `thres_tau()` rounded differently. In the
comparison above `cum_cloglog_gr` and `acat_gr` differ only by m7's
prior rows (their summaries and `vcov()` are identical).

**m3.** The comment in `test-brms-likelihood.R` now says the cratio
cloglog NaN depends on the data (a row with `disc * (mu - thres_k)`
above about 6.6 below its category), and a gated cratio cloglog row
with disc runs on seed 20260931, whose largest such value is 4.70:
identity to 5.7e-14 absolute (1.3e-16 relative), brms's gradient 1.6e-4
(`dev/ordinal-p1-cloglog.R`, `dev/ordinal-p1-log-cloglog.txt`; seed
20260930 has 8.75).

**m4.** brms defect 4 is wider than round 0 recorded: a multivariate
model with ONE equidistant response also fails ("delta_y not in
scope"), and so does an equidistant mixture (review log
`dev/ordinal-rev-log-defects.txt`). The hurdle logit `nthres + 2`
defect the review found is brms-1 of `dev/upstream-bugs.md` (on main
and in this tree); round 0's defects 1 to 5 are not in that file yet.

**m5.** NEWS says what changed on every output (the hidden disc, the
hurdle's hidden disc, the per-level prior rows) and replaces the
"unchanged to the bit" claim with the 570 of 581 comparison. The
multivariate `coef()` change the review saw is gone with B2.

**m6.** The disc-intercept warning also fires when the design has no
intercept column but spans one (`disc ~ 0 + h`), naming the columns
and advising to leave one out or hold them with a class `"b"` prior;
a prior on them silences it. Tested both ways.

**m7.** `default_prior()` shows `lb = 0` on `delta` for cumulative and
the hurdle (not for the unordered families, where brms has none), and
lists one `Intercept` row per level of `thres(gr = )`. Tested.

**New tests seen to fail on the round-0 lane build** (the review's
scratch install of it, `mut:` arm of `dev/ordinal-runtest.R`, logs
`dev/ordinal-p1-log-prepunch-*.txt`): core
`test-ordinal-disc-thres.R` pass=143 fail=10 err=2 (B2, B3, m6, m7),
`test-ordinal.R` fail=1 err=1, `test-thres-refit.R` fail=5,
`test-xbeta-zibb-hurdle-cum.R` fail=3, frmtmb.sample
`test-ordinal-disc-thres-draws.R` fail=3 err=1 (B1).

**Every suite again**, `dev/ordinal-run-par.sh` over all 329 files with
every gate set, logs `dev/ordinal-suite-p1all/`, totals
`dev/ordinal-p1-log-suite-sum.txt` (every log loads frmtmb from the
lane library):

    pkg files  pass fail err skip warn
    frmtmb   194 16058    8   0    0    0
    frmtmb.coupling    11   542    0   0    5    0
    frmtmb.eam    29  1743    0   0    3    0
    frmtmb.latent    10   359    0   0    2    0
    frmtmb.learn    15   500    0   0    2    0
    frmtmb.ode    11   547    0   0    1    0
    frmtmb.sample    44  2472    1   0    1    0
    frmtmb.spline    15   553    0   0    1    0

The 9 failures are the same 9 stale "now HOLDS" rows; every skip is a
`test-scale.R`. The gated `test-brms-likelihood.R`: pass=543 fail=0.

**R CMD check --as-cran** (`dev/ordinal-check.ps1`, built with vignettes
and manual from the final source): frmtmb `Status: 1 NOTE`, the
environmental V8 math-rendering note; frmtmb.sample `Status: OK`.

## What brms does

Read from brms 2.23.0, not recalled.

- **disc.** `cumulative()`, `sratio()`, `cratio()`, `acat()` and
  `hurdle_cumulative()` all list `dpars = c("mu", "disc")` (plus `hu`)
  and take `link_disc = "log"`; `links_dpars("disc")` is `log`,
  `identity`, `softplus`, `squareplus`. `brmsterms.brmsformula()` sets
  `pfix$disc <- 1` for every disc parameter nobody wrote, and Stan data
  then carries `disc = 1` (`real disc = 1` in transformed parameters).
  The density reads `disc * (thres - mu)` for cumulative and sratio and
  `disc * (mu - thres)` for cratio and acat
  (`brms:::stan_ordinal_lpmf()`, `dcumulative()` and kin), with `cs()`
  subtracted from the thresholds first (`Intercept - transpose(mucs)`).
  With `disc ~ z`, the Stan program declares `Intercept_disc` and
  `b_disc`, centers `X_disc` and reports `b_disc_Intercept`; its
  default prior is `normal(0, 1)` on `Intercept_disc` under every link
  but the identity, where it is `lognormal(0, 1)`
  (`brms:::def_dpar_prior()`). disc works with `cs()`, with
  `thres(gr = )` (`*_merged_lpmf` passes it) and with mixtures
  (`disc1`, `disc2`).
- **equidistant.** Parameters `first_Intercept` and `delta`;
  `Intercept[k] = first_Intercept + (k - 1) * delta`. `delta` has
  `lb = 0` for the ordered families (`prior_thres()`, `lb <-
  str_if(has_ordered_thres(bframe), "0")`) and no bound otherwise. The
  design is centered as for flexible thresholds. Class `Intercept`'s
  prior is on `first_Intercept`; no per-coef Intercept rows exist
  (`validate_prior(prior(..., class = Intercept, coef = 1))` fails:
  "do not correspond to any model parameter"). Default: `delta` flat.
  Under `thres(gr = )`: `first_Intercept_<k>` and `delta_<k>` per level
  (numbered, not named), prior rows `delta` plus one per level.
- **sum_to_zero.** `Intercept` declared as usual (ordered for the two
  ordered families), `Intercept_stz = Intercept - mean(Intercept)` used
  in the likelihood, `b_Intercept = Intercept_stz`, and the design is
  NOT centered (`stan_center_X()` is FALSE under sum_to_zero). The
  common location of `Intercept` is a direction the likelihood does not
  see; brms's default `student_t(3, 0, 2.5)` on class `Intercept`
  places it.
- **Refusals.** The constructor `match.arg`s `threshold` against the
  three names. `Cannot use equidistant and fixed thresholds at the
  same time` applies to mixtures only. Every structure's Stan data
  declares `int<lower=2> nthres` ungrouped and `int<lower=1>` per group
  under `thres(gr = )`, so brms cannot run an ungrouped two-category
  ordinal model at all.
- **acat off the logit.** `brms:::inv_link_acat()` and the generated
  `acat_<link>_lpmf` use `prod_{j<k} F(x_j) * prod_{j>=k} (1 - F(x_j))`,
  normalized, with `x_j = disc * (mu - thres_j)`.
- **fitted(scale = "linear").** `posterior_epred()` returns `mu` (plus
  `cs()` per threshold); `incl_thres = TRUE` of `posterior_linpred()`
  is the only per-threshold `disc * (thres - mu)` output. frmtmb's
  `fitted(scale = "linear")` already returned `eta`, or `eta + cs_k`
  per threshold under `cs()`, and disc does not enter either.

## What changed

Core (`R/families.R`, `R/thres.R`, `R/priors.R`, `R/confint.R`,
`R/methods-fit.R`, `R/frame.R`, `R/links.R`, `R/sampling-api.R`):

1. `cumulative()`, `sratio()`, `cratio()`, `acat()` take
   `link_disc = "log"` and `threshold = "flexible"` in brms's order and
   have dpars `mu`, `disc`, with `fixed_dpars = list(disc = 1)`, the
   route `hurdle_cumulative()` already used. Every density branch reads
   disc: the data path, the one-step (OSA) path, the grouped
   `thres(gr = )` densities, `ord_cat_probs()` and the simulators.
   `fam_sratio()` and `fam_cratio()` became one `fam_sequential()`.
2. The disc-intercept warning (`hurdle_cum_fit_check()`, now
   `ord_fit_check()`) runs for every ordinal family and names it.
3. `threshold = "equidistant"` and `"sum_to_zero"` on all five
   families. `thres_finalizer()` builds the structure once the count is
   known, from two factories every constructor leaves on the family
   (`ord_lpdf_make`, `ord_sim_make`, set by `ord_family_tail()`).
   `thres_layout()` carries the internal-parameter slice of each vector
   (`rstart`, `rend`) beside the threshold slice; `thres_tau()` maps
   the internal vector to the thresholds (taped and plain),
   `thres_raw_from_tau()` is its inverse, `thres_start()` gives start
   values. Internal parameters: equidistant `(tau_1, log delta)` for
   the ordered families and `(tau_1, delta)` otherwise; sum_to_zero the
   `K - 2` log increments of a centered ordered vector, or the first
   `K - 2` thresholds with the last minus their sum.
4. Consumers that read the threshold count off the internal vector now
   read `fam$thres$nthres`: `thres_ncat()`, the `cs()` coefficient
   length in `frame.R`, the refit pin `thres_pin_of_fit()`.
5. Priors: class `"delta"` (frmtmb's class list, brms's direct
   classes, resp-keyed), placed on delta with the log-Jacobian for the
   ordered families (the `"sd"` placement) and on the internal value
   otherwise; `group` selects a level of `thres(gr = )`. Under
   equidistant, class `"Intercept"` is the first threshold only, with
   the centering offset. Under sum_to_zero, class `"Intercept"` is
   refused with the reason (below). `default_prior()` lists `delta`
   (once, and once per level under `thres(gr = )`) and no `Intercept`
   under sum_to_zero.
6. Names: `variables()`, `hypothesis()` and `summary()$spec_pars`
   report `delta` (`delta_<k>`, `delta_<resp>`) through
   `ord_delta_info()`; `summary()` gives its Wald interval on the log
   scale, mapped back. `b_Intercept[k]` comes from the family's map.
7. `acat()` takes brms's six links; off the logit the density is
   `acat_general_E()` (data, OSA, grouped and probability paths).
8. `hurdle_cumulative()` declares its inverse threshold map, so its
   draws carry `b_Intercept[k]`.

frmtmb.sample (`R/draws-brms.R`, `R/sample.R`, `R/methods-draws.R`):

9. An ordinal draws block may report more names than it has internal
   parameters. The first `R` replace the internal columns in place and
   the rest (the thresholds past `R`, and `delta`) are appended before
   `lp__`, the mixture-weight pattern; the inverse drops them.
10. Default prior `normal(0, 1)` on the intercept of `disc` under every
    link but the identity; the identity case is named in the
    disclosure.
11. `posterior_linpred(incl_thres = TRUE)` is still refused, with a
    reason that no longer misstates brms (below).

Tests: `test-ordinal-disc-thres.R` (new), rows 12e to 12h of
`test-brms-likelihood.R` (gated), frmtmb.sample's
`test-ordinal-disc-thres-draws.R` (new); `helper-brms.R` translates
`first_Intercept`, `delta`, `Intercept_<k>` and the uncentered
sum-to-zero `Intercept`; `helper-brms-suite.R` (and its frmtmb.sample
copy) views a constant dpar as brms's data slot, the family's threshold
count as `nthres`, and brms's `Z_<id>[_<dpar>]_<k>`. Adjusted:
`test-ordinal.R` (the link string, acat links), `test-brms-families.R`
(acat links), `test-mv-gaps.R` and frmtmb.sample's
`test-ordinal-draws-names.R` (hand-built parameter vectors that must
now drop the mapped disc coefficient), `test-thres-refit.R`
(`fixef(flatten = TRUE)` carries `disc_(Intercept)`),
`test-xbeta-zibb-hurdle-cum.R` (the lifted refusal).

## Measurements

### The flexible model is unchanged to the bit

`dev/ordinal-bitwise.R` (seed 20260930, n = 300), base against lane,
`identical()` on `opt$par`, `logLik`, `fitted()` and the OSA residuals
(`dev/ordinal-log-bitwise-compare.txt`):

    cumulative                 par TRUE  logLik TRUE  fitted TRUE  osa TRUE
    cumulative_probit          par TRUE  logLik TRUE  fitted TRUE  osa TRUE
    sratio_cs                  par TRUE  logLik TRUE  fitted TRUE  osa TRUE
    cratio_cloglog             par TRUE  logLik TRUE  fitted TRUE  osa TRUE
    acat_cs                    par TRUE  logLik TRUE  fitted TRUE  osa TRUE
    cumulative_thres_gr        par TRUE  logLik TRUE  fitted TRUE  osa TRUE
    hurdle                     par TRUE  logLik TRUE  fitted TRUE  osa TRUE

### Against brms's compiled Stan program

`dev/ordinal-lpcheck.R` runs the suite's own `brms_lp_check()`
(flat priors; check A value, check B brms's gradient at frmtmb's
optimum) at seed 20260930, n = 300, fits to `grad_tol = 1e-6`. Output
`dev/ordinal-log-lpcheck.txt`, summarized by
`dev/ordinal-lpcheck-sum.R` into `dev/ordinal-log-lpcheck-sum.txt`:

    rows run: 31
    identities held (LPCHECK lines): 28
    rows that stopped with an error: 3
    largest |const| / max(1, |logLik|): 2.6e-16
    largest brms gradient at frmtmb's optimum: 0.000229

The 28 cover disc on all four families (cumulative under the logit
and the probit, sratio and acat under the logit, cratio under the
probit), disc beside `cs()`, `link_disc = "softplus"`, equidistant and
sum_to_zero on all four families and the hurdle, both structures and
disc under `thres(gr = )`, and acat under probit, cloglog, cauchit and
probit_approx. The same shapes are rows 12e to 12h of
`test-brms-likelihood.R`. The three errors are brms's (next section).

### brms 2.23.0 defects met (recorded, not filed)

`dev/ordinal-brms-defects.R`, output `dev/ordinal-log-brms-defects.txt`:

1. `cratio("cloglog")` with a modeled disc: brms's `log_prob` at
   frmtmb's optimum equals frmtmb's logLik (-441.1943809182 both), but
   its gradient is NaN in four of six entries, while frmtmb's own
   gradient there is 1.9e-4. At other points (disc slope 0; x slope
   halved) brms's gradient is finite. brms's line is `q[k] =
   log1m_exp(-exp(disc * (mu - thres[k])))`. The Stan identity row uses
   cratio probit instead; cratio cloglog without disc is row 12.
2. `cumulative(threshold = "sum_to_zero")`, logit, disc at 1: brms
   emits `ordered_logistic_glm_lpmf(Y | X, b, 0)`, the literal 0 as the
   thresholds, and stanc refuses it. The probit link and a modeled disc
   take other paths and compile; both are rows.
3. `acat("softit")`: brms's `softit()` helper divides a vector by a
   vector with `/`, which the installed stanc refuses, so no softit
   model compiles. acat softit is checked on brms's R side
   (`brms:::dacat()`) instead.
4. A multivariate or mixture model with equidistant thresholds declares
   `delta` twice (`stan_thres()` suffixes the prior with the group only)
   while the transformed parameters use `delta_y`, `delta_mu1`; stanc
   refuses it (`dev/ordinal-mv.R`, `dev/ordinal-log-brms-code.txt`).
   frmtmb names them as brms's transformed-parameter code does,
   `delta_y` and `delta_y2`.
5. brms's R side reads `probit_approx` as the exact `pnorm()`
   (`brms:::inv_link()`), its Stan program as `Phi_approx`. frmtmb
   follows the likelihood; `test-ordinal-disc-thres.R` checks its
   probit_approx acat against brms's formula with `Phi_approx`.

### Against brms's R-side densities

`test-ordinal-disc-thres.R` compares frmtmb's `fitted()` probabilities,
its simulator's `ord_cat_probs()` and its logLik with
`brms:::dcumulative()`, `dsratio()`, `dcratio()`, `dacat()` at frmtmb's
estimates, for 4 families x 3 structures with `disc ~ 0 + z`, grouped
acat, the hurdle and the five non-logit acat links: probabilities to
`1e3 * .Machine$double.eps`, logLik to 1e-12 relative. The OSA branch
of every density against its data branch, all rows kept
(`dev/ordinal-osa2.R`, `dev/ordinal-log-osa2.txt`): at most 1.43e-15
relative (cumulative equidistant).

### The base build refuses what the lane fits

`dev/ordinal-base-behavior.R` (seed 20260930), outputs
`dev/ordinal-log-base-behavior-{base,lane}.txt`: on the base,
`bf(y ~ x, disc ~ 0 + z)` with `cumulative()` stops "dpar(s) not
available", `hurdle_cumulative(threshold = "equidistant")` and
`"sum_to_zero"` stop "is not implemented", `acat("probit")` stops, and
`brmsfamily("sratio", threshold = "equidistant")` stops "has no
argument"; the lane fits all five. The new tests on the base:
`dev/ordinal-log-t-disc-thres-base.txt` (pass=1 err=10: nine blocks stop
at the missing `threshold` argument and the disc block behaviorally,
"dpar(s) not available") and
`dev/ordinal-log-ts-disc-thres-base.txt` (pass=0 fail=3 err=3: the
hurdle's draws are `tau_raw_k`, and the disc intercept has no default).

### Other

- Multivariate (`dev/ordinal-mv.R`, `dev/ordinal-log-mv.txt`): two
  ordinal responses, equidistant each, disc on one: `variables()`
  gives `delta_y`, `delta_y2`, and the joint logLik is the sum of the
  univariate fits' to 7.54e-12 relative. That is a measurement at two
  separately converged optima, not an identity at a shared point.
- frmtmb.sample (`dev/ordinal-sample-smoke.R`,
  `dev/ordinal-log-sample-smoke.txt`, seed 20260930, sampler seed 3):
  sum-to-zero draws sum to zero to 1.11e-16; equidistant draws carry
  `b_Intercept[1..4]` and `delta`. The round-trip lines in that log
  compare the wrong pair of matrices and are superseded by
  `test-ordinal-disc-thres-draws.R`, which maps the internal matrix
  forward and compares with the stored columns.
- With `disc ~ z` and the thresholds flat (frmtmb.sample's standing
  policy), the smoke's equidistant posterior was wide (thresholds from
  about -60 to 79 in the 95% intervals, `Intercept_disc` -2.3): only the
  `normal(0, 1)` on the disc intercept places the common scale, and
  flat thresholds reward large ones through their volume. brms's
  `student_t(3, 0, 2.5)` on the thresholds is what keeps brms's version
  sensible. `disc ~ 0 + z` is the spelling to use.

## Every suite, one file per process

`dev/ordinal-run-par.sh` over `dev/ordinal-jobs-all.txt` (all 329 test
files of frmtmb and its seven extensions, lane library first, the base
build behind it for the six extensions this lane did not change;
`NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_DRMTMB_FIT_TESTS` and
`FRMTMB_FUZZ` all set), then the 11 files changed or fixed after that
run again (`dev/ordinal-jobs-rerun.txt`). Logs in
`dev/ordinal-suite-all/` and `dev/ordinal-suite-rerun/`, totals from
`dev/ordinal-suite-sum.R` (`dev/ordinal-log-suite-sum.txt`):

    first run: 329 files; rerun: 11 files
                 pkg files  pass fail err skip warn
              frmtmb   194 15993    8   0    0    0
     frmtmb.coupling    11   542    0   0    5    0
          frmtmb.eam    29  1743    0   0    3    0
       frmtmb.latent    10   359    0   0    2    0
        frmtmb.learn    15   500    0   0    2    0
          frmtmb.ode    11   547    0   0    1    0
       frmtmb.sample    44  2443    1   0    1    0
       frmtmb.spline    15   553    0   0    1    0

The 9 failures are the 9 ported rows below, failing as "now HOLDS"
against their recorded verdicts. Every skip is a `test-scale.R`.
Four files failed on the first run because of this change and were
fixed in the tests, not the package: `test-mv-gaps.R` and
`test-ordinal-draws-names.R` hand-build a parameter vector from the
template and now drop the mapped disc coefficient;
`test-thres-refit.R` expects `fixef(flatten = TRUE)` to carry
`disc_(Intercept)`; `test-xbeta-zibb-hurdle-cum.R` asserted the old
refusal of `threshold = "equidistant"`.

## R CMD check --as-cran

`dev/ordinal-check.ps1`, built with vignettes and the manual inside
`dev/ordinal-check/<pkg>/`, `R_LIBS` = lane library, base build, user
library, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`:

- frmtmb: `Status: 1 NOTE`, the environmental "Skipping checking math
  rendering: package 'V8' unavailable" on the HTML manual. The tarball
  predates one comment-only edit (the log-space note moved into
  `ord_cumulative_logpmf()`) and the helper edit below; neither is
  code a check runs.
- frmtmb.sample: first run `Status: 1 WARNING`, "'::' or ':::'
  import not declared from: 'Matrix'": the `Z` view in the shared
  brms-suite helper called `Matrix::rowSums()`, and frmtmb.sample runs
  a copy of that helper without declaring Matrix. The view now uses a
  matrix product; the second run is `Status: OK`, and the six standata
  rows still hold after the change (`dev/ordinal-log-g-test-brms-suite-standata.R.txt`).

## Ported brms-suite rows that now hold

Gated runs on the lane build (`dev/ordinal-log-g-*.txt`, and the full
suite above). The verdict file is not edited; these rows fail as stale
"now HOLDS" until consolidation:

| row | package | recorded | now |
|---|---|---|---|
| priors:14 | frmtmb | cannot transfer | holds |
| standata:724 | frmtmb | cannot transfer | holds (constant dpar as a data slot) |
| standata:725 | frmtmb | cannot transfer | holds |
| standata:738 | frmtmb | cannot transfer | holds |
| standata:739 | frmtmb | cannot transfer | holds (Z view) |
| standata:740 | frmtmb | cannot transfer | holds (Z view present, so not hollow) |
| standata:745 | frmtmb | cannot transfer | holds (Z view) |
| families:37 | frmtmb | divergence | holds |
| families:37 | frmtmb.sample | divergence | holds |

standata:724 and :745 are outside this lane's families: they flip
because the standata view now carries constant dpars and brms's `Z`
columns, which :725, :739 and :740 needed.

## Decided, and why

- **sum_to_zero estimates `K - 2` directions**, not brms's `K - 1`
  parameters with a flat direction: under maximum likelihood a flat
  direction is a singular Hessian. The reported thresholds are the
  same function of the likelihood's parameters as brms's
  `b_Intercept`. The price is brms's class `Intercept` prior, which is
  about the uncentered vector; frmtmb has no such parameter and refuses
  the class there, naming `prior = list(tau_raw = )` for the free
  directions. A prior on the centered thresholds would be a different
  prior under the same name.
- **Two-threshold minimum, for equidistant only.** Equidistant on a
  single threshold is refused: `delta` has nothing to measure. WITHDRAWN
  in punch round 1 for sum_to_zero: round 0 also refused a response
  whose every threshold vector has one threshold, and said brms refuses
  it, which is false under `thres(gr = )` (brms declares the grouped
  count `int<lower=1>`). A single zero-sum threshold is 0, a vector with
  no parameter; it is fitted now, grouped and ungrouped (see Punch
  round 1, B3).
- **`Fixed dpar: disc = 1` printed on every ordinal summary** in round
  0, with `disc = log` on the `Links:` line and `disc_(Intercept)` in
  `fixef(flatten = TRUE)`. WITHDRAWN in punch round 1: the project rule
  is to match brms, which shows none of it, so a disc held at its
  family default is now hidden everywhere it showed, the hurdle's
  included (see Punch round 1, B2).
- **The disc constant is a mapped coefficient** (the existing
  `fixed_dpars` route), so every ordinal fit's `betad` gains one mapped
  entry. `opt$par` is unchanged; code that hand-builds a parameter
  vector from `parList()` must drop `frame$betad_fixed_idx`, which is
  what `test-mv-gaps.R` needed.

## Not done, and why

- **Ordinal mixtures** stay refused (unchanged): brms allows them with
  `disc1`, `disc2` and fixed thresholds across components, and its
  equidistant mixture does not compile (defect 4).
- **`thres(gr = )` and `cs()` for `hurdle_cumulative()`** stay refused
  (unchanged); brms fits both.
- **`posterior_linpred(incl_thres = TRUE)`** in frmtmb.sample stays
  refused; brms returns `disc * (thres - mu)` per threshold for any
  ordinal family. Its refusal text said brms supports it "for
  cumulative families alone", which brms 2.23.0 does not
  (`posterior_epred.brmsprep()` tests `is_ordinal()`); the text and the
  argument's documentation now say what brms does.
- **The scale tier** (`test-scale.R` of the extensions,
  `FRMTMB_SCALE_TESTS`) was not run: it times the extensions' own
  objectives, which no ordinal family reaches. Those files report
  skips in the run below.
- **A fixed disc in `variables()`.** brms lists `disc` (a transformed
  parameter, 1) on every ordinal fit; frmtmb lists no fixed dpar
  anywhere (pre-existing, all families).
- **Class `"Intercept"` with `coef` on ordinal thresholds** stays
  unaddressable (pre-existing); brms lists per-threshold rows under
  flexible and sum_to_zero.
- **`confint()`** names ordinal internal parameters `tau_raw_k`
  (pre-existing); under equidistant `tau_raw_2` is `log(delta)`.
- **The compat table** (`R/compat.R`) gained no rows for the
  structures or disc; frmtmb.eam's tests read that table.
