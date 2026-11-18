# Lane squant: a smoothed quantile family

Round of 2026-11-16, base b255350e (frmtmb 0.83.0), branch
`wt-r84-squant`. Brief: `round-briefs/r084-squant.md`. Backlog items
5449 (HIGH) and 4982 (MEDIUM) are the background; the user decided on
2026-11-13 that a NEW family with curvature is the answer and that
`asym_laplace()` stays as brms has it. This file is the record; logs
are local (`dev/squant-out/`, `dev/squant-log/`,
`dev/release/squant-*`, gitignored). Library:
`C:/Users/adf44/source/r/wt-squant-lib` (core only; the extensions are
read from rellib-r22). qgam 2.0.0, quantreg 6.1 and mgcv 1.9-5 were in
the user library; nothing was installed but frmtmb.

## Summary

- New family `quantile_elf(link = "identity", link_sigma = "log",
  link_quantile = "logit", smoothing = 0.25)` in a new file
  `R/families-quantile.R`. It is the extended log-F (ELF) density of
  Fasiolo et al. (2021), the density qgam fits. Parameters `mu`,
  `sigma`, `quantile` as `asym_laplace()`, and a FIXED smoothing
  `lambda` relative to `sigma`.
- It lives in CORE (argument below; the rule decides it).
- Validation (a) to (e) below, every number from a log with its script.
- Hessian SEs exist but are model-based; the sandwich
  (`vcov_cluster(fit, ~ row, "CR0")`) is the calibrated one. The
  parametric `frm_bootstrap()` is NOT calibrated (for this family and,
  measured, for `asym_laplace()` too: filed).
- With `(1 | g)` the Laplace objective has a maximum on every design
  of 5449 (all fits code 0, no warning); its accuracy needs a larger
  smoothing (`smoothing = 1`) or quadrature, which fitted all 200 such
  models (5449's quadrature broke down on 24 of 40 for asym_laplace).
- Bump: MINOR (a new exported family). Floors: none change (no
  extension calls anything new; the extensions' tests read the family
  registry and pass with this core, see "Suites").

## Core or extension

Core. Rewritten in fix round 1 on the points of `dev/core-boundary.md`
that decide it (the review, item 8, found that the first version
leaned on "would most users meet it", which is a stretch for quantile
regression):

1. "What core holds" lists "the general families: the GLM families,
   the distributional families, ordinal families, mixtures, ...". The
   ELF is a distributional family with a location, a scale and a shape
   (`quantile`), of no domain, and its unsmoothed sibling
   `asym_laplace()` is a core family. "How to decide" question 3 (a
   model class or workflow for one domain) is "no".
2. Guardrail 2 ("interface before move"): the family reads core
   internals that no seam exports: `log_ibeta_two()` (the tape-safe
   incomplete beta pair), `lbeta_ad()`, `log1m_exp_ad()`, `ad_min()`,
   `ad_pos()`, `sign_pos()`, `cdf_rows()`, `dpar_link()`, `mu_link()`,
   `link_arg_value()` and `frm_stop()` (`dpar_complement()` is
   exported). In an extension it would first need those as new public
   seams, and no extension needs them now.
3. Guardrail 3 ("measure before move"): 435 lines in
   `R/families-quantile.R` (116 of code, the rest roxygen and
   comments), 353 lines of tests that run in about 4 s, no dependency of
   its own (the tests use quantreg, already in Suggests; qgam only in
   `dev/squant-qgam.R`), no ported brms row. By the guardrail's own
   words that is a weak candidate for an extension: small, stable, no
   dependency.

The rule decides it, so no user question.

## The density, and how the smoothing is chosen

Read from qgam 2.0.0's source (`qgam::elf`, `qgam::qgam`,
`qgam:::.getErrParam`; dumped by `Rscript -e` into
`dev/squant-out/qgam-src1.txt`, `qgam-src2.txt`):

- qgam's `elf()` `aic()` is `-2 sum log f` with
  `log f = (1 - tau) r / sig - lam log1pexp(r / lam) / sig -
  log(lam beta(lam (1 - tau) / sig, tau lam / sig))`, `lam = co` and
  `sig = exp(theta) * lam / mean(lam)`.
- `qgam()` sets `co <- err * sqrt(2 pi varHat) / (2 log 2)`, with
  `varHat` from a gaussian pilot fit and `err` from `.getErrParam()`:
  a SHASH fit to the pilot residuals, then
  `h = (d 9 / (n pi^4))^(1/3) exp(lf0 / 3 - 2 lf1 / 3)` (an MSE-optimal
  kernel width, rate `n^(-1/3)`, `d` the model df), `err = min(h 2 log
  2 / sqrt(2 pi), 1)`. Near the residual mode it moves the quantile by
  0.05 because `f'` is 0 there.
- The learning rate `sigma = exp(lsig)` is chosen by `tuneLearnFast()`
  to calibrate posterior credible intervals (or by pinball loss), NOT
  by likelihood.

In the paper's notation the absolute kernel width is `lambda sigma`,
and `lambda` is relative to `sigma`. `quantile_elf()` takes the
paper's `lambda` as `smoothing`, FIXED, and estimates `sigma` by
maximum likelihood (it is a proper density's scale). With `lambda`
fixed the family is a location-scale family: `y = mu + lambda sigma
logit(Z)`, `Z ~ Beta(lambda (1 - p), lambda p)`. That gives exact CDF
(incomplete beta), quantile (qbeta), draws (two log-gammas), mean
`mu + lambda sigma (digamma(a) - digamma(b))` and variance. As `lambda
-> 0` the density is `asym_laplace()`'s with the same `mu`, `sigma`,
`quantile` (normalizer `lambda sigma B(a, b) -> sigma / (p (1 - p))`).
A predicted `sigma` scales the width row by row, which is what qgam's
`elflss` does (`sig <- sig * lam / mean(lam)` keeps `lam / sig`
constant).

Estimating `lambda` is not sound for quantile regression: the score in
`mu` sets `mean(plogis(r / (lambda sigma))) = 1 - p`, so `mu` is the
`p` quantile of `y - lambda sigma L` (L standard logistic); ML would
pick the `lambda` that fits the residual shape, and a large `lambda`
moves `mu` away from the quantile (bias about `-(pi^2 / 6) (lambda
sigma)^2 f'(q) / f(q)`, table below). qgam does not estimate it either.

The default 0.25 is fixed (memory rule "no special-case defaults"; a
data rule like qgam's was considered and not built, see "Not done").
It was chosen from the tables below: at 400 rows it keeps the sandwich
coverage at 0.89 to 0.95 on six designs, while 1.0 gives bias up to
1.7 sd (lognormal median intercept) and 0.1 a noisier kernel.

## (a) bias against the quantile, and rq() as the smoothing goes to 0

Population (n -> Inf) fit by numerical integration, no randomness,
`dev/squant-pop.R`, log `dev/squant-out/pop.txt`. Bias of `mu*` in
residual sd, and in sandwich SEs at n = 1000 (selected rows; the log has
84):

| residuals | quantile | h 0.1 | h 0.25 | h 0.5 | h 1 | h 1, SE at n 1000 | h 0.25, SE at n 1000 |
|---|---|---|---|---|---|---|---|
| gaussian | 0.5 | 0 | 0 | 0 | 0 | 0 | 0 |
| gaussian | 0.75 | 0.00111 | 0.00679 | 0.02511 | 0.08025 | 2.122 | 0.167 |
| gaussian | 0.9 | 0.00065 | 0.00402 | 0.01569 | 0.05740 | 1.176 | 0.077 |
| t3 | 0.75 | 0.00171 | 0.01019 | 0.03562 | 0.10313 | 3.401 | 0.351 |
| t3 | 0.9 | 0.00093 | 0.00576 | 0.02252 | 0.08174 | 1.578 | 0.112 |
| lognormal | 0.5 | 0.00239 | 0.01374 | 0.04000 | 0.08271 | 4.504 | 0.811 |
| lognormal | 0.9 | 0.00110 | 0.00685 | 0.02746 | 0.10887 | 1.236 | 0.078 |

Finite samples against `quantreg::rq()` (same tau): the test
"a fit reaches the quantile regression as the smoothing shrinks"
(seed 41, 300 rows, sd 0.5 + x, tau 0.75): max |coef - rq| / rq's nid
SE is 0.986 at h = 1, 0.194 at 0.25, 0.018 at 0.02; the same on all 21
ulp perturbations of y (`dev/squant-ulps.R`, `ulps.txt`, "(b) 21 of
21"). Monte Carlo bias of the slope in MC sd is in the SE table
(column bias1/sd of `se-summary.txt`); at h = 0.25 it is within 0.09
of rq's on every design.

## (b) against qgam on its own examples, same smoothing

`dev/squant-qgam.R`, log `dev/squant-out/qgam.txt`. qgam called with
`lsig` given (no tuning); `quantile_elf(smoothing = co / exp(lsig))`
with `sigma` held at `exp(lsig)`, which is the same density.

Parametric (`gamSim(1)`, n 400, `y ~ x0 + x1 + x2 + x3`):

| qu | lsig | smoothing | logLik qgam - frm | max rel coef diff | max score qgam | max score frm |
|---|---|---|---|---|---|---|
| 0.1 | -2 | 3.65720 | 6.4e-12 | 2.24e-07 | 1.2e-07 | 4e-05 |
| 0.1 | 0 | 0.49495 | 9.3e-12 | 1.52e-06 | 1.7e-08 | 5.7e-06 |
| 0.5 | -2 | 10.71550 | 5.1e-09 | 2.74e-06 | 1.6e-08 | 9.1e-04 |
| 0.5 | 0 | 1.45019 | 3.3e-12 | 4.46e-07 | 2.1e-09 | 7.6e-06 |
| 0.8 | -2 | 3.58772 | 5.2e-09 | 2.56e-05 | 7.2e-10 | 7.2e-04 |
| 0.8 | 0 | 0.48555 | 0 | 3.47e-08 | 9.6e-11 | 5.2e-08 |

The log-likelihoods agree to 5e-9 (2.6e-12 relative), so the models
are the same; the coefficients agree to 1e-6 relative on 4 of 6 and to
2.6e-5 on the other two, where frm's stop has a score of 7e-4 to 9e-4
(qgam's 1e-8): the difference is frm's stopping point, not the model.
`frmtmb_control(optCtrl = list(rel.tol = 1e-14))` returned the same
estimates bit for bit with code 1 (filed LOW below).

qgam's examples with smooths (`?qgam`): `gamSim(1)` with
`s(x0) + s(x1) + s(x2) + s(x3)` at qu 0.5, and `MASS::mcycle` with
`s(times, k = 20)` at qu 0.8 (`bs = "ad"` of the example has no frmtmb
equivalent, so the default basis). qgam (mgcv REML/LAML) against
`frm(..., REML = TRUE)`, sigma held: max |mu difference| relative to
the range of mu 3.5e-6 (gamSim, edf 14.3) and 1.9e-8 (mcycle, edf
17.9); both frm fits code 0, no warning.

## (c) Hessian and sandwich SEs against the Monte Carlo spread

`dev/squant-se.R <design> <nrep> <n>` (seeds 20261116 + rep), driver
`dev/squant-se-run.sh`, summary `dev/squant-se-sum.R`, log
`dev/squant-out/se-summary.txt`. Three designs (`homo`: gaussian errors;
`hetero`: sd 0.5 + x; `skew`: lognormal errors), 300 data sets each at
n 400, plus homo at n 100 (300) and n 2000 (122; stopped early for
time). Slope: ratio of mean SE to MC sd, and coverage of 95% Wald
intervals for the TRUE conditional quantile slope. All 10571 ELF fits
stopped at code 0 with no warning (counted from the tsv files).

| design | n | tau | h | MC sd | Hessian/sd | sandwich/sd | cover Hess | cover sand | rq nid cover |
|---|---|---|---|---|---|---|---|---|---|
| homo | 400 | 0.5 | 0.25 | 0.2085 | 0.829 | 0.977 | 0.893 | 0.943 | 0.947 |
| homo | 400 | 0.9 | 0.25 | 0.2871 | 0.592 | 1.002 | 0.777 | 0.913 | 0.947 |
| hetero | 400 | 0.5 | 0.25 | 0.2074 | 0.808 | 0.966 | 0.887 | 0.933 | 0.937 |
| hetero | 400 | 0.9 | 0.25 | 0.2773 | 0.591 | 0.995 | 0.767 | 0.907 | 0.953 |
| skew | 400 | 0.5 | 0.25 | 0.1999 | 1.028 | 0.985 | 0.953 | 0.950 | 0.940 |
| skew | 400 | 0.9 | 0.25 | 1.0451 | 0.505 | 1.012 | 0.700 | 0.913 | 0.957 |
| homo | 100 | 0.9 | 0.25 | 0.5552 | 0.595 | 1.123 | 0.707 | 0.837 | 0.963 |
| homo | 100 | 0.9 | 1 | 0.4970 | 0.692 | 1.052 | 0.827 | 0.910 | 0.963 |
| homo | 2000 | 0.9 | 0.25 | 0.1370 | 0.566 | 0.939 | 0.760 | 0.950 | 0.950 |

The ELF estimate is more efficient than rq's (MC sd 0.2085 against
0.2263 at homo 400 median). The intercept at h = 1 is biased (skew
median: 1.74 MC sd, coverage 0.64; homo 2000 tau 0.9: 0.81 sd,
coverage 0.88), which is the population bias above.

Parametric bootstrap (`dev/squant-boot.R`, `boot.txt`, homo 400 tau 0.9,
seeds 20261117 to 20261119, 200 draws): slope SE 0.1112, 0.1249,
0.1198 against Hessian 0.1693, 0.2004, 0.1835, sandwich 0.2867,
0.3668, 0.3371 and MC sd 0.2871. It draws from the fitted ELF, so it
is model-based. The same on `asym_laplace()` with rellib-r22
(`dev/squant-boot-ald.R`, `boot-ald.txt`): slope SE 0.1062, 0.1254,
0.1116 (38, 47, 39 of 200 refits converged) against rq's MC sd 0.2977;
`?asym_laplace` ("Quantile regression inference") recommends
`frm_bootstrap()` "for intervals you can defend". Filed (MEDIUM).

## (d) with (1 | g): the Laplace fit has a maximum

`dev/squant-re.R <design> <quantile>`, 5449's designs and seeds
(`dev/reviews/wcheck2-review-al.R`'s `mk_data`, seeds 20261114 to
20261123), logs `dev/squant-out/re-*.txt`. For each fit: Laplace fit,
the Laplace objective profiled in log sd at offsets -3 to +8 with every
other parameter re-optimized, and `quadrature = TRUE`.

| design | q | h | conv 0 | no warning | fit at profile max | quad fits | max abs log sd - quad | logLik - quad, median |
|---|---|---|---|---|---|---|---|---|
| g4 | 0.5 | 0.1 | 10 | 10 | 9 | 10 | 0.5174 | -1.2364 |
| g4 | 0.5 | 0.25 | 10 | 10 | 10 | 10 | 0.2969 | -0.9050 |
| g4 | 0.5 | 1 | 10 | 10 | 10 | 10 | 0.0347 | -0.2542 |
| g10 | 0.5 | 0.25 | 10 | 10 | 10 | 10 | 0.1357 | -0.1043 |
| g10 | 0.5 | 1 | 10 | 10 | 10 | 10 | 0.0091 | -0.0660 |
| g2 | 0.5 | 0.25 | 10 | 10 | 9 | 10 | 0.5908 | -2.7455 |
| g2 | 0.5 | 1 | 10 | 10 | 10 | 10 | 0.1818 | -0.9021 |
| al4 | 0.5 | 0.1 | 10 | 10 | 7 | 10 | 0.7312 | -0.7648 |
| al4 | 0.5 | 0.25 | 10 | 10 | 10 | 10 | 0.5465 | -0.5105 |
| al4 | 0.5 | 1 | 10 | 10 | 10 | 10 | 0.0809 | -0.1417 |
| g4 | 0.9 | 0.25 | 10 | 10 | 9 | 10 | 1.2352 | -4.8403 |
| g4 | 0.9 | 1 | 10 | 10 | 10 | 10 | 0.2793 | -1.9403 |

(All 20 SUMMARY lines are in the logs.) Every profile falls at both
ends (0.14 to 24.3 below at -3 and 52.6 to 289 below at +8, 200 fits):
the maximum exists, against asym_laplace's profile that rose in 30 of
30 (5449). But the fit is not always at the profile's maximum: 7, 2, 1
and 0 of 50 fits (h 0.1, 0.25, 0.5, 1) stop at code 0, with no warning,
below a higher profile point, by up to 1.76 at h 0.25 (g4 q 0.9 seed
10, where quadrature's logLik is 11.05 higher). Filed (MEDIUM) and
documented in `?quantile_elf` ("Group-level terms"): use smoothing 1
and quadrature. Quadrature fitted all 200 models; for asym_laplace on
the same designs "the Gauss-Kronrod objective broke down" on 24 of 40
(5449). `frm_sample()` was not run: quadrature is exact for one scalar
grouping factor, which these are.

The test "with (1 | g) in mu the Laplace approximation has a maximum"
(g4, seed 20261114, h 0.5) holds on 21 of 21 ulp perturbations,
profile 3.945 and 21.632 below at -2 and +2 (`ulps.txt`).

## (e) d/p/q/r consistency

- `dev/squant-cdf.R` (`cdf.txt`): the log CDF pair against
  `stats::pbeta(log.p = TRUE)` (each tail from the side below one half)
  over smoothing 0.02 to 20 and quantile 0.01 to 0.99 (shapes down to
  2e-4), `t` from -680 to 680: worst relative error 4.13e-14.
  `exp(lF) + exp(lS) - 1` at most 6.7e-15. Past |t| = 690 each tail is
  continued by its exact slope and equals the leading term of the
  incomplete beta to the last digit at t = -5000 and 5000. The taped
  pair equals the doubles' (-10.712557094055816 both) and its gradient
  the finite difference to 10 digits.
- `tests/testthat/test-families-quantile.R`: density against the
  logit-beta law built from `dbeta()`; integrates to 1; limit to
  `asym_laplace()`; CDF pair against `pbeta()`; far tails; quantile
  inverts the CDF; tape value, score and Hessian at the kink; fitted()
  against the integrated mean; osa residuals against `qnorm(pbeta())`;
  `cens()` likelihood against the hand-rolled one; the estimating
  equation at the fit; the weights refusal.
- `tests/testthat/test-simulate-density.R`: a `quantile_elf` row in
  the agreement tier (simulator against the family's own density),
  which its "every registry family is in the tier" test requires.

## Other changes

- `R/weights-unbounded.R`: `quantile_elf = "sigma"` in
  `unbounded_spike_dpar` (own line). Measured first
  (`dev/squant-weights.R`, `weights.txt`): with `(1 | obs)` and weights
  1.2 the Laplace objective rises without end in -log sigma (seed 1:
  -95.90 at the fit, +30.59 at -20) and the fit stopped at code 0 with
  no warning; at 0.9 it falls (-73.39 to -408.34). After
  (`weights-after.txt`): 1.2 refused with "The weighted likelihood of
  this model has no maximum", 0.9 unchanged, three seeds each.
- `R/start-ls.R`: a `quantile_elf` rule (own entry), the same as
  asym_laplace's (`laplace_sigma`).
- `R/links-brms.R` regenerated by `dev/famlink-gen-links.R`, which
  gains `quantile_elf = "asym_laplace"` as an analog (like huber's
  gaussian). The regenerated file differs from the base only in those
  lines and the header comment.
- `R/families.R`: one registry line (`quantile_elf = fam_quantile_elf`,
  after a comment line), nothing else.
- The `lcdf` slot: the family fills `cdf_log` (the log pair), from
  which `frmtmb_family()` derives `lcdf` and `lccdf`; no line of this
  lane names `lcdf`, so lane cdfname's rename script has nothing of
  this lane's to rename.

## The asym_laplace refusal's alternative (lane alapse)

Lane alapse's tree (read, not edited) has `ald_se_routes` in its new
`R/ald-fixed.R` ("A family with curvature that estimates the same
quantile goes first when there is one") and `ald_laplace_message()` in
`R/ald-laplace.R`. For the consolidation, the suggested first route:

    "quantile_elf(), a smoothed asymmetric Laplace with curvature
    (?quantile_elf), with vcov_cluster(fit, ~ row, \"CR0\") for the
    quantile estimate's standard errors"

and, in `ald_laplace_message()`, before frm_sample(): "quantile_elf()
fits the same quantile with a log density that has curvature, so its
Laplace approximation has a maximum (use smoothing = 1; see
?quantile_elf)." Note also that alapse's first route today is
`frm_bootstrap()`, which this lane measured as not calibrated for
asym_laplace (above).

## Not done, and why

- A data-driven default for the smoothing (qgam's err rule, rate
  `n^(-1/3)`): it needs a pilot fit and a residual density fit, and a
  default that switches with the data; the fixed 0.25 is documented
  with its bias and coverage tables.
- qgam's calibration of sigma: sigma is the density's ML scale here,
  and the sandwich gives calibrated SEs instead.
- `predictor_unit_families` (R/autoscale.R,
  lane slopetrig's file) and gp-start's family list: like huber, the
  family is not in them; filed LOW. (Compat rows: added in fix round 1.)
- No export of d/p/q/r functions (frmtmb exports none for any family);
  `elf_quantile()` is internal and tested.
- `frm_sample()` was not run on the family; quadrature was the
  reference for (d).

## Suites

One file per process, two at a time, NOT_CRAN=true, every log's `lib:`
line the lane core (`wt-squant-lib`; 550 of 550), extensions from
rellib-r22: `bash dev/release/run-par.sh squant-all dev/squant-all.jobs 2`
with `FRMTMB_RELLIB` = the lane library, log
`dev/release/squant-all.log` ("SQUANT-ALL ran 550 of 550"):

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb | 350 | 25788 | 0 | 0 | 243 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1851 | 0 | 0 | 5 | 0 |
| frmtmb.latent | 13 | 700 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 474 | 0 | 0 | 13 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 92 | 4286 | 0 | 0 | 18 | 0 |
| frmtmb.spatial | 8 | 551 | 0 | 0 | 6 | 0 |
| frmtmb.survival | 16 | 575 | 0 | 0 | 0 | 0 |

New and changed test files: `test-families-quantile.R` pass=109 fail=0
err=0 skip=0 warn=0; `test-simulate-density.R` pass=630 (fail 0);
`test-brms-families.R` pass=1280 (the new family's 5 links construct and
the roster's other links are refused).

Seen failing on rellib-r22 (`dev/squant-out/base-test.txt`):
`test-families-quantile.R` pass=0 err=15, each "could not find function
"quantile_elf"" (the weak form; there is no older behavior to pin: the
family is new). The behavior the family exists for is the contrast with
`asym_laplace()` on 5449's design, refused by `check_ald_laplace()` on
rellib-r22 and fitted here with an interior Laplace maximum (section d).

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`):
`test-families-quantile.R` pass=109 fail=0 err=0
(`openblas-test.txt`); `test-simulate-density.R` pass=630 fail=0
(`openblas-simdens.txt`); `dev/squant-ulps.R` "HELD of 21: (a) 21 (b) 21
(c) 21" (`ulps-openblas.txt`), as on reference BLAS (`ulps.txt`).

The contrast on the test's own data (`dev/squant-ald-contrast.R`,
`dev/squant-out/ald-contrast-*.txt`): `asym_laplace()` with `(1 | g)`
is refused on rellib-r22 and on the lane build ("frm() cannot fit
`(1 | g)` under the Laplace approximation ..."); `quantile_elf(smoothing
= 0.5)` fits, code 0, logLik -127.2105.

Single-session runs (`bash dev/release/run-sessions.sh squant-sess 2`,
`FRMTMB_RELLIB` = the lane library, log `dev/release/squant-sess.log`,
"SQUANT-SESS ran 27 of 27"): every session fail=0 err=0 warn=0; no
untagged LEAK line. frmtmb forward/reverse/preload pass 25788 / 25953 /
26057 (files 337 of 350: the 13 others are gated brms files with no
test outside the gate, which read pass=0 in the one-file run too);
frmtmb.sample 4286 in each order (87 of 92 files, likewise); every
other extension the same count in all three orders as in the one-file
run.

Gated: `test-brms-suite-families.R` with `FRMTMB_BRMS_FIT_TESTS=true`
pass=84 fail=0 err=0 skip=0 (`gated-brms-suite-families.txt`).

`test-simulate-density.R`'s header count regenerated by
`dev/rel071-simcount.R` (`simcount.txt`): 389 statistical assertions
(115 chi-square, 127 mean, 127 variance, 18 P(Y = k), 2 correlation),
P(any false alarm) 0.0004 at 1e-6, as before.

R CMD check --as-cran of core (`bash dev/squant-check.sh`, built with
vignettes, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, R_LIBS = lane library,
rellib-r22, user library; log kept as `dev/squant-out/check.log`, the
check tree deleted): "Status: 1 NOTE", the expected V8 math-rendering
NOTE on the HTML manual; tests 47 min, OK. No extension changed, so no
extension was checked.

## Files touched

New: `R/families-quantile.R`, `man/quantile_elf.Rd`,
`tests/testthat/test-families-quantile.R`, and the lane's `dev/squant-*`
scripts and this file. Changed: `R/families.R` (one registry entry and
its comment), `R/start-ls.R` (one rule), `R/weights-unbounded.R` (one
entry and its comment), `R/links-brms.R` (regenerated),
`dev/famlink-gen-links.R` (the analog), `NAMESPACE` (roxygen),
`NEWS.md`, `_pkgdown.yml` (one reference entry),
`tests/testthat/test-simulate-density.R` (one tier row, header count),
`dev/test-backlog.md` (section "Filed in the round of 2026-11-16", four
items, updated in fix round 1 with a fifth). Fix round 1 also changed `R/compat.R` (the rows listed under "7. frm_compat() rows"). No ledger row changed.

## Fix round 1 (review `dev/reviews/2026-11-17-squant.md`, verdict FIX)

### 1. The taped derivative of log F (MEDIUM, fixed)

The review found `d lF / d mu` of `elf_cdf_log()` off by about
`eps / F` for `-690 < t < 0`, and 0, doubled or garbage once `F <
1e-16`. Cause: the continuation past `|t| = 690` was added as a
difference `log1m_exp(ls_high) - log1m_exp(ls0)`, 0 in value inside
but with slope `S / F` in each term; the tape's two copies did not
cancel at that size. Fix (`R/families-quantile.R`): the continued
complement is switched in by `sign_pos()` weights, whose derivative is
0, and the clamp and the slopes use a new `elf_pos0(u) = u +
ad_pos(-u)`, which is `max(u, 0)` with derivative 0 at `u = 0` (RTMB's
`abs()` has derivative 1 there), so at `|t| = 690` itself the
derivatives are those of the pair read inside. Values are unchanged
(`dev/squant-cdf.R` after: worst relative error 4.13e-14, the far tails
equal the leading term to the last digit, `cdf-after.txt`).

`dev/squant-tapegrad.R` (the new test's grid: t from -700 to 700 with
+-689, +-690, p 0.1, 0.5, 0.9, smoothing 0.25 and 1, d/dmu of lF and lS
against -f/F and f/S), on the build before the fix
(`tapegrad-before.txt`): worst relative error 2.9e+67 to 5.1e+269 at t
= -690 in every setting, and 1 (a wrong derivative) from t = -200 to
-500 at p 0.1, smoothing 0.25. After (`tapegrad-after.txt`,
`tapegrad2.txt`): max |tape - analytic| / max(p, 1 - p) is 4.9e-15
over every setting and side. The relative error is below 1e-8 wherever
the derivative is above eps of that slope; where it is not (a log tail
within 1e-31 of 0, e.g. lS at t = -689 at p 0.1, |f/S| = 3.8e-69 of the
slope) the tape gives 0. The new test asserts both. Under OpenBLAS
0.3.26 with LAPACK the same (`tapegrad2-openblas.txt`, 4.9e-15).

The test "the taped log CDF pair has the derivatives -f / F and f / S"
fails on the pre-fix build by the same computation (the before log;
its relative errors of 1 at t = -200 are on derivatives far above the
eps bound).

The review's left-censored fit (`dev/reviews/squant-review-tape4.R`,
rerun into `dev/squant-out/review-tape4-after.txt`, rows censored at 60
in data near 100, p 0.1): max |tape gradient| against the finite
difference was 1.39e-05 against 0.0185 (smoothing 0.25) and 9.32e-05
against 0.0186 (smoothing 1) before; after, 1.33e-05 against 1.28e-05
and 4.23e-06 against 4.09e-06. Every other case of that script agrees
too (e.g. L 80, p 0.1, smoothing 0.25: 6.72e-06 against 6.76e-06,
before 6.37e-05 against 3.29e-04).

### 2 and 3. Standard errors in the docs

`?quantile_elf`, "Standard errors": one sentence that the default
`summary()`, `fixef()`, `confint()` and `hypothesis()` show the ELF
model's Hessian standard errors, which cover about 0.7 to 0.8 away
from the median; the steps to report the sandwich
(`V <- vcov_cluster(fit, ~ row, "CR0", full = TRUE)`,
`summary(fit, vcov = V)`, `confint(fit, vcov = V)`, and
`vcov(fit, cluster = ~ row)`); the example does all four. "Group-level
terms": which standard error each route has. The Laplace fit takes
`vcov_cluster(fit, ~ g, "CR0", full = TRUE)` (robust; coverage not
measured for this family); the quadrature fit is refused by
`vcov_cluster()` (`dev/reviews/squant-review-re.txt`), so only its ELF
Hessian standard errors exist, not calibrated. The default is not
changed (memory rule: no special-case defaults).

### 4. The (1 | g) MEDIUM

Backlog item and "Group-level terms" rewritten: the profile at a small
smoothing is a sawtooth with several local maxima (the log determinant
jumps as a group's mode crosses a data point); the counts 7, 2, 1, 0 of
50 are lower bounds from a 10-point profile; on the review's case the
finer profile is 1.764 higher at smoothing 0.25 and 3.420 at 0.1. The
doc's "by up to 1.76" is replaced.

### 5. elf_quantile() (LOW, fixed)

`qbeta()` is now only the start; Newton steps on `log F` (or `log S`
above one half), both concave in `t` because the density is
log-concave, solve to 1e-13. The review's cases are in the test: 229.34003,
68.313613, 68.309686 at (p 0.01; smoothing 0.25, 0.05, 0.001; pr 0.9,
0.5, 0.5), and the p 0.99 mirror images, to 1e-7 (pass).

### 6. The bias formula (LOW, docs)

"What mu estimates" now says the formula holds for a kernel width small
against the scale on which `f` changes, can be off by a factor of 2 at
smoothing 1, and can have the wrong sign near an end of the support;
the table is from integration, not from the formula.

### 7. frm_compat() rows (LOW, done)

`quantile_elf` is in the compat vocabulary. Every row changed, each on
its own line of `R/compat.R`:

- `compat_features_build()`: `"quantile_elf"` in `fams`.
- groups `cdf`, `cdf_continuous`, `mi_families`: `"quantile_elf"`.
  So `trunc()` and `cens()` read "works" from `group:cdf` /
  `group:cdf_continuous`, and `mi()` "works" from `group:mi_families`.
- the loop of continuous families with no truncated mean: `trunc()`
  x `quantile_elf` is "conditional" ("fitted() and residuals(type =
  "response") refuse a truncated response"), as measured:
  `fitted()` on a `trunc(lb = 0)` fit stops with "Family 'quantile_elf'
  has no truncated mean" (`dev/squant-compat.R`, `compat.txt`).
- `r("nl", "quantile_elf", "works", ...)`: an exponential body fits,
  code 0, no warning (`dev/squant-compat3.R`, `compat3.txt`).
- `r("mixture", "quantile_elf", "works", ...)`: two components with
  the quantiles held, code 0, no warning, locations -1.965 and 1.934
  for -2 and 2 (`compat3.txt`).

With all eight extensions loaded (`dev/squant-compat2.R`,
`compat2.txt`) `frm_compat("quantile_elf")` has 104 rows, the same
features as `asym_laplace`'s: 38 works, 40 conditional, 20 refused, 6
untested. Rows that differ from `asym_laplace`'s, and why: `us`,
`diag`, `homdiag`, `cs`, `s()`, `t2()`, `gp_pred()` and `mi()` are
"works" (asym_laplace's "conditional" rows are its Laplace refusal);
`ar()`, `ma()`, `arma()`, `cosy()`, `unstr()` are "conditional" through
the latent-residual rule (asym_laplace refuses them by name); an `ar()`
fit runs at code 0 with the latent residual's own bias warning
(`compat3.txt`); `nl` is "works" (asym_laplace "conditional").
`frm_sample` is "conditional" for both from frmtmb.sample's own
contribution (the review sampled a quantile_elf fit: posterior means
1.715, 3.067 against ML 1.710, 3.085). `mi_pred()` is "works" for both;
a predictor `mi(xm)` modelled by `quantile_elf()` fits at code 0 with
the Laplace-of-missing-values warning (`compat.txt`).

### 8. Core or extension

Rewritten above.

### 9. OSA on right-censored fits (filed, pre-existing)

Filed as MEDIUM in this lane's backlog section, owner OSA residuals,
citing `dev/reviews/squant-review-osa-base.txt` and `-lane.txt`.

### Sentences

The family page's sentences over 25 words were split where cheap (the
description, "The density", "What mu estimates", "Choosing the
smoothing", "Standard errors", "Group-level terms").

### Fix round 1: runs

- `test-families-quantile.R`: pass=136 fail=0 err=0 skip=0 warn=0 on
  reference BLAS (`fix1-test.txt`) and under OpenBLAS 0.3.26 with LAPACK
  (`fix1-openblas-test.txt`; the OpenBLAS copy deleted after).
- The suite files the change reaches (every core file that names the
  compat table, the family, the registry, the link sets, a CDF or
  weights, and every extension file that names the compat table),
  one per process, two at a time: `bash dev/release/run-par.sh
  squant-fix1 dev/squant-fix1.jobs 2`, log
  `dev/release/squant-fix1.log`, "SQUANT-FIX1 ran 121 of 121", every
  `lib:` line the lane core:

  | package | files | pass | fail | error | skip | warn |
  |---|---|---|---|---|---|---|
  | frmtmb | 100 | 11540 | 0 | 0 | 74 | 0 |
  | frmtmb.coupling | 1 | 58 | 0 | 0 | 0 | 0 |
  | frmtmb.eam | 3 | 113 | 0 | 0 | 0 | 0 |
  | frmtmb.latent | 1 | 311 | 0 | 0 | 0 | 0 |
  | frmtmb.learn | 6 | 302 | 0 | 0 | 0 | 0 |
  | frmtmb.ode | 1 | 15 | 0 | 0 | 0 | 0 |
  | frmtmb.sample | 2 | 293 | 0 | 0 | 1 | 0 |
  | frmtmb.spatial | 3 | 359 | 0 | 0 | 4 | 0 |
  | frmtmb.survival | 4 | 307 | 0 | 0 | 0 | 0 |

  Among them `test-compat.R` 876 (core), `test-compat-extensions.R` 2,
  `test-guardrail1.R` 42, `test-simulate-density.R` 630, all fail 0.
- One-session core runs after fix round 1 (`bash
  dev/release/run-sessions.sh squant-sess2 2 frmtmb`, log
  `dev/release/squant-sess2.log`, "SQUANT-SESS2 ran 3 of 3"): forward
  pass=25823, reverse 25991, preload 26097, each fail=0 err=0 warn=0
  (files 337 of 350, the 13 gated brms files as before); no untagged
  LEAK line.
- R CMD check --as-cran of core after fix round 1 (`bash
  dev/squant-check.sh`, log kept as `dev/squant-out/check2.log`, tree
  deleted): "Status: 1 NOTE", the V8 math-rendering NOTE on the HTML
  manual; tests OK; examples (with the new sandwich calls) OK.

## Fix round 2 (short check `dev/reviews/squant-short.md`, one MEDIUM)

### 1. gp() and the group-level rows of frm_compat() (MEDIUM, fixed)

The short check measured `ys ~ gp(x) + x2` at the default smoothing
0.25: a "Laplace approximation is not reliable" warning and 129 s at 40
rows, 501 s at 80, over 19 min at 150, against 0.6 s and 3.1 s at
smoothing 1 (`dev/reviews/squant-short-gp.txt`; 99 % in the inner
Newton search, `squant-short-gpprof.txt`). I did not rerun it.

Rows changed in `R/compat.R`, in one block after the asym_laplace rows
(so they win over the generic covstruct-by-family rule), with a shared
note `elf_note` that points at `?quantile_elf`, "Group-level terms"
(local maxima at a small smoothing; use smoothing 1 and quadrature):

- `us`, `diag`, `homdiag`, `cs`, `toep`, `homtoep`, `homcs`, `us_t`,
  `diag_t`, `s()`, `t2()` x `quantile_elf`: "works" -> "conditional".
  "works" came from the generic rule ("covariance structures act on the
  linear predictor"); I have no measurement that the local-maximum
  problem of `(1 | g)` (sections d and fix round 1, item 4) is absent
  for them, and the short check found clean fits only on well
  identified data, so "conditional" is the honest status.
- `gp_pred()` x `quantile_elf`: "works" -> "conditional", with the
  timings above and "Use smoothing = 1".
- The rows already "conditional" for other reasons (`ar1`, `hetar1`,
  `ou`, `exp`, `gau`, `mat`, `rr`, `equalto`, `gr_cov`, `gr_prec`,
  `smooth`, `gp`, `hsgp`) keep their generic notes; overriding them
  would drop those notes.

After (`dev/squant-compat4.R`, core only; `compat4-before.txt`,
`compat4-after.txt`): quantile_elf has 84 rows, works 38 -> 26,
conditional 32 -> 44, refused 13, untested 1.

`?quantile_elf`, "Group-level terms": a paragraph on the exact `gp()`
term (slow and warning at 0.25, fast at 1). Its first sentence now says
the approximation falls as the group sd grows without end, so it has a
highest point, which can be the boundary sd = 0.

### 3. cs() at a limiting correlation (filed LOW, pre-existing)

In this lane's backlog section, from `squant-short-compat2.txt` (a).

### Fix round 2: runs

Every test file that names the compat table (core and extensions) plus
`test-families-quantile.R`, one per process
(`bash dev/release/run-par.sh squant-fix2 dev/squant-fix2.jobs 2`, log
`dev/release/squant-fix2.log`, "SQUANT-FIX2 ran 55 of 55", every `lib:`
line the lane core): frmtmb 34 files pass 5123, coupling 1/58, eam
3/113, latent 1/311, learn 6/302, ode 1/15, sample 2/293, spatial 3/359,
survival 4/307; fail 0, error 0, warn 0 throughout. `test-compat.R`
878, `test-compat-extensions.R` 2, `test-families-quantile.R` 136.

R CMD check --as-cran of core after fix round 2 (log kept as
`dev/squant-out/check3.log`, tree deleted): "Status: 1 NOTE", the V8
math-rendering NOTE on the HTML manual; tests OK.
