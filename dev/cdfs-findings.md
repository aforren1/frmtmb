# Lane cdfs: CDFs for the core families that lacked one

Round opened 2026-10-20, base 4854e425 (frmtmb 0.72.0). Fix round 1
(section 11) followed the review. Worktree
`frmtmb-wt-cdfs`, branch `wt-cdfs`, private library `wt-cdfs-lib`.
Nothing committed. Blocks marked GENERATED are pasted from the logs
named beside them by `dev/cdfs-findings-build.R`; the logs live under
`dev/cdfs-out/` and `dev/release/cdfs-*` (gitignored).

## 1. What changed

- `R/families-cdf.R` (new, about 1050 lines): tape-safe log-scale CDF
  pairs `(lF, lS)` and the helpers under them.
  - `log_ibeta_two()`: both tails of the incomplete beta from
    `log_ibeta_half()` (fams2's), each direct below one half and the
    complement of the other above 0.6, with a second weight in the
    log so that a tail below 1e-300 is never taken as a complement.
  - `log_igamma_two()`: both tails of the incomplete gamma without
    `RTMB::pgamma()`: the power series of `P` (170 terms), the Legendre
    continued fraction of `Q` (60 terms, bottom up), and for shape 60
    and above 40-point Gauss-Legendre quadrature of the far tail near
    the mode, with the prefactor around the mode from Binet's remainder.
  - `student_cdf_log()` (incomplete beta, plus a quadrature near
    `z = 0`, where every `z^2` form has a zero slope),
    `skew_normal_cdf_log()` (Owen's T: absolute form where the tail is
    shallow, the positive angle integral by tanh-sinh where it is deep),
    `exgaussian_cdf_log()`, `asym_laplace_cdf_log()`,
    `von_mises_cdf_log()` (the tail as a tanh-sinh quadrature of the
    density, wrapped with `sign()`), `inverse_gaussian_cdf_log()` (both
    terms in logs), `nbinom_cdf_log()`, `binomial_cdf_log()`,
    `poisson_cdf_log()`, `beta_cdf_log()`, `gamma_cdf_log()`,
    `lognormal_cdf_log()`, `shifted_lognormal_cdf_log()`,
    `beta_binomial_cdf_log()` and `compois_cdf_log()` (log-space sums of
    the mass function), `zi_cdf_log()` and `hurdle_cdf_log()` (each base
    has one implementation), and the per-family `*_cdf_pair()` the
    slots call.
  - The truncated means brms has: `student_trunc_mean()`,
    `gamma_trunc_mean()`, `binomial_trunc_mean()`,
    `nbinom_trunc_mean()`; the count ones take the inclusive lower
    bound the likelihood reads.
- `R/families.R`: `lcdf`/`lccdf` slots and `cens`, `trunc` in
  `accepts_aterms` for student, skew_normal, exgaussian,
  shifted_lognormal, asym_laplace, zero_inflated_asym_laplace,
  von_mises, Gamma, Beta, binomial, beta_binomial, negbinomial,
  geometric, compois, zero_inflated_poisson, zero_inflated_negbinomial,
  zero_inflated_binomial, zero_inflated_beta_binomial,
  zero_inflated_beta, hurdle_poisson, hurdle_negbinomial, hurdle_gamma,
  hurdle_lognormal; an `lccdf` for poisson (its `lcdf` is unchanged,
  `RTMB::ppois()`); a log-scale `lcdf` and an `lccdf` for
  inverse.gaussian; `trunc_mean_fn` for student, Gamma, binomial,
  negbinomial, geometric; the `frmtmb_family()` documentation of which
  families declare `lccdf`.
- **Where negbinomial's CDF lives (for lane nb2dm).** In
  `negbinomial_cdf_pair()` and `nbinom_cdf_log()` in
  `R/families-cdf.R`, called from `fam_negbinomial()`'s own `lcdf` and
  `lccdf` slots with the family's links; the log shape is read with
  `dpar_log(dpars, "shape", lk_shape)` and `rate()` is folded in as
  the density folds it. nb2dm's `negbinomial2()` (`fam_reparam()`)
  wraps every slot of negbinomial() and sets `.eta_shape`, so it
  inherits the CDF and `trunc_mean_fn` at the merge with no change;
  `accepts_aterms` comes with it. Test the merged pair at consolidation.
- `R/frame.R`: the CDF refusal names the family and points at
  `frm_compat("cens()")`; brms's two truncation checks, "Truncation
  bounds are invalid: lb >= ub" and "Some responses are outside of the
  truncation bounds" (BREAKING: such a fit used to run).
- `R/compat.R`: `compat_cdf_discrete`; the families join `group:cdf`
  and `group:cdf_continuous`; the discrete ones get `cens()` and
  `trunc()` rows of their own, in a loop (a group rule would tie with
  `group:cdf`); bernoulli, xbeta and zero_one_inflated_beta are refused
  with brms's message; the continuous families without a truncated mean
  say so; `cens() x group:discrete` is gone, the kind rule covers the
  discrete families still without a CDF (nbinom1). zero_inflated_asym
  _laplace is not a compat feature, so it has no row.
- Tests: `tests/testthat/test-cdfs.R` (new, ungated) and
  `test-brms-lik-cdfs.R` (new, gated: Stan). Updated because the
  behavior they pinned changed: `test-cens-lccdf.R` (poisson and
  inverse.gaussian now declare `lccdf`), `test-cens-trunc.R`
  (negbinomial now censors; nbinom1 is the refused one),
  `test-compat.R` (beta works; the bernoulli refusal row named
  binomial), `test-hurdle-nb-zoib.R`, `test-xbeta-zibb-hurdle-cum.R`,
  `test-bcm-data-analysis.R` (the core binomial censors ChaSaSoon and
  gives `bcm_binomial_cdf()`'s likelihood).
- Docs: `vignettes/bayesian-cognitive-modeling.Rmd` and
  `inst/bcm/binomial-extras.R` (ChaSaSoon runs on the core binomial);
  NEWS.
- dev: `dev/brmsport-bin2-adapters.R` and `-spec.R` (the 15 brms
  p-function shims); `dev/upstream-bugs.md` (RTMB-8, RTMB-9,
  RTMBdist-2, RTMBdist-3, brms-31, brms-32; RTMB-7 and RTMBdist-1
  workarounds updated).

## 2. brms 2.23.0: which families it censors and truncates

`dev/cdfs-brms-survey.R`, `stancode()` with `cens(cc, y2)` (all four
codes), `trunc(lb, ub)` and `trunc(lb)`. GENERATED from
`dev/cdfs-out/survey.txt`; the Stan function named is the censored
rows' or the normalizer's:

<!-- GEN:survey -->
| family and term | brms 2.23.0 |
|---|---|
| student, cens(cc, y2) | student_t_lccdf, student_t_lcdf |
| student, trunc(lb = -1, ub = 30) | student_t_lcdf |
| student, trunc(lb = -1) | student_t_lccdf |
| skew_normal, cens(cc, y2) | skew_normal_lccdf, skew_normal_lcdf |
| skew_normal, trunc(lb = -1, ub = 30) | skew_normal_lcdf |
| skew_normal, trunc(lb = -1) | skew_normal_lccdf |
| exgaussian, cens(cc, y2) | exp_mod_normal_lccdf, exp_mod_normal_lcdf |
| exgaussian, trunc(lb = 1, ub = 30) | exp_mod_normal_lcdf |
| exgaussian, trunc(lb = 1) | exp_mod_normal_lccdf |
| shifted_lognormal, cens(cc, y2) | lognormal_lccdf, lognormal_lcdf |
| shifted_lognormal, trunc(lb = 1, ub = 30) | lognormal_lcdf |
| shifted_lognormal, trunc(lb = 1) | lognormal_lccdf |
| asym_laplace, cens(cc, y2) | asym_laplace_lccdf, asym_laplace_lcdf |
| asym_laplace, trunc(lb = -1, ub = 30) | asym_laplace_lcdf, asym_laplace_lccdf |
| asym_laplace, trunc(lb = -1) | asym_laplace_lccdf, asym_laplace_lcdf |
| zero_inflated_asym_laplace, cens(cc, y2) | zero_inflated_asym_laplace_lccdf, zero_inflated_asym_laplace_lcdf, asym_laplace_lcdf, asym_laplace_lccdf |
| zero_inflated_asym_laplace, trunc(lb = -1, ub = 30) | zero_inflated_asym_laplace_lcdf, asym_laplace_lcdf, asym_laplace_lccdf, zero_inflated_asym_laplace_lccdf |
| zero_inflated_asym_laplace, trunc(lb = -1) | zero_inflated_asym_laplace_lccdf, asym_laplace_lcdf, asym_laplace_lccdf, zero_inflated_asym_laplace_lcdf |
| von_mises, cens(cc, y2) | von_mises_lccdf, von_mises_lcdf |
| von_mises, trunc(lb = -2, ub = 2) | von_mises_lcdf |
| von_mises, trunc(lb = -2) | von_mises_lccdf |
| Gamma, cens(cc, y2) | gamma_lccdf, gamma_lcdf |
| Gamma, trunc(lb = 1, ub = 30) | gamma_lcdf |
| Gamma, trunc(lb = 1) | gamma_lccdf |
| Beta, cens(cc, y2) | beta_lccdf, beta_lcdf |
| Beta, trunc(lb = 0.05, ub = 0.95) | beta_lcdf |
| Beta, trunc(lb = 0.05) | beta_lccdf |
| xbeta, cens(cc, y2) | refused: Argument 'cens' is not supported for family 'xbeta(logit)'. |
| xbeta, trunc(lb = 0.05, ub = 0.95) | refused: Argument 'trunc' is not supported for family 'xbeta(logit)'. |
| xbeta, trunc(lb = 0.05) | refused: Argument 'trunc' is not supported for family 'xbeta(logit)'. |
| binomial, cens(cc, y2) | binomial_lccdf, binomial_lcdf |
| binomial, trunc(lb = 1, ub = 30) | binomial_lcdf |
| binomial, trunc(lb = 1) | binomial_lccdf |
| bernoulli, cens(cc, y2) | refused: Argument 'cens' is not supported for family 'bernoulli(logit)'. |
| bernoulli, trunc(lb = 1, ub = 30) | refused: Argument 'trunc' is not supported for family 'bernoulli(logit)'. |
| bernoulli, trunc(lb = 1) | refused: Argument 'trunc' is not supported for family 'bernoulli(logit)'. |
| beta_binomial, cens(cc, y2) | beta_binomial_lccdf, beta_binomial_lcdf |
| beta_binomial, trunc(lb = 1, ub = 30) | beta_binomial_lcdf |
| beta_binomial, trunc(lb = 1) | beta_binomial_lccdf |
| negbinomial, cens(cc, y2) | neg_binomial_2_lccdf, neg_binomial_2_lcdf |
| negbinomial, trunc(lb = 1, ub = 30) | neg_binomial_2_lcdf |
| negbinomial, trunc(lb = 1) | neg_binomial_2_lccdf |
| negbinomial2, cens(cc, y2) | neg_binomial_2_lccdf, neg_binomial_2_lcdf |
| negbinomial2, trunc(lb = 1, ub = 30) | neg_binomial_2_lcdf |
| negbinomial2, trunc(lb = 1) | neg_binomial_2_lccdf |
| geometric, cens(cc, y2) | neg_binomial_2_lccdf, neg_binomial_2_lcdf |
| geometric, trunc(lb = 1, ub = 30) | neg_binomial_2_lcdf |
| geometric, trunc(lb = 1) | neg_binomial_2_lccdf |
| com_poisson, cens(cc, y2) | com_poisson_lccdf, com_poisson_lcdf |
| com_poisson, trunc(lb = 1, ub = 30) | com_poisson_lcdf, com_poisson_lccdf |
| com_poisson, trunc(lb = 1) | com_poisson_lccdf, com_poisson_lcdf |
| poisson, cens(cc, y2) | poisson_lccdf, poisson_lcdf |
| poisson, trunc(lb = 1, ub = 30) | poisson_lcdf |
| poisson, trunc(lb = 1) | poisson_lccdf |
| inverse.gaussian, cens(cc, y2) | inv_gaussian_lccdf, inv_gaussian_lcdf |
| inverse.gaussian, trunc(lb = 1, ub = 30) | inv_gaussian_lcdf, inv_gaussian_lccdf |
| inverse.gaussian, trunc(lb = 1) | inv_gaussian_lccdf, inv_gaussian_lcdf |
| zero_inflated_poisson, cens(cc, y2) | zero_inflated_poisson_lccdf, zero_inflated_poisson_lcdf |
| zero_inflated_poisson, trunc(lb = 1, ub = 30) | zero_inflated_poisson_lcdf, zero_inflated_poisson_lccdf |
| zero_inflated_poisson, trunc(lb = 1) | zero_inflated_poisson_lccdf, zero_inflated_poisson_lcdf |
| zero_inflated_negbinomial, cens(cc, y2) | zero_inflated_neg_binomial_lccdf, zero_inflated_neg_binomial_lcdf |
| zero_inflated_negbinomial, trunc(lb = 1, ub = 30) | zero_inflated_neg_binomial_lcdf, zero_inflated_neg_binomial_lccdf |
| zero_inflated_negbinomial, trunc(lb = 1) | zero_inflated_neg_binomial_lccdf, zero_inflated_neg_binomial_lcdf |
| zero_inflated_binomial, cens(cc, y2) | zero_inflated_binomial_lccdf, zero_inflated_binomial_lcdf |
| zero_inflated_binomial, trunc(lb = 1, ub = 30) | zero_inflated_binomial_lcdf, zero_inflated_binomial_lccdf |
| zero_inflated_binomial, trunc(lb = 1) | zero_inflated_binomial_lccdf, zero_inflated_binomial_lcdf |
| zero_inflated_beta_binomial, cens(cc, y2) | zero_inflated_beta_binomial_lccdf, zero_inflated_beta_binomial_lcdf |
| zero_inflated_beta_binomial, trunc(lb = 1, ub = 30) | zero_inflated_beta_binomial_lcdf, zero_inflated_beta_binomial_lccdf |
| zero_inflated_beta_binomial, trunc(lb = 1) | zero_inflated_beta_binomial_lccdf, zero_inflated_beta_binomial_lcdf |
| zero_inflated_beta, cens(cc, y2) | zero_inflated_beta_lccdf, zero_inflated_beta_lcdf |
| zero_inflated_beta, trunc(lb = 0.05, ub = 0.95) | zero_inflated_beta_lcdf, zero_inflated_beta_lccdf |
| zero_inflated_beta, trunc(lb = 0.05) | zero_inflated_beta_lccdf, zero_inflated_beta_lcdf |
| zero_one_inflated_beta, cens(cc, y2) | refused: Argument 'cens' is not supported for family 'zero_one_inflated_beta(logit)'. |
| zero_one_inflated_beta, trunc(lb = 0.05, ub = 0.95) | refused: Argument 'trunc' is not supported for family 'zero_one_inflated_beta(logit)'. |
| zero_one_inflated_beta, trunc(lb = 0.05) | refused: Argument 'trunc' is not supported for family 'zero_one_inflated_beta(logit)'. |
| hurdle_poisson, cens(cc, y2) | hurdle_poisson_lccdf, hurdle_poisson_lcdf |
| hurdle_poisson, trunc(lb = 1, ub = 30) | hurdle_poisson_lcdf, hurdle_poisson_lccdf |
| hurdle_poisson, trunc(lb = 1) | hurdle_poisson_lccdf, hurdle_poisson_lcdf |
| hurdle_negbinomial, cens(cc, y2) | hurdle_neg_binomial_lccdf, hurdle_neg_binomial_lcdf |
| hurdle_negbinomial, trunc(lb = 1, ub = 30) | hurdle_neg_binomial_lcdf, hurdle_neg_binomial_lccdf |
| hurdle_negbinomial, trunc(lb = 1) | hurdle_neg_binomial_lccdf, hurdle_neg_binomial_lcdf |
| hurdle_gamma, cens(cc, y2) | hurdle_gamma_lccdf, hurdle_gamma_lcdf |
| hurdle_gamma, trunc(lb = 1, ub = 30) | hurdle_gamma_lcdf, hurdle_gamma_lccdf |
| hurdle_gamma, trunc(lb = 1) | hurdle_gamma_lccdf, hurdle_gamma_lcdf |
| hurdle_lognormal, cens(cc, y2) | hurdle_lognormal_lccdf, hurdle_lognormal_lcdf |
| hurdle_lognormal, trunc(lb = 1, ub = 30) | hurdle_lognormal_lcdf, hurdle_lognormal_lccdf |
| hurdle_lognormal, trunc(lb = 1) | hurdle_lognormal_lccdf, hurdle_lognormal_lcdf |
<!-- /GEN:survey -->

Read: brms censors and truncates every family above except xbeta,
bernoulli and zero_one_inflated_beta, which it refuses ("Argument
'cens' is not supported for family"). tweedie, huber and nbinom1 are
not brms families; frmtmb refuses them by the kind rule (no CDF). A
discrete truncation normalizer is `lcdf(lb - 1)` in brms as in frmtmb.
A discrete right or interval bound is exclusive in brms; it was
inclusive in frmtmb through the first build and is brms's since fix
round 1 (section 11; `?frmtmb_family`, "Censoring a discrete
response").

## 3. RTMB primitives (why most of this is hand-written)

GENERATED from `dev/cdfs-out/probe1.txt` (`dev/cdfs-probe1.R`) and
`probe3.txt` (`dev/cdfs-probe3.R`):

<!-- GEN:probe -->
```
-- dev/cdfs-out/probe1.txt (which RTMB CDFs tape; gradient vs FD)
pgamma lower                             value 0.69216221519748  ref 0.69216221519748  grad -0.61968106,-0.64038948  fd -0.61968106,-0.64038948
pgamma upper log                         ERROR: Non-numeric argument to mathematical function
pbeta                                    value 0.579825  ref 0.579825  grad 0,0.49722447  fd -0.50483156,0.49722447
pbeta upper                              ERROR: Non-numeric argument to mathematical function
pbinom prob                              value 0.378994383004009  ref 0.378994383004009  grad -0.60063618  fd -0.60063618
pnbinom mu                               ERROR: Non-numeric argument to mathematical function
pnbinom prob                             value 0.891843332177153  ref 0.891843332177153  grad 0,0.21644924  fd -0.18873799,0.21644925
ppois                                    value 0.433470120366709  ref 0.433470120366709  grad -0.78146726  fd -0.78146726
pnorm log upper                          value -2.33510327866244  ref -2.33510327866244  grad -1.7703278  fd -1.7703278
pnorm log lower                          value -2.33510327866244  ref -2.33510327866244  grad 1.7703278  fd 1.7703278
sign                                     value -0.49  ref -0.49  grad 1.4  fd 1.4
besselI nu=3                             value 0.297627709532404  ref 0.297627709532404  grad 1.0630419  fd 1.0630419
besselI expon                            value 0.0329780902323376  ref 0.0329780902323376  grad 0.045236598  fd 0.045236598
lgamma                                   value 0.0969474667906389  ref 0.0969474667906389  grad 1.1974456  fd 1.1974456
log1p                                    value 0.262364264467491  ref 0.262364264467491  grad 0.76923077  fd 0.76923077
expm1                                    value 0.349858807576003  ref 0.349858807576003  grad 1.3498588  fd 1.3498588
RTMBdist::pt                             value 0.769992966548131  ref 0.769992966548131  grad 0.26448836,0.017257752  fd 0.26448836,0.017257752
RTMBdist::pskewnorm                      ERROR: Can only determine finite status of constant expressions
RTMBdist::pexgauss                       value 0.412084042516506  ref 0.412084042516506  grad 0.30940875,0.047290809,0.17105431  fd 0.30940875,0.047290809,0.17105431
RTMBdist::pvm                            ERROR: Comparison is generally unsafe for AD types
pgamma upper (ppois tail)                value 0.566529879633291  ref 0.566529879633291  grad 0.78146726  fd 0.78146726
-- dev/cdfs-out/probe3.txt (tails; first lines)
== RTMB::ppois left tail, log, against stats::ppois(log.p = TRUE)
lambda    5 k   0  rtmb -5  ref -5  rel 0.00e+00
lambda    5 k   2  rtmb -2.08222926791572  ref -2.08222926791572  rel 0.00e+00
lambda    5 k  10  rtmb -0.0137899139115043  ref -0.0137899139115043  rel 5.38e-17
lambda   50 k   0  rtmb -50  ref -50  rel 0.00e+00
lambda   50 k   2  rtmb -42.8291115214875  ref -42.8291115214875  rel 0.00e+00
lambda   50 k  10  rtmb -25.7669172771077  ref -25.7669172771077  rel 0.00e+00
lambda  200 k   0  rtmb -200  ref -200  rel 0.00e+00
lambda  200 k   2  rtmb -190.086512612886  ref -190.086512612886  rel 0.00e+00
lambda  200 k  10  rtmb -162.070219566164  ref -162.070219566164  rel 1.75e-16
lambda  700 k   0  rtmb -700  ref -700  rel 0.00e+00
lambda  700 k   2  rtmb -687.588129371495  ref -687.588129371495  rel 0.00e+00
lambda  700 k  10  rtmb -649.579241427471  ref -649.579241427471  rel 0.00e+00
== RTMB::pgamma lower tail (poisson upper via P(Y > k) = P(G(k+1) < lam))
lambda  0.5 k  10  rtmb -25.5845108115866  ref -25.5845108115866  rel 0.00e+00
lambda  0.5 k  50  rtmb -188.25043865934  ref -188.25043865934  rel 0.00e+00
lambda  0.5 k 200  rtmb -Inf  ref -1008.35539710733  rel Inf
lambda  0.5 k 700  rtmb -Inf  ref -4382.89880123452  rel Inf
lambda    5 k  10  rtmb -4.29070486358006  ref -4.29070486358006  rel 0.00e+00
lambda    5 k  50  rtmb -75.2273746987172  ref -75.2273746987172  rel 0.00e+00
lambda    5 k 200  rtmb -550.013210904247  ref -550.013210904247  rel 0.00e+00
lambda    5 k 700  rtmb -Inf  ref -2773.28021562938  rel Inf
== RTMB::pgamma lower tail at small x
a   0.3 x   0.0003 rtmb -2.32541284278079 ref -2.32541284278079 rel 1.91e-16
a   0.3 x     0.03 rtmb -0.950681018757343 ref -0.950681018757343 rel 0.00e+00
a     3 x    0.003 rtmb -19.221438271392 ref -19.221438271392 rel 0.00e+00
a     3 x      0.3 rtmb -5.626962225906 ref -5.626962225906 rel 1.58e-16
a    30 x     0.03 rtmb -179.884005087249 ref -179.884005087249 rel 0.00e+00
a    30 x        3 rtmb -44.5984388730115 ref -44.5984388730115 rel 1.59e-16
a   300 x      0.3 rtmb -Inf ref -1776.39669407139 rel Inf
a   300 x       30 rtmb -424.441684517911 ref -424.441684517911 rel 0.00e+00
== log_ibeta_half over (0, 1): worst relative log error by shape pair
a    0.05 b   0.05  worst 1.37e-01 at x = 0.999999  nonfinite 0
a    0.05 b    0.5  worst 5.55e-16 at x = 0.04  nonfinite 0
a    0.05 b      1  worst 4.44e-16 at x = 0.4  nonfinite 0
a    0.05 b    3.1  worst 5.45e-16 at x = 0.24  nonfinite 0
a    0.05 b     25  worst 1.53e-16 at x = 0.01  nonfinite 0
a    0.05 b    900  worst 4.38e-16 at x = 1e-12  nonfinite 0
a     0.5 b   0.05  worst 1.59e-15 at x = 0.9  nonfinite 0
a     0.5 b    0.5  worst 1.26e-02 at x = 0.999999  nonfinite 0
```
<!-- /GEN:probe -->

- `RTMB::pbeta()`'s derivative in a shape is exactly 0 at a whole shape
  (`dev/cdfs-probe2.R`: 0 against -0.2524 at `pbeta(0.3, 2, 5)`, fine
  at 2.0001). `nu / 2` is whole at the student-t start `nu = 4`.
- `RTMB::pgamma()`'s third derivatives are not finite at 960 of 2000
  random points (`dev/cdfs-dbg9.R`); `RTMB::ppois()`'s are finite at
  500 of 500 (`dev/cdfs-dbg13.R`), so poisson's `lcdf` keeps it.
- `RTMBdist::pt()` is `pbeta()` underneath; `pskewnorm()` and `pvm()`
  do not tape; `pexgauss()` tapes but forms `1 - F` for the upper tail.
- `log_ibeta_half()` holds to 1.1e-14 on `x <= 0.6`
  (`dev/cdfs-probe4.R`), and above it wherever its result is not near
  one (`dev/cdfs-dbg4.R`, against 300-bit `Rmpfr::pbetaI()`: 5e-16 at
  `log I = -10991`, shapes (50000, 25), where `stats::pbeta()` is
  1.3e-3 off or `-Inf`); it fails near `x = 1` with both shapes small
  (0.21 at (0.01, 0.01), `dev/cdfs-dbg5.R`), which `log_ibeta_two()`
  avoids.

## 4. Accuracy against independent references

`dev/cdfs-acc.R` (seed 20261020), GENERATED from `dev/cdfs-out/acc.txt`.
Relative error of the log, floored at one. The `log_ibeta_two` grid's
8.2e-5 is `stats::pbeta()`'s error, not the pair's: the referee line
under it (300-bit `pbetaI()` at whole shapes beside the worst point)
gives the pair's value to all ten digits and `stats` off in the sixth.

<!-- GEN:acc -->
```
== log_ibeta_two against stats::pbeta(log.p = TRUE)
shapes 9 x 9, x 105 points: worst lF 8.16e-05  worst lS 8.16e-05
x 0.99 shapes (1e5, 7): referee -970.16014318  log_ibeta_two -970.16014318  stats -970.160912419
x 0.995 shapes (1e5, 7): referee -470.533587487  log_ibeta_two -470.533587487  stats -470.533587487
== log_igamma_two against stats::pgamma(log.p = TRUE)
worst lP 1.02e-14 at a = 10000 x = 6000; worst lQ 8.87e-14 at a = 30 x = 44.6931
== student-t against stats::pt(log.p = TRUE)
t nu = 1.01  lF                                n    17  worst 2.78e-16  median 4.7e-18  nonfinite 0
t nu = 1.01  lS                                n    17  worst 2.78e-16  median 4.7e-18  nonfinite 0
t nu = 1.5  lF                                 n    17  worst 1.64e-16  median 5.0e-17  nonfinite 0
t nu = 1.5  lS                                 n    17  worst 1.64e-16  median 5.0e-17  nonfinite 0
t nu = 2  lF                                   n    17  worst 1.46e-16  median 0.0e+00  nonfinite 0
t nu = 2  lS                                   n    17  worst 1.46e-16  median 0.0e+00  nonfinite 0
t nu = 3.7  lF                                 n    17  worst 2.67e-16  median 1.1e-16  nonfinite 0
t nu = 3.7  lS                                 n    17  worst 2.67e-16  median 1.1e-16  nonfinite 0
t nu = 4  lF                                   n    17  worst 2.65e-16  median 9.6e-18  nonfinite 0
t nu = 4  lS                                   n    17  worst 2.65e-16  median 9.6e-18  nonfinite 0
t nu = 10  lF                                  n    17  worst 3.55e-16  median 1.0e-17  nonfinite 0
t nu = 10  lS                                  n    17  worst 3.55e-16  median 1.0e-17  nonfinite 0
t nu = 30  lF                                  n    17  worst 3.00e-16  median 9.4e-31  nonfinite 0
t nu = 30  lS                                  n    17  worst 3.00e-16  median 9.4e-31  nonfinite 0
t nu = 1000  lF                                n    17  worst 9.65e-16  median 4.8e-17  nonfinite 0
t nu = 1000  lS                                n    17  worst 9.65e-16  median 4.8e-17  nonfinite 0
t nu = 400000  lF                              n    17  worst 2.00e-15  median 1.1e-16  nonfinite 0
t nu = 400000  lS                              n    17  worst 2.00e-15  median 1.1e-16  nonfinite 0
== skew-normal: alpha = +-1 is exact, Phi(z)^2 and 1 - Phi(-z)^2
skew_normal alpha = 1  lF                      n    17  worst 5.29e-14  median 1.1e-16  nonfinite 0
skew_normal alpha = 1  lS                      n    13  worst 1.18e-12  median 1.3e-17  nonfinite 0
skew_normal alpha = 1  lS (exact form)         n    17  worst 8.15e-14  median 5.6e-17  nonfinite 0
skew_normal alpha = -1 lS                      n    17  worst 5.28e-14  median 5.6e-17  nonfinite 0
skew_normal alpha = -1 lF                      n    17  worst 8.15e-14  median 1.7e-17  nonfinite 0
== skew-normal against sn::psn (body) and angle quadrature (tails)
sn::psn body (both tails above 1e-6): worst 6.12e-13
angle-quadrature tails, 49 points per side (z to +-38, alpha to +-300), 12 without a reference: worst lF 5.47e-14  lS 1.29e-14
== ex-gaussian against the closed form in 256-bit arithmetic
exgaussian 9 (sigma, beta) x 10 q: worst lF 1.55e-13  lS 9.80e-14
== asymmetric Laplace: brms::pasym_laplace() lower tail, and the
   closed form in 200-bit arithmetic for both (brms forms the
   upper tail as 1 - F)
asym_laplace p = 0.01 lF vs brms               n     9  worst 3.12e-16  median 0.0e+00  nonfinite 0
asym_laplace p = 0.01 lF                       n     9  worst 2.11e-16  median 0.0e+00  nonfinite 0
asym_laplace p = 0.01 lS                       n     9  worst 4.51e-17  median 8.7e-18  nonfinite 0
asym_laplace p = 0.3 lF vs brms                n     9  worst 3.33e-16  median 1.4e-16  nonfinite 0
asym_laplace p = 0.3 lF                        n     9  worst 2.22e-16  median 3.0e-17  nonfinite 0
asym_laplace p = 0.3 lS                        n     9  worst 1.11e-16  median 5.6e-17  nonfinite 0
asym_laplace p = 0.5 lF vs brms                n     9  worst 0.00e+00  median 0.0e+00  nonfinite 0
asym_laplace p = 0.5 lF                        n     9  worst 1.11e-16  median 0.0e+00  nonfinite 0
asym_laplace p = 0.5 lS                        n     9  worst 1.11e-16  median 0.0e+00  nonfinite 0
asym_laplace p = 0.9 lF vs brms                n     9  worst 1.11e-16  median 0.0e+00  nonfinite 0
asym_laplace p = 0.9 lF                        n     9  worst 4.08e-17  median 0.0e+00  nonfinite 0
asym_laplace p = 0.9 lS                        n     9  worst 1.50e-16  median 0.0e+00  nonfinite 0
== von Mises against adaptive quadrature of the density, its
   largest value on the interval divided out so that the far side
   of a concentrated distribution does not underflow
von_mises: worst lF 1.74e-13 at k=0.01 mu=-3 y=-3.141 ref=-8.736
von_mises: worst lS 1.26e-13 at k=0.01 mu=-3 y=3.141 ref=-8.736
== inverse gaussian against statmod::pinvgauss(log.p = TRUE)
inv.gaussian mu 0.1 shape 0.1 lF               n     8  worst 4.41e-16  median 1.1e-16  nonfinite 0
inv.gaussian mu 0.1 shape 0.1 lS               n     8  worst 1.91e-14  median 8.3e-17  nonfinite 0
inv.gaussian mu 0.1 shape 1 lF                 n     8  worst 7.77e-16  median 1.7e-56  nonfinite 0
inv.gaussian mu 0.1 shape 1 lS                 n     8  worst 9.99e-16  median 1.3e-16  nonfinite 0
inv.gaussian mu 0.1 shape 50 lF                n     8  worst 8.85e-16  median 0.0e+00  nonfinite 0
inv.gaussian mu 0.1 shape 50 lS                n     8  worst 1.54e-15  median 0.0e+00  nonfinite 0
inv.gaussian mu 0.1 shape 1000 lF              n     8  worst 5.39e-16  median 0.0e+00  nonfinite 0
inv.gaussian mu 0.1 shape 1000 lS              n     8  worst 2.73e-15  median 0.0e+00  nonfinite 0
inv.gaussian mu 1 shape 0.1 lF                 n     8  worst 0.00e+00  median 0.0e+00  nonfinite 0
inv.gaussian mu 1 shape 0.1 lS                 n     8  worst 1.71e-16  median 0.0e+00  nonfinite 0
inv.gaussian mu 1 shape 1 lF                   n     8  worst 4.41e-16  median 8.1e-17  nonfinite 0
inv.gaussian mu 1 shape 1 lS                   n     8  worst 1.91e-14  median 8.4e-17  nonfinite 0
inv.gaussian mu 1 shape 50 lF                  n     8  worst 2.92e-16  median 7.5e-212  nonfinite 0
inv.gaussian mu 1 shape 50 lS                  n     8  worst 2.61e-15  median 2.1e-17  nonfinite 0
inv.gaussian mu 1 shape 1000 lF                n     8  worst 5.66e-15  median 1.1e-123  nonfinite 0
inv.gaussian mu 1 shape 1000 lS                n     8  worst 5.77e-15  median 0.0e+00  nonfinite 0
inv.gaussian mu 10 shape 0.1 lF                n     8  worst 1.38e-16  median 0.0e+00  nonfinite 0
inv.gaussian mu 10 shape 0.1 lS                n     8  worst 2.45e-16  median 1.5e-16  nonfinite 0
inv.gaussian mu 10 shape 1 lF                  n     8  worst 0.00e+00  median 0.0e+00  nonfinite 0
inv.gaussian mu 10 shape 1 lS                  n     8  worst 1.71e-16  median 0.0e+00  nonfinite 0
inv.gaussian mu 10 shape 50 lF                 n     8  worst 6.33e-16  median 0.0e+00  nonfinite 0
inv.gaussian mu 10 shape 50 lS                 n     8  worst 5.24e-16  median 0.0e+00  nonfinite 0
inv.gaussian mu 10 shape 1000 lF               n     8  worst 2.92e-16  median 0.0e+00  nonfinite 0
inv.gaussian mu 10 shape 1000 lS               n     8  worst 3.87e-15  median 0.0e+00  nonfinite 0
== negative binomial against stats::pnbinom(log.p = TRUE)
nbinom mu 1e-06 shape 0.05 lF                  n     7  worst 4.90e-17  median 1.1e-54  nonfinite 0
nbinom mu 1e-06 shape 0.05 lS                  n     8  worst 1.41e-16  median 0.0e+00  nonfinite 0
nbinom mu 1e-06 shape 1 lF                     n     7  worst 4.99e-17  median 1.0e-66  nonfinite 0
nbinom mu 1e-06 shape 1 lS                     n     8  worst 1.29e-16  median 0.0e+00  nonfinite 0
nbinom mu 1e-06 shape 2 lF                     n     7  worst 4.33e-17  median 5.9e-69  nonfinite 0
nbinom mu 1e-06 shape 2 lS                     n     8  worst 3.11e-16  median 1.3e-16  nonfinite 0
nbinom mu 1e-06 shape 7.5 lF                   n     7  worst 3.74e-17  median 4.8e-72  nonfinite 0
nbinom mu 1e-06 shape 7.5 lS                   n     8  worst 1.23e-16  median 0.0e+00  nonfinite 0
nbinom mu 1e-06 shape 1000 lF                  n     7  worst 4.00e-17  median 2.6e-74  nonfinite 0
nbinom mu 1e-06 shape 1000 lS                  n     8  worst 1.68e-16  median 5.6e-17  nonfinite 0
nbinom mu 0.3 shape 0.05 lF                    n     7  worst 2.78e-17  median 3.5e-18  nonfinite 0
nbinom mu 0.3 shape 0.05 lS                    n     8  worst 5.26e-16  median 1.8e-16  nonfinite 0
nbinom mu 0.3 shape 1 lF                       n     7  worst 6.94e-17  median 1.8e-20  nonfinite 0
nbinom mu 0.3 shape 1 lS                       n     8  worst 1.51e-16  median 0.0e+00  nonfinite 0
nbinom mu 0.3 shape 2 lF                       n     7  worst 1.17e-15  median 1.2e-17  nonfinite 0
nbinom mu 0.3 shape 2 lS                       n     8  worst 2.52e-15  median 1.2e-16  nonfinite 0
nbinom mu 0.3 shape 7.5 lF                     n     7  worst 2.78e-16  median 2.8e-17  nonfinite 0
nbinom mu 0.3 shape 7.5 lS                     n     8  worst 6.89e-16  median 2.6e-16  nonfinite 0
nbinom mu 0.3 shape 1000 lF                    n     7  worst 1.67e-16  median 1.5e-17  nonfinite 0
nbinom mu 0.3 shape 1000 lS                    n     8  worst 4.59e-16  median 1.3e-16  nonfinite 0
nbinom mu 4 shape 0.05 lF                      n     7  worst 3.75e-16  median 1.4e-16  nonfinite 0
nbinom mu 4 shape 0.05 lS                      n     8  worst 2.98e-15  median 6.7e-16  nonfinite 0
nbinom mu 4 shape 1 lF                         n     7  worst 1.54e-15  median 7.8e-26  nonfinite 0
nbinom mu 4 shape 1 lS                         n     8  worst 6.88e-15  median 2.8e-17  nonfinite 0
nbinom mu 4 shape 2 lF                         n     7  worst 1.11e-15  median 1.5e-16  nonfinite 0
nbinom mu 4 shape 2 lS                         n     8  worst 1.33e-15  median 6.4e-16  nonfinite 0
nbinom mu 4 shape 7.5 lF                       n     7  worst 5.08e-15  median 2.1e-17  nonfinite 0
nbinom mu 4 shape 7.5 lS                       n     8  worst 1.64e-15  median 2.6e-16  nonfinite 0
nbinom mu 4 shape 1000 lF                      n     7  worst 5.96e-15  median 2.1e-17  nonfinite 0
nbinom mu 4 shape 1000 lS                      n     8  worst 2.58e-15  median 1.5e-16  nonfinite 0
nbinom mu 100 shape 0.05 lF                    n     7  worst 3.33e-16  median 1.9e-16  nonfinite 0
nbinom mu 100 shape 0.05 lS                    n     8  worst 2.00e-15  median 3.8e-16  nonfinite 0
nbinom mu 100 shape 1 lF                       n     7  worst 3.33e-16  median 1.1e-16  nonfinite 0
nbinom mu 100 shape 1 lS                       n     8  worst 6.63e-16  median 5.7e-17  nonfinite 0
nbinom mu 100 shape 2 lF                       n     7  worst 1.33e-15  median 4.8e-16  nonfinite 0
nbinom mu 100 shape 2 lS                       n     8  worst 1.89e-15  median 5.7e-17  nonfinite 0
nbinom mu 100 shape 7.5 lF                     n     7  worst 5.00e-15  median 1.8e-16  nonfinite 0
nbinom mu 100 shape 7.5 lS                     n     8  worst 6.11e-15  median 2.1e-17  nonfinite 0
nbinom mu 100 shape 1000 lF                    n     7  worst 2.10e-15  median 9.4e-16  nonfinite 0
nbinom mu 100 shape 1000 lS                    n     8  worst 1.67e-15  median 2.3e-42  nonfinite 0
nbinom mu 10000 shape 0.05 lF                  n     7  worst 3.33e-16  median 2.2e-16  nonfinite 0
nbinom mu 10000 shape 0.05 lS                  n     8  worst 7.77e-16  median 4.4e-16  nonfinite 0
nbinom mu 10000 shape 1 lF                     n     7  worst 7.55e-16  median 2.6e-16  nonfinite 0
nbinom mu 10000 shape 1 lS                     n     8  worst 1.39e-16  median 1.8e-17  nonfinite 0
nbinom mu 10000 shape 2 lF                     n     7  worst 3.28e-16  median 1.1e-16  nonfinite 0
nbinom mu 10000 shape 2 lS                     n     8  worst 4.86e-17  median 1.7e-17  nonfinite 0
nbinom mu 10000 shape 7.5 lF                   n     7  worst 1.92e-16  median 0.0e+00  nonfinite 0
nbinom mu 10000 shape 7.5 lS                   n     8  worst 4.64e-17  median 5.3e-20  nonfinite 0
nbinom mu 10000 shape 1000 lF                  n     7  worst 1.02e-15  median 5.7e-16  nonfinite 0
nbinom mu 10000 shape 1000 lS                  n     8  worst 0.00e+00  median 0.0e+00  nonfinite 0
== binomial against stats::pbinom(log.p = TRUE)
binomial size 1 p 1e-08 lF                     n     6  worst 5.02e-17  median 2.5e-17  nonfinite 0
binomial size 1 p 1e-08 lS                     n     4  worst 0.00e+00  median 0.0e+00  nonfinite 0
binomial size 1 p 0.05 lF                      n     6  worst 3.47e-17  median 1.7e-17  nonfinite 0
binomial size 1 p 0.05 lS                      n     4  worst 0.00e+00  median 0.0e+00  nonfinite 0
binomial size 1 p 0.5 lF                       n     6  worst 0.00e+00  median 0.0e+00  nonfinite 0
binomial size 1 p 0.5 lS                       n     4  worst 0.00e+00  median 0.0e+00  nonfinite 0
binomial size 1 p 0.97 lF                      n     6  worst 0.00e+00  median 0.0e+00  nonfinite 0
binomial size 1 p 0.97 lS                      n     4  worst 3.47e-18  median 3.5e-18  nonfinite 0
binomial size 10 p 1e-08 lF                    n     6  worst 5.19e-17  median 1.0e-30  nonfinite 0
binomial size 10 p 1e-08 lS                    n     5  worst 0.00e+00  median 0.0e+00  nonfinite 0
binomial size 10 p 0.05 lF                     n     6  worst 2.50e-16  median 4.8e-18  nonfinite 0
binomial size 10 p 0.05 lS                     n     5  worst 9.06e-16  median 0.0e+00  nonfinite 0
binomial size 10 p 0.5 lF                      n     6  worst 7.57e-16  median 0.0e+00  nonfinite 0
binomial size 10 p 0.5 lS                      n     5  worst 1.44e-15  median 0.0e+00  nonfinite 0
binomial size 10 p 0.97 lF                     n     6  worst 0.00e+00  median 0.0e+00  nonfinite 0
binomial size 10 p 0.97 lS                     n     5  worst 5.55e-17  median 3.5e-17  nonfinite 0
binomial size 1000 p 1e-08 lF                  n     6  worst 4.34e-17  median 0.0e+00  nonfinite 0
binomial size 1000 p 1e-08 lS                  n     5  worst 1.65e-16  median 0.0e+00  nonfinite 0
binomial size 1000 p 0.05 lF                   n     6  worst 1.50e-16  median 2.7e-175  nonfinite 0
binomial size 1000 p 0.05 lS                   n     5  worst 1.42e-16  median 3.4e-37  nonfinite 0
binomial size 1000 p 0.5 lF                    n     6  worst 1.66e-16  median 0.0e+00  nonfinite 0
binomial size 1000 p 0.5 lS                    n     5  worst 7.61e-41  median 0.0e+00  nonfinite 0
binomial size 1000 p 0.97 lF                   n     6  worst 1.26e-28  median 0.0e+00  nonfinite 0
binomial size 1000 p 0.97 lS                   n     5  worst 1.17e-16  median 0.0e+00  nonfinite 0
== poisson against stats::ppois(log.p = TRUE)
poisson mu 1e-06 lF                            n     8  worst 4.48e-17  median 1.3e-74  nonfinite 0
poisson mu 1e-06 lS                            n     9  worst 3.65e-16  median 0.0e+00  nonfinite 0
poisson mu 0.5 lF                              n     8  worst 3.27e-17  median 2.6e-103  nonfinite 0
poisson mu 0.5 lS                              n     9  worst 2.78e-16  median 0.0e+00  nonfinite 0
poisson mu 5 lF                                n     8  worst 1.80e-15  median 3.1e-44  nonfinite 0
poisson mu 5 lS                                n     9  worst 1.24e-15  median 8.7e-18  nonfinite 0
poisson mu 80 lF                               n     8  worst 6.86e-14  median 6.9e-17  nonfinite 0
poisson mu 80 lS                               n     9  worst 3.67e-15  median 6.1e-23  nonfinite 0
poisson mu 3000 lF                             n     8  worst 1.77e-16  median 0.0e+00  nonfinite 0
poisson mu 3000 lS                             n     9  worst 7.73e-15  median 0.0e+00  nonfinite 0
== beta-binomial against a 200-bit sum of its mass function
beta_binomial size 1 a 0.3 b 0.3 lF            n     5  worst 3.33e-16  median 3.3e-16  nonfinite 0
beta_binomial size 1 a 0.3 b 0.3 lS            n     4  worst 3.33e-16  median 3.3e-16  nonfinite 0
beta_binomial size 1 a 2 b 9 lF                n     5  worst 1.28e-15  median 1.3e-15  nonfinite 0
beta_binomial size 1 a 2 b 9 lS                n     4  worst 4.36e-60  median 0.0e+00  nonfinite 0
beta_binomial size 1 a 400 b 30 lF             n     5  worst 5.00e-14  median 5.0e-14  nonfinite 0
beta_binomial size 1 a 400 b 30 lS             n     4  worst 2.29e-13  median 2.3e-13  nonfinite 0
beta_binomial size 7 a 0.3 b 0.3 lF            n     5  worst 6.53e-16  median 5.6e-17  nonfinite 0
beta_binomial size 7 a 0.3 b 0.3 lS            n     5  worst 6.53e-16  median 1.7e-16  nonfinite 0
beta_binomial size 7 a 2 b 9 lF                n     5  worst 3.61e-15  median 2.2e-15  nonfinite 0
beta_binomial size 7 a 2 b 9 lS                n     5  worst 1.89e-15  median 8.7e-16  nonfinite 0
beta_binomial size 7 a 400 b 30 lF             n     5  worst 8.40e-14  median 1.3e-14  nonfinite 0
beta_binomial size 7 a 400 b 30 lS             n     5  worst 4.57e-14  median 5.5e-15  nonfinite 0
beta_binomial size 60 a 0.3 b 0.3 lF           n     5  worst 2.63e-14  median 2.3e-14  nonfinite 0
beta_binomial size 60 a 0.3 b 0.3 lS           n     5  worst 2.61e-14  median 2.3e-14  nonfinite 0
beta_binomial size 60 a 2 b 9 lF               n     5  worst 9.27e-14  median 3.2e-14  nonfinite 0
beta_binomial size 60 a 2 b 9 lS               n     5  worst 9.32e-14  median 2.6e-14  nonfinite 0
beta_binomial size 60 a 400 b 30 lF            n     5  worst 7.50e-13  median 3.6e-15  nonfinite 0
beta_binomial size 60 a 400 b 30 lS            n     5  worst 7.42e-13  median 7.4e-13  nonfinite 0
```
<!-- /GEN:acc -->

## 5. Gradients and third derivatives

`dev/cdfs-grad.R` (seed 7), GENERATED from `dev/cdfs-out/grad.txt`: the
AD gradient against Richardson finite differences (numDeriv) of an
independent reference where one exists, else of the same function on
doubles; and how many of 40 random points have a non-finite second or
third derivative.

<!-- GEN:grad -->
```
student lF, lS vs stats::pt                  grad rel err 7.2e-09  nonfinite 2nd 0/40  3rd 0/40
student at nu = 4 (start) vs stats::pt       grad rel err 7.2e-09  nonfinite 2nd 0/0  3rd 0/0
Gamma lF, lS vs stats::pgamma                grad rel err 1.3e-10  nonfinite 2nd 0/40  3rd 0/40
Beta lF, lS vs stats::pbeta                  grad rel err 5.9e-11  nonfinite 2nd 0/40  3rd 0/40
negbinomial vs stats::pnbinom                grad rel err 1.8e-10  nonfinite 2nd 0/40  3rd 0/40
negbinomial at shape = 2 (whole) vs stats    grad rel err 4.1e-10  nonfinite 2nd 0/0  3rd 0/0
binomial vs stats::pbinom                    grad rel err 1.0e-10  nonfinite 2nd 0/40  3rd 0/40
poisson vs stats::ppois                      grad rel err 4.7e-12  nonfinite 2nd 0/40  3rd 0/40
inverse.gaussian vs statmod::pinvgauss       grad rel err 1.1e-10  nonfinite 2nd 0/40  3rd 0/40
skew_normal (mu, log sigma, alpha)           grad rel err 3.1e-09  nonfinite 2nd 0/40  3rd 0/40
skew_normal at alpha = 0 (the start's stationary point) grad rel err 1.3e-07  nonfinite 2nd 0/40  3rd 0/40
exgaussian (mu, log sigma, log beta)         grad rel err 8.0e-09  nonfinite 2nd 0/40  3rd 0/40
asym_laplace (mu, log sigma, logit quantile) grad rel err 6.1e-09  nonfinite 2nd 0/40  3rd 0/40
von_mises (2 atan(mu), log kappa)            grad rel err 4.7e-08  nonfinite 2nd 0/40  3rd 0/40
shifted_lognormal (mu, log sigma, log ndt)   grad rel err 1.6e-10  nonfinite 2nd 0/40  3rd 0/40
beta_binomial (logit mu, log phi), 15 trials grad rel err 2.7e-08  nonfinite 2nd 0/40  3rd 0/40
compois (log mu, log nu)                     grad rel err 1.2e-06  nonfinite 2nd 0/40  3rd 0/40
```
<!-- /GEN:grad -->

## 6. Fits against brms: the censored rows' likelihood

`tests/testthat/test-brms-lik-cdfs.R` through `brms_lp_check()`: brms's
Stan program under flat priors, at frmtmb's estimates, against
`logLik(fit)` (difference below 1e-6 relative) and brms's gradient
there (below 1e-3), for every family with all four censoring codes and
for seven truncated fits. Since fix round 1 brms is given the same data
frame frmtmb is; through the first build a discrete right or interval
bound was handed to brms less one. GENERATED from
`dev/cdfs-out/brms-lik-report.txt`:

<!-- GEN:brmslik -->
```
lib: C:/Users/adf44/source/r/wt-cdfs-lib 
CASE order: student skew_normal exgaussian shifted_lognormal von_mises Gamma inverse.gaussian Beta | asym_laplace | negbinomial geometric binomial beta_binomial poisson | zero_inflated_poisson zero_inflated_negbinomial zero_inflated_binomial zero_inflated_beta_binomial hurdle_poisson hurdle_negbinomial hurdle_gamma hurdle_lognormal zero_inflated_beta | trunc: negbinomial student von_mises hurdle_gamma zero_inflated_poisson poisson+cens poisson-lb0
LPCHECK const 0 grad 1.69e-06 ours -227.6963549
LPCHECK const -1.051603249e-12 grad 0.000136 ours -201.3458129
LPCHECK const 0 grad 0.000122 ours -211.7995168
LPCHECK const 0 grad 4.29e-06 ours -125.291135
LPCHECK const 8.526512829e-14 grad 0.000642 ours -183.6775147
LPCHECK const -2.842170943e-13 grad 0.000217 ours -214.6374429
LPCHECK const -4.263256415e-13 grad 0.000281 ours -224.0987816
LPCHECK const 1.421085472e-14 grad 0.000754 ours -69.69797914
LPCHECK const 0 grad 0.666 ours -319.7001378
LPCHECK const 5.684341886e-14 grad 0.000136 ours -283.7450305
LPCHECK const 0 grad 4.45e-05 ours -292.7566341
LPCHECK const 3.410605132e-13 grad 5.12e-05 ours -243.1870764
LPCHECK const 1.250555215e-12 grad 0.000117 ours -283.6865441
LPCHECK const 0 grad 0.000208 ours -241.5133563
LPCHECK const 0 grad 2.43e-05 ours -301.5135401
LPCHECK const 0 grad 1.32e-05 ours -325.4245926
LPCHECK const 1.705302566e-13 grad 5.49e-05 ours -319.5067948
LPCHECK const 9.094947018e-13 grad 1.7e-05 ours -338.9213104
LPCHECK const -5.684341886e-14 grad 9.36e-06 ours -306.3750221
LPCHECK const 0 grad 0.000429 ours -329.6507124
LPCHECK const -2.842170943e-13 grad 0.182 ours -336.8666052
LPCHECK const 0 grad 7.73e-05 ours -266.807381
LPCHECK const -2.842170943e-14 grad 0.000471 ours -172.7121731
LPCHECK const 2.273736754e-13 grad 1.33e-05 ours -775.4797725
LPCHECK const -5.684341886e-14 grad 9.84e-05 ours -452.586882
LPCHECK const -3.183231456e-12 grad 1.27e-05 ours -452.9433707
LPCHECK const 2.273736754e-13 grad 0.0281 ours -610.8269786
LPCHECK const -1.136868377e-13 grad 6.07e-07 ours -695.8308777
LPCHECK const -1.136868377e-13 grad 0.000706 ours -646.9991728
LPCHECK const 0 grad 0.000428 ours -683.9881053
RESULT test-brms-lik-cdfs.R pass=70 fail=0 err=0 skip=0
```
<!-- /GEN:brmslik -->

compois is not in it: frmtmb's `compois()` takes the mean and brms's
`com_poisson()` the rate, so the two are different models of the same
data; its CDF is checked against the sum of frmtmb's own mass function
in `test-cdfs.R`.

## 7. Ported rows that flip

The distributions and log_lik blocks regenerated into
`dev/cdfs-out/bin2/` (`dev/cdfs-bin2-gen.R`, the tracked
`dev/brmsport-bin2/` left alone) and recorded on rellib-r10 and on this
lane's library; GENERATED from `dev/cdfs-out/bin2-flips.txt`:

<!-- GEN:flips -->
```
== distributions: 175 records, 43 moved
assert distributions:125 1   held FALSE -> TRUE | lane: 
assert distributions:150 1   held FALSE -> TRUE | lane: 
assert distributions:164 1   held FALSE -> TRUE | lane: 
assert distributions:169 1   held FALSE -> TRUE | lane: 
assert distributions:174 1   held FALSE -> TRUE | lane: 
assert distributions:179 1   held FALSE -> TRUE | lane: 
assert distributions:185 1   held FALSE -> TRUE | lane: 
assert distributions:195 1   held FALSE -> TRUE | lane: 
assert distributions:200 1   held FALSE -> TRUE | lane: 
assert distributions:205 1   held FALSE -> TRUE | lane: 
assert distributions:210 1   held FALSE -> TRUE | lane: 
assert distributions:52 1    held FALSE -> TRUE | lane: 
assert distributions:6 1     held FALSE -> TRUE | lane: 
assert distributions:64 1    held FALSE -> TRUE | lane: 
assert distributions:83 1    held FALSE -> TRUE | lane: 
setup distributions:124 1    held FALSE -> TRUE | lane: 
setup distributions:149 1    held FALSE -> TRUE | lane: 
setup distributions:163 1    held FALSE -> TRUE | lane: 
setup distributions:168 1    held FALSE -> TRUE | lane: 
setup distributions:173 1    held FALSE -> TRUE | lane: 
setup distributions:178 1    held FALSE -> TRUE | lane: 
setup distributions:184 1    held FALSE -> TRUE | lane: 
setup distributions:194 1    held FALSE -> TRUE | lane: 
setup distributions:199 1    held FALSE -> TRUE | lane: 
setup distributions:204 1    held FALSE -> TRUE | lane: 
setup distributions:209 1    held FALSE -> TRUE | lane: 
setup distributions:51 1     held FALSE -> TRUE | lane: 
setup distributions:62 1     held FALSE -> TRUE | lane: 
setup distributions:81 1     held FALSE -> TRUE | lane: 
value distributions:125 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 2.22e-16, max rel 5.93e-16
value distributions:150 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 5.55e-17, max rel 1.17e-16
value distributions:164 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 1.11e-16, max rel 1.46e-16
value distributions:169 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 1.11e-16, max rel 2.08e-16
value distributions:174 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 1.11e-16, max rel 1.05e-15
value distributions:179 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 2.22e-16, max rel 5.08e-16
value distributions:185 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 4.44e-16, max rel 8.25e-16
value distributions:195 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 2.22e-16, max rel 4.16e-16
value distributions:200 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 1.11e-16, max rel 4.16e-16
value distributions:205 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 4.16e-17, max rel 4.16e-16
value distributions:210 1    held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 4.16e-17, max rel 4.16e-16
value distributions:52 1     held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 1.5e-14, max rel 2.39e-11
value distributions:64 1     held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 9.06e-14, max rel 4.84e-11
value distributions:83 1     held FALSE -> TRUE | lane: HOLDS: 10 values (0 NA, 0 infinite), max abs diff 1.25e-16, max rel 5.8e-15
== log_lik: 261 records, 3 moved
assert log_lik:392 1         held FALSE -> TRUE | lane: 
setup log_lik:391 1          held FALSE -> TRUE | lane: 
value log_lik:392 1          held FALSE -> TRUE | lane: HOLDS: 15 values (0 NA, 0 infinite), max abs diff 3.5e-12, max rel 2.87e-13
```
<!-- /GEN:flips -->

All 15 `no-cdf` rows of `tests.distributions.R` (pstudent_t 6,
pvon_mises 52, pskew_normal 64, pexgaussian 83, pbeta_binomial 125,
pasym_laplace 150, pzero_inflated_poisson 164, _negbinomial 169,
_binomial 174, _beta_binomial 179, _beta 185, phurdle_poisson 195,
_negbinomial 200, _gamma 205, _lognormal 210) and `tests.log_lik.R:392`
(censored von_mises) now hold, each with a value check against brms's
own function that holds. Nothing else moved in either block.

## 8. Tests

GENERATED:

<!-- GEN:tests -->
```
ungated (dev/release/cdfs-fx-all-files): files 413: pass 25669 fail 0 err 0 skip 231 warn 0
gated_reference (dev/ciharden-log/cdfs-fx-ref-suite): files 413: pass 28632 fail 1 err 0 skip 0 warn 0
RESULT test-fuzz.R pass=1 fail=1 err=0 skip=0 warn=0
openblas326 (dev/ciharden-log/cdfs-fx-ob326-suite): files 18: pass 2939 fail 0 err 0 skip 0 warn 0
openblas332 (dev/ciharden-log/cdfs-fx-ob332-suite): files 18: pass 2939 fail 0 err 0 skip 0 warn 0
cens_trunc_ob326 (dev/ciharden-log/cdfs-fx-ct-ob0326-suite): files 1: pass 66 fail 0 err 0 skip 0 warn 0
cens_trunc_ob332 (dev/ciharden-log/cdfs-fx-ct-ob0332-suite): files 1: pass 66 fail 0 err 0 skip 0 warn 0
cens_trunc_ref (dev/ciharden-log/cdfs-fx-ct-ref-suite): files 1: pass 66 fail 0 err 0 skip 0 warn 0
base_r10_brmslik (dev/release/cdfs-base-brmslik-files): files 1: pass 0 fail 0 err 5 skip 0 warn 0
RESULT test-brms-lik-cdfs.R pass=0 fail=0 err=5 skip=0 warn=0
test-cdfs.R on rellib-r10 (dev/cdfs-out/base-test-cdfs.txt):
RESULT test-cdfs.R pass=2 fail=1 err=14 skip=0 warn=0
```
<!-- /GEN:tests -->

One session (`dev/release/run-sessions.sh`), GENERATED:

<!-- GEN:sessions -->
```
SESSION frmtmb forward files=246/256 pass=18490 fail=0 err=0 skip=196 warn=0
LEAK test-arg-refusal.R compat features: 1 -> 2 entries [with load of frmtmb.sample]
LEAK test-arg-refusal.R compat rules: 0 -> 1 entries [with load of frmtmb.sample]
LEAK test-arg-refusal.R prior defaults: 0 -> 1 entries [with load of frmtmb.sample]
LEAK test-families-extreme.R compat features: 2 -> 4 entries [with load of frmtmb.eam]
LEAK test-families-extreme.R compat rules: 1 -> 3 entries [with load of frmtmb.eam]
LEAK test-families-extreme.R aterm registry: dec ndt_group [with load of frmtmb.eam]
LEAK test-families-extreme.R frame checks: 0 -> 2 entries [with load of frmtmb.eam]
LEAK test-tmb-examples.R envvars: GDAL_DATA [with load of class classInt e1071 fmesher KernSmooth proxy sf units]
SESSION frmtmb reverse files=246/256 pass=18588 fail=0 err=0 skip=196 warn=0
LEAK test-tmb-examples.R envvars: GDAL_DATA [with load of class classInt e1071 fmesher KernSmooth proxy sf units]
LEAK test-families-extreme.R compat features: 1 -> 3 entries [with load of frmtmb.eam]
LEAK test-families-extreme.R compat rules: 0 -> 2 entries [with load of frmtmb.eam]
LEAK test-families-extreme.R aterm registry: dec ndt_group [with load of frmtmb.eam]
LEAK test-families-extreme.R frame checks: 0 -> 2 entries [with load of frmtmb.eam]
LEAK test-arg-refusal.R compat features: 3 -> 4 entries [with load of frmtmb.sample]
LEAK test-arg-refusal.R compat rules: 2 -> 3 entries [with load of frmtmb.sample]
LEAK test-arg-refusal.R prior defaults: 0 -> 1 entries [with load of frmtmb.sample]
SESSION frmtmb preload files=246/256 pass=18684 fail=0 err=0 skip=197 warn=0
LEAK test-tmb-examples.R envvars: GDAL_DATA [with load of class classInt e1071 fmesher KernSmooth proxy sf units]
SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
LEAK test-newdata-dec.R compat features: 2 -> 3 entries [with load of frmtmb.sample]
LEAK test-newdata-dec.R compat rules: 2 -> 3 entries [with load of frmtmb.sample]
LEAK test-newdata-dec.R prior defaults: 0 -> 1 entries [with load of frmtmb.sample]
SESSION frmtmb.eam reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
LEAK test-sampling.R compat features: 2 -> 3 entries [with load of frmtmb.sample]
LEAK test-sampling.R compat rules: 2 -> 3 entries [with load of frmtmb.sample]
LEAK test-sampling.R prior defaults: 0 -> 1 entries [with load of frmtmb.sample]
SESSION frmtmb.eam preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent forward files=11/11 pass=368 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent reverse files=11/11 pass=368 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent preload files=11/11 pass=368 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
LEAK test-records.R envvars: rxTempDir [with load of backports cachem checkmate data.table dparser dplyr farver fastmap ggplot2 gtable lotri memoise pillar pkgconfig PreciseSums RColorBrewer RcppParallel rxode2 S7 scales sys tibble tidyselect vctrs]
SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
LEAK test-records.R envvars: rxTempDir [with load of backports cachem checkmate codetools data.table dparser dplyr farver fastmap ggplot2 glue gtable lotri MASS memoise pillar pkgconfig PreciseSums RColorBrewer Rcpp RcppParallel RTMB rxode2 S7 scales sys tibble tidyselect vctrs]
SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
LEAK test-records.R envvars: rxTempDir [with load of backports cachem checkmate data.table dparser dplyr farver fastmap ggplot2 gtable lotri memoise pillar pkgconfig PreciseSums RColorBrewer RcppParallel rxode2 S7 scales sys tibble tidyselect vctrs]
SESSION frmtmb.sample forward files=56/61 pass=2849 fail=0 err=0 skip=8 warn=0
LEAK test-compat-preflight.R compat features: 1 -> 6 entries [with load of frmtmb.eam frmtmb.latent frmtmb.ode]
LEAK test-compat-preflight.R compat rules: 1 -> 6 entries [with load of frmtmb.eam frmtmb.latent frmtmb.ode]
LEAK test-compat-preflight.R aterm registry: 0 -> 2 entries [with load of frmtmb.eam frmtmb.latent frmtmb.ode]
LEAK test-compat-preflight.R frame checks: 0 -> 3 entries [with load of frmtmb.eam frmtmb.latent frmtmb.ode]
SESSION frmtmb.sample reverse files=56/61 pass=2849 fail=0 err=0 skip=8 warn=0
LEAK test-simulators.R compat features: 1 -> 3 entries [with load of frmtmb.latent]
LEAK test-simulators.R compat rules: 1 -> 3 entries [with load of frmtmb.latent]
LEAK test-compat-preflight.R compat features: 3 -> 6 entries [with load of frmtmb.eam frmtmb.ode]
LEAK test-compat-preflight.R compat rules: 3 -> 6 entries [with load of frmtmb.eam frmtmb.ode]
LEAK test-compat-preflight.R aterm registry: 0 -> 2 entries [with load of frmtmb.eam frmtmb.ode]
LEAK test-compat-preflight.R frame checks: 0 -> 3 entries [with load of frmtmb.eam frmtmb.ode]
SESSION frmtmb.sample preload files=56/61 pass=2849 fail=0 err=0 skip=8 warn=0
SESSION frmtmb.spline forward files=16/16 pass=594 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.spline reverse files=16/16 pass=594 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.spline preload files=16/16 pass=594 fail=0 err=0 skip=1 warn=0
CDFS-FX-SESS ran 24 of 24
```
<!-- /GEN:sessions -->

The sessions ran with NOT_CRAN only. Every LEAK line above names a
registry entry an extension's load adds, or an environment variable
that sf or rxode2 sets; none names a cdfs file.

The gated fuzz tier is red, with three REAL-NEW findings. None is a
CDF. The five families that gained a CDF joined the cens and trunc
columns of the fuzz grid, which reshuffles the pairwise plan, and the
new plan draws three specs the old one did not. Each fails the same
way on rellib-r10 (`dev/cdfs-fuzz-base.R`), GENERATED:

<!-- GEN:fuzzbase -->
```
-- dev/cdfs-out/fuzz-base-r10.txt (rellib-r10)
lib: 0.72.0 C:/Users/adf44/source/r/rellib-r10/frmtmb 
specs: 3 
real_new | [REAL-NEW] summary_prints
    spec: family=cumulative aterm=weights re=nested special=none dpar=none mode=quadrature op=vcov (seed 20322452, cover)
    summary() errored: no 'dimnames' attribute for array
real_new | [REAL-NEW] fit_error
    spec: family=weibull aterm=trunc re=diagbar special=none dpar=none mode=ml op=confint (seed 20446531, pad)
    frm() failed: The optimizer failed on this model (weibull, ML, nlminb): NA/NaN gradient evaluation. The likelihood was undefined or unbounded somewhere the optimizer stepped. Refit with verbose = TRUE to see which stage broke, try another optimizer (frmtmb_control(optimizer = "optim")), or start the fit nearer the optimum with the `start` argument of frm()
real_new | [REAL-NEW] fit_error
    spec: family=Gamma aterm=weights re=ar1 special=mo_int dpar=dpar_x mode=profile op=vcov (seed 20502220, pad)
    frm() failed: The optimizer failed on this model (Gamma, ML, profile, nlminb): NA/NaN gradient evaluation. The likelihood was undefined or unbounded somewhere the optimizer stepped. Refit with verbose = TRUE to see which stage broke, try another optimizer (frmtmb_control(optimizer = "optim")), or start the fit nearer the optimum with the `start` argument of frm()
-- dev/cdfs-out/fuzz-base-mine.txt (wt-cdfs-lib)
lib: 0.72.0 C:/Users/adf44/source/r/wt-cdfs-lib/frmtmb 
specs: 3 
real_new | [REAL-NEW] summary_prints
    spec: family=cumulative aterm=weights re=nested special=none dpar=none mode=quadrature op=vcov (seed 20322452, cover)
    summary() errored: no 'dimnames' attribute for array
real_new | [REAL-NEW] fit_error
    spec: family=Gamma aterm=weights re=ar1 special=mo_int dpar=dpar_x mode=profile op=vcov (seed 20502220, pad)
    frm() failed: The optimizer failed on this model (Gamma, ML, profile, nlminb): NA/NaN gradient evaluation. The likelihood was undefined or unbounded somewhere the optimizer stepped. Refit with verbose = TRUE to see which stage broke, try another optimizer (frmtmb_control(optimizer = "optim")), or start the fit nearer the optimum with the `start` argument of frm()
```
<!-- /GEN:fuzzbase -->

The first fuzz rerun also had 13 generator findings: a zero-heavy
count response has its 15% quantile at its minimum, so no row was
left-censored. `fuzz_data()` now moves a lattice bound one step in
from the minimum and the maximum; the second rerun has none.

R CMD check --as-cran, core (`dev/cdfs-check.sh`), GENERATED:

<!-- GEN:check -->
```
* checking HTML version of manual ... [21s] NOTE
Status: 1 NOTE
```
<!-- /GEN:check -->

## 9. Decisions

- The discrete censoring convention is brms 2.23.0's, by the user's
  decision in fix round 1 (section 11). Through the first build it was
  frmtmb's inclusive one, as for poisson and discrete_weibull since
  0.54.0.
- `$lcdf` stays a probability (core's contract). Core reads both log
  tails through the new `cdf_log` slot where a family declares it
  (fix round 1), so the floor at 1e-308 remains only for the families
  without one.
- compois's survivor is the complement of its summed CDF, as brms's
  `com_poisson_lccdf()` is, until that is below about 1e-3, and an
  upper sum over a window whose length is data beyond (fix round 1).
- brms's `zero_inflated_asym_laplace_lccdf()` at `y = 0` includes the
  point mass (brms-31); frmtmb's is `P(Y > y)` everywhere.
- Truncated means only where brms has one; elsewhere `fitted()` refuses
  a truncated response by name, as brms does.

## 10. Defects found and not fixed

- `asym_laplace()` with a free quantile stops at nlminb's false
  convergence on 160 rows drawn from it with no censoring, on 0.72.0 and
  on this build alike (`dev/cdfs-dbg14.R`); the likelihood has a kink at
  every observation.
- `trunc()` did not check that the responses lie in the window
  (pre-existing; fixed here because a hurdle fit then failed with a NaN
  gradient instead of a refusal).
- Upstream: `dev/upstream-bugs.md` RTMB-8, RTMB-9, RTMBdist-2,
  RTMBdist-3, RTMBdist-4, brms-31, brms-32, brms-33, rstan-3.
- Three fuzz specs failed on 0.72.0 as on the first build (section 8).
  The weibull `trunc()` one was finding #4 of the review and fits now;
  the other two are in section 12.

## 11. Fix round 1

The review (`dev/reviews/2026-10-20-cdfs.md`, MERGEABLE AFTER FIXES)
and the coordinator's fix list. Findings are numbered as in the review
file. Finding 6 was decided by the user: match brms 2.23.0's discrete
conventions exactly. The pre-fix lane build is the reviewer's
`wt-cdfs-rev-lib` (same source, installed after the last `R/` edit), and
every new or changed test was run on it first.

### What changed, by finding

- **#1 `log_igamma_two()` NaN past shape 12,600.** Both quadrature sums
  in `igamma_gl_tails()` cap the integrand's exponent at 600, active
  only where the arm is weighted zero. Measured over shapes 0.01 to 1e7
  against `stats::pgamma()` (`dev/cdfs-fix1-igamma.R`), and the two
  fits reach the independent maximum (GEN:fix1igamma below). Tests:
  `test-cdfs.R` gamma grid to 1e6, the censored Gamma at shape 1e5, and
  poisson right censoring at 75000.
- **#2 zi and hurdle row sizing.** `zi_cdf_log()` and `hurdle_cdf_log()`
  size rows by the longest of `q`, the base pair and the dpars, and
  recycle `q` and `base_zero`; `row_lpdf()` also recycles every edge to
  the rows. Test: every zi and hurdle family, the bound as a constant,
  as a column and per row, against the hand-built likelihood at the
  same parameters (`test-cdfs.R`).
- **#3 interval cancellation; and #4 the `trunc(lb = )` normalizer.**
  A family can now declare `cdf_log` (both tails on the log scale; new
  argument of `frmtmb_family()`), and for such a family `row_lpdf()`
  reads every censored row and every truncation window through one
  function, `log_window()`: the log mass of `(lo, hi]`, open where an
  edge is absent. Both edges present: from
  `S` where `S(lo) >= 1/2`, from `F` below (`log_window_two()`); one
  edge: `lS(lo)` or `lF(hi)`. So `trunc(lb = )` alone is `lccdf(lb')`,
  as brms. Point-mass families carry their base pair and weight, and a
  window above zero is the base's window times the weight (exact for a
  hurdle `(0, y2]`); von Mises carries its own arc mass (`lmass`),
  because an arc on the far side of the mean is below both tails.
  `log1m_exp_ad()` switches to `log1p(-exp(l))` below `-log 2`
  (Maechler): `log(-expm1(l))` there had kept only absolute accuracy,
  which was the 20 to 45 percent gradient error of the review's grad5.
  A family WITHOUT `cdf_log` is read by `row_lpdf_prob()`, the 0.72.0
  probability-scale consumer with each slot called where it was. The
  first version of this fix read every family through log pairs built
  from `lcdf` and `lccdf`, and the full run caught what that broke:
  frmtmb.eam's Wiener `lcdf` is a boundary's defective distribution
  function on left- and interval-censored rows, not the complement of
  its `lccdf`, and its slots read the censoring codes to know which
  rows a call is for, so `test-wiener-cdf.R` "all four censoring
  codes" was 85 relative off the hand-written likelihood. The legacy
  path gains #4 too: `trunc(lb = )` with no upper bound divides by
  `lccdf(lb')`, called without the censoring codes (a Wiener slot fills
  the uncensored rows with a neutral value when it sees them).
  exponential, weibull, cox, frechet, gen_extreme_value and
  discrete_weibull stay on it.
  Tests: `test-cens-windows.R` (new): upper-tail intervals against
  300-bit Rmpfr quadrature, the von Mises far side, the hurdle base
  mass, AD against Richardson differences, `trunc(lb = )` against
  stats' own log tails, the weibull fuzz spec 20446531, and a custom
  family's `cdf_log`.
- **#5 exgaussian.** `exgauss_lx()` takes `log phi(z) + log M(u - z)`
  with Mills's ratio (`log_mills()`: direct below 6, Laplace's
  continued fraction above 7) where `u - z >= 1`, and the written form
  below `-1`, blended; the CDF's `lF` takes `log M(u - z) - log M(-z)`
  in the Mills region. The family's density is now the closed form
  (`exgauss_lx() - log(beta)`), not `RTMBdist::dexgauss()` (filed as
  RTMBdist-4). frmtmb's own likelihood WAS affected (GEN:fix1exg):
  0.63 off on 400 rows at `sigma / beta = 50`, and three of four
  near-gaussian fits stopped at `sigma / beta = 20` below the maximum.
  Tests: the density against 300-bit Rmpfr at `u` 5 to 1e3, and the
  CDF's derivatives at `u` 1e3 and 1e4 against 400-bit differences.
- **#6 the discrete convention (decided: brms).** `row_lpdf()` applies
  no censoring shift; only `trunc_edges()` shifts `lb` to `lb - 1`, and
  `lb <= 0` on a count is absent. The once-per-session notice is gone;
  `frame.R` refuses a fractional `cens()` or `trunc()` bound on a count
  and an interval with `y2 <= y` on any family. The size of the change
  (GEN:fix1conv): 30.2 in log-likelihood and 5.5 percent in `lambda` on
  200 poisson draws recorded as "6 or more". Changed in step: the tests
  (`test-cdfs.R` `hand_ll()`, `test-cens-trunc.R`, the shim removed
  from `test-brms-lik-cdfs.R`, `test-families-extreme.R`,
  `test-review-v25.R`, `test-v07.R`, `test-bcm-data-analysis.R`), the
  compat rows (discrete `cens()` and the poisson and discrete_weibull
  `trunc()` rows are "works"), `?frmtmb_family`, `?discrete_weibull`,
  the OSA refusal text, `vignette("brms-migration")`,
  `vignette("compatibility")`, the BCM vignette and `inst/bcm` (the
  ChaSaSoon band is now the row `z = 14`, `(14, 25]`), and the bin-2
  shim `brms_g1_pd()` (P(Y > q) is now the row right-censored at q).
- **#7 compois survivor.** Where the complement is below about 1e-3,
  the survivor is the upper sum over 100 terms past `q` plus a
  geometric remainder at the last term's ratio.
- **#8 hurdle `F(0)`.** `F` is formed directly where `S >= 1/2`, and at
  `q = 0` (data) `F = hu` and `S = 1 - hu` are set exactly.
- **#9 below 1e-308.** The new `cdf_log` slot of `frmtmb_family()`
  carries both tails on the log scale; all the lane's families and
  gaussian, lognormal and poisson declare it, so a left-censored row
  keeps its log probability. A beta-binomial log CDF is held at or
  below zero (it passed by 3e-10 at shapes 1e5).
- **#10 von Mises near +-pi.** A short arc (log density changing by
  less than about 4 across it) is a 20-point quadrature of the density
  over the arc; the comment is corrected. Student at `nu = 1e12` is
  unchanged (3.6e-9 at `z = 50`, an input that does not occur).
- **#11 the BREAKING check.** An infinite bound is absent (no
  `Inf * 0`); `trunc(lb = 0)` on a count is accepted; tests for the
  absent case (a response on either bound fits) and for `lb >= ub`.
- **Cost (found by the one-session run, not by the review).** With
  `cdf_log` on poisson, a zero-truncated poisson's normalizer went
  through the incomplete gamma on every row, about 170 times
  `RTMB::ppois()` per row, and the reverse one-session run of core hit
  `std::bad_alloc` in `test-importance.R` ("trunc() and cens() are
  corrected"). Measured on that fit with 500 importance draws: 428
  seconds and a 27 GB peak working set (transcribed: the log was
  overwritten by the rerun on the fixed build), against 5.4 seconds and 0.75
  GB on rellib-r10 (`dev/cdfs-fix1-imp.R`). `poisson_cdf_log()` now
  takes each row's cheapest exact form by its count (data): closed
  form at 0, finite sums of the masses to a count of 30, the incomplete
  gamma above; `nbinom_cdf_log()` takes the closed form at 0. The same
  fit is now 10.6 seconds and 1.3 GB (GEN:fix1cost below).
- **#12 NEWS.** The claims are rewritten from the logs; the rstan-3
  sentence says where the gradient was measured.
- **#13 NITs.** "Neighboring", "centered"; the refusal names
  `frm_compat()` for the term used; RTMB-8 says 0 or NaN by case
  (`dev/cdfs-rev-out/upstream.txt`); RTMB-9 has the review's second
  derivatives (52 of 500); RTMBdist-4 is filed; rstan-3 says the review
  did not reproduce it.

### Incident: another lane's R CMD check stopped

At about 19:51 on 2026-10-20 I stopped this lane's R CMD check of core
(restarted after adding Suggests) by matching process command lines on
`frmtmb_0.72.0.tar.gz`. That also matched lane fuzzfix's check of core
in `frmtmb-wt-fuzzfix/dev/fuzzfix-check/`, which was in its tests step:
its `check-frmtmb.txt` stops at "Running 'testthat.R'" with no Status
line, and its script went on to frmtmb.sample. That check needs a
rerun. Processes are now matched on this worktree's own paths only.

### Files touched in fix round 1

- R: `R/objective.R` (`fam_log_pair()`, `log_window()`,
  `log_window_two()`, `trunc_edges()`, `row_lpdf()`, `row_lpdf_prob()`),
  `R/families-cdf.R`, `R/families.R` (`frmtmb_family(cdf_log =)`, the
  families' `cdf_log` slots, the exgaussian density, the documentation
  sections), `R/frame.R`, `R/compat.R`, `R/predict.R` (the OSA refusal
  text), `R/families-extreme.R` (the discrete_weibull text); `man/`
  regenerated (`frmtmb_family.Rd`, `frmtmb-extreme-families.Rd`,
  `residuals.frmtmb_fit.Rd`).
- Tests: `test-cens-windows.R` (new); `test-cdfs.R`,
  `test-cens-trunc.R`, `test-brms-lik-cdfs.R`, `test-cens-lccdf.R`,
  `test-trunc-postfit.R` (a relative tolerance on the mean residual,
  `dev/cdfs-fix1-tpois.R`), `test-v07.R`, `test-review-v25.R`,
  `test-families-extreme.R`, `test-bcm-data-analysis.R`,
  `test-compat.R`, `test-hurdle-nb-zoib.R`,
  `test-xbeta-zibb-hurdle-cum.R`.
- `DESCRIPTION`: extraDistr and Rmpfr in Suggests, which the new
  tests reach through `::` (R CMD check warned).
- Docs: NEWS, `vignettes/brms-migration.Rmd`,
  `vignettes/compatibility.Rmd`,
  `vignettes/bayesian-cognitive-modeling.Rmd`,
  `inst/bcm/binomial-extras.R`.
- dev: `dev/upstream-bugs.md`, `dev/brmsport-bin2-adapters.R`
  (`brms_g1_pd()`), `dev/cdfs-findings-build.R`,
  `dev/cdfs-brms-lik-report.R`, `dev/cdfs-fix1-*.R` (new).

### Evidence

`dev/cdfs-fix1-igamma.R`, GENERATED:

<!-- GEN:fix1igamma -->
```
lib: C:/Users/adf44/source/r/wt-cdfs-lib/frmtmb 
== accuracy: relative error of the log, floored at one
a 0.01     points  25 NaN 0 worst 4.78e-14
a 0.1      points  25 NaN 0 worst 9.19e-15
a 0.5      points  25 NaN 0 worst 2.38e-15
a 1        points  24 NaN 0 worst 5.33e-15
a 2        points  26 NaN 0 worst 7.24e-16
a 5        points  27 NaN 0 worst 2.18e-14
a 12       points  27 NaN 0 worst 5.72e-14
a 30       points  28 NaN 0 worst 7.17e-14
a 60       points  28 NaN 0 worst 1.12e-14
a 80       points  28 NaN 0 worst 2.96e-14
a 100      points  23 NaN 0 worst 1.54e-14
a 300      points  29 NaN 0 worst 3.26e-15
a 1000     points  30 NaN 0 worst 4.46e-15
a 3000     points  31 NaN 0 worst 1.13e-14
a 10000    points  24 NaN 0 worst 1.02e-14
a 12500    points  31 NaN 0 worst 6.71e-15
a 13000    points  31 NaN 0 worst 1.18e-14
a 30000    points  31 NaN 0 worst 7.04e-15
a 100000   points  31 NaN 0 worst 7.31e-15
a 300000   points  31 NaN 0 worst 1.29e-14
a 1e+06    points  29 NaN 0 worst 2.47e-14
a 1e+07    points  31 NaN 0 worst 8.04e-14
ALL: NaN 0 worst 8.04e-14
== gradient finite at a = 1e5, x = 0.7 a and 0.3 a
x/a 0.3: -120398.21 70000.429
x/a 0.7: -35670.327 30002.333
== censored Gamma at shape 1e5 (review's data)
frmtmb logLik 539.69911610 conv 0 | independent 539.69911605 | diff 5.24e-08
== poisson right censored at 75000 (counts near 1e5), brms's P(Y > k)
censored rows: 63 
frmtmb logLik -1644.674816 conv 0 | independent -1644.674816 | diff 6.02e-09
q 100000 mu 70000 lS: -5673.323 ref -5673.323 
END
```
<!-- /GEN:fix1igamma -->

The cost of the pairs per row (`dev/cdfs-fix1-cost.R`) and the
importance-sampled zero-truncated poisson (`dev/cdfs-fix1-imp.R`) on
rellib-r10 and on this build, GENERATED:

<!-- GEN:fix1cost -->
```
lib: C:/Users/adf44/source/r/wt-cdfs-lib/frmtmb 
RTMB::ppois                                  3.00 ms per sweep
poisson_cdf_log lS                         207.00 ms per sweep
poisson_cdf_log lF                         244.50 ms per sweep
log_igamma_two at a = q + 1                422.00 ms per sweep
pnorm both tails                             7.00 ms per sweep
nbinom_cdf_log lS                          677.00 ms per sweep
gamma_cdf_log lS                           419.00 ms per sweep
END
-- rellib-r10
lib: C:/Users/adf44/source/r/rellib-r10/frmtmb 
importance 500: objective 448.721461 seconds 5.4
laplace: objective 448.856654 seconds 0.1
peak memory (MB, memory.size): 748.3984375 
END
-- this build
lib: C:/Users/adf44/source/r/wt-cdfs-lib/frmtmb 
importance 500: objective 448.721461 seconds 4.3
laplace: objective 448.856654 seconds 0.0
peak memory (MB, memory.size): 764.84375 
END
```
<!-- /GEN:fix1cost -->

`dev/cdfs-fix1-conv.R` (finding 6), GENERATED:

<!-- GEN:fix1conv -->
```
lib: C:/Users/adf44/source/r/wt-cdfs-lib/frmtmb 
censored rows 50
this build: lambda 4.4005 logLik -379.3214
inclusive reading: lambda 4.1722 logLik -349.0793
logLik difference 30.24 (0.605 per censored row), lambda 5.47 percent
6 or more written as cens at 5: lambda 4.1722 logLik -349.0793
END
```
<!-- /GEN:fix1conv -->

`dev/cdfs-fix1-exgfit.R` on rellib-r10 and on this build (finding 5),
GENERATED:

<!-- GEN:fix1exg -->
```
-- rellib-r10
lib: C:/Users/adf44/source/r/rellib-r10/frmtmb 
u    10: -objective -574.563430 closed form -574.563430 diff 6.821e-13
u    20: -objective -574.621752 closed form -574.621752 diff 2.046e-12
u    50: -objective -575.265651 closed form -574.638056 diff -0.6276
u   100: -objective -574.933330 closed form -574.640286 diff -0.293
u  1000: -objective -574.668442 closed form -574.641002 diff -0.02744
seed 1: logLik -590.764386 closed form -590.764386 diff 0 | u 20 | b 1.00948 0.44768 | conv 1
seed 2: logLik -562.727235 closed form -562.727235 diff 4.55e-13 | u 20 | b 0.98137 0.58878 | conv 1
seed 3: logLik -569.787859 closed form -569.787859 diff 3.41e-12 | u 20 | b 1.04982 0.49075 | conv 1
seed 4: logLik -584.055852 closed form -584.055852 diff -3.41e-13 | u 2.514 | b 1.01710 0.48452 | conv 0
-- this build
lib: C:/Users/adf44/source/r/wt-cdfs-lib/frmtmb 
u    10: -objective -574.563430 closed form -574.563430 diff 1.478e-12
u    20: -objective -574.621752 closed form -574.621752 diff 6.253e-12
u    50: -objective -574.638056 closed form -574.638056 diff -5.457e-12
u   100: -objective -574.640286 closed form -574.640286 diff 1.357e-10
u  1000: -objective -574.641002 closed form -574.641002 diff 4.218e-09
seed 1: logLik -590.763594 closed form -590.763594 diff 1.64e-09 | u 421.7 | b 1.00948 0.44768 | conv 0
seed 2: logLik -562.725616 closed form -562.725616 diff 1.63e-09 | u 589.8 | b 0.98137 0.58878 | conv 0
seed 3: logLik -569.785480 closed form -569.785480 diff 6.35e-09 | u 695.1 | b 1.04982 0.49075 | conv 0
seed 4: logLik -584.055852 closed form -584.055852 diff 1.14e-13 | u 2.514 | b 1.01710 0.48452 | conv 0
```
<!-- /GEN:fix1exg -->

The review's own scripts rerun on this build (`dev/cdfs-out/rerev/`,
the library path changed and nothing else), GENERATED:

<!-- GEN:fix1rerev -->
```
-- acc.txt (the review's tail grid): comparisons, and at most 1e-12
582 comparisons, 554 at most 1e-12; the rest:
student nu=1e+08 lF (ref)                    n   10 worst 1.22e-12 at z=-50 (ref -1254.82 got -1254.82) nonfinite 0
student nu=1e+08 lS (ref)                    n   10 worst 1.22e-12 at z=50 (ref -1254.82 got -1254.82) nonfinite 0
student nu=1e+12 lF (ref)                    n   10 worst 3.61e-09 at z=-50 (ref -1254.83 got -1254.83) nonfinite 0
student nu=1e+12 lS (ref)                    n   10 worst 3.61e-09 at z=50 (ref -1254.83 got -1254.83) nonfinite 0
Gamma shape=0.0001 lS                        n   23 worst 2.29e-12 at x=0.4001 (ref -9.56382 got -9.56382) nonfinite 0
nbinom shape=7.5 mu=100 lS                   n   10 worst 2.04e-10 at q=10500 (ref -724.085 got -724.085) nonfinite 0
nbinom shape=7.5 mu=10000 lS                 n   10 worst 2.25e-02 at q=1000500 (ref -698.841 got -714.591) nonfinite 0
nbinom shape=7.5 mu=1e+06 lS                 n   10 worst 6.85e-02 at q=100000500 (ref -668.697 got -714.496) nonfinite 0
binom size=1e+06 p=0.05 lF                   n    9 worst 8.41e-12 at q=48910 (ref -15.1462 got -15.1462) nonfinite 0
binom size=1e+06 p=0.95 lS                   n    9 worst 8.35e-12 at q=951089 (ref -15.1462 got -15.1462) nonfinite 0
skew_normal alpha=-1000 lF                   n   11 worst 1.18e-12 at q=-1 (ref -1.82446 got -1.82446) nonfinite 0
skew_normal alpha=-50 lF                     n   14 worst 1.18e-12 at q=-1 (ref -1.82455 got -1.82455) nonfinite 0
skew_normal alpha=50 lS                      n   14 worst 1.18e-12 at q=1 (ref -1.82455 got -1.82455) nonfinite 0
skew_normal alpha=1000 lS                    n   11 worst 1.18e-12 at q=1 (ref -1.82446 got -1.82446) nonfinite 0
exgauss sigma=0.001 beta=1000 lF             n    7 worst 3.90e-12 at q=-1000 (ref -14.7345 got -14.7345) nonfinite 0
bbinom n=1 (100000, 100000) lF               n    1 worst 2.04e-10 at q=0 (ref -0.693147 got -0.693147) nonfinite 0
bbinom n=1 (100000, 100000) lS               n    1 worst 2.04e-10 at q=0 (ref -0.693147 got -0.693147) nonfinite 0
bbinom n=2 (100000, 100000) lF               n    2 worst 1.98e-10 at q=0 (ref -1.38629 got -1.38629) nonfinite 0
bbinom n=2 (100000, 100000) lS               n    2 worst 1.98e-10 at q=1 (ref -1.38629 got -1.38629) nonfinite 0
bbinom n=50 (100000, 100000) lF              n    5 worst 3.59e-10 at q=49 (ref -1.8447e-11 got -3.77079e-10) nonfinite 0
bbinom n=50 (100000, 100000) lS              n    5 worst 3.79e-10 at q=25 (ref -0.812225 got -0.812225) nonfinite 0
bbinom n=1000 (2, 9) lF                      n    5 worst 1.23e-12 at q=500 (ref -0.0110223 got -0.0110223) nonfinite 0
bbinom n=1000 (2, 9) lS                      n    5 worst 1.23e-12 at q=0 (ref -8.8318e-05 got -8.8318e-05) nonfinite 0
bbinom n=1000 (100000, 100000) lF            n    5 worst 3.00e-10 at q=500 (ref -0.668296 got -0.668296) nonfinite 0
bbinom n=1000 (100000, 100000) lS            n    5 worst 3.21e-10 at q=500 (ref -0.718631 got -0.718631) nonfinite 0
compois mu=30 nu=0.3 lF                      n    5 worst 1.47e-12 at q=100 (ref -8.43099e-09 got -8.43246e-09) nonfinite 0
compois mu=30 nu=0.3 lS                      n    5 worst 3.20e-12 at q=30 (ref -0.782485 got -0.782485) nonfinite 0
compois mu=30 nu=1 lS                        n    5 worst 1.23e-12 at q=30 (ref -0.794851 got -0.794851) nonfinite 0
-- zitrunc4.txt (finding 2)
zero_inflated_poisson trunc(ub = 30)           -fn: constant bound -127.88216946 | column bound -127.88216946 | diff 0
zero_inflated_poisson trunc(ub = 10)           -fn: constant bound -127.83716240 | column bound -127.83716240 | diff 0
zero_inflated_negbinomial trunc(ub = 30)       -fn: constant bound -137.87544592 | column bound -137.87544592 | diff 0
zero_inflated_binomial trunc(ub = 9)           -fn: constant bound -114.07365959 | column bound -114.07365959 | diff 0
zero_inflated_beta_binomial trunc(ub = 9)      -fn: constant bound -128.88956275 | column bound -128.88956275 | diff 0
zero_inflated_poisson trunc(lb = 1)            -fn: constant bound -119.97539844 | column bound -119.97539844 | diff 0
zero_inflated_beta trunc(ub = 0.95)            -fn: constant bound -21.47648529 | column bound -21.47648529 | diff 0
hurdle_poisson trunc(ub = 10)                  -fn: constant bound -127.71752123 | column bound -127.71752123 | diff 0
poisson trunc(ub = 10)                         -fn: constant bound -129.46141839 | column bound -129.46141839 | diff 0
== hand-rolled ZIP trunc(ub = 10) at the column fit's start
hand -127.83716240  column -127.83716240
-- trlb.txt (finding 4)
gaussian mean far below lb: fn 1491.728  gradient  -154.6309 -3259.9632
student  mean far below lb: fn 63.72416  gradient  -3.2816408   0.4739358 -26.7694548
Gamma    mean far below lb: fn 16366.16  gradient -18364.84  16617.70
-- fuzzweib.txt (finding 4, the unpatched spec)
as is      seed 20446531 family=weibull aterm=trunc re=diagbar special=none dpar=none mode=ml op=confint
   fits: logLik -145.969682 conv 0
-- grad5.txt (finding 3)
skew_normal    codes kept all  (n  71): max rel |AD - FD| 3.68e-11  fn 557.2365
skew_normal    codes kept -1   (n  21): max rel |AD - FD| 3.02e-11  fn 531.6824
skew_normal    codes kept 1    (n  31): max rel |AD - FD| 4.91e-11  fn 547.26442
skew_normal    codes kept 2    (n  19): max rel |AD - FD| 3.11e-11  fn 544.10947
hurdle_gamma   codes kept all  (n  85): max rel |AD - FD| 7.18e-09  fn 982.24152
hurdle_gamma   codes kept -1   (n  19): max rel |AD - FD| 2.11e-09  fn 926.01277
hurdle_gamma   codes kept 1    (n  33): max rel |AD - FD| 4.28e-09  fn 949.67885
hurdle_gamma   codes kept 2    (n  33): max rel |AD - FD| 1.18e-09  fn 1274.3034
von_mises      codes kept all  (n  83): max rel |AD - FD| 1.82e-08  fn 1876.5121
von_mises      codes kept -1   (n  19): max rel |AD - FD| 1.07e-08  fn 2000.0035
von_mises      codes kept 1    (n  22): max rel |AD - FD| 1.15e-08  fn 2062.3439
von_mises      codes kept 2    (n  42): max rel |AD - FD| 1.29e-08  fn 2020.2911
-- grad4.txt (finding 5)
mu -2 log sigma 0 log beta -10 (u = 2.2e+04): worst rel err 8.27e-08
mu 2 log sigma 0 log beta -10 (u = 2.2e+04): worst rel err 8.27e-08
mu -2 log sigma 3 log beta -10 (u = 4.42e+05): worst rel err NaN
mu -2 log sigma 6 log beta -10 (u = 8.89e+06): worst rel err NaN
mu -2 log sigma 10 log beta 0 (u = 2.2e+04): worst rel err 5.39e-07
mu 2 log sigma 0 log beta -5 (u = 148): worst rel err 6.08e-09
mu -2 log sigma 0 log beta -7 (u = 1.1e+03): worst rel err 6.08e-09
mu 0 log sigma 0 log beta -3 (u = 20.1): worst rel err 6.08e-09
mu -2 log sigma 2 log beta -2 (u = 54.6): worst rel err 6.08e-09
alpha -10: worst rel err 3.38e-10 at q = 1.5
alpha -3: worst rel err 2.54e-10 at q = 1.7
alpha 3: worst rel err 1.18e-10 at q = 5.7
alpha 10: worst rel err 1.69e-10 at q = -5.9
-- exgdens2.txt (finding 5)
sigma / beta 1      max |lpdf - brms| 0.000e+00  at y = -3
sigma / beta 2      max |lpdf - brms| 4.441e-16  at y = 2.5
sigma / beta 5      max |lpdf - brms| 4.219e-15  at y = -1
sigma / beta 8      max |lpdf - brms| 3.997e-15  at y = -1.5
sigma / beta 10     max |lpdf - brms| 5.773e-15  at y = 0
sigma / beta 12     max |lpdf - brms| 1.554e-14  at y = -2
sigma / beta 15     max |lpdf - brms| 3.508e-14  at y = -1.5
sigma / beta 20     max |lpdf - brms| 5.862e-14  at y = -3
sigma / beta 40     max |lpdf - brms| 1.488e-13  at y = 1.5
sigma / beta 100    max |lpdf - brms| 1.443e-12  at y = -3
sigma / beta 1000   max |lpdf - brms| 7.155e-11  at y = -1.5
RTMBdist::dexgauss directly, mu 0 sigma 1 lambda 20, y = 0: -0.9189385  closed form: -0.9214231 
```
<!-- /GEN:fix1rerev -->

The new and changed tests on the pre-fix lane build
(`wt-cdfs-rev-lib`), GENERATED:

<!-- GEN:fix1pre -->
```
RESULT test-cdfs.R pass=151 fail=53 err=2 skip=0 warn=3
RESULT test-cens-windows.R pass=4 fail=8 err=3 skip=0 warn=1
RESULT test-cens-trunc.R pass=51 fail=9 err=3 skip=0 warn=0
RESULT test-v07.R pass=24 fail=0 err=1 skip=0 warn=0
lib: frmtmb 0.72.0 C:/Users/adf44/source/r/wt-cdfs-rev-lib; frmtmb 0.72.0 C:/Users/adf44/source/r/wt-cdfs-rev-lib
30 blocks fail or error there:
a censored Gamma at shape 1e5 reaches the independent maximum
a constant trunc() bound gives every row the point mass
a count's cheap forms keep both tails exact
a custom family's cdf_log is read for every censored row
a fractional discrete bound is refused by name
a hurdle's F(0) is hu to the last digit at small hu
a hurdle interval (0, y2
a left-censored row below 1e-308 keeps its log probability
a one-point discrete interval is an event of probability zero
a short von Mises arc next to +-pi keeps its digits
a von Mises interval on the far side keeps its digits
all four discrete censoring codes match a hand-rolled likelihood
an infinite truncation bound is absent and fits
an interval in the upper tail keeps its digits
compois's survivor is the upper sum where the complement is lost
discrete cens() composes with trunc() as brms composes them
discrete censoring gives no notice: it is brms's reading
each family's censored fit is its hand-rolled likelihood
interval rows' AD gradient is the finite difference's
poisson right censoring at counts near 1e5 is exact
the discrete convention is brms's, not the inclusive one
the ex-gaussian CDF's derivatives hold at sigma / beta of 1e3 and 1e4
the ex-gaussian density is the closed form past sigma / beta = 20
the incomplete gamma pair matches stats::pgamma() in both tails
the other continuous fits are their hand-rolled likelihoods
the weibull trunc() fuzz spec that stopped on a NaN gradient fits
the zero-inflated and hurdle fits are hand-rolled likelihoods
trunc(lb <= 0) on a count is no truncation, as in brms
trunc(lb = ) is divided by the log survivor at lb
zero-truncated poisson matches glmmTMB truncated_poisson
```
<!-- /GEN:fix1pre -->

The gated fuzz tier on this build is red with the two section-12 specs
only. The weibull spec now fits, and its `trunc_support` invariant then
meets the entry `trunc-postfit` already on `FUZZ_KNOWN_PENDING`
(`simulate()` found that its rejection sampler did not fill 1 of 120
rows), which triages it as known. GENERATED:

<!-- GEN:fix1fuzz -->
```
RESULT test-fuzz.R pass=1 fail=1 err=0 skip=0 warn=0
-- dev/cdfs-fuzz-base.R on this build, the three specs of section 8
specs: 3 
real_new | [REAL-NEW] summary_prints
    spec: family=cumulative aterm=weights re=nested special=none dpar=none mode=quadrature op=vcov (seed 20322452, cover)
    summary() errored: no 'dimnames' attribute for array
real_new | [REAL-NEW] fit_error
    spec: family=Gamma aterm=weights re=ar1 special=mo_int dpar=dpar_x mode=profile op=vcov (seed 20502220, pad)
    frm() failed: The optimizer failed on this model (Gamma, ML, profile, nlminb): NA/NaN gradient evaluation. The likelihood was undefined or unbounded somewhere the optimizer stepped. Refit with verbose = TRUE to see which stage broke, try another optimizer (frmtmb_control(optimizer = "optim")), or start the fit nearer the optimum with the `start` argument of frm()
```
<!-- /GEN:fix1fuzz -->

## 12. To file

- The cumulative + weights + nested random effects + quadrature fuzz
  spec (seed 20322452): `summary()` errors with "no 'dimnames'
  attribute for array". No CDF is involved; it fails the same way on
  rellib-r10. Lane fuzzfix.
- The Gamma + weights + ar1 + `mo()` + distributional profile fuzz spec
  (seed 20502220): NaN gradient. No `cens()` or `trunc()`; fails the
  same way on rellib-r10. Lane fuzzfix.
- A right-censored count at its binomial's trials, or an interval that
  starts there, is an event of probability zero under brms's reading
  (`P(Y > n) = 0`), and scores `-Inf`, or `NaN` for an interval, as in
  brms. The review's `dev/cdfs-rev-grad.R` data, drawn for the
  inclusive reading, have such rows. A refusal by name would be
  friendlier. Done in fix round 2 (N5, section 13).
- The weibull fuzz spec 20446531, which fits now, meets
  `FUZZ_KNOWN_PENDING`'s `trunc-postfit` entry through its
  `trunc_support` invariant: `simulate()`'s rejection sampler did not
  fill 1 of 120 rows in 100 passes at an acceptance rate of 0.48. That
  entry says the fix is on an unmerged sibling branch; whoever owns it
  should confirm that this case is the one it covers. Done in fix round
  2: the entry was stale and is deleted, and simulate() is fixed (N3).

## 13. Fix round 2

The review's "Re-check after fix round 1" (logs in
`dev/cdfs-rev-out/r1/`) and the coordinator's list. The pre-fix build
for this round is `wt-cdfs-rev-lib`, which the reviewer reinstalled
after fix round 1's last `R/` edit.

- **N1 (summed CDFs).** `poisson_cdf_log()`'s summed form (counts 1 to
  30) and `compois_cdf_log()` now take `log F` as `log1m_exp(lS)` from
  the exact upper sum wherever they take the survivor from it (S below
  about 1e-3), the mirror of the survivor. Every summed-mass CDF was
  checked: `discrete_sum_cdf_log()` serves beta_binomial, which sums
  both tails directly (right in the review), and compois; the
  zero-inflated and hurdle poisson read the poisson pair; no other
  family sums. Tests (`test-cens-windows.R`, failing on the pre-fix
  build): poisson, ZIP and hurdle_poisson `(30, 31]` at `mu = 3`, the
  poisson `(20, 25]` and `(25, 30]`, compois `(20, 25]` and `(31, 33]`
  against mass sums; `trunc(lb = 31, ub = 60)` at intercepts 3.3 to 0.5
  and at its optimum against the mass-sum likelihood; compois interval
  rows' AD gradient against Richardson differences at `nu = e^3` and
  `e^4`.
- **N2 (cens inside a two-sided trunc).** frmtmb's composition kept.
  `stancode()` run here confirms brms also writes a LEFT-censored row
  under `trunc(lb = )` alone as `lcdf(Y) - lccdf(lb - 1)`, so the two
  differ for right censoring with an upper bound and left censoring
  with a lower one. NEWS ("reads its bounds", a sub-bullet with 6.16 and
  14.85), `?frmtmb_family`, `row_lpdf()`'s documentation and the compat
  rows (`cens() x trunc()`, poisson, the discrete loop, discrete_weibull)
  say so; filed as brms-34 in `dev/upstream-bugs.md`. Test: right- and
  left-censored rows under `trunc(lb = 2, ub = 9)` (poisson) and
  `trunc(lb = -1, ub = 3)` (gaussian) against hand-rolled likelihoods,
  and the poisson against brms's whole-line form (differs by more than
  one). `test-cens-trunc.R`'s block is renamed "a right-censored count
  under trunc(lb = ) alone is brms's row".
- **N3 (stale fuzz entry, simulate()).** `trunc-postfit` is deleted
  from `FUZZ_KNOWN_PENDING`, with a comment. `sim_response()` falls
  back to `sim_trunc_inverse()` for rows rejection has not filled after
  20 passes, for this package's families with a CDF and no extra
  parameters (cox reads its CDF from a design built at the recorded
  bounds; another package's family keeps rejection alone, as the
  frmtmb.eam Wiener `lcdf` is a boundary's). The draw is a uniform point
  of the window's mass in the tail the lower edge leaves more of, then
  bisection (counts: the smallest value whose CDF reaches the point).
  Measured (GEN:fix2simcost): rejection is 20 to 100 times cheaper
  where the window keeps most of the mass, and inversion is no worse
  where it keeps little, so inversion is the fallback, not the default.
  Tests: the weibull fuzz spec simulates every row inside `[0.5, Inf)`,
  and 4000 draws of its least-mass row are within the 0.1 percent KS
  critical value of the truncated CDF; a poisson at `mu = 2` in
  `[12, 15]` (1e-6 of its mass) draws inside the window with its mean.
  The deleted entry had also been absorbing 14 findings on the trunc
  specs of Beta and zero_inflated_poisson, which fix round 1 added to
  the grid's trunc column: `predict_eq_fitted` and `simulate_mean` meet
  `fitted()`'s deliberate refusal, "Family 'beta' has no truncated
  mean", the brms 2.23.0 behavior (`dev/cdfs-fix2-fuzztrunc.R`, all 45
  trunc specs of the plan). That is a refusal on purpose, so it is
  `FUZZ_KNOWN_REFUSAL`'s new `trunc-no-mean` entry, matched on the
  trunc aterm, those two invariants and the message.
- **N4 (importance cost).** Profiled (`dev/cdfs-fix2-imp-prof.R`): the
  time was the tape, `MakeADFun` and its optimizer, made longer by
  `log1m_exp_ad()`'s two arms on every row's normalizer. With every row
  at a count of 0, `poisson_cdf_log()` now returns `-mu` and
  `log(-expm1(-mu))` (the survivor to its absolute accuracy, which is
  what a normalizer adds) with no sub-assignment, and `log_window()`
  reads a window open on the same side on every row without one. The
  fit is now 4.1 to 4.4 seconds and 751 to 762 MB, against 4.7 s and
  720 MB on rellib-r10 (GEN:fix2imp).
- **N5 (impossible censored counts).** brms 2.23.0, run here
  (`dev/cdfs-fix2-brms-top.R`): a binomial or beta_binomial row
  right-censored at its trials compiles, its `log_prob` is `-Inf`, and
  it makes no draws; a binomial left-censored at 0 samples (it is
  `P(Y = 0)`). frmtmb refuses a right- or interval-censored row at or
  above its trials by name, and fits the left-censored zero. A count
  below 0 is already refused as a response.
- **N6.** NEWS says 8.04e-14. The censored Gamma at shape 3000 stopping
  at nlminb's false convergence at the maximum (review `g3000.txt`:
  logLik equal to the independent one, AD and FD gradients agreeing at
  -0.0116) is left as it is: nlminb's tolerance at a CV of 1.8 percent,
  not a defect, but its warning now reaches a user.

GENERATED:

<!-- GEN:fix2 -->
```
-- the round's new tests on the pre-fix build (wt-cdfs-rev-lib)
RESULT test-cens-windows.R pass=23 fail=9 err=2 skip=0 warn=0
RESULT test-cens-trunc.R pass=69 fail=0 err=0 skip=0 warn=0
RESULT test-cdfs.R pass=234 fail=0 err=0 skip=0 warn=0
5 blocks fail or error there:
a summed CDF keeps its upper tail in a window
compois interval rows have a smooth objective at large nu
inverse draws of a count stay in the window with its mean
simulate() fills every row of a valid truncated fit, from the window
trunc(lb = 31, ub = 60) on a poisson is exact away from its maximum
-- N1: the review's sumtail.R and sumfit.R on this build
poisson mu 3 interval (10, 12]                           got -8.19442912    ref -8.19442912    rel 3.0e-14
poisson mu 3 trunc(lb = 11, ub = 12), row at lb          got -0.223143551   ref -0.223143551   rel 2.5e-13
zero_inflated_poisson mu 3 interval (10, 12]             got -8.55110406    ref -8.55110406    rel 2.9e-14
hurdle_poisson mu 3 interval (10, 12]                    got -8.50003488    ref -8.50003488    rel 2.9e-14
compois mu 3 nu 1.5 interval (10, 12]                    got -10.8804185    ref -10.8804185    rel 3.3e-16
poisson mu 3 interval (20, 25]                           got -25.1637613    ref -25.1637613    rel 4.2e-16
poisson mu 3 trunc(lb = 21, ub = 25), row at lb          got -0.145519524   ref -0.145519524   rel 1.1e-14
zero_inflated_poisson mu 3 interval (20, 25]             got -25.5204363    ref -25.5204363    rel 4.2e-16
hurdle_poisson mu 3 interval (20, 25]                    got -25.4693671    ref -25.4693671    rel 4.2e-16
compois mu 3 nu 1.5 interval (20, 25]                    got -35.430419     ref -35.430419     rel 2.0e-16
poisson mu 3 interval (25, 30]                           got -35.5805589    ref -35.5805589    rel 0.0e+00
poisson mu 3 trunc(lb = 26, ub = 30), row at lb          got -0.11722337    ref -0.11722337    rel 0.0e+00
zero_inflated_poisson mu 3 interval (25, 30]             got -35.9372338    ref -35.9372338    rel 0.0e+00
hurdle_poisson mu 3 interval (25, 30]                    got -35.8861646    ref -35.8861646    rel 0.0e+00
compois mu 3 nu 1.5 interval (25, 30]                    got -50.598768     ref -50.598768     rel 1.4e-16
poisson mu 3 interval (30, 31]                           got -47.0352426    ref -47.0352426    rel 1.5e-16
poisson mu 3 trunc(lb = 31, ub = 31), row at lb          got -7.10542736e-15 ref 0              rel 7.1e-15
zero_inflated_poisson mu 3 interval (30, 31]             got -47.3919175    ref -47.3919175    rel 1.5e-16
hurdle_poisson mu 3 interval (30, 31]                    got -47.3408484    ref -47.3408484    rel 1.5e-16
compois mu 3 nu 1.5 interval (30, 31]                    got -67.216541     ref -67.216541     rel 2.1e-16
poisson mu 3 interval (31, 33]                           got -49.3153548    ref -49.3153548    rel 1.4e-16
poisson mu 3 trunc(lb = 32, ub = 33), row at lb          got -0.087011377   ref -0.087011377   rel 7.1e-15
zero_inflated_poisson mu 3 interval (31, 33]             got -49.6720298    ref -49.6720298    rel 1.4e-16
hurdle_poisson mu 3 interval (31, 33]                    got -49.6209606    ref -49.6209606    rel 1.4e-16
compois mu 3 nu 1.5 interval (31, 33]                    got -70.6518908    ref -70.6518908    rel 6.0e-16
beta_binomial(50) mu .2 phi 5 interval (30, 35]          got -3.99954073    ref -3.99954073    rel 7.1e-15
binomial(50) p .2 interval (30, 35]                      got -22.92836      ref -22.92836      rel 4.6e-16
beta_binomial(50) mu .2 phi 5 interval (40, 45]          got -6.19504116    ref -6.19504116    rel 3.9e-15
binomial(50) p .2 interval (40, 45]                      got -46.2988738    ref -46.2988738    rel 4.6e-16
beta_binomial(50) mu .2 phi 5 interval (45, 50]          got -8.41579624    ref -8.41579624    rel 1.2e-11
binomial(50) p .2 interval (45, 50]                      got -62.5581979    ref -62.5581979    rel 1.1e-16
independent: 3.2865436 0.40089486 logLik -4137.953281
frmtmb:      3.2865436 0.40089478 logLik -4137.953281 conv 0
frmtmb objective at the independent optimum -4137.953281 vs -4137.953281
== rows whose normalizer is wrong at the independent optimum
== the objective along the intercept (slope 0.4): frmtmb vs independent
intercept 3.3: frmtmb -fn -4140.6813  independent -4140.6813
intercept 2.5: frmtmb -fn -10277.0817  independent -10277.0817
intercept 2.0: frmtmb -fn -16284.7749  independent -16284.7520
intercept 1.5: frmtmb -fn -22582.3035  independent Inf
intercept 1.2: frmtmb -fn -26423.2395  independent Inf
intercept 1.0: frmtmb -fn -28999.5407  independent Inf
intercept 0.5: frmtmb -fn -35476.2464  independent Inf
intercept 2.0: frmtmb -fn -16284.7749  exact -16284.7749
intercept 1.5: frmtmb -fn -22582.3035  exact -22582.3035
intercept 1.0: frmtmb -fn -28999.5407  exact -28999.5407
intercept 0.5: frmtmb -fn -35476.2464  exact -35476.2464
-- N1: the review's compois.R on this build
codes none nu eta 3: fn 1992.6081  d/d eta_nu AD 1938.74 FD(1e-3) 1938.74 FD(1e-5) 1938.74
codes none nu eta 4: fn 5415.0237  d/d eta_nu AD 5450.47 FD(1e-3) 5450.47 FD(1e-5) 5450.47
codes all  nu eta 3: fn 1689.9039  d/d eta_nu AD 1615.24 FD(1e-3) 1615.24 FD(1e-5) 1615.24
codes all  nu eta 4: fn 4340.3497  d/d eta_nu AD 4114.16 FD(1e-3) 4114.16 FD(1e-5) 4114.16
codes -1   nu eta 3: fn 1919.9626  d/d eta_nu AD 1875.03 FD(1e-3) 1875.03 FD(1e-5) 1875.03
codes -1   nu eta 4: fn 5269.238  d/d eta_nu AD 5355.19 FD(1e-3) 5355.19 FD(1e-5) 5355.19
codes 1    nu eta 3: fn 1915.5373  d/d eta_nu AD 1851.96 FD(1e-3) 1851.96 FD(1e-5) 1851.96
codes 1    nu eta 4: fn 5094.7596  d/d eta_nu AD 5015.51 FD(1e-3) 5015.51 FD(1e-5) 5015.51
codes 2    nu eta 3: fn 1839.6201  d/d eta_nu AD 1765.72 FD(1e-3) 1765.72 FD(1e-5) 1765.72
codes 2    nu eta 4: fn 4806.3994  d/d eta_nu AD 4644.41 FD(1e-3) 4644.41 FD(1e-5) 4644.41
-- N3: dev/cdfs-fix2-simcost.R
gaussian, trunc(lb = -1), keeps 84 percent     inv_after 20:     1.5 ms per draw of 10000 rows
poisson mu 3, trunc(lb = 1), keeps 95 percent  inv_after 20:     0.5 ms per draw of 10000 rows
weibull, trunc(lb = 0.5), keeps 10 percent     inv_after 20:   163.5 ms per draw of 10000 rows
Gamma shape 3, trunc(lb = 6), keeps 6e-3       inv_after 20:  2886.0 ms per draw of 10000 rows
gaussian, trunc(lb = -1), keeps 84 percent     inv_after  0:    23.5 ms per draw of 10000 rows
poisson mu 3, trunc(lb = 1), keeps 95 percent  inv_after  0:    56.0 ms per draw of 10000 rows
weibull, trunc(lb = 0.5), keeps 10 percent     inv_after  0:   244.5 ms per draw of 10000 rows
Gamma shape 3, trunc(lb = 6), keeps 6e-3       inv_after  0:  2865.0 ms per draw of 10000 rows
-- N4: dev/cdfs-fix1-imp.R, this build and rellib-r10
importance 500: objective 448.721461 seconds 4.3
peak memory (MB, memory.size): 764.84375 
importance 500: objective 448.721461 seconds 5.4
peak memory (MB, memory.size): 748.3984375 
-- N5: dev/cdfs-fix2-brms-top.R (brms 2.23.0)
binomial, right-censored at its trials           brm ran; log_prob at 0: -Inf; draws 0
beta_binomial, right-censored at its trials      brm ran; log_prob at 0: -Inf; draws 0
binomial, left-censored at 0                     brm ran; log_prob at 0: -64.58502; draws 20
```
<!-- /GEN:fix2 -->

The gated fuzz tier after the refusal entry, and the plan's trunc specs
before it, GENERATED:

<!-- GEN:fix2fuzz -->
```
RESULT test-fuzz.R pass=1 fail=1 err=0 skip=0 warn=0
[REAL-NEW] summary_prints
    spec: family=cumulative aterm=weights re=nested special=none dpar=none mode=quadrature op=vcov (seed 20322452, cover)
[REAL-NEW] fit_error
    spec: family=Gamma aterm=weights re=ar1 special=mo_int dpar=dpar_x mode=profile op=vcov (seed 20502220, pad)
-- dev/cdfs-fix2-fuzztrunc.R, before the entry: classes, and the
   REAL-NEW records' messages
cls
known_pending      real_new   unconverged 
            5            14             2 
    fitted/predict errored: Family 'beta' has no truncated mean; fitted(), residuals() and frm_linpred(type = "response") would report the untruncated mean
    simulate/predict errored: Family 'beta' has no truncated mean; fitted(), residuals() and frm_linpred(type = "response") would report the untruncated mean
    fitted/predict errored: Family 'zero_inflated_poisson' has no truncated mean; fitted(), residuals() and frm_linpred(type = "response") would report the untruncated mean
    simulate/predict errored: Family 'zero_inflated_poisson' has no truncated mean; fitted(), residuals() and frm_linpred(type = "response") would report the untruncated mean
```
<!-- /GEN:fix2fuzz -->
