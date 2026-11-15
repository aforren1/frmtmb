# Short check of lane predraw2's fix round (round of 2026-11-13)

Base ec8178d3, lane library `wt-predraw2-lib` (installed 21:11 UTC,
after the last edit of `R/predict-brms.R`; no file under R/, man/ is
newer), rellib-r21 behind it. Only the fixes were checked. Scripts
`dev/reviews/predraw2-short-*.R`, logs with the same names (`.txt`).
No lane file was edited, no git writes. Reference BLAS throughout.

## Verdict: FIX

Two small items: file the dilution limit of the ratio rule with the
numbers below (and say it in `?predict.frmtmb_fit`), and make the
no-loo allowance in `test-predict-draw-edge.R` conditional. FIX 1's
rule does what it claims on its own designs, and FIX 2 is right.

## 1. FIX 1, `draw_check_verdict()`

(a) Dilution. Yes: a block of well-identified parameters hides a
real non-normal direction. The ratio is median(e) / (median(m) / 2).
A separable block of K coordinates adds about K to the median of m and
nothing to e. Script `predraw2-short-a.R`, part S. Each of the 18
mixture fits that fire on their own (gaussian seeds 20 to 40, poisson
seeds 5 to 39, the `mc42` fits) uses its own e and m at the check's test
draws, with the same 20/40 sequential rule:

- Exactly quadratic block (e 0, m chi-square K): the check stops firing
  at K 8 on 2 fits, K 16 on 4, K 32 on 5, K 64 on 3, K 128 on 2. Only
  gaussian seeds 23 and 25 (excess 23.5 and 47.3) still fire at K 256.
  Example: gaussian seed 20 alone fires at excess 2.38, ratio 0.708.
  With K 8 it reads ratio 0.122 and does not fire. The weak direction
  and its excess are unchanged.
- The brief's design: the review's 200-level poisson `y ~ 0 + f` block
  (seed 1, alone excess 2.23, ratio 0.0237) added to each mixture.
  16 of 18 no longer fire (ratio 0.037 to 0.119). Only gaussian seeds
  23 and 25 still fire.

This is inherent to a summed check. The K200 block alone reads excess
2.23. Gaussian seed 29's one weak direction reads 2.26. No summed
statistic can separate the two. Before the fix round, the rule was
excess > 1 only, and it fired on both. The lane did not try a
per-direction reading (findings section 11).

Real fits, part R: the same data with K = 5, 10 or 20 noise covariates
in both components (54 fits). 11 fire. 36 have a summed excess of 1 or
less, so neither the old rule nor the new rule sees them. In 7 fits
the ratio decides: all are poisson, with excess 1.07 to 1.62 and
ratio 0.049 to 0.111 (seed 13 K 20, 20 K 10, 25 K 20, 26 K 20, 36 K 10,
36 K 20, 39 K 20). In those 7 fits the mean of predict() is 18.3 to
32.6 Monte Carlo SE from fitted(), with no message (the 4882
symptom). With the check forced on, every one warns
`frmtmb_draw_pareto` (k > 0.7). The error against the truth improves
on 3 and gets worse on 4 (for example 0.4665 to 0.4062, 0.5508 to
0.5957). So the hidden fits lose the warning, not better draws.

Not filed: `dev/test-backlog.md` "From lane predraw2" has four
items. None covers this, and `?predict.frmtmb_fit` says only that
"a model with many well-identified parameters does not" fire. FIX:
file it (suggested MEDIUM: a silent predict() on a weakly identified
mixture with many covariates), with the numbers above. Add one
sentence to the docs: a weak direction among many well-identified
parameters can go unseen.

(b) The review's K reproducer, `predraw2-review-d.R pois_K200_m5_mu2
pois_K100_m5_mu1.5 pois_K100_m5_mu1`, REVD_SEEDS 1:10 (log
`predraw2-short-d.txt`): 17 fits with a verdict, 0 fired, 1000
distinct draws each. The coefficient and predictive 90% coverage are
equal to the mocked 0.82.0 draw on all 17. The first four columns are
identical to the lane's `dev/predraw2-review-d-after.txt`. Holds.

(c) Mixture designs, seeds 1 to 40 (`predraw2-short-c.R`): g fires on
6, p on 12, c on 26. The same seeds as `dev/predraw2-mc42-*-lane.tsv`
on 40 of 40 for each design. Excess equal to the tsv within 3.4e-13.
Holds.

## 2. FIX 2, `lw` with log g'(x)

Code: `lw = -drop + m/2 + la` with `la` the log slope itself. The
anchor `log g'_side(0)` is only in `e`. The review's 1-D test
(`predraw2-review-a.R jac`, log `predraw2-short-a1jac.txt`) reads
0.2354 (exact 0.2308) and 0.0975 (exact 0.0909), both with k 0.00. A
single seed does not tell noise from bias, so the same test was run
over 40 seeds (`predraw2-short-e.R`):

    a 0.30 exact 0.2308 | IS mean 0.23101 (z 0.62) | resampled 0.23067 (z -0.19)
    a 0.10 exact 0.0909 | IS mean 0.09103 (z 0.63) | resampled 0.09152 (z 1.61)

No bias. The single-seed gap is Monte Carlo error. Holds.

k = Inf to 0 (`draw_resample()`), `predraw2-short-b.R`, S = 1000
(tail size 95), loo 2.10.1:

- All weights equal: loo gives Inf, the lane gives 0, 1000 distinct
  draws kept, ESS 1000. Right: no reweighting needed.
- One finite weight, the rest -Inf: the lane gives Inf (warning), 1
  draw kept. One weight at 0, the rest at -1000: also Inf. The
  `all(is.finite(tl))` and range guards keep the degenerate case
  visible.
- 95 equal top weights, the rest at -50: the lane gives 0, ESS 95. With
  94 equal top weights loo itself gives k -5.45 (ESS 94). So the
  override gives the same verdict as loo's neighboring case, and it
  hides nothing that loo would flag.
- Two values (0 and log 0.3): the lane gives 0, ESS 775.
- 24 equal weights: Inf is kept (the 25-draw floor).

Holds.

## 3. LOW, 3277 test allowance

`tests/testthat/test-predict-draw-edge.R:348-349` allows "They are NOT
corrected; install loo" UNCONDITIONALLY, beside "the reweighted draws
are not reliable", with no `require`. Only the `!requireNamespace("loo")`
branch of `draw_reweight()` prints that text. If a regression sends a
fit with loo present down that branch, the draws are returned
uncorrected and this test still passes. FIX, by the brief's rule.
Allow it only where loo is absent, for example
`c("the reweighted draws are not reliable", if (!requireNamespace("loo",
quietly = TRUE)) "They are NOT corrected; install loo")`.

## 4. Docs and NEWS

- `?predict.frmtmb_fit` (`R/predict-brms.R:144-163`, Rd regenerated)
  and NEWS lines 5-35 state both conditions (median excess above 1 AND
  above 0.15 times the law's median drop), the weighting "by the
  likelihood over the density it was drawn from" (now true, item 2),
  and "a parameter held at its estimate in part of the draws keeps that
  part" (true: `draw_resample()` resamples within each stratum and
  keeps its count). NEWS "stayed under 0.12 times the normal law's
  drop": the fa-sum table's largest unfired ratio outside the mixtures
  is 0.1163. True.
- Missing: the limit of item 1(a), as above.
- Minor, internal: the roxygen of `draw_check_tol`
  (`R/predict-brms.R:2183`) still says "the median ... above which
  draw_check() fires" and quotes the 430-fit corpus max 0.29. The K
  designs reach 7.4, and the rule needs both conditions. Not blocking.

## 5. Reruns on the lane library, reference BLAS

- `test-predict-draw-check.R`: 8 tests, pass 31, fail 0, err 0, skip 0,
  warn 0 (`predraw2-short-t-check.txt`).
- `test-predict-draw-edge.R`: 10 tests, pass 45, fail 0, err 0, skip
  0, warn 0 (`predraw2-short-t-edge.txt`). The printed lines are
  boundary-fit messages, not warnings.
- `dev/predraw2-ulp2.R` (21 roundings, `predraw2-short-ulp2.txt`):
  every fit keeps one verdict on all 21. The output is line for line
  identical to the lane's `dev/predraw2-ulp2-ref.txt`. The two fits
  that fire read 8.9851 (ratio 3.9137) and 52.29 to 52.32 (11.84). The
  unfired verdict nearest a boundary is E10_lam1, ratio 0.0551 to
  0.0592, 0.237 log units from the half point 0.075.

## Process

At most 2 R processes of mine at any time, all small (1-D, GLM,
300-row mixtures). Free space on C: 145.9 GB at start, 146.5 GB at
end. No temp output outside the `predraw2-short-*` scripts and logs.
