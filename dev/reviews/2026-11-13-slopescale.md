# Review of lane slopescale (round of 2026-11-13, base ec8178d3)

Adversarial review of item 1018 (autoscale and random-slope columns).
Lane tree `frmtmb-wt-slopescale`, library `wt-slopescale-lib`, base
`rellib-r21`. Scripts and logs: `dev/reviews/slopescale-review-*`. Where
base and lane are compared bitwise, both run in ONE R process (frmtmb
has no compiled code, so the namespace is unloaded and loaded again from
the other library; rtmb-pitfalls item 21), with base run twice to show
the process is deterministic.

## Verdict: FIX

One regression against 0.82.0 in the lane's own scope (F1). The rest of
the change holds up: the reparameterization is exact, ordinary fits are
bit-identical to 0.82.0, the far-point hunk gives no false verdict, and
every merge is clean.

## FIX items

### F1. The default no longer engages on an uncentered slope column of small sd (regression)

`R/autoscale.R:113` measures a Z column by its root mean square, and
`R/autoscale.R:338` (`autoscale_decide()`) engages the default when that
value is below `autoscale_slope_sd`. 0.82.0 used the fixed twin's sd. For
an uncentered column, x = c0 + s z, the RMS is near c0, so the default
no longer engages where 0.82.0 did.

Reproducer: `dev/reviews/slopescale-review-uncentered.R`. Data: 20
groups of 15, gaussian, slope sd 0.5 on z, x = 5 + 0.03 z, seed 2,
`y ~ x + (1 + x | g)`.

- Expected (0.82.0 and lane `autoscale = TRUE`): the maximum.
- Actual (lane default): 9.84 log-likelihood units short, no warning.

Over 80 fits (c0 in {5, 7}, s in {0.03, 0.01}, seeds 1 to 10, both
forms, `slopescale-review-uncentered.log`):

```
fits: 80 ; short > 0.001 : base 58 (silent 34 ) ; lane 44 (silent 22 ) ; lane TRUE 0
lane short where base was not: 3 ; base short where lane is not: 17 ; max lane short 41.1 ; max TRUE short 5.94e-09
```

The 3 regressions (all `y ~ x + (1 + x | g)`, all engaged on 0.82.0):
c5 s0.03 seed 2 (9.84 short, silent), c5 s0.03 seed 4 (6.05, silent),
c7 s0.03 seed 4 (6.09, one warning).

The lane's plan, when it runs, fixes all 80 (`autoscale = TRUE`: max
short 5.9e-9). So the plan and its scale are right; only the engagement
rule is wrong. `slopescale-review-uncentered2.log` shows the mechanism
on 4 cases. The plain fit and a restart from its own stop end at the
same point. The data with the column divided by its RMS (about 5),
fitted alone, reaches the maximum of the centered fit. Example, seed 2,
s 0.03, `y ~ 1 + (1 + x | g)`: plain -491.6586, restart -491.6586,
pre-fit alone -450.6303, centered -450.6303.

Fix: engage on the column's spread about its mean (sd over the rows
where it is nonzero), and keep the RMS as the scale. At minimum, keep
0.82.0's engagement where it engaged. Add the 3 regression cases to
test-autoscale.R, and show them failing on the lane as it is now.

### F1b. Same mechanism, wider, not a regression: the default never engages on a poorly centered slope column (needs a user decision)

`slopescale-review-uncentered3.log`, seeds 1 to 5, short against the
centered, standardized fit, which is the same model:

```
    c0    s                form short silent shortT
 5e+00 0.10 y ~ 1 + (1 + x | g)     4      4      0
 5e+00 1.00 y ~ 1 + (1 + x | g)     0      0      0
 2e+03 3.00 y ~ 1 + (1 + x | g)     5      5      0
 5e+00 0.10 y ~ x + (1 + x | g)     4      3      0
 5e+00 1.00 y ~ x + (1 + x | g)     0      0      0
 2e+03 3.00 y ~ x + (1 + x | g)     5      2      0
```

A calendar year as a random slope (x = 2000 + 3 z, sd 3) is short on
10 of 10 fits, 7 of them silent, by 3.3 to 41.1 units. The rows for
the year are in the log. `autoscale = TRUE` is within 1.1e-8 on all 30
cases. What fails is the ratio of RMS to sd, not a small sd: x = 5 + z
is fine, and x = 5 + 0.1 z (sd above 0.05) is not. 0.82.0 does the
same, so this is not a regression. It is a common design, and the
silent shortfalls are as large as those of item 1018. Recommendation:
fix it the same way as F1 (engage when the RMS-to-sd ratio of a slope
column is large, calibrated as `autoscale_slope_sd` was), or file it as
HIGH. Changing the rule for when the default engages is the user's
call. It extends the existing default uniformly; it is not a special
case.

### F2 (LOW). The default-mode note fires on fits that are not short

`R/autoscale.R:267` with `R/fit.R:1368`: under the default, with nothing
planned, a slope that is not rescaled is named when its RMS is below
0.05. In `slopescale-review-rms.log` (b) the note fired on 4 of 21
small-unit designs. All 4 are `homdiag(1 + t | g)` with t a time in
years (sleepstudy Days/250 and Days/365) or a proportion (0 to 0.08).
`autoscale = FALSE` reaches the maximum on every one of them (short
3.7e-10, 4.5e-10, 0, 6.3e-12). A homdiag block shares its sd with the
intercept, so the stall mechanism does not apply. It is a message, not
a warning, and its text is conditional ("If the fit stops short").
Still, it fires 4 of 4 times where its condition is absent, and the
lane's test (`names every structure`, the default arm) pins it on such
a case without showing a short fit. Either drop the default-mode arm
(keep the note for `TRUE` and for engaged plans, as 0.82.0 did), or show
a fit that the note's case leaves short.

## Checked and holding

### 1. When the default engages; ordinary fits

`slopescale-review-ordinary.log`: 21 ordinary fits at spread near 1.
They are sleepstudy (3 forms); gaussian, poisson and bernoulli
`(1 + x | g)` and `x + (1 + x | g)`; cs, ar1 and homdiag; uniform(0, 1)
x; x = 5 + z; and two singular fits.

```
base deterministic: 21 of 21 ; lane identical to base (all fields): 21 of 21
```

The fields are opt$par, logLik, vcov, VarCorr with Est.Error, and every
condition with its text. This is an existing default extended
uniformly (engage below 0.05 on a slope column, now with or without a
fixed twin and in more structures), not a special case, except for the
change of measure in F1.

### 2. Back-transform: exact at a fixed point

`slopescale-review-exact.log` (seeds 1 and 3, 13 structures, plus
`-xtwin.log`). Both parameterizations are evaluated at the lane's
optimum and at its `autoscale_map(.., "scaled")` image, with
`optCtrl = list(iter.max = 0)`. Largest relative differences:

| quantity | non-nested | nested |
|---|---|---|
| objective | 4.8e-16 | 2.0e-13 |
| VarCorr sd | 8.0e-16 | |
| cor (absolute) | 1.1e-16 | |
| ranef, coef | 2.2e-16 | |
| fitted | 4.3e-13 | |
| vcov | 8.1e-13 | 6.7e-9 |
| predict draws, seed for seed (of sd(y)) | 6.0e-15 | 2.2e-8 |
| VarCorr Est.Error | 1.4e-9 | 1.5e-6 |

The remaining Est.Error values come from finite differences, at a flat
ridge in the nested case.

The lane's 1e-5 to 1.2e-4 is the gap between two optima, not a
back-transform error. The spread-1 optimum mapped onto the 1e-4 data
differs from the lane's optimum by at most 7.4e-5 in a parameter. In
objective it differs by at most 6.2e-9, which matches the quadratic
model 0.5 d'Hd (6.5e-9, idlink seed 1). No accessor reads the plan. The
reported fit is the unscaled model warm-started from the template, so
a map error could only change the start, not an accessor.

### 3. autoscale_block_exact() is a spot check, not a proof (LOW, doc only)

The check tests one invariance (density of b*c at theta + log c equals
density less the Jacobian) at one point, theta = start + 0.25 cos(k).
`slopescale-review-blockexact.log` registers a valid one-to-one sd
parameterization, sd = exp(theta) * (1 + exp(5 (theta - 7))). It is
scale equivariant to 1e-15 near theta 0.25 and not near theta 8.5,
where a slope on x * 1e-4 lands. It passes the check on 6 of 6 seeds,
so the slope is rescaled. The outcome is still right on 6 of 6 (short
at most 1.3e-9, no warning), because the reported fit re-optimizes the
user's model from the template. A false pass costs a worse start, not a
wrong answer, as far as measured. The roxygen of
`autoscale_block_exact()` ("Whether scaling ... is an exact
reparameterization") overstates it. Say "at one point". The
`mybump` test case (a constant added to the variance) is correctly
refused.

### 4. The far-point hunk (`R/se-check.R:2890-2898`)

`slopescale-review-singular.log`: 260 fits with true slope sd 0
(gaussian, poisson and bernoulli, both forms, seeds 1 to 20, spread 1
under the default and 1e-4 under `autoscale = FALSE`), plus 20
bernoulli `(1 | g)` fits with true sd 0 and groups of 3.

```
fits: 260 ; lane 'stopped short': 0 ; lane 'is at zero': 54 ; lane probe gain > 1e-6: 4 ; max lane gain 19.5
base 'stopped short': 0 ; base 'is at zero': 54 ; conditions identical: 260 of 260 ; logLik identical: 260
```

The 4 probe gains are poisson `(0 + x | g)` at 1e-4 under
`autoscale = FALSE`. There the slope term absorbs the group effect that
the model leaves out, so its sd is not at zero. Base and lane give the
same single warning. Under OpenBLAS 0.3.32 with LAPACK
(`slopescale-review-singular-ob0332.log`, lane only), the counts are the
same: 0 short, 54 at zero, 4 gains, max 19.5.

The hunk also feeds `sd_up_escape()` (`R/fit.R:4432`) and
`headroom_drop_rows()` (`R/fit.R:6177`). The findings do not say so.
With all 260 fits and 21 ordinary fits identical to 0.82.0, no change
was seen through them. The findings should name these two other
callers.

selimit2: a three-way merge of R/se-check.R (CRLF to LF) is clean, and
selimit2 keeps "the fit stopped short of the maximum: " (its line
3914). The merged tree (lane plus selimit2's se-check.R, fit.R,
joint-sparse.R, predict.R and se-reparam.R) was loaded from source with
`pkgload::load_all()`. Nothing was installed.
`RESULT merged test-autoscale.R pass= 176 fail= 0 err= 0` and
`RESULT merged test-se-limit.R pass= 362 fail= 0 err= 0`.

### 5. Note threshold on ordinary small-unit designs

See F2. The default engaged on 6 of 21 small-unit designs (RMS 0.015
to 0.047). Each one ended within 2.2e-11 of the best fit, and none
warned.

### 6. Nested seed 3, logLik only

Justified. At a fixed point the nested objective agrees to 2.0e-13 and
VarCorr sd to 7.6e-16. Between the two optima, nested seed 1 moves a
parameter by 1.02 for an objective change of 8.3e-9, and seed 3 by
0.263 for 1.7e-8. That is a flat ridge, and the spread-1 fit is the
better one on seed 3. No back-transform defect hides there.

### 7. Extensions

`frm_sample()` builds its objective with `dry_run = "objective"`
(`extensions/frmtmb.sample/R/sample.R:1369`), and `fit_assembled()`
skips autoscale under `objective_only` (`R/fit.R:1363`). The lane does
not reach the sampler. In `slopescale-review-sample.log`, laplace draws
on the autoscaled fit (x * 1e-4) are on the original scale: sd_g__x
median 3680 raw, 0.368 times 1e-4. The spread-1 fit gives 0.445. The
gap comes from the brms default sd prior (student_t(3, 0, 2.5)), which
is not invariant to the units of x. 0.82.0 gives the same, and brms
does the same. frmtmb.spatial's `car`/`spde` and frmtmb.latent register
no structure that carries a slope column (`car` has `sd_idx = 1` over
indicator columns), so nothing is planned there. The lane's
all-suites log has 529 RESULT lines, all fail 0 and err 0. Every `lib:`
line names the lane core.

### 8. Tests

The rewritten "names every structure" test is not a loosening. It
still asserts a registered non-scaling structure (`mybump`) and
`homdiag(1 + x | g)` are named, and its absent cases are explicit.
OpenBLAS 0.3.32 with LAPACK (`dev/ciharden-openblas.sh 0.3.32 lapack`,
built and then deleted):

- `RESULT test-autoscale.R pass= 176 fail= 0 err= 0 skip= 0 warn= 0`
- `RESULT test-se-limit.R pass= 240 fail= 0 err= 0 skip= 0 warn= 0`

The `expect_same_model()` tolerance (1e-4 relative) is loose against
the fixed-point agreement of 1e-15. It is acceptable because the
reported fit never reads the map. A fixed-point check with
`iter.max = 0`, as in `slopescale-review-exact.R`, would pin the map
itself.

### 9. NEWS and merges

NEWS has a new top `# frmtmb (development version)` section, with
BREAKING marked, no lane names, no `dev/` paths and no version numbers.
Three-way merges (CRLF to LF) against base, pairwise and cumulative,
are all clean:

- R/fit.R with priors2, quadseam2, selimit2, gpfield, wcheck2, refits
  and famend.
- R/se-check.R with quadseam2 and selimit2.

nlstop2 and arspat2 leave R/fit.R at base in their current trees.

## Needs the user

- F1b: whether the default should engage on a poorly centered slope
  column (calendar year, x = 5 + 0.1 z). `autoscale = TRUE` already
  fixes 30 of 30. Recommend yes, calibrated. Otherwise file it as HIGH.
