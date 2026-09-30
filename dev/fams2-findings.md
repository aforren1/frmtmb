# Lane fams2: xbeta(), zero_inflated_beta_binomial(), hurdle_cumulative()

Lane fams2, 2026-09-29, base commit 1f40800d (frmtmb 0.65.0). The three
brms 2.23.0 families that blocked 32 bin-1 assertions of the ported
brms suite now exist, fit, and reach the post-fit methods their
siblings reach. brms's `cse()` is accepted as `cs()`. The template was
lane fams (`dev/fams-findings.md`).

The brms facts below were READ from brms 2.23.0, not recalled:
`dev/fams2-brms-src.R` (output `.txt`) prints the constructors,
`.family_<name>()`, `posterior_epred_*`, `posterior_predict_*`,
`log_lik_*` and the Stan chunk files; `dev/fams2-brms-stancode.R`
prints the generated Stan programs, `standata()` and `default_prior()`;
`dev/fams2-brms-src2.R` prints brms's `extra_cat` handling of the
hurdle response.

## What changed

### The families (`R/families.R`)

- `xbeta(link = "logit", link_phi = "log", link_kappa = "log")`, dpars
  `mu`, `phi`, `kappa`, response in `[0, 1]`. The density is brms's
  `xbeta_lpdf()`: `P(Y = 0) = I_q(a, b)`, `P(Y = 1) = I_q(b, a)` with
  `q = kappa / (1 + 2 kappa)` (the second because
  `1 - (1 + kappa) / (1 + 2 kappa)` is that same `q`), and
  `dbeta((y + kappa) / d, a, b) / d` inside, `d = 1 + 2 kappa`.
  `post$mean_fn` is `posterior_epred_xbeta()` term for term;
  `post$var_fn` the exact second moment (`xbeta_moments()`); `sim` is
  the censored stretched beta. `post$fit_check` (`xbeta_fit_check()`)
  warns when the response has no exact 0 or 1 AND `kappa` ran to 0
  (punch round 1, M1). The interior density is written out by hand
  (punch round 1, B1).
- `log_ibeta_half()` / `log_ibeta_cf()`: the incomplete beta the two
  ends need. It replaces `RTMB::pbeta()`, whose third derivatives are
  `NaN` at ordinary points (next section). ROUND 1 used Lentz's
  continued fraction with 50 steps and a hard switch to the complement
  at `m`; the review falsified its accuracy near `m` (B2). The
  punch-round form is below, under "Punch round 1". What survives from
  round 1: each branch is evaluated at `x` clamped to its own side, so
  the discarded one is still a converging fraction and cannot put
  `NaN` through its zero weight, and the clamp is `x - max(x - m, 0)`,
  which returns `x` bit for bit on its own side; the pitfalls-file form
  `0.5 (x + m - |m - x|)` rounds `x + m` and cost 11 digits of `x` at
  `x = 6.5e-6` against `m = 0.69` (seen on the first sweep: error 7e-12
  where the fixed form gives 4.6e-15).
- `zero_inflated_beta_binomial(link = "logit", link_phi = "log",
  link_zi = "logit")`, dpars `mu`, `phi`, `zi`, with `trials()`:
  `P(0) = zi + (1 - zi) BB(0)`, `P(y) = (1 - zi) BB(y)`. The shapes come
  off the log-odds as in `beta_binomial()`; `mean_fn` is
  `posterior_epred_zero_inflated_beta_binomial()`, `var_fn` exact.
- `beta_binomial()` gains `post$var_fn` (`beta_binomial_var()`), so
  pearson residuals work. PIN: on rellib-r3 the test fails with "Family
  'beta_binomial' has no variance function; pearson residuals are
  unavailable" (`dev/fams2-pins.txt`).
- `hurdle_cumulative(link = "logit", link_hu = "logit", link_disc =
  "log", threshold = "flexible")`, dpars `mu`, `hu`, `disc`, `type =
  "ordinal"`, response `0..K`. `P(0) = hu`; `P(k) = (1 - hu) [F(disc
  (tau_k - eta)) - F(disc (tau_{k-1} - eta))]`, with brms's links
  (logit, probit, probit_approx, cloglog, cauchit, softit) through
  `ord_cumulative_logpmf()`, the data path of `cumulative()`'s density
  with `disc` in it (log-space differences for the logit-type links).
  Thresholds are ordered, `(tau_1, log increments)`, as for
  `cumulative()`. `disc` is held at 1 through the existing
  `fixed_dpars` route unless the formula models it, as in brms. A new
  family field `extra_cat = TRUE` (brms's own name for the special)
  says the codes start at 0; `ord_code0()` reads it. An ordered factor
  codes its FIRST level as 0, as `brms:::data_response` does.
  `thres(x = )` works, counted over the rows above the hurdle;
  `thres(gr = )` is refused by name (the grouped densities in
  `R/thres.R` read `mu` alone and would drop the hurdle silently), as
  are `cs()` (as for `cumulative()`) and `threshold != "flexible"`.
- `cse()`: `R/parse.R` renames a barless top-level `cse(...)` term to
  `cs(...)` after the whole-term refusals, so `treat * cse(carry)` is
  still refused in its own words. `rating ~ treat + period + cse(carry)`
  gives the identical logLik and `fixef()` as `cs(carry)`. PIN: on
  rellib-r3, "could not find function "cse"" (`dev/fams2-pins.txt`; the
  weak form, since the name is what is missing).

### The ordinal machinery, generalized for a family with a category 0
and dpars beside `mu`

Found by a read-only survey of every ordinal consumer (the agent's
table is summarized here). Each change is a no-op for the four
existing ordinal families, whose `ord_code0()` is 1 and whose dpars are
`mu` alone.

- `R/predict.R`: `ordinal_ncat()` counts the category 0;
  `ordinal_codes()` new; `ord_probs_from_eta()` evaluates the density
  at the codes and takes the other dpars' row values (`more`);
  `ord_probs()` evaluates them at the rows through a new
  `ord_dpars_more()`; `ord_prob_se()` adds each estimated extra dpar's
  design and `dP / d eta` to the delta method (so a band on a hurdle
  category includes `hu`'s uncertainty) and takes `allow_new_levels`;
  `ord_cat_moments()` scores by the codes; `sim_restore_type()` maps
  code 0 to the first level; `hurdle_cumulative` joins
  `osa_point_mass_families`.
- `R/frame.R`: ordered-factor coding, the `non-negative`/`positive`
  wording of brms's refusal, the `thres(x = )` against factor levels
  check, the `cs()` gate (now `ord_ordered_families`).
- `R/thres.R`: `thres_pin_recode()` shifts by the code.
- `R/conditional-effects.R`: the `categorical = FALSE` expected category
  is weighted by the codes; `pp_check()` and newdata `y` shift factor
  positions to codes. `R/diagnostics.R` the same for DHARMa's draws.
  `R/predict-brms.R` `predict_category_props()` tabulates codes and
  labels `P(Y = 0)`.
- `R/priors.R`: the threshold prior's `ordthres` scale is chosen by
  `ord_ordered_families`, not by the name "cumulative".
  `R/methods-fit.R` `coef_shift_thresholds()`: hurdle_cumulative is a
  `thres_minus_eta` family (brms's specials say so).
- `family_link_str()`: an ordinal family names its other dpars' links,
  `cdf = logit; hu = logit; disc = log` (unchanged `cdf = <link>` for
  the four).

DIVERGENCE FROM BRMS, on purpose: brms's `conditional_effects()` scores
the expected category by COLUMN POSITION
(`brms:::ordinal_probs_continuous()` multiplies column `k` by `k`), so
for a hurdle model its curve is `E[Y] + 1`, off the response's own
scale. frmtmb scores by the codes. Documented in `?frmtmb-families`
and the migration vignette.

### Registration

- `R/links-brms.R` regenerated by `dev/famlink-gen-links.R` (three rows
  added to its map); brms bars none of the three as a mixture
  component. `R/links.R` tables. `R/compat.R`: vocabulary, `discrete`,
  `point_mass`, `trials_families`, `ordinal`, `categorical_probs`, and
  two named rules for `hurdle_cumulative` (`thres()` conditional,
  `residuals_osa` refused). `R/frame.R` `trials_families`.
- frmtmb.sample: `ce_cat_mean()` takes the code of the first category
  (read off the family field, so it needs no new core export); the
  default-prior disclosure names `kappa` beside `shape`, `phi`, `nu`
  (brms gives `kappa` a gamma prior for `xbeta` and `von_mises`).

## RTMB::pbeta() and the Laplace approximation

`dev/fams2-pbeta-nan-rate.R` (seed 20260929, 4000 points, q = kappa /
(1 + 2 kappa), log kappa ~ U(-6, 2), log a and log b ~ U(-4, 6)):
`RTMB::pbeta()`'s value, gradient and Hessian are finite at 3839 of
3840 points with log P > -700, and its THIRD derivative at 3780 (60
non-finite, mostly where P is near 1). `dev/fams2-pbeta-third.R` shows
one such point (q = 0.2689, a = 0.2231, b = 2.718), where every third
partial is `NaN`, and a neighbour where all are finite: the failure is
sporadic in the inputs.

The Laplace approximation needs third derivatives, so the first xbeta
build (on `RTMB::pbeta()`) stopped at "NA/NaN gradient evaluation" on
a `(1 | g)` fit on the first try. `dev/fams2-validate.R` section 2 fits
the SAME density on `RTMB::pbeta()` as a custom family beside
`xbeta()` on ten seeds with a random intercept:

```
  seed 1001: xbeta() ok; the same density on RTMB::pbeta() NA/NaN gradient
  ... (all ten the same)
  seed 1010: xbeta() ok; the same density on RTMB::pbeta() NA/NaN gradient
```

ROUND-1 MEASUREMENT, SUPERSEDED by "Punch round 1" below. The round-1
fraction on the tape (`dev/fams2-ibeta-tape.R`, output `.txt`; 3000
points, x from e^-12 to 0.5, shapes e^-4 to e^16; 1604 with
-600 < log P). This sweep drew `x` and the shapes independently, so
it seldom landed near `m`, where the round-1 fraction was wrong (B2):

```
value rel err   ours   max 1.19e-12 q99 4.20e-14 q999 4.49e-13
value rel err   RTMB   max 9.99e-16 q99 2.22e-16 q999 8.19e-16
third   finite  ours 1604  RTMB 1588
gradient vs RTMB::pbeta AD: max 4.79e-10 q99 1.73e-11 q999 2.71e-10
value rel err by max(a, b):  (0,100] 4.6e-15  (100,1e4] 3.8e-14
                             (1e4,1e6] 1.2e-12  (1e6,Inf) 6.7e-14
gradient vs RTMB AD by max(a, b): 1.3e-13  1.4e-11  1.4e-10  4.8e-10
```

Round 1 claimed from this that "at the shapes an xbeta fit visits
(`phi` up to 1e4)" the value is good to 4e-14. That was WRONG twice:
near `m` the error reached 1.7e-6 at shape sum 1e4 and 1.0 at 1e6
(`dev/fams2-rev-thresh.txt` of the review), and an xbeta fit did not
reach `phi` 1e4 at all, because of `RTMB::dbeta()` (B1). Round 1 also
said "the value error does not improve from 50 to 100 steps"; that
holds only away from `m`, where the fraction converges in a few
steps. Near `m` it needs O(sqrt(max(a, b))) steps (the review: 1041
to 1088 at shape 1e6).

The unit test pins the NaN point with an expectation to FLIP, not
delete, if RTMB's `pbeta()` gains finite third derivatives there;
`log_ibeta_half()` could then go.

## Punch round 1

The review (`dev/reviews/2026-09-29-fams2.md`) found two blockers in
`xbeta()`, one major and four minor items. The reviewer's scripts were
rerun on the final build, output `dev/fams2-p1-rev-<name>.txt`.

### B1. xbeta did not fit precise data

`RTMB::dbeta()` with its first argument on the tape has a `NaN`
gradient once the shape sum passes about 1e3
(`dev/fams2-p1-rev-thresh.txt`: finite at 1046, `NaN` from 1047).
The interior density is now written out,
`(a - 1) log z + (b - 1) log1p(-z) - lbeta_ad(a, b) - log(d)`, where
`lbeta_ad()` forms `log B(a, b)` from `lgamma_shift_diff()` based at
the larger shape (a smooth blend of the two orderings, so the
derivative in `a` is right at `a == b`; a `min`/`max` form gave 0
there, because RTMB's `abs()` has derivative 1 at 0).
`RTMB::lbeta()` cancels two large `lgamma()` values for a small shape
beside a large one, and cost the incomplete beta's gradient seven
digits at `a = 0.01`, `b = 5e6`.

`dev/fams2-p1-dbeta.R` (output `.txt`), value and gradient in
`(z, a, b)` against a 200-bit Rmpfr reference, `phi` from 1e0 to 1e6,
five means, five points each: `RTMB::dbeta()`'s gradient is finite at
96 of 168 points (non-finite somewhere at every `phi` from 1e2 up);
the package form's at 168 of 168, value within 1.7e-10 relative
(floored at one)
and gradient within 1.3e-11 at `phi` 1e6, within 2.8e-12 and 5e-14 at
1e4.

The review's scripts now: `dev/fams2-p1-rev-phi5k.txt` fits at `phi`
3e3, 5e3 and 1e4 (round 1: `NA/NaN gradient evaluation` at all
three); `dev/fams2-p1-rev-kfix.txt` reaches the reference optimum to
1.2e-8 in logLik at every row, `phi` 2e3 to 2e5, free and held `kappa`.

### B2. log_ibeta_half() near its switch

The new `log_ibeta_half()` (its `@noRd` help says why each piece is
there):

- The continued fraction is evaluated BOTTOM-UP from a fixed 2N + 1
  terms (N = 50), not by Lentz's method. Lentz's first denominator
  `1 - (a + b) x / (a + 1)` is 0 at `x = (a + 1) / (a + b)`, inside
  the range the direct side reads once it runs past `m`; bottom-up
  has no such pole (seen: `NaN` at `x = 0.041`, `a = 40`, `b = 960`).
- The direct fraction and the complement are blended with a C2
  quintic smoothstep over `[m + sd/2, m + 3 sd/2]`, not switched at
  `m`. `dev/fams2-p1-sides.R` (output `.txt`) measured each side by
  `(x - m) / sd`: the complement loses digits below `m` for unbalanced
  shapes (2.7e-4 at `a = 0.3`, `b = 7`, two sd below), the direct side
  past `m + 1.5 sd` at large shapes (1.3e-8 at `a = 400`, two sd
  above). The complement takes `log(1 - x)` and `log(x)` from `x`
  itself, not from a rounded `1 - x`.
- Where the fraction needs O(sqrt(max(a, b))) steps, near the mean at
  large shapes, the value is `log(RTMB::pbeta())`: once
  `a b / (a + b)` passes 150, within 6 sd of the mean. There the
  probability is moderate, so `log()` loses nothing, and
  `RTMB::pbeta()`'s third derivatives were finite at 20000 of 20000
  random points (`dev/fams2-p1-pthird.R`, seed 11, shapes 150 to 1e7,
  within 6.5 sd). Punch 1 said "its NaN points are at small shapes,
  which this branch never reads"; that is FALSE at the mean itself,
  where every derivative is `NaN` (the re-check's n1, fixed in punch
  2 below). The blends are smoothsteps in `|x - m| / sd`
  from 5 to 6 and in `log(a b / (a + b))` from `log(150)` to
  `log(450)`; where the branch has no weight it is read at clamped
  shapes and `x`, in the same safe region.
- The exact-tie `x == m` gradient blow-up is gone: no switch weight
  divides by `|x - m|` any more. The clamps are one-sided at a tie
  (RTMB's `abs()'(0) = 1`), which is harmless, because each clamp sits
  where its branch has zero weight.

Measurements on the punch-1 build (punch 2 below supersedes the
accuracy and cost figures):

- `dev/fams2-p1-sweep.R pkg` (output `dev/fams2-p1-sweep-pkg.txt`,
  seed 20260930): 3167 points, shapes log-uniform on 1e-3 to 1e7,
  half within 5 sd of `m`, a tenth at `m` exactly. Value error
  against `stats::pbeta(log.p = TRUE)`, relative and floored at one:
  max 8.49e-12, q99 3.10e-13; at `x == m` (297 points) max 4.97e-13;
  below 5e-13 while both shapes are below 1e4. Non-finite value,
  gradient, Hessian or third derivative: 0. Gradient against
  `RTMB::pbeta()`'s AD at the ten worst stencil points: max 7.33e-10.
- The review's `dev/fams2-rev-ibeta.R` (`dev/fams2-p1-rev-ibeta.txt`):
  grid A max 2.35e-12 (round 1: 4.4e-9); grid B at the switch, worst
  4.02e-11 at shape sum 1e6, where the value is -44.09 (relative
  9e-13; the Rmpfr check at that point: ours -44.0916348894, Rmpfr
  -44.0916348876); section C (xbeta-shaped, `phi` to 1e6) max 2.0e-12;
  D, the one-sided derivatives across `m` agree to the last printed
  digit at `eps = 1e-12` at every shape, including `a = 3e4`,
  `b = 7e4` (round 1: 613 and 809 against 549); E, no spike at
  `x == m`; F, 0 non-finite third derivatives of 2061 (`RTMB::pbeta()`
  35). `dev/fams2-p1-rev-thresh.txt`: worst error within 1e-2 of `m`
  5.8e-13 at every shape sum from 1e3 to 1e5 (round 1: 1.7e-6 at 1e4).
- `dev/fams2-p1-joins.R` (output `.txt`): at every blend edge and at
  `m`, the value jump over a step of 1e-9 sd matches the predicted
  `g dx` to 6e-13 (except 1.0e-10 at `a = 2e5`, `b = 8e5`, 6 sd below,
  where the value is large), and the Hessian and third derivatives
  are finite at each edge.
- Cost, CORRECTED in punch 2: `dev/fams2-p1-cost.R` reported 21 us a
  row and 1.18 times "the round-1 fraction", but its "round 1" was the
  round-1 switch over the NEW bottom-up fraction, not round 1. The
  re-check's `dev/fams2-rev2-mean.R` (section 4; `MakeADFun` gradient,
  500 rows, seed 5, true round 1 rebuilt from the reviewed diff) gives
  10.3 us a row for the punch-1 build, 0.74 times true round 1
  (13.9 us) and 6.9 times `log(RTMB::pbeta())` (1.5 us). Only the rows
  at 0 or 1 pay it.

`dev/fams2-p1-pins.R` (output `.txt`) replays the round-1 density
(copied verbatim) as a custom family on the new tests' designs:
at `phi` 1e4 and 2e4 the round-1 form stops at `NA/NaN gradient
evaluation`, xbeta fits and matches the reference; at `phi` 2e5 with
`kappa` held, the round-1 form's logLik is 1.68 below the reference at
its own estimates (the review's figure), xbeta's 1.2e-8 above; on the
new test's incomplete-beta grid the round-1 form's worst error is
6.1e7 tolerance units, the new one's 0.36. These are the pins seen
failing: the round-1 build is not installed anywhere to rerun the
test file against.

### M1. The kappa warning

`xbeta_fit_check()` now warns only when the response has no exact 0
or 1 AND the fitted `kappa` is below 1e-6 at every row, and says what
the fit shows: "kappa ran to 0 (at most X), where the model is
Beta(); its standard error is not usable. Without 0s and 1s only the
interior shape places kappa." It no longer says the likelihood rises
as `kappa` falls or recommends `Beta()` outright. The review's A1
(`dev/fams2-p1-rev-a1.txt`): the fit reaches the reference optimum,
objective -1666.6332, `kappa` 65.8, and no `kappa` warning. It does
warn "false convergence" and has `NaN` standard errors: the optimum is
a flat ridge in (`phi`, `kappa`) (`phi` 4.0e5 here, 2.0e5 at the
reference's), which frmtmb's existing convergence warning already
reports. The A2-type case (zeros, no ones) stays quiet, and a
`kappa` held with `bf(kappa = 0.1)` stays quiet (tests).

### Minor

- m1 (the `x == m` gradient): fixed with B2.
- m2 (`hurdle_cumulative()` with `disc ~ 1 + z` fit silently). brms
  2.23.0 ACCEPTS it: `dev/fams2-p1-brmsdisc.R` (output `.txt`) shows
  the generated Stan code puts `normal(0, 1)` on `Intercept_disc` and
  says nothing. Refusing would break a brms formula that brms fits, so
  frmtmb WARNS, through a new `post$fit_check`
  (`hurdle_cum_fit_check()`), when `disc`'s design has an intercept
  and no prior covers it: "disc has an intercept, which the likelihood
  cannot tell apart from the scale of the thresholds ... Write
  disc ~ 0 + ..., or hold the intercept with a prior, as brms does with
  its default set_prior("normal(0, 1)", class = "Intercept",
  dpar = "disc")". `dev/fams2-p1-rev-disc1.txt`: 1 warning. With a
  prior on that intercept the fit is quiet and its standard errors are
  finite (test).
- m3: the findings errors are corrected above ("eight of the sixteen";
  `zi_(Intercept)` -2.61 and the hurdle `x` +2.09 added; the 50-to-100
  steps statement limited to points away from `m`; the `phi` 1e4 claim
  withdrawn).
- m4: the incomplete-beta test now covers shapes 0.3 to 8e5, `x` at
  `m`, at `m (1 +- 1e-3)`, 2 and 5.5 sd either side, and the blends at
  `a = 3e4`, `b = 7e4`, with a tolerance in ulps of the prefactor's
  scale. The absolute tolerances are replaced: the moment checks by
  `integrate()`'s own `abs.error` and ulps times the count, the
  probability sums by ulps relative to the value, the expected
  category by ulps times its own size.

Not changed: `RTMB::pbeta()` (with `pbinom()`, `pnbinom()`) and
`RTMB::dbeta()` are recorded for the user under "Defects found and not
fixed", not filed.

## Punch round 2

The re-check (`dev/reviews/2026-09-29-fams2.md`, "Re-check after punch
round 1") found the lane MERGEABLE, with the round-1 build rebuilt and
the new tests seen failing there on exactly the four punch items, and
left five minors, n1 to n5. Scripts are `dev/fams2-p2-*.R` and the
re-check's own `dev/fams2-rev2-*.R`, rerun on the punch-2 build with
output in `dev/fams2-p2-rev2-*.txt`.

### n1. RTMB::pbeta() at the mean

`dev/fams2-p2-tie.R` (output `.txt`): `RTMB::pbeta()`'s gradient,
Hessian and third derivative are all `NaN` at `x == a / (a + b)`
exactly, at every pair tried from (150, 150) to (1e7, 1e7), and
finite a relative 1e-15 away; at `(a + 1) / (a + b + 1)` they are
finite. The mean is where the pbeta branch of `log_ibeta_half()` is
read hardest, so `log_ibeta_half(0.3, 300, 700)` had a `NaN` gradient,
and so did an xbeta pair of rows (0, 0.3) at `mu == q == 1/6`.

The branch is now `log_pbeta_ad()`: below the mean it uses
`I_x(a, b) = I_x(a + 1, b) + x^a (1 - x)^b / (a B(a, b))`, whose
`pbeta()` ties at `(a + 1) / (a + b + 1)` instead. The two forms are
blended over the middle half of the gap between the two ties, each
read at `x` clamped out of its own tie. Both are exact, so the value
does not move. Routing the tie to the fraction, as suggested, would
not have worked: the fraction needs O(sqrt(a + b)) steps there, which
is why the pbeta branch exists.

- `dev/fams2-p2-pthird.R` (output `.txt`): `log_pbeta_ad()` over the
  20000 random points of `dev/fams2-p1-pthird.R` (shapes 150 to 1e7,
  within 6.5 sd), 0 non-finite derivatives to third order, value within
  3.5e-13 of `stats::pbeta()`; at the mean, at the second tie and at
  both blend edges for seven shape pairs, all finite, within 2.6e-15
  (4.3e-13 at (1e7, 1e7)).
- `dev/fams2-rev2-mean.R` sections 2 and 3 (output
  `dev/fams2-p2-rev2-mean.txt`): `log_ibeta_half()` at the mean has
  finite Hessian and third derivatives at all six pairs, and the xbeta
  rows at `mu == q == 1/6` have finite third derivatives.
- Tests: the gradient and third derivative at (0.3, 300, 700) and at
  the xbeta case, and the mean `a / (a + b)` added to the grid.
- The findings sentence "Its NaN points are at small shapes" is
  corrected above. For RTMB upstream: the gradient fails at the mean,
  not only the third derivatives.

### n2. Accuracy at the largest shapes

The re-check's grid put the error at 3.2e-13 below shapes of 1e4,
5.5e-12 to 1e5, 4.2e-11 to 1e6 and 6.5e-10 to 1e7, with a 4e-10 value
step at the 6 sd edge at (1e6, 3e6). Both came from rounding, not from
the method, and the fix was cheap, so it is fixed rather than stated:

- The prefactor `a log x + b log(1 - x) - log B(a, b)` cancels terms of
  size `a + b` down to a few units. Once both shapes pass 12 it is now
  formed around the mean `mu = a / (a + b)` as
  `a log1p((x - mu) / mu) + b log1p(-(x - mu) / (1 - mu)) + C(a, b)`,
  with `C = log(a b / (2 pi (a + b))) / 2` plus Binet remainders from
  `lgamma_binet()`, exact from Stirling's formula; the sum is stationary
  in `mu`, so the rounding of `mu` drops out (`log_ibeta_pre()`).
- The complement's odd terms `1 + e_odd(k)` cancelled at a large first
  shape; they are now formed from `xc = 1 - y` directly, with the large
  terms removed algebraically (`log_ibeta_cf(xc = )`). This was the
  last 1e-11, at a small shape beside one of 1e7, one sd above the mean.
- The "step" at 6 sd was not a step: past 6 sd the value is the
  fraction alone, and its prefactor rounding showed as noise of that
  size between two points 1e-9 sd apart.

On the punch-2 build (`dev/fams2-rev2-ibeta.R`, output
`dev/fams2-p2-rev2-ibeta.txt`; 4867 points, shapes 1e-3 to 1e7), the
error against `stats::pbeta(log.p = TRUE)`, relative and floored at
one, by the larger shape: 1.4e-15 to 1, 3.1e-14 to 1e2, 3.4e-14 to
1e4, 4.9e-14 to 1e5, 1.7e-13 to 1e6 and 3.0e-13 to 1e7; at `x == m`
2.6e-14. Against Rmpfr at the eight worst points (all at shapes of
1e6 to 1e7, within 6 sd of the mean) ours is within 3.9e-13 and
`stats::pbeta()` within 1.0e-13. Derivatives to third order finite at
all 4867. The
edges (section C): the value jump matches `g dx` to 5.2e-13 at every
edge, 6 sd at (1e6, 3e6) included (was 4e-10).
`dev/fams2-p1-sweep.R pkg` (output `dev/fams2-p2-sweep-pkg.txt`): max
5.45e-13, q99 1.8e-14, 3.3e-13 at `x == m`; gradient within 6.7e-14 of
`RTMB::pbeta()`'s AD at the ten worst stencil points.
`dev/fams2-p1-joins.R pkg` (output `dev/fams2-p2-joins.txt`): worst
1.6e-13, all edges finite. The test's value tolerance is now 1.4e-12
relative, floored at one, in place of ulps of the prefactor's scale.

Cost (`dev/fams2-rev2-mean.R` section 4, output
`dev/fams2-p2-rev2-mean.txt`): 14.5 us a row, 1.04 times true round 1
(13.9 us) and 9.1 times `log(RTMB::pbeta())` (1.6 us); the punch-1
build was 10.3 us. The extra 4 us buys 3e-13 at shapes of 1e7 in place
of 6.5e-10, and finite derivatives at the mean.

`dev/fams2-p1-cost.R` and `dev/fams2-p1-sides.R` were written against
the punch-1 `log_ibeta_cf()` signature and no longer run; their output
is kept.

### n3. The kappa check

SUPERSEDED IN PART by the final check's K1 (below): the some-rows rule
described here was dropped. `xbeta_fit_check()` now also warns, with
no row at 0 or 1, when `kappa` is below 1e-6 at SOME row ("kappa ran
to 0 at some rows (down to X)"), and, for the log link, when a
coefficient of `kappa` has a
standard error above 10 on the log scale or none that is finite ("the
standard error of kappa_(Intercept) is 96 on the log scale"). The
standard errors come from `sdr_of()`, which the fit caches for
`summary()`, and are read only when no row is at 0 or 1
(`dev/fams2-p2-kappa.R`, output `.txt`: a fit of 400 rows took 50 ms
with the check's `sdreport()` and 60 ms with a row at 0 and none).

On the re-check's cases (`dev/fams2-p2-kappa.txt`):

- the seed-7 miss (n 2000, kappa 3.4e-5, SE 96): warns;
- `kappa ~ x` at seed 31 (row kappa 1.3e-10 to 0.024): warns;
- the eight Beta() data sets whose `kappa` the interior shape placed
  (SE 1.35 to 6.91): quiet;
- 3b, `kappa` truly 0.3 and 0.5 (SE 2.06 and 3.15): quiet;
- A1 (kappa 65.8 on the ridge, SE `NaN`): NOW WARNS, "the standard
  error of kappa_(Intercept) is not finite on the log scale", beside
  the optimizer's own "false convergence". Both statements are true
  there; the message no longer claims `kappa` ran to 0 and recommends
  only a comparison with `Beta()`.

The threshold of 10 sits between the largest placed case (6.91) and
the smallest unplaced one the SE rule has to catch (22, the intercept
of `kappa ~ x` at seed 31).

### After the final check

- K1 (blocker). The some-rows rule warned on well-identified fits: a
  real slope in log kappa (-1 + 3x, x in (-4.7, 0), phi 200, n 2000,
  seeds 402 and 403) runs some rows to 2e-10 while the kappa standard
  errors are 0.26 and 1.38 (and 0.24 and 0.82), 28 to 31 log-likelihood
  units above `Beta()` (the reviewer's `dev/fams2-rev3-kappa.txt`,
  3b). The rule is dropped; the check keeps `max(kappa) < 1e-6` and
  the standard-error rule, which alone catches seed 31's `kappa ~ x`
  (intercept SE 22) and seed 7 (96). Seeds 402 and 403 are now
  must-stay-quiet cases in the test; seed 31 now warns through its
  standard error.
- m1, recorded, not fixed. Inside the tie blend of `log_pbeta_ad()`
  the third derivatives are off by up to 0.12 relative at (1e6, 3e6)
  (1.6e-2 at `(a + 1) / (a + b + 2)`, 7.7e-5 at (3e4, 7e4)), and the
  second by up to 1.3e-6, against 3e-9 for plain `log(RTMB::pbeta())`
  (the reviewer's `dev/fams2-rev3-near.txt`). The blend weight's k-th
  derivative scales as `gap^-k`, with `gap = b / ((a + b)(a + b + 1))`
  about `1 / (a + b)`, and it multiplies the rounding difference of
  the two exact forms. Widening does not fix it: the blend cannot
  leave the gap between the two ties, so the most it can gain is the
  middle four-fifths in place of the middle half, `(0.5 / 0.8)^3`, a
  factor of 4 (0.12 to 0.03). A shift of `a` by `j` moves the second
  tie `j` gaps away but costs `j` terms, and a width of one sd needs
  `j` near `sqrt(a + b)`. The band is about `1 / sqrt(a)` sd wide, so
  few rows land in it; values and first derivatives are unaffected. It
  is stated in `log_pbeta_ad()`'s comment.
- m2. The value bound at shapes of 1e7 is 7.2e-13 against Rmpfr, not
  3.0e-13: the reviewer's scan across both ties
  (`dev/fams2-rev3-clamp.txt`) finds 7.2e-13 at (1e7, 1e7) just below
  the tie, where `stats::pbeta()` itself is 2.9e-13 out. The 3.0e-13
  above was against `stats::pbeta()` on the grid, which does not
  include that point. Corrected in the code comment.

### n4. The disc warning's wording

With a prior on the thresholds (class `Intercept`) or on `mu`'s
coefficients and none on the disc intercept, the warning now says
"only the priors on the thresholds or on the coefficients of mu place
it" instead of that nothing has a usable standard error. The re-check's
section 4 on the punch-2 build (`dev/fams2-p2-rev2-guards.txt`): every
row warns or stays quiet as before; the wording differs only in the
threshold-prior row. Test added.

### n5. Recorded: xbeta at phi 2e4 stops short without a warning

The re-check's `dens` B at `phi` 2e4 (`dev/fams2-rev2-rerun-dens.txt`):
the fit ends 2.4e-5 below the reference optimum in logLik, on the flat
mu, kappa and phi ridge that few zeros leave (mu 0.152 against 0.158),
without a warning; the density at its estimates agrees with the
reference to 2.4e-8. It is the model's geometry and the optimizer's
tolerance, not the density. Not changed.

### Seen failing

The punch-2 pins were seen failing on the punch-1 build by the
re-check itself, which is the build they describe: the `NaN` gradient
at (0.3, 300, 700) and at `mu == q == 1/6`
(`dev/fams2-rev2-mean.txt`), no warning for the seed-7 case or
`kappa ~ x` (`dev/fams2-rev2-kappa9.txt`, `dev/fams2-rev2-guards.txt`),
the "no usable standard error" wording under a threshold prior
(`dev/fams2-rev2-guards.txt`), and 6.5e-10 at shapes of 1e7, above the
test's new 1.4e-12 (`dev/fams2-rev2-ibeta.txt`). The punch-1 build was
not kept, so the test file was not rerun against it.

## Validation

`dev/fams2-validate.R`, output `dev/fams2-validate.txt` (seeds in the
script). References are written there from `stats::dbeta()`,
`stats::pbeta()`, `lbeta()`, `lchoose()` and the link CDFs, after
brms's Stan functions, separately from the package code.

```
== 1. xbeta log-likelihood vs dbeta/pbeta ==           (n = 800: 8 zeros, 90 ones)
  optimum      rel diff 1.84e-15;  perturbed 1  2.53e-15;  perturbed 2  1.81e-16
  link = probit / cloglog / cauchit at their optima: 5.52e-15 / 5.33e-15 / 4.23e-15
  grid of 45 (eta, log phi, log kappa) x 5 responses: off-tape 7.57e-14, on-tape 7.59e-14
== 3. zero_inflated_beta_binomial vs lbeta/lchoose ==
  optimum 1.37e-16;  perturbed 1  4.02e-16;  perturbed 2  2.27e-15
== 4. zero_inflated_beta_binomial vs glmmTMB(betabinomial, ziformula) ==
  fixed:  logLik diff 3.26e-09;  max |coef diff| mu 9.09e-07  zi 5.65e-06  log phi 1.56e-06
  (1|g):  logLik diff -1.14e-08; max |coef diff| mu 1.33e-09  zi 1.03e-08
== 5. hurdle_cumulative vs hand-written density ==
  logit, probit, cauchit, cloglog at the optimum and a perturbed point: all <= 9.94e-16
  disc ~ 0 + z at its optimum: 1.99e-16
== 6. hurdle_cumulative factorization identity (glm + MASS::polr) ==
  logit    residual  1.32e-09  max |zeta diff| 4.61e-06
  probit   residual -1.04e-10  max |zeta diff| 1.28e-06
  cloglog  residual  1.31e-09  max |zeta diff| 2.60e-06
  cauchit  residual  2.94e+00  (polr did not converge; 6b)
== 6b. cauchit: MASS::polr against ordinal::clm and frmtmb ==
  reference ordinal log-likelihood at frmtmb's estimates -653.9340558458, at clm's -653.9340558557
  max |threshold diff| frmtmb vs clm 2.28e-05, |b_x diff| 1.38e-06
== 7. brms:::log_lik_hurdle_cumulative() at the frmtmb optimum ==
  logit diff 6.82e-13;  probit diff 1.14e-12
== 8. fitted() against brms's posterior_epred formulas ==
  xbeta 9.50e-16 (rel);  zero_inflated_beta_binomial 2.60e-16 (rel);
  hurdle_cumulative 2.50e-16 (abs, probabilities)
== 9. simulate() against the fitted distribution (seed 109, 400 draws per row) ==
  xbeta  P(Y = 0) z -1.49;  P(Y = 1) z -1.31;  sum of Y z -0.89;  mean (Y - m)^2 z -0.57
  zibb   counts 0..20: chi-square 14.2 on 20 df, p 0.822
  hurdle categories 0..4: chi-square 0.41 on 4 df, p 0.982
```

How to read it:

- Sections 1, 3, 5, 7 are MEASUREMENTS of the density at shared
  parameter points. Section 4 is agreement at the ML optimum with an
  independent fitter (glmmTMB's `betabinomial` is the same
  mean-precision parameterization).
- Section 6 is an IDENTITY: with separate predictors the hurdle
  likelihood is a bernoulli on `1{y = 0}` plus a cumulative model on
  `y > 0`, sharing no parameter, so the joint ML fit is `glm()` plus
  `polr()`. The residuals are optimizer precision. For cauchit
  `MASS::polr()` stops 2.9 nats short; `ordinal::clm()` agrees with
  frmtmb to 2e-5 in the thresholds, and the reference density at
  frmtmb's estimates is 1e-8 above its value at clm's. clm's REPORTED
  logLik there (-653.93698) is 2.9e-3 below the reference density at
  its own estimates (-653.93406); that is clm's, not this lane's, and
  not investigated.
- brms's R-side `log_lik_xbeta()` needs betareg and
  `log_lik_zero_inflated_beta_binomial()` needs extraDistr; neither is
  installed, so section 7 covers the hurdle only. The Stan identity
  below covers all three.

### Against brms's compiled Stan program

`brms_lp_check()` with flat priors, four new rows in
`tests/testthat/test-brms-likelihood.R` (gated), fitted to
`grad_tol = 1e-6` (see "Decided"). Check A is `log_prob()` at frmtmb's
estimates against `logLik(fit)`; check B brms's gradient there. From
the full gated run of the file (`dev/fams2-suite-gated/`, all 44 blocks,
450 expectations pass):

    row                                             lp - logLik   max |grad|
    16d xbeta, kappa ~ x (rows at 0 and at 1)       -4.26e-14     1.14e-05
    16e zero_inflated_beta_binomial, zi ~ x          0            1.83e-04
    16f hurdle_cumulative, hu ~ x (logit)            0            1.33e-04
    16g hurdle_cumulative, probit, disc ~ 0 + x      0            1.13e-04

16f runs brms's `ordered_logistic` path, 16g its generic
`hurdle_cumulative_probit_lpmf` with a modeled `disc`.

### brms defect: hurdle_cumulative, logit link, modeled disc

brms's generic `hurdle_cumulative_logit_lpmf()` (`brms:::
stan_hurdle_ordinal_lpmf()`, the `inv_logit` branch) tests
`y == nthres + 2` for the top category, which never occurs; the top
category `y = nthres + 1` falls to the interior branch and reads
`thres[nthres + 1]`. The non-logit branch tests `nthres + 1` and is
right. brms takes the generic logit function whenever `disc` is
modeled or `cs()` is present (`use_ordered_builtin()`).
`dev/fams2-brms-hc-bug.R` (output `.txt`, seed 18, 43 rows in the top
category):

    == link logit ==
    lpmf used: target += hurdle_cumulative_logit_lpmf
    log_prob at a legal point: ERROR: ... vector[uni] indexing: accessing
      element out of range. index 4 out of range; expecting index to be
      between 1 and 3
    == link probit ==
    lpmf used: target += hurdle_cumulative_probit_lpmf
    log_prob at a legal point: -466.1419

So brms 2.23.0 cannot sample that model at all. frmtmb fits it (the
logit `disc ~ 0 + z` fit in validate section 5 agrees with the
reference density to 2e-16). RECORDED FOR THE USER, NOT FILED
upstream (coordinator's instruction, punch round 1). brms's R-side
`log_lik_hurdle_cumulative()` is correct.

### Recovery at a known truth

`dev/fams2-coverage.R` (one RDS per replicate in `dev/fams2-cov/`),
summarized by `dev/fams2-coverage-sum.R` into `dev/fams2-coverage.txt`:
400 replicates per family, seeds 1..400, n = 500, 95% Wald intervals
from `confint()` on the link scale. All 1200 fits ran, none errored or
warned. Every coverage's Wilson interval contains 0.95:

```
hurdle_cumulative  x 0.945 [0.918, 0.963]  hu_(Intercept) 0.953 [0.927, 0.969]
                   hu_z 0.945 [0.918, 0.963]  tau_raw_1 0.945 [0.918, 0.963]
                   tau_raw_2 0.950 [0.924, 0.967]  tau_raw_3 0.968 [0.945, 0.981]
xbeta              (Intercept) 0.958 [0.933, 0.973]  x 0.963 [0.939, 0.977]
                   phi_(Intercept) 0.940 [0.912, 0.959]
                   kappa_(Intercept) 0.963 [0.939, 0.977]  kappa_x 0.948 [0.921, 0.965]
zibb               (Intercept) 0.945 [0.918, 0.963]  x 0.942 [0.915, 0.961]
                   phi_(Intercept) 0.955 [0.930, 0.971]
                   zi_(Intercept) 0.968 [0.945, 0.981]  zi_x 0.950 [0.924, 0.967]
```

The mean estimates are not all within Monte Carlo error of the truth
(`bias/mcse` in the output, the distance in Monte Carlo standard errors
of the mean over 400): eight of the sixteen are within 1 (corrected
in punch round 1; round 1 said nine), and the largest are xbeta
`kappa_x` +3.27, hurdle `tau_raw_1` -3.16 and `tau_raw_2` +2.74, zibb
`zi_(Intercept)` -2.61, the two `log phi` +2.45 (zibb) and +2.29
(xbeta), and the hurdle's `x` +2.09 (`dev/fams2-rev-cov-cmp.txt`
of the review). Every one is below
0.17 of the estimator's own standard deviation (largest: `kappa_x`,
0.026 against sd 0.159), and the directions are the ML small-sample
ones: an ordinal slope and the threshold spread inflated
(`x` 0.8126 against 0.8), a precision overestimated. Read as O(1/n)
bias at n = 500, which the coverage already absorbs; not investigated
further.

### frmtmb.sample

`dev/fams2-sample-smoke.R` (output `.txt`, seed 20260929, n = 200, 2
chains of 600): `log_lik()` against the reference density at 20 draws,
max relative difference 4.66e-15 (xbeta), 1.31e-14 (zibb), 9.30e-16
(hurdle); `posterior_epred()` is `ndraws x n x 5` for the hurdle,
`posterior_predict()` codes 0..4; `conditional_effects()` and `loo()`
run. The zibb run reported divergences and an R-hat of 1.88 from Stan
under frmtmb.sample's flat default on `zi` and `phi`; that is the
existing policy (lane fams, "Decided not to do"), not a family defect,
and a short chain.

## The brms-suite ports

The gated `test-brms-suite-families.R` on the lane build, in both
packages (`dev/fams2-port-flips.R`, uncapped, output
`dev/fams2-port-flips-<pkg>.txt`): 27 assertions now HOLD where the
committed verdict says `cannot transfer`, and nothing else fails:

    families:52 53 54 55 56 57 58 59   zero_inflated_beta_binomial (8 of 8)
    families:60 61 62 63 64 65 66 67   hurdle_cumulative (8 of 8)
    families:68 69 70 71 72 73 74 75 76 81 113   xbeta (11 of 16)

65, 66 and 67 (and 57, 58, 59) were VACUOUS before: their
`expect_error()` passed on "could not find function"; they now hold on
frmtmb's own refusal. The five xbeta rows that still do not hold read
brms-only fields: `families:107` `$closed`, `108` `$ybounds`, `109`
`$type == "real"` (frmtmb says `"continuous"`), `117`/`118` `$prior`
(a Stan prior string) and `122` `$include` (a Stan file name). No
frmtmb family carries those fields, and giving one family brms's
internals would not be a port. `priors:14` still fails: it needs
`sratio(threshold = "equidistant")` (not done, below); `cse()` alone
would not flip it.

The generated files are regenerated from `dev/brmsport-verdicts.tsv`,
which this lane does not edit, so both gated `test-brms-suite-families.R`
runs report `fail=27` until the consolidating session sets those 27
rows to `pass` and regenerates:

    RESULT frmtmb test-brms-suite-families.R pass=57 fail=27 err=0 skip=0
    RESULT frmtmb.sample test-brms-suite-families.R pass=57 fail=27 err=0 skip=0

## Decided not to do

- **`threshold = "equidistant"` and `"sum_to_zero"`** on any ordinal
  family. Not small: equidistant thresholds make the threshold COUNT a
  fact of the data that every consumer of `tau_raw` (five densities,
  the simulators, `ord_tau_from_raw()`, `thres_ncat()`, `ordinal_ncat()`,
  `thres(x = , gr = )`, the naming in `variables()`/`fixef()`, the
  threshold priors) reads today from the raw vector's length, and brms
  gives the spacing its own prior class `delta`, which frmtmb's prior
  machinery and frmtmb.sample do not have. `sum_to_zero` has a flat
  direction under maximum likelihood (brms declares `nthres` free
  thresholds for `nthres - 1` degrees of freedom and relies on the
  prior). `hurdle_cumulative()` takes `threshold` for brms's signature
  and refuses the other two by name; the other four constructors take
  no `threshold`, as before. `priors:14` stays blocked on this.
- **`disc` on the other four ordinal families** (lane links2's SEAM).
  `hurdle_cumulative()` has it because brms's constructor does and the
  generalized probability path reads every dpar; the four would need
  their densities, the grouped `thres()` densities and `ord_cat_probs()`
  rewritten. The route is now clear: `fixed_dpars = list(disc = 1)`
  plus `disc` in the density.
- **`thres(gr = )` for the hurdle.** brms fits it; the grouped
  densities would need the hurdle term. Refused by name.
- **`cens()`/`trunc()` for zibb.** brms allows both; frmtmb's
  `beta_binomial()` has no CDF either, so both are refused with the
  existing message, as for the sibling. brms allows neither for xbeta.
- **brms's `kappa`, `phi`, `zi`, `hu` default priors on the sample
  route.** frmtmb.sample's policy for dispersion and gate dpars (lane
  fams) is unchanged; only the disclosure now names `kappa`.
- **The simulate-density header count** (`test-simulate-density.R`
  line 55) was already stale per lane fams and was not re-derived; this
  lane adds 3 rows to that tier (xbeta, zibb, and a hurdle block).
- **Mixtures of hurdle_cumulative** are refused by the existing
  "mixture() does not support component family" message, as for the
  four ordinal families; brms allows ordinal mixtures.
- **Speed** was not benchmarked. No existing objective's arithmetic
  changed (the ordinal generalization is off the tape).

## Defects found and not fixed

1. brms 2.23.0 `hurdle_cumulative_logit_lpmf()` off-by-one (above):
   it tests `y == nthres + 2` where the top category is
   `nthres + 1`. Recorded for the user; not filed.
2. RTMB 2.0 (recorded for the user; not filed): `pbeta()`, `pbinom()`
   and `pnbinom()` have non-finite third derivatives at ordinary
   points (the review's `dev/fams2-rev-pthird.R`, seed 7: 47 of 2000,
   96 of 472 and 196 of 2000), and `dbeta()` gives a `NaN` gradient
   with its first argument on the tape once the shape sum passes
   about 1e3 (`dev/fams2-p1-rev-thresh.txt`: finite at 1046, `NaN`
   from 1047; `dev/fams2-p1-dbeta.txt`: finite at 96 of 168 points).
   No frmtmb family calls any of them now; the exposure is nonlinear
   formula bodies (through `nl_rtmb_shadow`) and custom families.
3. brms scores the expected ordinal category by column position, which
   is one off for a hurdle model (above); frmtmb departs.
4. `ordinal::clm(link = "cauchit")`'s reported logLik differs from the
   cauchit likelihood at its own estimates by 2.9e-3 (section 6b);
   `MASS::polr(method = "cauchit")` stops 2.9 nats short with
   convergence code 0. Not ours; noted because a reader of section 6
   would otherwise suspect frmtmb.
5. frmtmb.sample: zi/hu/phi flat on the sample route produced a
   poorly mixed zibb posterior on a short run (above). Policy, lane
   `sampfix` territory.

## Tests

New files:

- `tests/testthat/test-xbeta-zibb-hurdle-cum.R` (22 blocks, 406
  expectations after punch round 1, 602 after punch 2 (604 after the
  final check's K1 pins), which added
  derivative checks to the incomplete-beta grid and the n1, n3 and n4
  pins): constructors and brms's link
  sets; refusals (links, responses, `threshold`, `cens()`, `trunc()`,
  `cs()`, `thres(gr = )`, osa); the three densities against the
  references on and off the tape; `log_ibeta_half()` against
  `stats::pbeta()` on a grid of shapes 0.3 to 8e5 with `x` at, near
  and away from `m`, its `x`-derivative against `dbeta / I`, and its
  gradient and third derivatives across the blends at `a = 3e4`,
  `b = 7e4`; finite third derivatives at a point where `RTMB::pbeta()`
  has `NaN` ones (an expectation to FLIP when RTMB is fixed); xbeta
  fits of a precise response (`phi` 1e4 and 2e4 free, 2e5 with `kappa`
  held) against the reference logLik; glmmTMB agreement for zibb,
  fixed and `(1 | g)`; the glm + polr identity for the hurdle; xbeta
  and zibb with `(1 | g)`; `fitted()` against brms's epred formulas;
  the variance functions against summed and integrated moments;
  beta_binomial pearson residuals; the hurdle's 0..K coding through
  factor responses, `simulate()`, `predict()`, `conditional_effects()`
  (the expected category by codes, and `P(Y = 0)`'s band equal to
  `hu`'s own); the `kappa` warning and its no-false-alarm cases,
  including the review's A1 data; the `disc` intercept warning and its
  two quiet cases; `thres(x = )`; `disc ~ 0 + z`; `cse()`; the
  registry.
- `extensions/frmtmb.sample/tests/testthat/test-xbeta-zibb-hurdle-draws.R`
  (3 blocks, 12 expectations): `log_lik()` against the references at
  draws, the hurdle's expected category by codes on the draws surface,
  the `kappa` disclosure.

Changed: `test-brms-likelihood.R` (rows 16d to 16g),
`helper-brms.R` (`brms_ord_thresholds()` reads hurdle_cumulative as
ordered), `test-simulate-density.R` (xbeta and zibb rows, a hurdle
block, the coverage list).

Pins seen failing on rellib-r3 (`dev/fams2-pins.txt`,
`dev/fams2-t-new-base.log`): the beta_binomial pearson residual
(behavioural: "has no variance function"), `cse()` ("could not find
function"), and the whole new file (pass=1 fail=1 err=19; the families
do not exist there). The round-1 `kappa` warning test was seen failing
on the lane build before the check existed. The punch-round B1 and B2
pins are seen failing on the round-1 density in
`dev/fams2-p1-pins.txt` (under "Punch round 1"). The `disc` intercept
test was seen failing before `hurdle_cum_fit_check()` existed (the
review's `dev/fams2-rev-disc1.txt`: 0 warnings).

Punch 2 reran, on the punch-2 build: the core ungated suite, 183
files, pass 12969, fail 0, err 0, skip 156, warn 0
(`dev/fams2-suite-core/RESULTS.txt`; `test-xbeta-zibb-hurdle-cum.R`
602 and `test-bracket-access.R` 33 rerun after a `[[` fix that changed
no value); gated `test-brms-families.R` 1102 and
`test-brms-likelihood.R` 450, fail 0 (`dev/fams2-suite-p2-gated/`);
and frmtmb.sample's `test-xbeta-zibb-hurdle-draws.R` 12, fail 0. The
table below is the punch-1 build's, unchanged in punch 2 except the
core ungated row.

Runs on the punch-1 build, one file per process, `NOT_CRAN=true`
(`dev/fams2-suite-*/RESULTS.txt` by `dev/fams2-suite.sh`, summed by
`dev/fams2-sum.sh`):

    core, ungated:        183 files, pass 12773, fail 0, err 0, skip 156, warn 0
    core, gated test-brms-*.R: 21 files, pass 3743, fail 27, err 0, skip 0, warn 0
    frmtmb.sample, ungated: 37 files, pass 1958, fail 0, err 0, skip 4, warn 0
    frmtmb.sample, gated test-brms-*.R: 8 files, pass 587, fail 27, err 0, skip 0, warn 0
    frmtmb.coupling, ungated: 11 files, pass 542, fail 0, err 0, skip 5, warn 0
    frmtmb.eam, ungated:    29 files, pass 1733, fail 0, err 0, skip 3, warn 0
    frmtmb.latent, ungated: 10 files, pass 359, fail 0, err 0, skip 2, warn 0
    frmtmb.learn, ungated:  15 files, pass 429, fail 0, err 0, skip 13, warn 0
    frmtmb.ode, ungated:    11 files, pass 547, fail 0, err 0, skip 1, warn 0
    frmtmb.spline, ungated: 15 files, pass 553, fail 0, err 0, skip 1, warn 0

The five extensions beside frmtmb.sample were installed from this
worktree into the lane library, so each loaded the lane's core
(`lib:` line of each log, `dev/fams2-suite-ext-<name>/`). The two
`fail 27` are the 27 ported rows that now hold (above), in
`test-brms-suite-families.R` of each package; `dev/fams2-port-flips.R`
rerun on the final build finds 27 "now HOLDS" and no other failure in
each (`dev/fams2-port-flips-<pkg>.txt`). The 156 ungated core skips
are the gated tiers (`test-brms-*`, `test-bcm-*`,
`test-drmtmb-agreement.R` and others behind `FRMTMB_BRMS_FIT_TESTS`),
none in a file this lane touched. The ungated sample skips:
`test-loo.R` 2, `test-sampling-ported.R` 1, `test-scale.R` 1, their
own gates.

`?frm` has no family list to update; the families are listed in
`?frmtmb-families` (two new sections and parameters) and
`?frmtmb-links` (the tables). `vignette("brms-migration")` gains a row
in its port table and a bullet under "What changes".

## Version

A minor bump for frmtmb (three exported families and `cse()`) and a
patch bump for frmtmb.sample (display and disclosure only). frmtmb.sample
needs no new core floor: `ce_cat_mean()` reads the family field
`extra_cat`, which older cores simply do not set.

## R CMD check --as-cran

frmtmb rerun on the punch-2 build (frmtmb.sample, unchanged since, on
the punch-1 build). Built and checked in
`dev/fams2-check/` (logs kept, tarballs and
`.Rcheck` removed), `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and
TinyTeX on PATH, `R_LIBS` = lane library, rellib-r3, user library:

    frmtmb         Status: 1 NOTE  (HTML manual: V8 unavailable, the expected one)
    frmtmb.sample  Status: OK
