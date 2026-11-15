# Lane wcheck2: weights check residues (round of 2026-11-13)

Base: ec8178d3 (frmtmb 0.82.0), branch wt-r83-wcheck2, library
`C:/Users/adf44/source/r/wt-wcheck2-lib`. "Base" or "before" is
rellib-r21. Brief: `round-briefs/r083-wcheck2.md`. Nothing committed.

Status: done. This file is the record.

## The brief's item numbers

The brief's item 2 is labeled 4977 but describes 4949 (asym_laplace,
"the other side costs 1/2"); 4977 at ec8178d3 is the treatment contrast
under an Intercept prior, and the brief's item 4 ("the pre-existing
MINOR") is 4985. This record uses the backlog's lines: 4949, 4968, 4977,
4985, 4993, 4998. All six are covered below.

## What changed

Core only. No extension changed. This section is the first round; the
section "Fix round" at the end supersedes it where they differ (the
pair rule's units and cost, and the new asym_laplace refusal in
`R/ald-laplace.R` with its compat rows).

1. `R/weights-unbounded.R`
   - 4968: `hurdle_lognormal` (sigma), `hurdle_gamma` (shape, rate 1/2),
     `zero_inflated_beta` and `zero_one_inflated_beta` (phi, rate 1/2)
     join `unbounded_spike_dpar` and `unbounded_spike_rate`. A new table
     `unbounded_point_mass` (0 for these four and for
     zero_inflated_asym_laplace, 0 and 1 for the zero-one-inflated beta)
     replaces the ZIAL special case: an uncensored row at the point mass
     is left out of the excess. The refusal says how many rows sit at
     the mass and why they do not count.
   - 4977: `unbounded_cell_pair()`, tried after the one-cell rule: with
     a centered intercept, cells of a scale on the log link, every
     moving column's class `b` prior free (none, flat, student_t, no
     bound), it refuses when `E_j / n_j > W_k / n_k` for two disjoint
     cells j and k (or one cell and the reference level, when the cells
     do not overlap). The message names both cells and the arithmetic.
   - 4985: `check_weights_unbounded()` resolves the prior
     (`resolve_prior_input()`) before it refuses, on the MAP route, so a
     prior that the fit would refuse gives its own error first.
2. `R/fit.R`, `check_re_structure()` only (4993): a new helper
   `re_check_term()` writes the term of a multivariate model as its
   formula does, followed by the response ("`1 | obs` of response y");
   both the observation-level and the one-level warning use it. The
   univariate text is unchanged.

## 4968: the hurdle and zero-inflated continuous families

The proof is the table's own with one more step. A row at the point
mass has density `hu^w` (`zi^w`, `(zoi coi)^w`): no mu and no spike
dpar, so its own effect sees only its prior, the Laplace step for it is
exact and constant, and the row neither grows nor curves. The other rows
are the lognormal, Gamma and Beta densities times `(1 - hu)^w`, a
constant in (0, 1) at any interior hu: the table's proofs for those
families with the excess summed over them alone. The mean's links (log,
logit, probit, cloglog, identity on the log scale) reach every positive
or (0, 1) response, so each row's own effect reaches its observation.
`hurdle_gamma` and `zero_one_inflated_beta` were not named in the item;
they are the same defect and the same proof, and both were measured to
run away, so they are in.

`Rscript dev/wcheck2-4968.R <base|lane> <hl|hg|zib|zoib> [quad]`, data
of `dev/reviews/wcheck-review-misc.R` part quad (seed 20261111, n = 80;
12 zeros for the hurdles, 10 for zib, 10 zeros and 6 ones for zoib).
Profile of the Laplace objective along the spike dpar's log, the rest
re-optimized, slope per unit toward the limit; "rate x E_nz" is the
prediction from the other rows' excess. Logs
`dev/wcheck2-out/4968-<fam>-<base|lane>.txt`.

| family, arm | E_nz | sum(w) - n | slopes | rate x E_nz | rellib-r21 | lane |
|---|---|---|---|---|---|---|
| hurdle_lognormal, all w 1.5 | 34 | 40 | 34 34 34 34 34 | 34 | logLik 968.909, log sigma -37.49, conv 1 | refused |
| hurdle_lognormal, other rows E +4, mass w 0.5 | 4 | -2 | 4 4 4 4 4 | 4 | 18.0545, -36.28, conv 1 | refused |
| hurdle_lognormal, other rows E -2, mass w 3 | -2 | 22 | -2 x5 | -2 | -141.305, conv 0 | same |
| hurdle_lognormal, other rows E -2, mass w 1 | -2 | -2 | -2 x5 | -2 | -108.569, conv 0 | same |
| hurdle_gamma, all w 1.5 | 34 | 40 | 17.3 17.1 17.0 17.0 17.0 | 17 | 7.03e24, log shape 84.6, conv 0 (gradient warning) | refused |
| hurdle_gamma, E +4, mass 0.5 | 4 | -2 | 2.08 2.03 2.01 2.00 2.00 | 2 | 6.03e10, conv 1 | refused |
| hurdle_gamma, E -2, mass 3 | -2 | 22 | -0.953 ... -1.000 | -1 | -141.939, conv 0 | same |
| hurdle_gamma, E -2, mass 1 | -2 | -2 | -0.953 ... -1.000 | -1 | -109.203, conv 0 | same |
| zero_inflated_beta, all w 1.5 | 35 | 40 | 17.7 17.9 17.6 17.5 17.5 | 17.5 | 2.79e10, log phi 52.87, conv 1 | refused |
| zero_inflated_beta, E +4, mass 0.5 | 4 | -1 | 1.83 2.48 2.11 2.02 2.00 | 2 | **0.544522, log phi 1.829, conv 0, no warning** | refused |
| zero_inflated_beta, E -2, mass 3 | -2 | 18 | -1.253 ... -0.997 | -1 | -42.8019, conv 0 | same |
| zero_inflated_beta, E -2, mass 1 | -2 | -2 | -1.253 ... -0.997 | -1 | -12.3085, conv 0 | same |
| zero_one_inflated_beta, all w 1.5 | 32 | 40 | 16.2 16.3 16.1 16.0 16.0 | 16 | 48973, conv 1 | refused |
| zero_one_inflated_beta, E +4, mass 0.5 | 4 | -4 | 1.83 2.40 2.09 2.01 2.00 | 2 | 597.236, conv 1 | refused |
| zero_one_inflated_beta, E -2, mass 3 | -2 | 30 | -1.248 ... -0.998 | -1 | -91.2092, conv 0 | same |
| zero_one_inflated_beta, E -2, mass 1 | -2 | -2 | -1.248 ... -0.998 | -1 | -34.2658, conv 0 | same |

The bold row was a silent wrong answer on 0.82.0: conv 0, no warning,
at a point the objective rises from by 2.00 per unit of log phi
(-fn 0.508 at log phi 2, 21.399 at 12). The "mass w 0.5" arms have
sum(w) < n and are unbounded; the "mass w 3" arms have sum(w) > n and are
bounded: the excess of all rows would decide both wrongly. Under
quadrature (`4968-<hl|zib>-quad-base.txt`) rellib-r21 ends the two
unbounded hurdle_lognormal arms at conv 1 with `frmtmb_quad_unsettled`
(479.136, -38.9559), zib the same (390.202, 35.4929); the bounded arms
fit as under Laplace. On the lane (`4968-<hl|zib>-quad-lane.txt`) the
two unbounded arms of each are refused under quadrature too (the exact
route counts censored rows at 0 and the mass rows not at all), and the
two bounded arms fit as on rellib-r21 (-141.305, -108.569; -42.8019,
-12.3085, conv 0, no warning).

Rounding: `Rscript dev/wcheck2-perturb.R lane` (the new test block's
data, x and the response times `1 + 2e-15 k`, k = 0..24; log
`dev/wcheck2-out/perturb-lane.txt`):

    hurdle_lognormal absent              code 0, 0 warnings x25
    hurdle_lognormal refused             REFUSED x25
    hurdle_gamma absent                  code 0, 0 warnings x25
    hurdle_gamma refused                 REFUSED x25
    zero_inflated_beta absent            code 0, 0 warnings x25
    zero_inflated_beta refused           REFUSED x25
    zero_one_inflated_beta absent        code 0, 0 warnings x25
    zero_one_inflated_beta refused       REFUSED x25

## 4949: censored rows cost 1/2 each; the 8.00 is the inner search

Misdiagnosed. The allowance's 1/2 per censored row that a constant mean
cannot bring to 1 is the Laplace approximation's own cost. The slope of
8.00 that the item measured is frmtmb's objective, which TMB evaluates
at inner points that are not the modes. No change to the rule.

Derivation. A right-censored row whose bound is above its mean m: its
own effect b goes past the bound, to where the log survivor
`log(1 - p exp((1 - p) u))`, `u = (c - m - b) / sigma`, has gradient
`b / sd^2`. Its curvature there is `(1 - p) / sigma` times that
gradient, so `sd^2 H = 1 + (1 - p) |b| / sigma`, and the row's Laplace
term is `-1/2 log(sd^2 H)`: 1/2 per unit of `-log sigma` while `|b|`
stays away from 0. Gaussian rows give the same (the 0.82.0 measurement
had 3.7 to 3.9 for 4 halves off 8).

Measurement: `Rscript dev/wcheck2-al-decomp.R base [both|right72]`
(logs `dev/wcheck2-out/al-decomp-base.txt`,
`al-decomp-right72-base.txt`). Each row's Laplace term is rebuilt from
TMB's inner values (`obj$env$last.par`) and matches TMB's objective to
the fourth decimal (row sum = TMB in every line), then recomputed with
each row's mode found by a 1-D maximization, and the exact marginal by
adaptive quadrature.

"both" (the review's construction: seed 20261111, 36 left, 36 right, 8
uncensored, weights 1), along TMB's own profile path:

    log sigma -14: TMB 95.5721 (row sum 95.5721) | exact inner mode -302.9174 | exact marginal -67.6054 | m 0.9274, quantile 0.8047, log sd 0.847
        left  rows: TMB mode u/sigma in [1.382e+04, 1.382e+04], curv*sd^2 max 0, term sum -0.0615 | exact mode u in [19.79, 19.79], curv*sd^2 in [1.21e+05, 1.21e+05], term sum -210.6848
        right rows: TMB mode u/sigma in [-6.352e+05, -6.352e+05], curv*sd^2 max 0, term sum -1.5549 | exact mode u in [-80.36, -80.36], curv*sd^2 in [3.7e+04, 3.7e+04], term sum -189.4213
      tmb   slope 33.23 16.61  8.19  8.00  8.00  8.00
      lap   slope   8.12  -4.87 -36.61 -26.90 -27.81 -27.96
      exact slope  0.32767  0.17872 -0.28418 -0.00155  0.00129  0.00051

"right72" (`dev/wcheck-cens-path.R`'s case: seed 20261021, 72 of 80
right censored, mean held 2 below the bound, quantile 0.5); the profile
of the Laplace approximation with exact inner modes, log sd re-optimized:

      tmb   slope 63.11  8.66  8.00  8.00  8.00  8.00
    own profile, Laplace with exact inner modes
      slope -19.4 -25.6 -27.5 -27.9 -28.0 -27.9
    own profile, exact marginal
      slope  2.48e-09 -1.75e-10 -1.35e-11 -1.90e-12 -3.00e-12 -1.39e-10

So with the modes found, each of the 72 censored rows costs 1/2: 8 -
72/2 = -28 against -27.96 and -28.0. The exact marginal is flat
(unweighted rows: bounded). In "both" the profile of the Laplace
approximation with exact modes is noisy (my outer optimizer over three
parameters around an inner `optimize()`): -23.56, -23.12, -22.92 at log
sigma -6, -8, -10, then -29.23 and -25.68; no growth beyond 0.22 per
unit, against TMB's +8.00.

Why TMB's objective rises: at the TMB values each censored row's effect
sits 1.4e4 to 6.4e5 sigmas past its bound, where
`p exp((1 - p) u)` underflows to 0, so the computed curvature is 0 and
the row costs nothing. That is not a stationary point. At log sigma -14
(arithmetic on the logged values) a right row's effect is `(c - m) - u
sigma` = (1.0849 - 0.9274) + 6.352e5 exp(-14) = 0.685, sd^2 = exp(2 x
0.847) = 5.44, and the prior's gradient alone is -0.685 / 5.44 = -0.126,
the data's exactly 0. A Newton step from the
linear side of the CDF, where the curvature is the prior's 1/sd^2 alone,
jumps there. The fits end at conv 1 with a gradient warning (39.17,
1.43e6). Filed (below) as a defect of the objective, not of the check.

WITHDRAWN in the fix round: "The check reads the Laplace approximation,
which has its maximum, so it is right to let these through." The
review (F3) showed that for a grouped effect the exact-mode Laplace
approximation of asym_laplace never falls in sd and has no maximum, so
it is no quantity with a maximum either. The weights check still fails
open on these models, which is right for its own question; the
asym_laplace refusal of the fix round (section "Fix round") decides
them where its proof holds. Censored rows, as here, are left open by
that refusal.

Is an ordinary asym_laplace fit hit? `Rscript dev/wcheck2-al-group.R
base` (log `dev/wcheck2-out/al-group-base.txt`; the right72 data with
`(1 | g4)`, at each fit's own outer parameters):

    72 of 80 right censored, (1 | g4)
      fit: logLik -20.7150 conv 1, log sigma -1.6808, quantile 0.1103, log sd -0.7926
      at the fit: TMB -21.0547 (rebuilt from TMB's modes -21.0547) | exact inner modes -21.0539 | exact marginal -23.1654 | max |TMB mode - exact mode| 0.016
    no row censored, (1 | g4)
      fit: logLik -82.7387 conv 1, log sigma -1.6823, quantile 0.3072, log sd 0.3484
      at the fit: TMB -83.5369 (rebuilt from TMB's modes -83.5369) | exact inner modes -82.7225 | exact marginal NA | max |TMB mode - exact mode| 0.639

With censoring and groups of 4 the modes are close (8e-4 in the
objective). With no row censored the inner modes are off by up to 0.639
and the objective by 0.81: asym_laplace's piecewise linear log density
gives Newton no curvature at its modes. That is 3389's ground ("no
asym_laplace fit passes the gradient check"); filed with these numbers.

## 4977: two cells and a held centered intercept

The linear program. frmtmb puts a class `Intercept` prior on
`Ic = b0 + means_X'b` (`lp_center_offset()`: colMeans of the parametric
columns over the design's rows). Hold Ic, so the prior is a constant.
Let the log scale of cell j move by `s_j` per unit t. Holding Ic is
`sum_j n_j s_j = 0` (n_j the cell's design rows). A cell that shrinks
gains `E_j |s_j|` (its excess, less the censored rows' cost); a cell
that grows loses `W_k s_k`, because each uncensored row off a point mass
loses its weight per unit (density at a fixed point like 1 / sigma, own
effect's data curvature to 0) and a censored or point-mass row loses
nothing. Normalized by `sum n_j s_j^- = sum n_k s_k^+ = 1` the growth
is `sum_j (E_j / n_j) a_j - sum_k (W_k / n_k) b_k` over two probability
vectors, so its maximum is at one cell each: unbounded exactly when
some `E_j / n_j > W_k / n_k`. For two cells this is the review's
condition; for more it is the pair. The directions: two treatment
columns move with b0 held; with the reference level, b0 moves and every
other cell's column takes it back (valid when the cells do not
overlap). Every column that moves must have a free prior. The loss of
weight per unit is proved for scales (sigma toward 0) on the log link
only; a shape going to 0 loses at a family's own rate (Gamma 1,
inverse.gaussian 1/2, von_mises 0, weibull faster than any power), not
derived here, so shapes stay "unknown" (filed).

`Rscript dev/wcheck2-4977.R <lib> <two|three>` (logs
`dev/wcheck2-out/4977-<two|three>-lane.txt`): the check's verdict, then
the Laplace objective along exactly that direction from the start, the
other parameters re-optimized (slope per unit t, where t moves log sigma
of cell j by -t / n_j and of cell k by +t / n_k). gaussian, `y |
weights(w) ~ x + (1 | obs)`, `sigma ~ 1 + fac`, `normal(0, 1)` on the
sigma Intercept. Two cells: the review's data (seed 20261111, n = 80);
three cells: seed 20261113, n = 90, 30 per level.

| case (down, up) | weights | G = E_j/n_j - W_k/n_k | slopes | lane |
|---|---|---|---|---|
| R1 (b, a) | a 0.8, b 1.05 | -0.75 | -0.706 -0.750 -0.750 -0.750 -0.750 | fit -94.2889, conv 0 (as 0.82.0) |
| R3 (b, a) | a 0.3, b 1.5 | 0.2 | 0.219 0.200 0.200 0.200 0.200 | refused (0.82.0: 612.75, conv 1) |
| edge + (b, a) | a 0.3, b 1.4 | 0.1 | 0.118 0.100 0.100 0.100 0.100 | refused |
| edge - (b, a) | a 0.3, b 1.2 | -0.1 | -0.0845 -0.1000 -0.1000 -0.1000 -0.1000 | fit 168.218, conv 1 |
| ref down + (a, b) | a 1.4, b 0.3 | 0.1 | 0.173 0.100 0.100 0.100 0.100 | refused |
| ref down - (a, b) | a 1.2, b 0.3 | -0.1 | -0.0274 -0.1000 ... | fit 192.043, conv 1 |
| 3 cells (c, b), b0 held | a 1, b 0.3, c 1.4 | 0.1 | 0.126 0.100 0.100 0.100 0.100 | refused |
| 3 cells (c, b) | a 1, b 0.3, c 1.2 | -0.1 | -0.077 -0.100 ... | fit 115.63, conv 1 |
| 3 cells (a, b), reference down | a 1.4, b 0.3, c 1 | 0.1 | 0.129 0.100 ... | refused |
| 3 cells (a, b) | a 1.2, b 0.3, c 1 | -0.1 | -0.0743 -0.1000 ... | fit 118.057, conv 1 |
| 3 cells (b, a), reference up | a 0.3, b 1.4, c 1 | 0.1 | 0.111 0.100 ... | refused |
| 3 cells (b, a) | a 0.3, b 1.2, c 1 | -0.1 | -0.0928 -0.1000 ... | fit 122.673, conv 1 |

The "-" arms are bounded along every pair direction and let through, as
they must be; they end at conv 1 far out (log sigma near -33, the
centered intercept near -16) because the bounded maximum is where the
normal prior on Ic and the cell's excess balance, the degenerate kind
the message describes. R3 with `student_t(3, 0, 2.5)` on b as well is
refused (the contrast is free); with `normal(0, 1)` on b it fits
(217.064, conv 0). R2 (no prior) is refused by the one-cell rule, as on
0.82.0.

## 4985: the prior's error first

`Rscript dev/wcheck2-4985.R <base|lane>` (logs
`dev/wcheck2-out/4985-<base|lane>.txt`; data of
`dev/reviews/wcheck-short-attacks.R` part couple, weights 1 + 8 / n):

| prior | rellib-r21 | lane |
|---|---|---|
| normal(0, 1), class b, coef Intercept, dpar sigma (invalid) | weights refusal | "names an intercept that this model centers" |
| the same, unweighted | the prior's error | the prior's error |
| class b, coef z, dpar sigma (invalid) | weights refusal | "Prior target not found" |
| student_t(3, 0, 2.5), class Intercept, dpar sigma (valid) | weights refusal | weights refusal |
| normal(0, 5), class b (valid) | weights refusal | weights refusal |
| no prior | weights refusal | weights refusal |

The resolution runs only on the refusal path, so a fit pays nothing for
it. frm_sample() resolves its priors before its own call of the check.

## 4993: the multivariate term

`Rscript dev/wcheck2-olre-label.R <base|lane>` (logs
`dev/wcheck2-out/olre-label-<base|lane>.txt`). The one-level warning had
the same defect, fixed in the same place:

| model | rellib-r21 | lane |
|---|---|---|
| univariate (1 \| obs) | `1 \| obs` gives every observation ... | same |
| mvbf, (1 \| obs) on y | `y 1 \| obs` gives every ... | `1 \| obs` of response y gives every ... |
| mvbf, (1 \| obs) on y2 | `y2 1 \| obs` gives every ... | `1 \| obs` of response y2 gives every ... |
| mvbf, (1 \| one) on y | Grouping factor 'one' in `y 1 \| one` has a single level | ... in `1 \| one` of response y has a single level |
| mvbf, sigma ~ (1 \| one) on y | ... in `y sigma: 1 \| one` ... | ... in `sigma: 1 \| one` of response y ... |

## 4998: mislabeled

`Rscript dev/wcheck2-4998.R base` (log `dev/wcheck2-out/4998-base.txt`).
In `dev/wcheck-cases.R` part cens the loop over k leaves `d$ca` and
`d$yac` at k = 80, so the arm "asym_laplace 72 right censored (1 | g4)"
fits 80 of 80 rows right censored at the minimum. With 72 censored the
same model ends at logLik -20.715, conv 1 (asym_laplace's gradient,
3389). The all-censored model has no maximum likelihood estimate (3392's
case): gaussian gets "Standard errors are not available for 4 of 4
parameters ... the maximum likelihood estimate does not exist",
asym_laplace only "Starting value 0 for sigma is -Inf through its log
link" and "false convergence (8)" (`Rscript dev/wcheck2-allcens.R
base`, log `dev/wcheck2-out/allcens-base.txt`). Filed below.

| model | k = 72 | k = 80 |
|---|---|---|
| asym_laplace (1 \| g4) | -20.715, conv 1 | -0, log sigma -38.61, conv 1 |
| gaussian (1 \| g4) | -21.6509, conv 0 | 3.55e-15, conv 0, SE warning |
| asym_laplace, no effect | -21.4738, conv 1 | -2.69e-215, log sigma -498.4, conv 1 |
| asym_laplace (1 \| obs) | refused | -1.42e-14, conv 1 |

Not a weights-check case: no weights, and `(1 | g4)` is not
observation-level.

## Acceptance

`bash dev/wcheck2-accept.sh base lane`, then `Rscript
dev/wcheck2-table.R` (log `dev/wcheck2-out/accept/table.txt`):

    case table: 58 of 58 rows identical on rellib-r21 and the lane
    common: 26 lines, 0 differ (apart from lib:)
    routes: 27 lines, 0 differ (apart from lib:)

wcheck's 58-case table (`dev/wcheck-cases.R`, three parts), the
review's 24 guard-absent models (`dev/reviews/wcheck-review-common.R`)
and the route table (`dev/reviews/wcheck-short2-routes.R`) are
unchanged. No row of them is one this lane fixes: the 4998 arm stays as
it was (it is the all-censored model), and the arm "1 + fac, normal on
Intercept only" is bounded under the pair rule (E_b / n_b 0.307 against
W_a / n_a 1.29).

## Tests

`tests/testthat/test-weights-unbounded.R`, three new blocks: "a hurdle
or zero-inflated point mass neither grows nor stops it" (the four
families: refused with the mass at 0.5 and sum(w) < n; fitted at conv 0
with no warning with the mass at 3 and sum(w) > n), "two cells of sigma
move with the centered intercept held" (R3 refused, with student_t on b
refused and normal on b passed; the edge "-" case passed; three cells,
non-reference and reference pairs, each with its "-" twin passed), "a
prior's own error comes before the refusal". `test-diagnostics-ux.R`:
"the structure checks name a multivariate term as written".

    rellib-r21:  RESULT test-weights-unbounded.R pass=161 fail=0 err=3 skip=0 warn=2
    lane:        RESULT test-weights-unbounded.R pass=201 fail=0 err=0 skip=0 warn=0
    rellib-r21:  RESULT test-diagnostics-ux.R pass=131 fail=3 err=0 skip=0 warn=0
    lane:        RESULT test-diagnostics-ux.R pass=134 fail=0 err=0 skip=0 warn=0

(`dev/wcheck2-out/tf-wu-<base|lane>.txt`, `tf-dux-<base|lane>.txt`.) On
rellib-r21 every new block fails as behavior: a `frmtmb_fit` where the
refusal is expected (the hurdle block), NULL where the pair refusal is
expected, the weights refusal where the prior's error is expected, and
the "`y 1 | obs`" / "`y2 1 | one`" labels.

OpenBLAS 0.3.26 with LAPACK (`bash dev/ciharden-openblas.sh 0.3.26
lapack`, its check printed LAPACK 3.12.0; log
`dev/wcheck2-out/openblas-build.txt`): `test-weights-unbounded.R`
pass=201 fail=0 err=0 skip=0 warn=0, `test-diagnostics-ux.R` pass=134
fail=0 err=0 skip=0 warn=0 (`dev/wcheck2-out/tf-ob-*.txt`).
`dev/ciharden-out/` deleted. 0.3.32 not run: no new assertion sits near
a rounding-decided branch (the pair margins are 0.1, the mass arms'
excesses 4 and -2).

## Suites

One file per process, every suite, 3 at a time (`bash
dev/release/run-par.sh wcheck2-suite-par dev/wcheck2-out/jobs-all.txt 3`
with `NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS=true`, `FRMTMB_RELLIB` = the
lane library; log `dev/release/wcheck2-suite-par.log`, "WCHECK2-SUITE-PAR
ran 529 of 529"). Summed from its RESULT lines:

    frmtmb files=330 pass=27487 fail=0 err=0 skip=15 warn=0
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=33 pass=1857 fail=0 err=0 skip=3 warn=0
    frmtmb.latent files=13 pass=697 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=16 pass=545 fail=0 err=0 skip=2 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=91 pass=4603 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial files=8 pass=612 fail=0 err=0 skip=0 warn=0
    frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0

All 529 `lib:` lines load frmtmb 0.82.0 from `wt-wcheck2-lib`. The 29
skips are gates not set here: 15 "set FRMTMB_SCALE_TESTS=true", 13 "set
FRMTMB_DRMTMB_FIT_TESTS=true", 1 "set FRMTMB_FUZZ=true".

One session per suite, three orders (`bash dev/release/run-sessions.sh
wcheck2-suite-sess 3`, same environment; log
`dev/release/wcheck2-suite-sess.log`, "WCHECK2-SUITE-SESS ran 27 of
27"): all 27 sessions fail=0 err=0 warn=0. Core: forward pass=27487
skip=15, reverse 27652 skip=15, preload 27755 skip=17. Every extension's
three sessions give its one-file-per-process counts above. 194 LEAK
lines, every one tagged "[with load of ...]"; none untagged.

R CMD check --as-cran, core (`bash dev/wcheck2-check.sh`, built with
vignettes in `dev/wcheck2-check/`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`,
pandoc and TinyTeX on PATH, `R_LIBS` = lane library, rellib-r21, user
library; log kept as `dev/wcheck2-out/check-frmtmb.txt`, check tree
deleted):

    * checking examples ... [94s] NOTE
    frm_hazard_reads 3.24   0.11    5.05
      Running 'testthat.R' [59m]
     [59m] OK
    * checking HTML version of manual ... [44s] NOTE
    Skipping checking math rendering: package 'V8' unavailable
    Status: 2 NOTEs

The V8 NOTE is the expected one. The examples NOTE is elapsed time under
load (3.35 s of CPU, 5.05 s elapsed) on an example this lane does not
touch, the same NOTE the 0.82.0 round's wcheck check had (lane-rules.md,
"An examples-timing NOTE on this box measures LOAD"). No extension
changed, so no extension check.

## For the merge

- Version: minor (two BREAKING behavior changes). No new dependency;
  the floors stay.
- lane priors2 has not touched `R/weights-unbounded.R`,
  `check_re_structure()` or the check's call in `frm()` (read on
  2026-11-13 near the end; its `R/fit.R` hunks are elsewhere). The new
  4985 step calls `resolve_prior_input()`; if priors2 changes what that
  resolver refuses, rerun `test-weights-unbounded.R` after the merge.
- `NEWS.md` and `dev/test-backlog.md`: the development-version heading
  and the "Filed in the round of 2026-11-13" section are ones every lane
  adds; consolidation conflicts, not hunks against another lane.

## Defects found, not fixed (filed under "From lane wcheck2")

- HIGH (MEDIUM in the first round, rewritten in the fix round):
  asym_laplace's Laplace approximation has no maximum with a random
  effect in mu, and TMB's inner search does not find the modes either
  (section 4949 and "Fix round"). The fix round ships a refusal where
  the proof holds and the profile confirms it; the rest is open. Near
  3389.
- LOW: asym_laplace with every row censored has no "estimate does not
  exist" warning (section 4998).
- LOW: the pair rule for a shape is not decided (section 4977).
- LOW: `dev/wcheck-cases.R` labels its all-censored arm "72 right
  censored" (section 4998).
- LOW (fix round): the pair rule's misses on a held column and on
  unsaturated designs; huber `(1 | g4)` false convergence with a logLik
  0.064 to 0.136 off a fresh evaluation; identity-link sigma stopping at
  0 under a normal Intercept prior.

## Not done, and why

- The pair rule for a shape (Gamma, Beta, ... `shape ~ 1 + fac`): the
  loss rate of a shape going to 0 is per family and not derived; those
  stay "unknown" (filed).
- 4949's inner-search defect: not this check's; the check reads the
  Laplace approximation, which has a maximum there.
- The all-censored asym_laplace stop (4998's real content): the SE
  check's limit rule (selimit2's area), filed.

## Files touched

- `R/weights-unbounded.R`; `R/fit.R` (`check_re_structure()` and the
  new `re_check_term()` beside it, and in the fix round one call of
  `check_ald_laplace()` after `check_weights_unbounded()` in `frm()`);
  new `R/ald-laplace.R` (fix round); `R/compat.R` (fix round: the
  asym_laplace rows, beside the `group:autocor` one).
- `tests/testthat/test-weights-unbounded.R` (three blocks, then in the
  fix round the saturated-design cases, a cost block and the ZIAL
  block's `(1 | g)` case), `tests/testthat/test-diagnostics-ux.R` (one
  block), new `tests/testthat/test-ald-laplace.R` (fix round).
- `NEWS.md` (a development-version section; after the fix round five
  bullets, three BREAKING).
- `dev/test-backlog.md` (CLOSED lines under 4949, 4968, 4977, 4985,
  4993, 4998; new items under "Filed in the round of 2026-11-13 / From
  lane wcheck2").
- New scripts: `dev/wcheck2-common.R`, `dev/wcheck2-4968.R`,
  `dev/wcheck2-4977.R`, `dev/wcheck2-4985.R`, `dev/wcheck2-4998.R`,
  `dev/wcheck2-al-decomp.R`, `dev/wcheck2-al-group.R`,
  `dev/wcheck2-olre-label.R`, `dev/wcheck2-perturb.R`,
  `dev/wcheck2-table.R`, `dev/wcheck2-accept.sh`, `dev/wcheck2-tf.sh`,
  `dev/wcheck2-install.sh`, `dev/wcheck2-allcens.R`,
  `dev/wcheck2-check.sh`; in the fix round `dev/wcheck2-time-pair.R`,
  `dev/wcheck2-ald-profile.R`, `dev/wcheck2-ald-blocks.R`,
  `dev/wcheck2-compat-ald.R`.

No compat row, ledger row or bin-2 row changed in the first round; the
fix round changes compat rows (section "Fix round"), and no ledger or
bin-2 row.

## Fix round (after `dev/reviews/2026-11-13-wcheck2.md`, verdict FIX)

The review's scripts are `dev/reviews/wcheck2-review-*.R`, its logs
`dev/reviews/wcheck2-review-out/`. The coordinator added a user
decision (2026-11-13): refuse asym_laplace() and
zero_inflated_asym_laplace() fits whose Laplace-integrated random
effects enter mu where the Laplace marginal provably has no maximum,
and point to frm_sample(). The line "No compat row, ledger row or bin-2
row changed" above holds for the first round only; this round changes
compat rows (below).

### F1: the pair rule's cost

`unbounded_cell_pair()` now builds its units once and reads each
unit's gain (E/n) and loss (W/n) in one pass over its rows (membership
by mask, no `intersect()`), returns at once when the largest gain is
not above the smallest loss, and walks the candidate pairs in order of
gain and loss, stopping at the first that is valid for each shrinking
unit. The caller skips the rule when no uncensored weight is above 1
(every unweighted model). `excess_of()` reads rows by mask too, and the
one-cell loop reads the design once (`prior_column_verdict()` had
converted it to a dense matrix once per cell, which is rellib-r21's own
cost).

Times of `check_weights_unbounded()` alone, in seconds. The review's
scripts: `dev/reviews/wcheck2-review-time2.R` (unweighted, seed 3) and
`dev/reviews/wcheck2-review-prior.R <arm> time` (weights 0.95); and
`dev/wcheck2-time-pair.R` (weights uniform in 0.9 to 1.1 with
`normal(0, 1)` on the sigma Intercept, where the pair rule runs in full
and refuses nothing; minimum of 3 calls). Logs
`dev/wcheck2-out/time2-*.txt`, `prior-time-lane-fr2.txt`,
`time-pair-*.txt`; the review's rellib-r21 and first-lane numbers are
from `dev/reviews/wcheck2-review-out/`.

| case | rellib-r21 | first lane (review) | fixed lane |
|---|---|---|---|
| unweighted, n 2000, 100 levels | 0.03 to 0.04 | 3.00 | 0.01 to 0.03 |
| unweighted, n 4000, 200 levels | 0.12 to 0.16 | 13.75 | 0.26 to 0.28 (one call); 0.06 to 0.10 (see below) |
| w 0.95, n 2000, 100 levels, no prior / Intercept prior | 0.04 / 0.03 | 2.69 / 2.69 | 0.07 / 0.01 |
| w 0.95, n 6000, 300 levels | 0.80 / 0.61 | 48.07 / 47.80 | 0.69 / 0.45 |
| w 0.95, n 10000, 1000 levels | 10.24 / 8.87 | stopped | 7.92 / 8.48 |
| w (0.9, 1.1) + prior, n 2000, 100 levels | 0.20 | not run | 0.03 |
| same, n 4000, 200 levels | 0.97 | not run | 0.22 |
| same, n 6000, 300 levels | 4.44 | not run | 0.70 |
| same, n 10000, 1000 levels | 190.93 | not run | 11.09 |

The unweighted 4000-row row is a single call per process, and three
alternating runs gave rellib-r21 0.14, 0.12, 0.13 and the lane 0.26,
0.27, 0.28 (`time2-repeat.txt`). Timed after `gc()` on a fresh frame
the same call takes 0.11 and 0.11 on rellib-r21 and 0.06 and 0.10 on
the lane, and timed without it 0.27, 0.29 against 0.25, 0.22; five
repeated calls in one process take 0.09 to 0.11 against 0.06 to 0.08.
The single-call difference is a garbage collection that falls inside
the timed window, not the check's work. The last four rows show
rellib-r21's own cost (the dense design per cell) that the fix removes
too.

### F2: saturated designs, and the claims

When the parametric dispersion design has as many distinct rows as its
rank (`unbounded_row_patterns()`: rows coded column by column, then
`qr()` of the distinct rows), each distinct row is a unit; a pair's
direction is `solve(D, s)`, and only its moving columns must be free.
Otherwise the units are treatment cells and the reference level, as
before. The review's grid on the lane (`Rscript
dev/reviews/wcheck2-review-pair.R lane grid`, log
`dev/wcheck2-out/pair-grid-lane-fr.txt`):

    three        cols (Intercept),facb,facc | agree_ref 42, agree_pass 108, false_refusal 0, miss 0, lp_disagree 0
    four_unbal   cols (Intercept),facb,facc,facd | agree_ref 77, agree_pass 73, false_refusal 0, miss 0, lp_disagree 1
    inter        cols (Intercept),f1B,f2D,f1B:f2D | agree_ref 73, agree_pass 77, false_refusal 0, miss 0, lp_disagree 5
    fac_x        cols (Intercept),facb,facc,z | agree_ref 47, agree_pass 103, false_refusal 0, miss 0, lp_disagree 2
    sum3         cols (Intercept),fac1,fac2 | agree_ref 49, agree_pass 101, false_refusal 0, miss 0, lp_disagree 0
    three_bnorm  cols (Intercept),facb,facc | agree_ref 0, agree_pass 143, false_refusal 0, miss 7, lp_disagree 0
    three_w0     cols (Intercept),facb,facc | agree_ref 48, agree_pass 102, false_refusal 0, miss 0, lp_disagree 1

Against the review's lane log: f1 * f2 from 73 misses to 0 (73
refused), contr.sum from 49 misses to 0 (49 refused); 0 false refusals
in all 1050. The 7 misses with `normal(0, 1)` on `facb` stay: every
pair there moves facb, and the best direction moves a and b together
against c, a linear program over the held column that the rule does not
solve (filed). NEWS bullet 2, the roxygen of `unbounded_cell_pair()`
and the 4977 CLOSED line now say exactly this: pairs of distinct rows
of a saturated design (exact when the pair's moving columns are free),
else treatment cells and the reference level (other directions not
searched).

### F3 and the user's decision: asym_laplace with random effects

The filed MEDIUM is rewritten as HIGH with the review's statement and
table (4 designs x 10 seeds; exact-mode Laplace sup minus fit logLik,
median 6.63, 2.90, 49.24, 5.35; quadrature with no logLik in 8, 5, 4
and 7 of 10), and now says what ships and what remains.

What ships: `R/ald-laplace.R`, `check_ald_laplace()`, called by
`frm()` after the weights check (`R/fit.R`, one call). It refuses, with
class `frmtmb_ald_no_laplace`, an asym_laplace() or
zero_inflated_asym_laplace() response whose Laplace-integrated block
enters mu alone, when mu has the identity link and no nonlinear body,
the response has no `cens()`/`trunc()`, no random effect of the
response sits in sigma, quantile or zi, the model has no `mi()`,
`me()` or residual correlation term, every prior entry is on a class
that leaves the scale free (b, Intercept, sigma, quantile, zi, cor, L,
lscale), and the fit is not `quadrature = TRUE` or `importance`. The
proof: the log density is piecewise linear in mu, so H = Sigma^-1 and
the approximation is `max_b (l(b) - b' Sigma^-1 b / 2)`, which does not
decrease as the block's covariance is scaled up.

Each structure was profiled before it went into the refused list
(`Rscript dev/wcheck2-ald-profile.R lane <designs>`, seed 20261114, 80
rows, groups of 4; TMB's Laplace objective along every log sd of the
block, the rest re-optimized; logs
`dev/wcheck2-out/ald-profile-lane-<1..5>.txt`):

| term | frm() on the lane | -fn at log sd + 0, 1, 2, 4, 8 | never falls |
|---|---|---|---|
| (1 \| g) | refused | -101.446 -94.363 -93.670 -93.567 -93.565 | TRUE |
| (1 + x \| g) | refused | -88.602 -87.199 -86.716 -86.596 -86.595 | TRUE |
| (1 + x \|\| g) | refused | -75.012 -67.710 -66.976 -66.711 -66.372 | TRUE |
| s(z) | refused | -130.691 -130.316 -129.576 -129.184 -129.183 | TRUE |
| t2(z, x) | refused | -130.349 -128.173 -126.032 -120.012 -119.906 | TRUE |
| gp(z, k = 8) (hsgp) | refused | -130.707 -129.880 -129.387 -129.320 -129.319 | TRUE |
| ar1(tt + 0 \| g) | refused (by the weights check first: each row owns its effect) | 533.10 552.00 557.07 558.46 558.51 | TRUE |
| cs(1 + x \| g) | refused | -90.414 -87.592 -87.205 -87.133 -86.934 | TRUE |
| homdiag(1 + x \| g) | refused | -90.640 -87.527 -87.122 -86.944 -86.938 | TRUE |
| (1 \| gr(g, cov = A)) | refused | -102.554 -95.538 -93.645 -93.523 -93.521 | TRUE |
| zero_inflated_asym_laplace (1 \| g) | refused | -94.786 -93.864 -93.750 -93.732 -93.732 | TRUE |
| exact gp(z) | open (taken out of the list) | -129.315 -128.134 -127.123 -127.528 -261.369 | FALSE |
| rr(1 + x \| g, d = 1) | open (taken out) | the profile moved the loadings by t, not the scale: no test | |
| log link, (1 \| g) | open | -86.29, then non-finite | FALSE |
| sigma ~ (1 \| g) | open | -131.406 ... -133.749 | FALSE |
| cens(), (1 \| g) | open | -98.461 -97.061 -105.620 -118.456 -138.404 | FALSE |
| me(x, sdx) + (1 \| g) | open | -200.714 -199.710 -199.585 -199.565 -199.480 | TRUE |

The exact gp()'s objective falls at log sd + 4 and + 8, where TMB's
inner search fails on its near-singular covariance (correlation plus a
1e-12 nugget, scaled by sd^2, so the proof holds on paper); it is left
open because the profile does not confirm it. The me() model's profile
rises too (its latent does not couple with b), but the check leaves
every `mi()`/`me()` model open: a latent in a slope would couple with
b, which is not read. hetar1, ou, toep, homtoep, homcs, the spatial
structures, rr and gr(prec = ) were not profiled along their scale and
stay open. With a log link the fit errors ("NA/NaN gradient
evaluation"); with sigma ~ (1 | g) or cens() it fits at conv 1 with the
kink's warnings, as before.

frm_compat(): explicit rows `<structure> x asym_laplace`, status
"conditional", for us, diag, homdiag, cs, ar1, gr_cov, smooth, hsgp,
s() and t2(), and a `gp_pred() x asym_laplace` row that separates the
approximate gp (refused) from the exact one (not refused). Before,
each was "works" (the kind default) or its structure's own
"conditional" note. zero_inflated_asym_laplace is not in the compat
vocabulary.

Equivalence survey: pair p072 (asym_laplace x gr_cov) changes from
"equivalent" (L1) in `dev/equiv-data/results.tsv` to "defect, frmtmb
refuses or errors where brms builds the model" on the lane (`Rscript
dev/equiv-measure.R p072 dev/wcheck2-equiv-out`, log
`dev/wcheck2-out/equiv-p072.txt`), by the user's decision; the record
is not edited (rerun at consolidation). No ported brms-suite ledger row
fits asym_laplace with a random effect (its rows are the distribution
functions and log_lik of drawn values).

Also filed (LOW): huber `(1 | g4)` false convergence and a logLik that
differs from a fresh evaluation by 0.064 to 0.136 in 5 of 10; an
identity-link `sigma ~ 1 + fac` under a normal Intercept prior that
stops at sigma = 0 (logLik 136.452).

Note on the 58-row table: in the review's rerun 4 rows of part mv
differ from rellib-r21 in their warning text only (the 4993 wording),
by design; `dev/wcheck2-table.R` compares verdicts and values, which
are identical.

### Tests and runs of the fix round

New and changed tests:
- `test-weights-unbounded.R`: the pair block gains contr.sum and f1 *
  f2 (refused, and passed on the other side of the boundary) and the
  held-column miss (asserted as let through, with the filed item named);
  a new block counts `base::intersect()` calls during one check on 100
  levels (fewer than 100; the first lane made about 3 per ordered pair,
  9900 pairs); the ZIAL block's `(1 | g)` fit becomes the weights
  check's NULL plus the new refusal.
- New `test-ald-laplace.R`: refused on `(1 | g)`, `(1 + x | g)`,
  `(1 + x || g)`, `s(z)`, ZIAL `(1 | g)` and under a prior on b; left
  open (the check returns NULL) under exact, a prior on sd, the log
  link, sigma ~ (1 | g) with and without a mu effect, cens(), an exact
  gp(), a gaussian family and no random effect.

    rellib-r21:  RESULT test-ald-laplace.R pass=0 fail=1 err=2 skip=0 warn=2
    lane:        RESULT test-ald-laplace.R pass=21 fail=0 err=0 skip=0 warn=0
    lane:        RESULT test-weights-unbounded.R pass=210 fail=0 err=0 skip=0 warn=0
    lane:        RESULT test-diagnostics-ux.R pass=134 fail=0 err=0 skip=0 warn=0

(`dev/wcheck2-out/tf-ald-<base|lane>.txt`, `tf-wu-lane-fr.txt`,
`tf-dux-lane-fr.txt`.) On rellib-r21 the refusal block fails as
behavior: `frm()` returns a `frmtmb_fit` (with the false-convergence and
gradient warnings) where the refusal is expected; the second block
errors on the missing function (the weak form). The contr.sum and f1 *
f2 refusals were seen as misses on the first lane build by the review
(`pair-grid-lane.txt`: 49 and 73 misses).

OpenBLAS 0.3.32 with LAPACK (`bash dev/ciharden-openblas.sh 0.3.32
lapack`: 617 LAPACK routines to OpenBLAS, LAPACK 3.12.0; OpenBLAS picks
its Haswell kernels on this CPU, as dev/reviews/2026-10-07-ciharden.md
recorded; log `dev/wcheck2-out/openblas032-build.txt`):
`test-weights-unbounded.R` pass=210, `test-diagnostics-ux.R` pass=134,
`test-ald-laplace.R` pass=21, each fail=0 err=0 skip=0 warn=0
(`tf-ob032-*.txt`). `dev/ciharden-out/` deleted.

Every suite, one file per process, 3 at a time (`bash
dev/release/run-par.sh wcheck2-suite-par2 dev/wcheck2-out/jobs-all.txt
3`, log `dev/release/wcheck2-suite-par2.log`, "WCHECK2-SUITE-PAR2 ran
530 of 530"; all 530 `lib:` lines name `wt-wcheck2-lib`):

    frmtmb files=331 pass=27517 fail=0 err=0 skip=15 warn=0
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=33 pass=1857 fail=0 err=0 skip=3 warn=0
    frmtmb.latent files=13 pass=697 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=16 pass=545 fail=0 err=0 skip=2 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=91 pass=4603 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial files=8 pass=612 fail=0 err=0 skip=0 warn=0
    frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0

Against the first round's run (`dev/wcheck2-out/suite-par-round1-results.txt`)
exactly two RESULT lines change: the new `test-ald-laplace.R` (21) and
core's `test-weights-unbounded.R` (201 to 210). The files that fit
asym_laplace or zero_inflated_asym_laplace (`grep -l asym_laplace` over
the tests: core's test-brms-lik-cdfs.R, test-brms-suite-brmsterms.R,
test-cdfs.R, test-families.R, test-latres.R, test-mi-families.R,
test-numerical-robustness.R, test-optstart.R, test-se-exact-hessian.R,
test-simulate-density.R, test-v14.R, test-compat*.R; frmtmb.sample's
test-brms-suite-brmsterms.R, test-default-priors-dpar.R,
test-mi-laplace-route.R, test-sds-prior.R; frmtmb.eam's test-family.R,
which reads the compat table) keep their counts: none fits an
asym_laplace random effect in mu under Laplace.

Acceptance on the final build (`bash dev/wcheck2-accept.sh base lane`,
`dev/wcheck2-out/accept-table2.txt`):

    case table: 58 of 58 rows identical on rellib-r21 and the lane
    common: 26 lines, 1 differ (apart from lib:)
      base: CASE asym_laplace (1|g) w                   ok logLik -286.2962004 | sum|est| 2.330028119 | conv 1 | warn 2
      lane: CASE asym_laplace (1|g) w                   ERROR frm() cannot fit `(1 | g)` under the Laplace approximation: ...
    routes: 27 lines, 0 differ (apart from lib:)

One of the review's 24 guard-absent models changes: weighted
asym_laplace with `(1 | g)`, which rellib-r21 fitted at conv 1 with two
warnings, is now refused by the user's decision; the weights check
itself still lets it through. Part mv of the case table differs only in
the warning text of 4 rows (the 4993 wording).

R CMD check --as-cran, core, on the final tree (`bash
dev/wcheck2-check.sh`, log kept as `dev/wcheck2-out/check-frmtmb-fr.txt`,
check tree deleted): examples OK, tests "[46m] OK", `Status: 1 NOTE`
(the expected V8 NOTE). No extension changed, so no extension check.

The machine-wide limit of 2 R processes (coordinator's notice) came
during this check, which was the only R run of the lane at the time; no
run failed for disk space or memory.

## Fix round 2 (S1 of `dev/reviews/2026-11-13-wcheck2-short.md`)

`ald_laplace_blocks()` refused a gaussian block beside a Student-t
block (`gr(h, dist = "student")`) or another package's block in the same
mu, where the t block's mode-dependent curvature breaks the proof (the
short check's TMB profile fell 10.4 at nu 2.5, the exact-mode Laplace
7.46). It now refuses nothing in a mu that holds such a block. The
roxygen, the compat note, NEWS and the HIGH backlog entry say so.

New "stays open" assertion in `test-ald-laplace.R` (`y ~ x + (1 | g) +
(1 | gr(h, dist = "student"))`), seen failing on the build before the
fix:

    before:  RESULT test-ald-laplace.R pass=19 fail=0 err=1 skip=0 warn=0
    after:   RESULT test-ald-laplace.R pass=22 fail=0 err=0 skip=0 warn=0
             RESULT test-weights-unbounded.R pass=210 fail=0 err=0 skip=0 warn=0
             RESULT test-compat.R pass=868 fail=0 err=0 skip=0 warn=0

(`dev/wcheck2-out/tf-ald-lane-s1before.txt`, `tf-ald-lane-s1.txt`,
`tf-wu-lane-s1.txt`, `tf-compat-lane-s1.txt`; at most 2 R processes.)
The "before" error is the refusal itself, class
`frmtmb_ald_no_laplace`, where NULL is expected.
