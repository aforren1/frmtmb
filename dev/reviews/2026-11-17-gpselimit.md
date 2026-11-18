# Review of lane gpselimit (round of 2026-11-16), 2026-11-17

Item 5571: an approximate `gp(k = )` at the white limit of its basis.
Reviewer scripts: `dev/reviews/gpselimit-review-*.R`; logs in
`dev/reviews/gpselimit-review-out/` (local). Arms: `base` = rellib-r22
(0.83.0), `lane` = wt-gpselimit-lib over rellib-r22. Every `lib:` line of
every log names the arm's frmtmb. Base and lane ran in separate
processes.

## Verdict: FIX

1. The certificate passes a stop inside the white plateau below a higher
   interior maximum, and the fit gets the boundary message where base
   warned. Add the way-back rule (finding 2).
2. The "misdiagnosed" claim is wrong for `gp(by = )`: base keeps finite
   SEs with no lost reason on white sub-GPs, and the lane's
   `gp_white_dirs()` fixes that. Correct the record (findings, NEWS,
   backlog CLOSED line), and add a test seen failing on rellib-r22
   (finding 3).
3. Small: do not certify a fit that has a convergence warning (finding
   4); the CLOSED line is four lines, the rule asks for one (finding 9).

Everything else checks out: the 30-fit table, bitwise equality on single
term fits, the hypothesis() string and its SE, `check_se` semantics, the
choice of "boundary" over "limit", the reproduced exact-gp MEDIUM, and the
80-column check.

## Findings

### 1. Misdiagnosis on the review's 30 fits: confirmed

`gpselimit-review-sweep.R` (same data as
`dev/reviews/gpfield-review-basis.R`, seeds 9600 + row), read through
`sdr_of()`, `vcov()` and `summary()`, not warning text
(`sweep-base.txt`, `sweep-lane.txt`):

- base: 13 of 30 white (rows 1-5, 9, 10, 14, 21, 22, 25-27), each with
  `se(theta)` NaN, NaN and `se_lost` `theta_1=boundary,
  theta_2=boundary` (kind "limit"); the 17 resolved rows lose nothing.
- The review's own metric (`grepl("Standard errors are not available",
  warnings)`) gives 0 on all 30 rows on base: that is where "se_lost 0
  on all 13" came from.
- lane: the same NaN and "boundary" (kind "gpwhite") on the 13; the
  condition is `m:frmtmb_gp_basis+frmtmb_boundary_fit` in place of the
  warning.

For single-term fields the lane is right. See finding 3 for `gp(by = )`.

### 2. The certificate does not look outward (FIX 1)

`gp_white_ridge()` probes only toward lscale 0 (t = 2, 12). At a point
deep in the white plateau the gradient is about 0 in every direction, so
the optimizer stops there whatever the likelihood does at larger length
scales, and the probe toward 0 is level by construction.
`gpselimit-review-cases.R` group A: resolved sweep rows started at
lscale = 3e-5 of the spacing (`start = list(theta = c(log(30),
log(3e-5 * sp)))`, `cases-lane.txt`, `cases-base.txt`):

| row | default fit ll | plateau stop ll | short by | lane | base |
|---|---|---|---|---|---|
| 6 | -35.399991 | -40.581497 | 5.18 | cert TRUE, message | warning |
| 12 | -20.543822 | -31.070729 | 10.53 | cert TRUE, message | warning |
| 18 | -37.172714 | -51.344372 | 14.17 | cert TRUE, message | warning |
| 24 | -73.639229 | -87.886060 | 14.25 | cert TRUE, message | warning |
| 30 | -16.543695 | -20.108328 | 3.56 | cert TRUE, message | warning |

All five have convergence code 0 and no other condition. The message says
"Boundary fit ... The fitted field and its predictions hold", which is
false 3.6 to 14.3 log-likelihood units below the maximum. This is the
case the lane rules name: a warning turned into a message on a short fit.

The objective along the ridge separates the two cases. Outward (t < 0,
toward larger length scales), f(t) - f0:

- the 13 true white fits (`sweep-lane.txt`): t = -4 from +4.71e-06 (row
  27) to +5.16e-05 (row 21), t = -2 from +8.47e-08 to +9.29e-07: a cost
  on the way back, every row;
- the 5 plateau stops: t = -4 from -1.52e-05 (row 30) to -1.09e-04 (row
  24), t = -6 from -8.31e-04 to -5.96e-03: a gain.

0.82.0's limit rule for coefficients already has this test
(`se_limit_dir()`: "one that gains on the way back is not at a maximum
that way either"). Fix: in `gp_white_ridge()`, also move outward along
each ridge (t = -2 * se_edge_step, or to the gate's edge, l = 0.1 of the
spacing, which is stronger for stops deeper than 3e-5) and refuse on a
gain above `se_edge_tol`. Test: row 6 started at 3e-5 of the spacing must
keep the warning (it gives the message on the lane build now).

How often without a user start: group B refitted each of the 13 white
rows from lscale = 0.3, 1 and 3 spacings. 36 of 39 came back to the
white limit within 8e-9 of its log-likelihood, 3 ended lower with
convergence warnings; none found a higher maximum. So the default start
did not produce a plateau stop here; `start` (and anything that refits
from a moved start) does.

Other cases of check 2, all on both builds (`cases-*.txt`):

- Priors (sweep row 3): lognormal on lscale at 0.01 of the spacing with
  sd 1 and 3, `student_t(3, 0, 2.5)` and `exponential(0.001)` on sdgp:
  certificate FALSE (ridge f(t) - f0 at t = 2 from 0.22 to 2.0), warning
  kept, SEs finite, identical to base. A proper prior cannot pass: the
  ridge moves log sdgp by D/2 per unit.
- Stop short toward 0: the lane's iter.max case, and group A's start at
  1e-3 of the spacing (ridge at t = 2: +1.4e-6 to +9.8e-6): certificate
  FALSE, warning, reason "flat" (base: "boundary").
- 2-D iso and aniso, exp_quad, matern32, matern52 (seed 5161 field):
  certified, message, all hyperparameters "gpwhite", SEs identical to
  base.
- Aniso with one white dimension (x1 at 5.03, x2 at 4.2e-4 of the
  spacing): not gated, `theta_1, theta_2` "flat" with the generic SE
  warning on both builds. Out of scope, as the lane says.
- `sigma ~ gp(x, k)`: certified where white (3 fits), names
  `sdgp_sigma_gpx`, identical to base.
- Two terms `gp(x1, k = 10) + gp(x2, k = 30 or 10)`: only the white one
  named; identical to base.
- REML (rows 3, 9, 6) and poisson (rows 3, 9, 21): certified where white,
  identical to base.
- Toward 0 the certificate is sound on every fit seen: the true white fits
  sit about 1e-8 short of the limit (t = 1 to 20: -1.4e-9 to -1.7e-8),
  two to three orders under the 1e-6 gate, and in the gate's range the
  objective is close to linear in lscale^2.

### 3. `gp(by = )`: base keeps finite SEs, the lane changes them (FIX 2)

`gpselimit-review-byjit.R`: `y ~ gp(x, by = f, k)`, f alternating a/b on
sweep rows 3 and 21, response moved by `y * (1 + 2e-15 * j)`, j = 0..20
(`byjit-base.txt`, `byjit-lane.txt`):

| row | build | white hyperparameters with a finite SE | lost |
|---|---|---|---|
| 21 | base | 4 of 4 on 21 of 21 draws | none |
| 3 | base | 2 of 4 on 21 of 21 draws | theta_3, theta_4 boundary |
| 21, 3 | lane | 0 of 4 on 42 of 42 draws | all four boundary |

Base row 21 j = 0: theta SEs 4916.05, 9832.09, 2030.81, 4061.62 (log
scale), so `summary()` shows sdgp 197.7 with an SE, and nothing in the SE
check says why. That is the item as filed. The lane fixes it through
`gp_white_dirs()` in `cov_from_hessian()`, which its findings call a
safety net that changed nothing. Here it is the fix, and it moves other
SEs: intercept 0.75213746542 to 0.75213709108, sigma 0.054298124 to
0.054298113 (largest relative change 5.0e-7, `cmp.txt` rows E-by-f-*).
That move is right: the ridge's near-zero curvature leaked into them.
Also seen on the `frmtmb.sample` test's data (`E-by-f-k5`: base theta_3,
theta_4 SE 3515, 7030).

Needed: findings section 1 and the backlog CLOSED line should say the
item was right for `gp(by = )` and is fixed; NEWS "no standard error, as
before" is false there (add it, BREAKING: SEs that were finite become
NaN, and other SEs move in the 7th digit); the `gp_white_dirs()` doc
comment should say it is load-bearing; a test on `gp(by = f)` seen
failing on rellib-r22, under ulp draws and OpenBLAS.

### 4. A fit with convergence warnings is certified (FIX 3, small)

`gpselimit-review-nl2.R` (`y ~ a * (1 - exp(-b * z))`, `a ~ gp(x, k =
20)`, b0 = -0.8), seeds 2 and 4: code 1, `frmtmb_conv_code` (false
convergence) and on seed 2 `frmtmb_conv_gradient` (0.00588 at
`b_(Intercept)`), and also the boundary message, which says "The fitted
field and its predictions hold". The ridge probe holds b at a stop that
is not a maximum. The warnings stay, so the user is told; but the message
asserts more than is known. Gate the message on no `frmtmb_conv_*`
verdict, else keep the warning.

### 5. Other SEs: identical where the Hessian removed the ridge

`gpselimit-review-cmp.R` (`cmp.txt`): estimates, point, log-likelihood,
every SE and the lost reasons are `identical()` on all 30 sweep rows.
Estimates, point and log-likelihood are identical on all 42 case fits,
the SEs on 39: the three that differ are the `gp(by = )` fits of finding
3. The lost reasons differ, as intended, on the five group-A starts at
1e-3 (boundary to flat).

### 6. The hypothesis() string: correct

On every certified fit the string evaluates, and its SE equals the SE of
psi = 2 log sdgp + sum of log length scales from a numerical Hessian
(`optimHess` on `obj$fn`, `obj$gr`) in (psi, ridge) coordinates: 13
sweep rows (largest relative difference 1.4e-6) and 27 case strings (1-D,
2-D iso `2 * log(lscale)` and aniso `log(l1) + log(l2)`, matern32 and 52,
`by = f` per level, `sigma`, REML, poisson; largest 2.5e-6). The ridge
directions are orthogonal to the gradient of psi ((2, D) . (D/2, -1) = 0;
(2, 1, 1) . (1/2, -1, 0) = 0), so psi is estimable and its SE does not
depend on which direction is removed. Sign: the string is the quantity
itself, positive orientation; the estimate matches `gp_brms_values()`.

One string fails: two terms on one covariate (`gp(x1, k = 10) + gp(x1, k
= 5)`) give "Internal error: two parameters of this model share the name
'sdgp_gpx1'" (`G-two-white`). That name clash predates the lane; the
message only repeats the name. Worth a LOW backlog item, not this lane's.

### 7. `check_se` semantics and class consumers: consistent

`gpselimit-review-misc.R` (`misc-*.txt`), a singular `(1 | g)` fit beside
white row 3:

| check_se | singular (both builds) | white, base | white, lane |
|---|---|---|---|
| warning (default) | boundary message | gp warning | boundary message |
| ignore | silent | gp warning | silent |
| stop | error at frm() | gp warning | error at frm() |

The lane's message now behaves as the SE check's own boundary message.
Consumers: nothing reads `frmtmb_boundary_fit` by class in R/ or the
extensions' R/ (only tests' `allow_boundary()`); `nl_flat_warning()` and
`nl_fit_off_saddle()` read `frmtmb_limit_fit`, which the message does not
carry; `diagnose()`, `confint()` and predict()'s draw law read `se_lost`
and V, which are unchanged on single-term fits.

Code-read risk, not reproduced: on a nonlinear body, `nl_fit_off_saddle()`
holds warnings and `frmtmb_limit_fit` messages until its fit is returned,
but not this message; `nl_fit_probed()` says a replaced first fit's
messages. A white first fit that a restart, sign probe or pole escape
replaces could so print a stale boundary message. 52 fits
(`gpselimit-review-nl.R`, `-nl2.R`; 16 with `frmtmb_nl_restart`) gave at
most one `frmtmb_gp_basis` condition each, always about the returned fit.
Holding `frmtmb_gp_basis` beside `frmtmb_limit_fit` would close it.

### 8. "boundary" against "limit": the lane is right, for a concrete reason

In `se_report()` a "limit" parameter is quiet only when
`se_limit_quiet()` holds, which needs a `se_limit_info()` kind for every
name and `se_limit_family_par()` for all of them. gp hyperparameters have
no `se_limit_info()` entry (it covers coefficients run off), so
`anyNA(kind)` makes the fit not quiet and `re_check_act()` gives the
"Standard errors are not available" warning: "limit" would turn the
message back into a warning. `se_limit_draw_held()` would also try to
hold their draws by the coefficient rule. "boundary" with its own kind for
the text is the right mapping.

### 9. Exact gp() MEDIUM: reproduces

Design of `dev/gpselimit-exact2.R`, rep 3, seeds 1 to 4
(`gpselimit-review-misc.R`, `-exact.R`): seed 2 ends at 0.182 of the
smallest gap, the objective -6.09e-9 at log lscale -2 and -12; in a fresh
process `summary()` shows `lscale(gpx) 0.00 0.12 0.00 Inf` with no
condition and no lost reason (both builds); after probing the objective
first, `sdr_of()` loses `theta_2` as "concave", no condition
(`exact-lane.txt`). Seed 4 (0.040 of the gap) gets "boundary". Seeds 1
and 3 are interior (2.1 and 7.3 gaps, objective rises 0.9 and 1.3 at -2).
The filed text is accurate.

Format: the CLOSED line under item 5571 runs four lines; the rule asks for
one line.

### 10. Test edits: nothing asserted was lost

Each site now wrapped in `allow_boundary()` used `allow_warnings(...,
"spacing its basis resolves")` with no `require`, so none asserted the
basis warning before. The one that did (test-gp-field.R) now counts one
`frmtmb_gp_basis` condition of either type, and test-gp-white-limit.R
asserts the message on row 3. `allow_boundary()` muffles any
`frmtmb_boundary_fit`, including an SE-check singular message, which is
wider than needed, but no wrapped site asserted one. Gaps: no test of
NEWS bullet 2 (a stop short now "flat"), none for `gp(by = )` (finding
3), none for the plateau (finding 2). The seen-to-fail log is real:
`new-base.log` `pass=13 fail=7 err=2 skip=0 warn=4` on rellib-r22.

### 11. The 80-column check

The lane's substitute (Grep for lines over 80 characters) is the same
check as the denied awk scan if it covered every changed file. Rerun with
awk on every added line of the 13 changed files (base from `git show
b255350e:<path>`, compared with `diff --strip-trailing-cr`) and on the new
test file: 0 lines over 80.

## Reviewer notes

- No R process of the review exceeded two at a time; no package code
  edited; no git operation beyond `git status` and `git show`.
