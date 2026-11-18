# Lane alapse: asym_laplace without curvature

Round of 2026-11-16, base b255350e (frmtmb 0.83.0), worktree
`frmtmb-wt-alapse`, library `wt-alapse-lib`. Brief:
`round-briefs/r084-alapse.md` (items 4982, 5449, 5490). User decision
2026-11-14: on a fit without random effects, refuse the Hessian
standard error by name; no sparsity or sandwich covariance here.

Logs are under `dev/alapse-out/` (gitignored); the numbers below are
pasted from them. Earlier lanes' logs this lane read were copied to
`dev/alapse-prev-out/` before their trees go.

## What changed

New file `R/ald-fixed.R`. Hooks, each a few lines, in shared files:

- `R/fit.R`: `fit_assembled()` calls `ald_polish()` after the escapes
  and mixture starts (one hunk); the `se_explained <- "convergence"`
  gate skips a fit with refused parameters (one condition);
  `check_convergence()` reads `fit$opt$ald_cert$ok` (`ald_ok`) to skip
  the boundary-stop test and the gradient warning, uses
  `ald_conv_clause()` in the code-8 message, calls `ald_cens_sep()`
  when `separation_check()` finds nothing, and returns
  `at_max = at_max || ald_ok`.
- `R/se-check.R`: `sdr_rescue()` ends with `ald_sdr_refuse()`;
  `se_fit_analysis0()` returns `ald_se_analysis()` first;
  `se_report()` hands reason "kink" to `ald_se_report()`;
  `se_lost_clauses()` has the "kink" text; `jc_nan_lost()` reads "kink"
  as "nocurv"; the "nocurv" prediction warning names asym_laplace's
  case; `sdr_sep_lost()` relabels a named parameter's "nocurv", "flat"
  or "limit" as "separation".
- `R/confint.R` (`diagnose()`): two printed lines for a proved fit and
  for refused parameters. `R/interop.R` (emmeans basis): the "nocurv"
  words when every lost parameter is refused by name.
- `R/ald-laplace.R`: `ald_scaled_covstructs` gains hetar1, ou, toep,
  homtoep, homcs, exp, gau, mat, gr_prec (one per line); the gr_prec
  label.
- `R/compat.R`: one new row, `emmeans x asym_laplace`; the loop over
  `ald_scaled_covstructs` adds nine rows (below).
- `R/families.R`: the "Quantile regression inference" section of the
  family help, rewritten (no `lcdf` line touched).
- `NEWS.md`: five bullets under a new `# frmtmb (development version)`.
- Tests: new `tests/testthat/test-ald-fixed.R`; new block in
  `test-ald-laplace.R`; allowed-warning lists in `test-mi-families.R`
  and `test-cdfs.R`; `test-se-limit.R`'s 5423 test reads "kink".

### 4982: the refusal by name, and where it sits

`ald_se_lost(fit)` names the refused parameters: every coefficient of an
asym_laplace() or zero_inflated_asym_laplace() response except `zi`'s,
on a fit with no random effect, no latent value and no integration
(`ald_fixed_resps()`), and only when the response has at least one
exact row (a response whose every row is censored reads only tail
probabilities, which are curved; 5490 below). The refusal sits in the
standard-error analysis itself, so every reader agrees:

- `se_fit_analysis0()` returns `ald_se_analysis()` (the fit-time check);
- `sdr_rescue()` applies `ald_sdr_refuse()` last, so `sdr_of()`, which
  `summary()`, `vcov()`, `confint()`, `fixef()`, predictions, the
  emmeans basis and marginaleffects (through `vcov()`) read, gets NaN
  and reason "kink" whatever the stop (a short stop's
  `sdr_rescue_hessian()` path included);
- `se_report()` gives the fit-time warning, class `frmtmb_ald_no_se`,
  through `check_se` ("ignore" silent, "stop" an error of that class);
- `diagnose()` prints the reason and, for a proved fit, that the
  gradient is a one-sided derivative at the kinks.

Which fits the rule covers, and why: asym_laplace() with any link, any
`sigma ~` and `quantile ~` model, a fixed quantile, weights, censored
rows beside exact ones; zero_inflated_asym_laplace() the same, with
`zi`'s coefficients kept: the zero rows read only `zi` and the other
rows only the asymmetric Laplace, so the exact Hessian's cross block is
0 and `zi`'s block keeps the analysis of its own (a logistic regression
on the zeros; test: ratio to glm()'s standard errors within 1e-4).
The covariance needs the inverse of the whole Hessian, whose location
block is 0 with the identity link (`dev/alapse-out/probe1-base.txt`:
the 2 x 2 location block of the AD Hessian is exactly 0), so sigma's
and the quantile's standard errors go too; on 0.83.0 a `sigma ~ z,
quantile ~ z` fit showed finite standard errors for those four
coefficients beside NaN for mu's (same probe). Another link is covered
because its curvature, `rho'(u) f''(eta)`, is not the sparsity either.
Not covered: fits with random effects (the outer Hessian is a Laplace
or quadrature marginal's; 5449).

The message names the routes from one vector, `ald_se_routes`, so the
smoothed quantile family (lane squant) can be put first there without
a rewrite. frmtmb's bootstrap reaches the family: `frm_bootstrap()` is
a parametric bootstrap that draws from the fitted asymmetric Laplace;
on the readers probe, 20 refits gave a bootstrap SE for x of 0.073
(`dev/alapse-out/readers-lane.txt`). `confint(method = "profile")`
still runs TMB's profile (0.464 to 0.773 for x on the same fit), and
TMB's `confint.tmbprofile()` lets a "collapsing to unique 'x' values"
warning escape, on 0.83.0 as on the lane (filed).

### 4982: the convergence warnings, and the proof

The optimizer stops on the kinks. `ald_polish()` (fit end, before the
sdreport and the convergence check) solves the location exactly: with
sigma and the quantile held, the location coefficients minimize the
weighted check function `sum(w_i / sigma_i * rho_{p_i}(y_i - x_i b))`,
a linear program. `ald_vertex()` starts at the basis of the k rows
nearest the stop and takes simplex steps of the L1 problem (the edge
along which a basic row's residual leaves zero, a weighted-median line
search); a vertex is optimal exactly when the subgradient condition
holds (Koenker 2005, Theorem 2.1, with weights and a quantile per row).
Then one Newton step in the other parameters (sigma, quantile, zi,
other responses), the location solved again, until the Newton gain is
below 1e-10. The proof at the end: the vertex condition (margin above
-1e-8), a positive definite Hessian in the other parameters, a Newton
gain within `min(conv_at_max_tol, grad_tol)`, and an objective not
above the stop's. The kinks lie on hyperplanes in the location alone,
so the derivative in the other parameters is continuous across them
and this is the first-order condition of the whole problem; with the
vertex condition strict the location is a sharp minimum, so the point
is a strict local maximum. Proved fits get code 0, a message saying
so, and `fit$opt$ald_cert` (the optimizer's code and message, the
gain, the margin, the gap closed, the move).

CORRECTED in fix round 1 (review finding 4, and finding 1): "strict"
was wrong, since `ald_vertex()` accepts a margin down to -1e-8 and 0
(a non-unique solution, the review's even-median case); the point is a
local maximum, strict only with a positive margin. And code 0 for a
local proof was wrong: see "Fix round 1" below; code 0 now needs a
global proof.

Not proved (the warnings stay, with `ald_conv_clause()` in place of
"the standard errors are those of the stop"): `cens()`/`trunc()` rows,
another link, a nonlinear body, `mo()`, `me()`, a residual correlation
term, a prior density or a bound on the location's coefficients, a
bound that holds a smooth parameter, a curvature that is not positive
definite.

The rq comparison (`Rscript dev/alapse-rq.R <base|lane> 1:25`, logs
`rq-base.txt`, `rq-lane.txt`; summaries by `dev/alapse-rqsum.R`,
`rqsum-base.txt`, `rqsum-lane.txt`). Designs (seed 1000 * design +
s, s = 1..25): A 160 rows, ALD noise p 0.3, quantile estimated; B 100,
gaussian; C 120, `y ~ x1 + x2`, t3 noise, quantile = 0.25 fixed; D 200,
`sigma ~ z`; E 200, `quantile ~ z`; F 200, zero_inflated, 20% zeros; G
300, `y ~ x + fac`, ALD p 0.7; plus the four 4982 cases of
`dev/selimit2-asym.R` (dbg14 seed 43, gauss5 seed 5, fixq, zi seed 9).
`gap_b` is the ALD log-likelihood at quantreg's coefficients (rq.wfit,
weights 1 / sigma_i, the fit's own quantile; a row-varying quantile as
a median regression with one pseudo-row) minus the fit's, at the fit's
own sigma and quantile; `gap_j` the joint maximum (the quantile process
for a constant sigma and quantile, else alternating rq and optim from
the fit) minus the fit. 179 fits per arm.

rellib-r22 (`rqsum-base.txt`):

    dbg14  n  1 | code!=0  1 | conv warn  1 | gap_b>2e-5  0 (max 4.21e-10) | gap_j>2e-5  1 (max 3.23e-01, median 3.23e-01) | short warned 0 of 0, at-rq silent 0 of 1
    gauss5 n  1 | code!=0  1 | conv warn  1 | gap_b>2e-5  0 (max 7.24e-09) | gap_j>2e-5  1 (max 4.08e-02, median 4.08e-02) | short warned 0 of 0, at-rq silent 0 of 1
    fixq   n  1 | code!=0  1 | conv warn  1 | gap_b>2e-5  1 (max 2.05e-05) | gap_j>2e-5  1 (max 2.27e-05, median 2.27e-05) | short warned 1 of 1, at-rq silent 0 of 0
    zi     n  1 | code!=0  0 | conv warn  1 | gap_b>2e-5  0 (max 3.41e-09) | gap_j>2e-5  0 (max 5.69e-09, median 5.69e-09) | short warned 0 of 0, at-rq silent 0 of 1
    A      n 25 | code!=0 24 | conv warn 25 | gap_b>2e-5  5 (max 9.63e-04) | gap_j>2e-5 15 (max 1.09e-01, median 8.52e-05) | short warned 5 of 5, at-rq silent 0 of 20
    B      n 25 | code!=0 20 | conv warn 25 | gap_b>2e-5  3 (max 1.66e-04) | gap_j>2e-5 19 (max 7.74e-01, median 2.33e-03) | short warned 3 of 3, at-rq silent 0 of 22
    C      n 25 | code!=0 20 | conv warn 25 | gap_b>2e-5  7 (max 4.21e-02) | gap_j>2e-5  8 (max 4.21e-02, median 2.40e-06) | short warned 7 of 7, at-rq silent 0 of 18
    D      n 25 | code!=0 24 | conv warn 25 | gap_b>2e-5  3 (max 1.09e-04) | gap_j>2e-5 14 (max 2.52e-03, median 2.32e-05) | short warned 3 of 3, at-rq silent 0 of 22
    E      n 25 | code!=0 22 | conv warn 25 | gap_b>2e-5  7 (max 1.16e-03) | gap_j>2e-5 15 (max 5.59e-03, median 8.70e-05) | short warned 7 of 7, at-rq silent 0 of 18
    F      n 25 | code!=0 20 | conv warn 25 | gap_b>2e-5  3 (max 4.66e-03) | gap_j>2e-5 12 (max 1.63e-02, median 1.78e-05) | short warned 3 of 3, at-rq silent 0 of 22
    G      n 25 | code!=0 25 | conv warn 25 | gap_b>2e-5 15 (max 7.30e-03) | gap_j>2e-5 17 (max 1.10e-01, median 7.93e-05) | short warned 15 of 15, at-rq silent 0 of 10

Lane (`rqsum-lane.txt`):

    dbg14  n  1 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 0.00e+00) | gap_j>2e-5  1 (max 3.23e-01, median 3.23e-01) | short warned 0 of 0, at-rq silent 1 of 1
    gauss5 n  1 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 0.00e+00) | gap_j>2e-5  1 (max 4.08e-02, median 4.08e-02) | short warned 0 of 0, at-rq silent 1 of 1
    fixq   n  1 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 0.00e+00) | gap_j>2e-5  0 (max 0.00e+00, median 0.00e+00) | short warned 0 of 0, at-rq silent 1 of 1
    zi     n  1 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 0.00e+00) | gap_j>2e-5  0 (max 5.68e-14, median 5.68e-14) | short warned 0 of 0, at-rq silent 1 of 1
    A      n 25 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 0.00e+00) | gap_j>2e-5  8 (max 1.09e-01, median 3.98e-13) | short warned 0 of 0, at-rq silent 25 of 25
    B      n 25 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 2.84e-14) | gap_j>2e-5 13 (max 7.74e-01, median 3.75e-04) | short warned 0 of 0, at-rq silent 25 of 25
    C      n 25 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 0.00e+00) | gap_j>2e-5  0 (max 4.54e-11, median 0.00e+00) | short warned 0 of 0, at-rq silent 25 of 25
    D      n 25 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 0.00e+00) | gap_j>2e-5  0 (max 9.09e-12, median 0.00e+00) | short warned 0 of 0, at-rq silent 25 of 25
    E      n 25 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 0.00e+00) | gap_j>2e-5  0 (max 1.92e-11, median 0.00e+00) | short warned 0 of 0, at-rq silent 25 of 25
    F      n 25 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 0.00e+00) | gap_j>2e-5  7 (max 1.63e-02, median 3.98e-13) | short warned 0 of 0, at-rq silent 25 of 25
    G      n 25 | code!=0  0 | conv warn  0 | gap_b>2e-5  0 (max 1.14e-13) | gap_j>2e-5  3 (max 1.10e-01, median 0.00e+00) | short warned 0 of 0, at-rq silent 25 of 25

On 0.83.0 every one of the 179 fits warned, 43 of them more than 2e-5
below rq at their own quantile (up to 0.042, design C), and 136 at
rq's objective warned too. On the lane all 179 are proved (179 lines
"code 0 the exact quantile regress..."), the largest `gap_b` is
1.14e-13 (G) and the largest coefficient difference 1e-15 order, and
none warns about convergence; each warns once, by name (`ald_no_se`).

WITHDRAWN in fix round 1 (review finding 1): the next sentence is
wrong. The D and E rows' `gap_j` were computed by `alt_max()` started
at the fit, which finds the fit again; the review's multi-start
(9 quantile starts) found 9 of 20 D and 9 of 20 E fits below a higher
maximum, by up to 0.313 and 1.885 (`dev/reviews/alapse-review-out/
multistart-lane.txt`). The table's D and E `gap_j` columns measure
nothing.

The joint gaps that remain are all in designs with a constant estimated
quantile (A, B, F, G, dbg14, gauss5): the stop is a strict local
maximum, and a higher one sits in another interval of the quantile
process. Largest (seed: gap_j, |dp| = distance of the quantile from the
joint maximum's):

    B 20 5.27e-02 5.51e-02 | A 10 9.50e-02 1.68e-02 | B 17 1.04e-01 4.02e-02
    A 12 1.09e-01 2.23e-02 | G 13 1.10e-01 2.02e-02 | B 2 1.52e-01 8.94e-02
    B 3 2.47e-01 1.32e-01  | dbg14 3.23e-01 3.63e-02 | B 24 3.58e-01 1.52e-01
    B 12 3.77e-01 9.41e-02 | B 4 4.32e-01 5.63e-02 | B 22 7.74e-01 3.34e-01

Why: with sigma profiled, the log-likelihood in p is `n log p(1-p) - n
log S(p)` with `S(p) = min_b sum rho_p` concave and piecewise linear,
so `-n log S(p)` has an upward kink at every breakpoint of the quantile
process and a local maximum can sit in many intervals near the top.
0.83.0 stopped at the same local maxima (its `gap_j` is no smaller).
Filed (MEDIUM, needs a decision: searching the quantile process costs
O(n) pivots per interval walked). The family help now says to fix the
quantile for a quantile regression.

A fit the proof does not cover still warns, measured three ways:
- `test-ald-fixed.R`, "a stop the proof does not cover still warns":
  with `optCtrl = list(iter.max = 3)` and `ald_vertex()` mocked to
  NULL, the stop is short (the unmocked fit's logLik is higher by more
  than 1e-6 relative), and the fit warns with `frmtmb_conv_code` or
  `frmtmb_conv_gradient` beside the refusal, the code-8 text ending
  "the proof that it is the quantile regression's optimum failed".
- censored rows (same file): `ald_cert$why` names the censoring and
  the warnings stay.
- 0.83.0's own short stops (`rqsum-base.txt`): 43 of 43 warned; on the
  lane none is short.

Rounding: 20 perturbations `y * (1 + 2e-15 k)` of a 60-row data set
(test "the proof holds on every rounding of the data") are proved on
every one.

### 5490: every row censored

`dev/alapse-allcens.R` (data of `dev/wcheck2-allcens.R`, seed
20261021, 80 rows right censored at their minimum), logs
`allcens-base.txt`, `allcens-lane.txt`:

    rellib-r22, asym_laplace (1 | g4): code 1, logLik -0, log sigma -38.61; warnings: start, false convergence (8); analysis lost (Intercept)=limit, x=limit, sigma_(Intercept)=flat, quantile_(Intercept)=limit, theta_1=flat; se_explained convergence
    rellib-r22, asym_laplace fixed:    code 1, logLik -2.687e-215, logit quantile -498.4; warnings: start, false convergence (8); lost (Intercept)=flat, x=flat, sigma_(Intercept)=flat, quantile_(Intercept)=limit
    rellib-r22, gaussian (1 | g4):     code 0, the limit warning ("the maximum likelihood estimate does not exist")
    rellib-r22, gaussian fixed:        code 1 (9), the limit warning (se_boundary_stop() made it a boundary stop)

Why the limit rule does not reach the asymmetric Laplace: the optimizer
runs until every probability is 1 in double precision (the quantile at
3.4e-217 without the group term, sigma at exp(-38.6) with it), where
every direction the edge probes move along is level, so the analysis
calls some parameters "flat" and `se_boundary_stop()` (all lost must be
"boundary" or "limit") does not fire; the stop stays a convergence
warning and `se_check()` never reports.

The right rule, `ald_cens_sep()`, a proof from the data alone: every
positive-weight row right censored (or every row left censored), an
intercept column in mu on the identity link free in that direction, no
prior density, no truncation, only gaussian frmtmb blocks in mu. Then
each row's probability is below 1 at every finite value and rises to 1
as the intercept grows, whatever sigma, the quantile and the effects
are; the Laplace approximation is below 1 too (each row's log
probability is concave in mu because the density is log-concave, so
the data only add curvature: `|Sigma|^-1/2 |Sigma^-1 + C|^-1/2 <= 1`)
and goes to 1 along the same path. It feeds `check_convergence()`'s
separation path (class `frmtmb_separation`), so every parameter is
lost as "separation". Lane log: both asym_laplace fits now warn "Every
row of the asym_laplace() response is right censored, so a larger
intercept of mu raises every row's probability toward 1, which no
finite value reaches: the likelihood has no maximum, and the maximum
likelihood estimates do not exist", with no code-8 warning; gaussian
unchanged. Guard-absent: 72 of 80 censored (8 exact rows) is not
named (`test-ald-fixed.R`). An objective-based rule (every row's
probability above max(p, 1 - p), then sigma to 0) was written first
and did not fire on the fixed fit, where the run-off is in the
quantile (`dev/alapse-dbg-cens.R`); it was replaced.

### 5449: the models the 0.83.0 refusal leaves open

The proof: with mu linear in the effects b on the identity link and
every row exact, the data add no curvature, H = Sigma^-1, and the
approximation is `max_b (l(b) - b' Sigma^-1 b / 2)`, which does not
fall as Sigma is scaled up. It needs (1) a gaussian block whose
covariance a free scale multiplies, (2) zero data curvature in b.

Extended (a proof, and profiled): hetar1, ou, toep, homtoep, homcs, exp,
gau, mat, gr(prec = ). For each, the structure's own `vcov()` at sds
times 3 is exactly 9 times its own at a random theta (premise (1)), and
TMB's Laplace objective profiled along every log sd (+0, 1, 2, 4, 8,
the rest re-optimized) never falls (`Rscript dev/alapse-ald-profile.R
lane`, seed 20261116, 20 groups x 4 times x 2 rows; logs
`ald-profile-1.txt`, `ald-profile-2.txt`):

    ou       -185.158 -113.301 -93.245 -92.084 -91.666   never falls TRUE
    toep     -170.747 -86.820 -83.672 -83.169 -83.089    never falls TRUE
    hetar1   -154.652 -102.409 -91.790 -90.680 -90.240   never falls TRUE
    homcs    -109.386 -96.038 -94.232 -77.934 -77.929    never falls TRUE
    homtoep  -90.576 -86.486 -84.960 -84.392 -84.392     never falls TRUE
    exp      -113.812 -101.835 -96.429 -94.555 -78.636   never falls TRUE
    gau      -89.324 -73.863 -71.771 -71.449 -71.443     never falls TRUE
    mat      -110.419 -82.941 -80.228 -79.528 -79.508    never falls TRUE
    grprec   -254.584 -253.347 -252.927 -251.575 -251.525 never falls TRUE

(vcov ratio 9 to 9 on all nine; homtoep 9 to 9.00000000000001.) On
0.83.0 each of these fits ended at code 1 with both convergence
warnings, at logLik values 30 to 60 below where the profile climbs
(e.g. homcs -135.42 against -77.93). Guard-absent cases in
`test-ald-laplace.R`: the same nine terms with gaussian() and with
`cens()` rows are not refused. Compat rows changed (by the loop over
`ald_scaled_covstructs`): `hetar1 x asym_laplace`, `ou x asym_laplace`,
`toep x asym_laplace`, `homtoep x asym_laplace`, `homcs x
asym_laplace`, `exp x asym_laplace`, `gau x asym_laplace`, `mat x
asym_laplace`, `gr_prec x asym_laplace`: now "conditional" with the
refusal note.

Left open, and why (no proof):
- another link: `rho'(u) f''(eta)` curves mu in b; wcheck2's profile
  of the log link went non-finite (`ald-profile-lane-*.txt` copies).
- censored rows: a tail probability is log-concave, so it adds
  curvature C, and `-1/2 log |I + Sigma C|` falls as Sigma grows;
  wcheck2's profile fell (-98.461 to -138.404).
- a random effect in sigma: smooth in sigma, and couples with mu's
  effects; its profile fell (-131.4 to -133.7).
- an exact gp(): the proof holds on paper (sdgp^2 times a correlation
  plus a nugget scaled the same way), but TMB's inner search fails on
  its near-singular covariance at large sd and the objective falls
  there (-129.315 ... -127.528, -261.369 at +4, +8, wcheck2); left
  open, as the profile does not confirm it.
- rr(): the covariance is Lambda Lambda', scaled by scaling the
  loadings; not profiled along a common scale of the loadings (the
  profile script moves log sds), so left open.
- a Student-t block in mu: its curvature depends on its mode (wcheck2
  S1, fell by 10.4); another package's block: no density is known here.
- mi()/me(): a latent with a gaussian measurement model entering mu
  linearly through a fixed slope keeps H block diagonal and the proof
  holds (wcheck2's me profile never fell), but a latent in a varying
  slope couples with b (`b * xi` is bilinear, so rho' times it curves);
  the check reads neither, so left open.

The quadrature breakdown ("the Gauss-Kronrod objective broke down"):
see the section below.

### 5449: the quadrature breakdown, diagnosed

`Rscript dev/alapse-quad.R lane` (log `quad-lane.txt`; the review's g4
design: 80 rows, 20 groups of 4, seeds 1 to 10, `y ~ x + (1 | g)`,
`quadrature = TRUE`). 9 of 10 fits stop with "the calibration expired
on '1 | g' ... at every point tried"; seed 10 fits with
`frmtmb_quad_unsettled` and both convergence warnings at logLik
-123.29. TMB's Laplace objective run by nlminb from the default start
(the warm start `quad_fit()` anchors on) ends at log sd 0.001 to 2.306,
and with log sd 2 higher and the rest re-optimized it is higher on 10
of 10, by 1.133 to 15.332: the anchor is not a maximum, as the proof of
`check_ald_laplace()` says.

Diagnosis. TMBad's marginal_gk transform freezes one width per effect,
read from the log-integrand's curvature at its mode by differences of
step 1 (`quad_fit()`'s roxygen; upstream TMB-7: the step does not reach
from an RTMB spec). With a log density linear in b on each side of
every kink, that curvature is the prior's 1 / sd^2 plus whatever kinks
fall inside the step, while the integrand's real width is of order
sigma over the rows per group; and the anchor's sd is wherever TMB's
inner search stopped. So the frozen widths are wrong at the anchor and
at every displaced point `quad_fit()` tries. This is the calibration's
premise (a smooth log-integrand with curvature at its mode) failing,
not an arithmetic defect of the integrator, and the fix (a rule whose
widths read the kinks, the quadseam2 item's route (a), a Gauss-Hermite
rule on the tape) is quadrature's, lane quadesc's file. Filed (MEDIUM).
No quadrature code was changed.

## Tests and runs

New and changed test files, one per process (`dev/release/run-tests.R`
with `FRMTMB_RELLIB` = the lane library; logs `tf-*.txt`,
`dev/release/alapse-h.log`, `dev/release/alapse-ald.log`):

    rellib-r22: RESULT test-ald-fixed.R pass=12 fail=69 err=3 skip=0 warn=3
    lane:       RESULT test-ald-fixed.R pass=102 fail=0 err=0 skip=0 warn=0
    rellib-r22: RESULT test-ald-laplace.R pass=22 fail=1 err=1 skip=0 warn=0
    lane:       RESULT test-ald-laplace.R pass=58 fail=0 err=0 skip=0 warn=0
    lane:       RESULT test-se-limit.R pass=380 fail=0 err=0 skip=0 warn=0
    lane:       RESULT test-cdfs.R pass=234 fail=0 err=0 skip=0 warn=0
    lane:       RESULT test-mi-families.R pass=194 fail=0 err=0 skip=0 warn=0
    lane:       RESULT test-bracket-access.R pass=36 fail=0 err=0 skip=0 warn=0
    lane:       RESULT test-compat.R pass=868 fail=0 err=0 skip=0 warn=0

On rellib-r22 test-ald-fixed.R fails as behavior: the fits give
frmtmb_conv_code and frmtmb_conv_gradient where only frmtmb_ald_no_se
is expected, the lost reasons are not "kink", the code is 1, rq's
objective is above the fit's by 1.6e-8 relative and the coefficients
differ by 1.6e-5 (case 1), the all-censored fits give no separation
warning; the mock of `ald_vertex()` errors (the weak form, the
function is new). test-ald-laplace.R's new block fails as behavior
(hetar1 not refused). Before the allowed-warning edits, the lane's
test-mi-families.R failed 3 (the mi() fit was not proved while the
complete-case one was, 8.6e-7 apart; fixed by reading mi()'s unread
rows as weight 0) and test-se-limit.R 2 (the 5423 test expected
"nocurv").

OpenBLAS 0.3.26 with LAPACK (`bash dev/ciharden-openblas.sh 0.3.26
lapack`: Rblas 74 routines to OpenBLAS, Rlapack 616, LAPACK 3.12.0; log
`openblas026-build.txt`; `bash dev/alapse-ob.sh`, logs `ob026-*.txt`):

    RESULT test-ald-fixed.R pass=102 fail=0 err=0 skip=0 warn=0
    RESULT test-ald-laplace.R pass=58 fail=0 err=0 skip=0 warn=0
    RESULT test-cdfs.R pass=234 fail=0 err=0 skip=0 warn=0
    RESULT test-mi-families.R pass=194 fail=0 err=0 skip=0 warn=0
    RESULT test-se-limit.R pass=380 fail=0 err=0 skip=0 warn=0

and the rq comparison, seeds 1 to 5 (`rq-ob026.txt`): 39 of 39 fits
proved (code 0), largest `gap_b` 5.68e-14, no convergence warning.
`dev/ciharden-out/` deleted. 0.3.32 was not run: no assertion here
sits near a rounding-decided branch (the vertex margin is a strict
inequality, and the 20-perturbation test passes on every rounding).

Every suite, one file per process, two at a time, with `NOT_CRAN=true`
and `FRMTMB_BRMS_FIT_TESTS=true` (`bash dev/release/run-par.sh
alapse-all dev/alapse-jobs-all.txt 2`, log `dev/release/alapse-all.log`,
"ALAPSE-ALL ran 550 of 550"; 550 of 550 `lib:` lines name
`wt-alapse-lib` for frmtmb; the extensions load from rellib-r22, as no
extension changed):

    frmtmb files=350 pass=28895 fail=0 err=0 skip=15 warn=1
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=33 pass=1857 fail=0 err=0 skip=3 warn=0
    frmtmb.latent files=13 pass=700 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=16 pass=545 fail=0 err=0 skip=2 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=92 pass=4651 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial files=8 pass=634 fail=0 err=0 skip=0 warn=0
    frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0

The one warning was the new refusal escaping the gated
`test-brms-lik-cdfs.R` (its allowed list named only the convergence
warnings); with "no Hessian standard error" added the file reads
`RESULT test-brms-lik-cdfs.R pass=70 fail=0 err=0 skip=0 warn=0`
(`tf-brmslikcdfs-lane.txt`).

## The routes the refusal names (after lane squant's measurement)

Lane squant (its `dev/squant-findings.md`, "(c)" and "The asym_laplace
refusal's alternative"; `dev/squant-boot-ald.R`) measured on rellib-r22
that `frm_bootstrap()` on `asym_laplace()` (gaussian data, 400 rows,
`quantile = 0.9`, seeds 20261117 to 20261119, 200 draws) gives slope
SEs 0.1062, 0.1254, 0.1116 against a Monte Carlo sd of rq()'s slope of
0.2977: the parametric bootstrap draws from the fitted asymmetric
Laplace, so it is model-based and not calibrated for quantile
regression. The first version of this lane named it as the first
route. Now `ald_se_routes` names `quantreg::rq()` with `se = "nid"` or a
bootstrap of the rows first, and `frm_bootstrap()` and
`frmtmb.sample::frm_sample()` only "where the data are asymmetric
Laplace", with the measured shortfall; the family help, NEWS and the
emmeans compat note say the same. The lane-build check of the same
script is below ("bootstrap on the lane").

Suggested wording for the consolidation, which adds quantile_elf() in
one place, as the first element of `ald_se_routes` in `R/ald-fixed.R`:

    paste0("quantile_elf(), a smoothed asymmetric Laplace that has ",
           "curvature, with vcov_cluster(fit, ~ row, \"CR0\") for its ",
           "standard errors")

and in `ald_laplace_message()` (`R/ald-laplace.R`), before the
frm_sample() sentence: "quantile_elf() fits the same quantile with a
log density that has curvature, so its Laplace approximation has a
maximum."

Every suite in one session, three orders (`bash
dev/release/run-sessions.sh alapse-sess 2`, `FRMTMB_RELLIB` = the lane
library, `NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS=true`; log
`dev/release/alapse-sess.log`, "ALAPSE-SESS ran 27 of 27"): every
SESSION line fail=0 err=0 warn=0 (frmtmb forward 350/350 pass=28895
skip=15, reverse pass=29060 skip=15, preload pass=29164 skip=17; the
extensions as in the per-file run, in all three orders). 218 LEAK
lines, every one tagged `[with load of ...]` or `[extension session
state]`; 0 untagged.

Bootstrap on the lane (`Rscript dev/alapse-boot.R lane`, log
`boot-lane.txt`, squant's design): slope bootstrap SE 0.1072, 0.1261,
0.1123 (200 of 200 refits converged on each, where rellib-r22 counted
38, 47 and 39), rq's nid SE 0.2451, 0.3623, 0.3589, against squant's
Monte Carlo sd 0.2977. The polish makes every refit converge and
leaves the bootstrap's shortfall as it was: 0.36 to 0.42 of the Monte
Carlo sd. `test-ald-fixed.R` after the routes change: pass=102 fail=0
err=0 skip=0 warn=0 (`tf-ald-fixed-lane2.txt`).

## What this lane did NOT do, and why

- No sparsity or sandwich covariance (the brief: a later decision).
- No search of the quantile process for an estimated quantile's highest
  maximum (filed MEDIUM; needs a decision).
- No quadrature change: the breakdown is the calibration's premise,
  in lane quadesc's files (filed MEDIUM).
- No refusal for an exact gp(), rr(), me()/mi(), another link, censored
  rows or a random effect in sigma under the Laplace approximation (no
  proof, or none confirmed; filed LOW).
- The proof of `ald_polish()` covers the identity link only; a log-link
  asym_laplace() fit without random effects keeps the convergence
  warnings beside the refusal.
- A degenerate vertex (ties in y with a discrete design) can cycle the
  simplex steps; the iteration cap then fails the proof and the fit
  keeps its warnings. Not measured on such data. (WITHDRAWN in fix round
  1: the review measured it, 17 of 18 integer-y fits kept their
  warnings; `ald_vertex_ties()` now proves all 18, see "Fix round 1".)
- `ald_lp_blocks()` makes the location design dense (`as.matrix()`); a
  very wide sparse design would cost memory there. Not measured.
- No equivalence-survey pair or ported-ledger row touches an
  asym_laplace() fit without random effects or the new structures
  (`dev/equiv-data/results.tsv` names asym_laplace only in p072,
  gr_cov, already a decided refusal; the bin-2 ledger rows are the
  distribution functions and log_lik), so none was remeasured and no
  manual ledger row changed.

## Bump and floors

Minor (0.84.0): behavior changes marked BREAKING (the refusal, the
moved estimates and code 0, nine more structures refused), a new
warning class `frmtmb_ald_no_se`, `fit$opt$ald_cert`. No extension
changed and none needs a new floor: every extension suite passed
against this core with the extensions of rellib-r22.

## Files touched

R/ald-fixed.R (new), R/ald-laplace.R, R/compat.R, R/confint.R,
R/families.R, R/fit.R, R/interop.R, R/se-check.R, man/frmtmb-families.Rd
(roxygen), NEWS.md, vignettes/diagnostics.Rmd,
tests/testthat/test-ald-fixed.R (new), tests/testthat/test-ald-laplace.R,
tests/testthat/test-brms-lik-cdfs.R, tests/testthat/test-cdfs.R,
tests/testthat/test-mi-families.R, tests/testthat/test-se-limit.R,
dev/test-backlog.md (3 CLOSED lines, section "Filed in the round of
2026-11-16"), dev/alapse-*.R, dev/alapse-*.sh, dev/alapse-findings.md,
dev/alapse-jobs-*.txt (gitignored), dev/alapse-prev-out/ (copies,
gitignored).

## R CMD check

`bash dev/alapse-check.sh frmtmb` (build and `--as-cran` inside
`dev/alapse-check/`, `R_LIBS` lane, rellib-r22, user; pandoc and
TinyTeX on PATH; `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`): tests "[45m]
OK", vignettes "[341s] OK", `Status: 1 NOTE` (the expected V8 math
rendering NOTE). Log kept as `dev/alapse-out/check-frmtmb.txt`; the
check tree was deleted. A first run under `timeout 3500` was killed
while rebuilding the vignettes (habit.Rmd), which printed a vignette
ERROR with no message (`check-frmtmb-timeout.txt`); that run is void.
No extension changed, so no extension was checked.

## Fix round 1 (review `dev/reviews/2026-11-17-alapse.md`, verdict FIX)

### 1 (HIGH): code 0 only where the maximum is proved the highest

The review showed the first build giving code 0 and no warning to fits
at a lower local maximum (design B: 24 of 46 below the global maximum,
seed 109 by 2.65 at quantile 0.420 against 0.725; `sigma ~ z` and
`quantile ~ z`: 9 of 20 each, up to 0.313 and 1.885), where 0.83.0
warned on every one. Now (`R/ald-fixed.R`):

- `ald_local()` is the old polish loop; `ald_polish()` runs it, then
  `ald_global_case()` decides which global proof applies:
  - fixed quantile, constant sigma (log link, intercept only), one
    response, no prior density, no other parameter but zi's: "convex".
    The proof: the negative log-likelihood of the asymmetric Laplace
    part is `W log sigma + sum w_i rho_p(y_i - x_i b) / sigma`; in
    `t = 1 / sigma`, `a = b / sigma` it is
    `-W log t + sum w_i rho_p(t y_i - x_i a)` (rho_p is positively
    homogeneous), a sum of convex functions, so every local maximum is
    the global one; a zi part is a logistic regression on the zeros,
    convex and apart.
  - estimated constant quantile (logit link, intercept only), same
    otherwise: "search", `ald_quantile_search()`: branch and bound over
    the quantile p on the profile `W (log p + log(1 - p) - log(S(p) / W)
    - 1)`, `S(p)` the exact LP minimum (concave, piecewise linear). On
    each interval S lies above its chord (through 0 at the ends), which
    gives an upper bound with a closed-form maximum (a root of `c1 p^2 +
    2 c0 p - c0`); it is exact where no breakpoint lies inside, so the
    search splits at the bound's maximum. When a higher maximum is
    found, the local run restarts there.
  - anything else (a `sigma ~` or `quantile ~` model, a fixed quantile
    with `sigma ~ z`, more than one response, a prior density): local
    only.
- Code 0 and no convergence warning only when proved global. A local
  proof keeps the location solved exactly, sets code 1 (the optimizer's
  own code in `ald_cert$code`), and warns once, class
  `frmtmb_ald_local` (`ald_local_message()`: why the proof is local, and
  to fix the quantile and drop the scale model). The refusal's text and
  `diagnose()` say "the highest maximum" or "a local maximum".

The fixed quantile with `sigma ~ z`, measured before deciding
(`Rscript dev/alapse-fixsig.R lane`, log `fixsig-lane.txt`; design D of
`dev/alapse-rq.R` and a stronger one, sd exp(1.2 z), quantile 0.5 and
0.2, seeds 4001 to 4020, 15-start alternating maximization):

    SUMMARY D p 0.2: n 20 | gap > 1e-6: 1 | max gap 1.656e-03
    SUMMARY H p 0.2: n 20 | gap > 1e-6: 3 | max gap 1.526e-01
    SUMMARY D p 0.5: n 20 | gap > 1e-6: 0 | max gap 9.283e-11
    SUMMARY H p 0.5: n 20 | gap > 1e-6: 0 | max gap 6.543e-11

4 of 40 below at 0.2, so it is not global in general: local, with the
warning.

On the fixed build:
- the review's `alapse-review-qglob.R lane "c(2,3,4,12,20,22,101:140)"`
  (log `fr1-qglob-lane.txt`): `B summary: n 46 | gap > 1e-6: 0 | proved
  and gap > 1e-6: 0 | ... | max gap 8.86e-11 | max dp 2.32e-07`; 46 of
  46 code 0. B109 now: quantile 0.7247, logLik -134.505740 (the global
  value), 34 linear programs, `moved_from` 0.4197 (`fix1-lane.txt`).
- the review's `alapse-review-multistart.R lane 1:20`
  (`fr1-multistart-lane.txt`): D and E fits are still below the
  multi-start maximum on 9 of 20 each (the same points), but now every
  one is code 1 with `frmtmb_ald_local` (the script's "silent" counts
  only the generic convergence classes). Three of them had nlminb's
  code 0 in a build between; the final build gives 40 of 40 code 1
  (`fr1-multistart2-lane.txt`).
- `dev/alapse-rq.R lane 1:25` (`fr1-rq-lane.txt`, summary
  `fr1-rqsum-lane.txt` by `dev/alapse-rqsum2.R`): A, B, C, F, G and the
  four cases are 125 + 4 of 129 code 0 with joint gaps at most 8.2e-11;
  D and E 25 + 25 of 50 code 1 with `ald_local`; no fit gives a generic
  convergence warning; largest `gap_b` 1.14e-13.
- cost (`dev/alapse-searchcost.R`, `searchcost-*.txt`, y ~ x, seed 7,
  minimum of 3): n 100, 1000, 5000 take 0.03, 0.04, 0.20 s on the lane
  (37, 38, 42 linear programs) against 0.02, 0.03, 0.06 s on
  rellib-r22, which stops short with a warning. `ald_search_evals` caps
  the search at 600 programs; none came near it.

### 2 (MEDIUM): tied data

`ald_vertex_ties()`: when `ald_vertex()` returns NULL, the search runs
on y perturbed by a fixed pattern of size 1e-9 (1 + max|y|), where no
vertex is degenerate; the basis then gives b on the data as they are,
and the proof is made there: every nonzero residual keeps its sign (else
NULL), each zero row off the basis takes the end of its interval its
perturbed residual's sign picks, and the basic rows' s must fall
inside. Those are subgradients, so the proof is exact; the perturbation
only finds the basis. The review's `alapse-review-rq.R lane`
(`fr1-review-rq-lane.txt`): 22 of 22 fits proved (`cert ok TRUE`),
the 18 integer-y fits among them (17 failed before); every LP equal to
`quantreg::rq.wfit(method = "br")` (rel 0 to -1.18e-16); the local
probes find no higher point (largest gain 1.78e-15, even-median's flat
face). Tied fits with a fixed quantile and constant sigma are global
(code 0); the `sigma ~ g` ones are local (code 1, `ald_local`).

### 3 (MEDIUM): profile intervals

`confint(method = "profile")` and `method = "uniroot"` stop, class
`frmtmb_ald_no_se`, for a refused parameter, with the refusal's words
(`R/confint.R`). Reason: both read the same working likelihood, whose
intervals the refusal says hold only where the data are asymmetric
Laplace; the first build gave 0.464 to 0.773 for x and let TMB's
warning escape. Tested.

### 4 (LOW)

- "strict" corrected above, in `ald_polish()`'s roxygen and here.
- The separation relabel's measured effect (review section 5,
  `alapse-review-relabel.R`): two wordings changed, both right:
  `bern_complete` (seed 4) from "nocurv" (whose text named
  asym_laplace() for a bernoulli) to "separation", and `probit_factor`
  (seed 7) from "limit (at 6.391, toward +Inf)" to "separation"; five
  other separated fits unchanged; predictions unchanged.
- The quadrature mechanism is marked inferred in the findings above and
  in the backlog item ("No calibrated width was read").
- The exact gp() reason restated (`R/ald-laplace.R` comment and the
  `gp_pred() x asym_laplace` compat note): its fall at a large sd is
  TMB's inner search failing, as it does for `(1 | g)`, so it is no
  evidence against the proof; it stays open because no profile confirms
  it.
- Backlog line 4917 wrapped; filed LOW: the generic "flat" words on
  `sigma ~ (1 | g)` fits, and the unbounded likelihood with one exact
  row above every censoring point. The HIGH item replaces the MEDIUM
  one and says the first build removed the warning.

### Fix round 1: tests and runs

New and changed test files (one per process; `dev/release/alapse-fr1.log`,
`alapse-fr1b.log` after the `<<-` fix below):

    rellib-r22: RESULT test-ald-fixed.R pass=15 fail=121 err=4 skip=0 warn=6
    lane:       RESULT test-ald-fixed.R pass=161 fail=0 err=0 skip=0 warn=0
    lane:       RESULT test-ald-laplace.R pass=58, test-se-limit.R 380,
                test-cdfs.R 234, test-mi-families.R 194,
                test-brms-lik-cdfs.R 70 (gated), test-bracket-access.R 36,
                test-compat.R 868, test-diagnostics-ux.R 134,
                test-se-check.R 134, test-no-superassign.R 16; all fail 0
                err 0 warn 0

On rellib-r22 the new blocks fail as behavior (`fr1-tf-ald-fixed-base.txt`):
B109 at quantile 0.42 with `frmtmb_conv_gradient`, 0.0193 relative below
the quantile process's maximum; D7 and E18 with the generic convergence
classes and no proof; the tied fits with a warning; the profile intervals
finite. The first lane build failed them too, by the review's logs
(`dev/reviews/alapse-review-out/qglob-lane.txt`: B 109 code 0, no
warning, 2.65 below; `multistart-lane.txt`: D 7 and E 18 silent;
`rq-lane.txt`: the tied fits unproved with warnings).

The rounding test runs 21 roundings (k = 0 to 20) of a continuous and an
integer-y data set; every one is proved the highest maximum. OpenBLAS
0.3.26 with LAPACK, on the final build (`ob026-*.txt`):
test-ald-fixed.R pass=161, test-ald-laplace.R pass=58 (and, on the build
before the `<<-` fix, test-se-limit.R 380, test-cdfs.R 234,
test-mi-families.R 194), all fail 0 err 0 warn 0. `dev/ciharden-out/`
deleted.

A first full run stopped at test-no-superassign.R (2 failures:
`ald_quantile_search()` used `<<-`); its driver and processes were
stopped (matched on this worktree's path), the state moved into an
environment, and the run restarted.

Every suite, one file per process, gated tier on (`bash
dev/release/run-par.sh alapse-all2 dev/alapse-jobs-all.txt 2`, log
`dev/release/alapse-all2.log`, "ALAPSE-ALL2 ran 550 of 550"; every
`lib:` line names `wt-alapse-lib`):

    frmtmb files=350 pass=28954 fail=0 err=0 skip=15 warn=0
    frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam files=33 pass=1857 fail=0 err=0 skip=3 warn=0
    frmtmb.latent files=13 pass=700 fail=0 err=0 skip=2 warn=0
    frmtmb.learn files=16 pass=545 fail=0 err=0 skip=2 warn=0
    frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample files=92 pass=4651 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial files=8 pass=634 fail=0 err=0 skip=0 warn=0
    frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0

Every suite in one session, three orders (`bash
dev/release/run-sessions.sh alapse-sess2 2`, log
`dev/release/alapse-sess2.log`, "ALAPSE-SESS2 ran 27 of 27"): every
SESSION line fail=0 err=0 warn=0 (frmtmb forward 350/350 pass=28954,
reverse 29119, preload 29223); 218 LEAK lines, 0 untagged.

### Fix round 1: files touched

R/ald-fixed.R (ald_vertex_ties(), ald_vertex_try(), ald_local(),
ald_global_case(), ald_quantile_search(), ald_search_evals,
ald_local_message(), the refusal's wording), R/fit.R
(check_convergence(): the `frmtmb_ald_local` warning and the code-8
branch gated by the proof), R/confint.R (the profile and uniroot
refusal; diagnose()'s wording), R/ald-laplace.R (exact gp() comment),
R/compat.R (the `gp_pred() x asym_laplace` note restated),
R/families.R (help section), man/frmtmb-families.Rd, NEWS.md,
vignettes/diagnostics.Rmd, tests/testthat/test-ald-fixed.R,
dev/test-backlog.md (HIGH item, quadrature item, two LOW items, the 4982
CLOSED line), dev/alapse-fixsig.R, dev/alapse-fix1.R,
dev/alapse-searchcost.R, dev/alapse-rqsum2.R, this file.

R CMD check --as-cran on the fix-round tree (`bash dev/alapse-check.sh
frmtmb`, log kept as `dev/alapse-out/check-frmtmb-fr1.txt`, check tree
deleted): tests "[40m] OK", vignettes "[272s] OK", `Status: 1 NOTE`
(the expected V8 NOTE).

## Fix round 2 (short review `dev/reviews/alapse-short.md`, verdict FIX)

The short review confirmed the convexity proof, the search (80 of 80 at
the exact maximum within 8.4e-11 by brute force), the ties proof (none
of 365 adversarial non-optimal bases certified), the local warning and
the confint refusal.

1. HIGH: the search ignored the optimizer's box. Design B seed 109 with
   a bound on the quantile's intercept ([-2, 0.3]) or sigma's
   ([-1.1, 2]) came back at logit 0.9681 or log sigma -1.2669, outside
   the bound, code 0, "global" (the review's
   `alapse-short-bounds-lane.txt`). Now `ald_global_case(fit, blocks,
   bounds)` allows "search" only when sigma's and the quantile's
   coefficients have no finite bound (else local, with its own reason
   and advice); `ald_polish()` accepts the restarted run only when
   `inside(run2$par)`; the concave case keeps a box on sigma (an
   interval in t is convex; the local proof needs an interior point,
   and the local run's line search stays inside). The review's script
   on the fixed build (`fr2-bounds-lane.txt`): all three cases inside
   the box, code 1, `frmtmb_ald_local`, at the local maximum p 0.4197
   that rellib-r22 also stops at (with its generic warnings there).
   Tests: the two bounded fits are inside their box, local and warned;
   a fixed quantile with a sigma box is inside, global, code 0.
2. MEDIUM: the constant's link. Any link in `dpar_links_positive` for
   sigma and `dpar_links_unit` for the quantile (each maps the
   coefficient one-to-one onto the range); the search restarts through
   the links' `linkfun`/`linkinv`. Test: softplus sigma, fixed
   quantile, code 0, global, logLik within 1e-10 relative of the
   log-link fit; an identity-link estimated quantile, global.
3. LOW: filed the two-response case; the local warning's advice now
   comes from the reason (`ald_cert$advice`): "Each response fitted
   alone ..." for several responses, "Without the prior density, ..."
   for a prior density, the bound, the quantile or the scale model
   otherwise. Test: two asym_laplace responses warn with "Each
   response fitted alone" and not "drop the scale model".

Seen failing: on rellib-r22 the new blocks fail as behavior
(`fr2-tf-base.txt`: RESULT test-ald-fixed.R pass=22 fail=132 err=5
warn=8; the bounded fits give the generic convergence classes, the
softplus fit code 1, the two-response fit no `frmtmb_ald_local`). On
the fix-round-1 build, by the short review's logs: the bounded fits
outside their box with code 0 (`alapse-short-bounds-lane.txt`), the
softplus fit `frmtmb_ald_local` with code 1 (`alapse-short-over.txt`).

Runs:

    lane: RESULT test-ald-fixed.R pass=181 fail=0 err=0 skip=0 warn=0
    lane: RESULT test-ald-laplace.R pass=58 fail=0 err=0 skip=0 warn=0
    OpenBLAS 0.3.26 + LAPACK (`ob026-*.txt`): test-ald-fixed.R 181,
    test-ald-laplace.R 58, both fail 0 err 0 warn 0; dev/ciharden-out
    deleted.
    core, one file per process, gated tier on (`dev/release/alapse-core3.log`,
    "ALAPSE-CORE3 ran 350 of 350"; every `lib:` line wt-alapse-lib):
    frmtmb files=350 pass=28974 fail=0 err=0 skip=15 warn=0
    core in one session, three orders (`dev/release/alapse-sess3.log`,
    "ALAPSE-SESS3 ran 3 of 3"): forward pass=28974, reverse 29139,
    preload 29243, each fail=0 err=0 warn=0; 135 LEAK lines, 0 untagged.

R CMD check --as-cran on the fix-round-2 tree (`bash dev/alapse-check.sh
frmtmb`, log `dev/alapse-out/check-frmtmb-fr2.txt`, check tree
deleted): tests "[38m] OK", vignettes "[266s] OK", `Status: 1 NOTE`
(the expected V8 NOTE).

Fix round 2 files touched: R/ald-fixed.R (`ald_global_case()` with
bounds, links and advice; `ald_polish()`'s restart through the links
and `inside(run2$par)`; `ald_local_message()`), R/families.R (help),
man/frmtmb-families.Rd, NEWS.md, tests/testthat/test-ald-fixed.R (three
blocks), dev/test-backlog.md (one LOW item), this file.
