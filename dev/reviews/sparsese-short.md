# Short check of lane sparsese, fix round 1 (2026-11-17)

Scope: F1 (`jcs_root_at()`, `lb_dense()`), F2 (`jcs_fd_root()`,
`fd_joint_root()`) and `jcs_psd_factor()` against an independent analytic
reference, not the package's dense route. Lane library `wt-sparsese-lib`
(installed after the last source edit), base `rellib-r22`. No package code
edited, no git operation. At most two R processes; the largest probe
peaked at 1801 MB. Scripts: `dev/reviews/sparsese-short-*.R`. Logs:
`dev/reviews/sparsese-short-out/`.

## Verdict: SHIP

F1 and F2 are fixed, and the algebra is exact. Given the package's own
outer covariance `Vf`, the sparse route agrees with the analytic reference
to 1.3e-8 in standard error and 3.4e-9 in covariance (F1), and 4e-8 (F2),
on every design, seed and size tested, lost-sd fits included. All remaining differences from the reference are in
`Vf` (the outer Hessian) and in the semantics of a held parameter. They
are the same on the dense route and on `rellib-r22`, and they are not in
this lane's scope. One latent LOW in new code (item 3b) is optional.

## 1. F1: frm_lp_basis() V_root and frm_curve() .se

Reference: `sparsese-review-indep.R`'s analytic joint covariance, copied
into `sparsese-short-f1int.R` with ONE line changed (see "Correction to the
review" below). Rows: 11 grid rows with no group effect, 11 rows at seen
levels, and every training row (100 groups only). "pkg Vf" puts the
package's `vcov(full = TRUE)` in place of the reference's own Richardson
`Vf`. This isolates the route from the outer Hessian.

| run (route) | lost | V_root recipe and frm_curve .se, max rel (own Vf / pkg Vf) | V_root V_root' + V_rest vs natural ref, scaled max (pkg Vf) |
|---|---|---|---|
| 100 s1 sd 0.5, limit 20 | none | 1.4e-7 / 2.1e-9 | 1.7e-9 |
| 100 s2 sd 0.5, limit 20 | none | 3.3e-9 / 3.0e-9 | 3.4e-9 |
| 100 s3 sd 0.5, limit 20 | none | 4.4e-9 / 1.3e-9 | 1.8e-9 |
| 100 s1 sd 0, limit 20 | sd, boundary | 2.5e-9 / 2.5e-9 | 4.1e-10 |
| 100 s2 sd 0, limit 20 | sd, boundary | 3.1e-9 / 3.1e-9 | 5.2e-10 |
| 100 s3 sd 0, limit 20 | none | 9.5e-6 / 8.8e-10 | 1.9e-9 |
| 4500 s1 sd 0, default | none | 1.1e-6 / 1.3e-8 | 1.7e-10 |
| 4500 s2 sd 0, default | sd, boundary | 8.0e-9 / 8.0e-9 | 3.9e-10 |
| 4500 s3 sd 0, default | sd, boundary | 8.7e-9 / 8.7e-9 | 4.2e-11 |

- `V_root` is present and `V_rest` is all zero in every run.
- `frm_curve()` returns in every run, with `cov_rel_error` from 7.6e-11
  to 9.7e-10.
- The natural `A V A'` stays 4.9e-4 to 1.5e-2 off. This is the loss the
  fix avoids.
- On lost fits, the b-row variances equal the reference: 6.45e-11 (s1),
  1.656e-10 (s2), 4.985e-10 and 3.479e-10 (4500 s2, s3). The dense and
  sparse factors agree to 4.6e-16 (`lostb-s1.txt`, `lostb-s2.txt`).

Random slope, `(1 + x1 | g)`, degree 5, 30 groups, limit 1e6 against 20
(`sparsese-short-slope.R`). The reference is written here: the exact
gaussian REML, my own whitener, and an analytic `d Sigma^-1 / d psi`.

| seed | lost | sparse V_root recipe / frm_curve .se / frm_linpred se.fit vs ref, pkg Vf | V_root V_root' + V_rest, scaled (pkg Vf) | own Vf |
|---|---|---|---|---|
| 1 | none | 7.1e-11 / 7.1e-11 / 1.9e-12 | 5.6e-11 | 2.4e-7 |
| 3 | none | 6.1e-11 / 6.1e-11 / 2.0e-12 | 7.0e-11 | 4.9e-7 |
| 2 | correlation, boundary (z = 307) | 3.5e-9 / 1.6e-9 / 1.6e-9 | 2.3e-9 | 5.1e-4 |

At seed 2, "pkg Vf" is `sdr$cov_fixed_prop`, the matrix the routes
propagate. The 5.1e-4 with the reference's own Vf comes from the held
direction: the package's repaired `Vf` against my axis restriction differ
by 4.5e-3 (scaled). The dense and sparse routes print the same digits.

## 2. F2: cumulative fitted() Est.Error

Reference (`sparsese-short-cum.R`), independent of the package's
covariance code. The script writes the logit cumulative likelihood
directly and finds the inner mode by Newton to 1e-13. The REML Laplace
objective uses block elimination, and the mode agrees with the package to
4e-13. `Vf` is the inverse of a Richardson Hessian of that objective in
standard-error-scaled coordinates. `J = -W H_r,psi` uses analytic mixed
partials, and the delta method runs over the joint (gamma, b, psi).

| design | route | new rows, max rel (own Vf / pkg Vf) | in-sample (own / pkg Vf) |
|---|---|---|---|
| s30 seed 1 (38 joint) | dense | 7.2e-5 / 1.6e-8 | 1.4e-4 / 3.4e-8 |
| s30 seed 1 | sparse | 7.2e-5 / 9.7e-9 | 1.4e-4 / 3.2e-8 |
| s30 seed 3 | dense | 4.1e-5 / 2.5e-8 | 6.9e-5 / 4.1e-8 |
| s30 seed 3 | sparse | 4.1e-5 / 1.3e-8 | 6.9e-5 / 3.3e-8 |
| s4500 (4508 joint), default | sparse | 3.7e-7 / 3.4e-9 | 2.3e-6 / 2.5e-8 |

- At 4508, `fitted()` (new rows plus in-sample) takes 4.76 s. The first
  three rows match the review's 62 s dense run in every printed digit.
- Seed 2 loses the sd (`theta_1` nonfinite). The root rows are absent,
  so the K branch is not taken: see section 4.

### Is Vrest + K K' exact? (derivation)

In whitened coordinates `x = (d, b)`, `d` the drawn parameters' loaded
positions `cd`, `V_gamma` positive semidefinite:

- Let `V_dd = Ld Ld'` with `Ld` of full column rank.
- `K_d = M_d Ld`, so `K_d K_d' = M_d V_dd M_d'`, the natural drawn block.
- `K_b = V_bd Ld'^+`.
- `K_d K_b' = M_d Ld Ld^+ V_db = M_d P V_db`, with `P` the projector onto
  `range(Ld) = range(V_dd)`. A PSD matrix has
  `range(V_db) in range(V_dd)`, so this equals `M_d V_db`, the natural
  cross block.
- `K_b K_b' = V_bd (Ld Ld')^+ V_db = V_bd V_dd^+ V_db`.
- `Vrest[b, b] = V_bb - V_bd V_dd^+ V_db`. This is the Schur complement:
  the covariance of the group effects given the drawn parameters.

So the factor is exact. The drawn parameters' whole contribution is in
`K`, and the rest is the conditional term. `M` is the identity on `b`, so
the group-effect rows need no transform.

Measured on s30 against the reference's natural joint covariance (scaled
max `|KK' + Vrest - Vref| / (s_i s_j)`):

- With the package's Vf: 3.9e-8 (seed 1) and 3.2e-8 (seed 3), sparse.
  The dense route gives 2.9e-8 and 2.2e-8.
- `Vrest[b, b]` against the reference's conditional covariance: 1.3e-8
  and 2.0e-8.
- `K` is 38 x 8 on the sparse route (one column per drawn direction) and
  38 x 38 on the dense route.

### Same drawn set as fit_fd_se()?

Yes.

- `fd_joint_root()` builds `pos_d` from `unique(map$comp)` of
  `fit_draw_space()$map`, as `fd_joint_cov()` does. It checks
  `length(pos_d) == length(map$names)`.
- `fit_fd_se()` perturbs `fit_outer_vector(fit, map)` through
  `fit_set_outer(..., map)` with `K$K[seq_len(p), k]`, `p = length(v0)`.
- Logged: drawn `x1 .. x5, theta_1, tau_raw_1, tau_raw_2`, positions
  `beta x5, theta, tau_raw x2`, same length, on both routes.

## 3. jcs_psd_factor() (`sparsese-short-psd.R`)

| input | cols | max \|LL' - V\| / max \|V\| | scaled max | invented variance (max diag excess / max \|V\|) |
|---|---|---|---|---|
| full rank, sds 1e-4 to 1e4 | 8 | 2.1e-20 | 2.3e-16 | +1.1e-20 |
| zero rows 3, 7 | 8 | 8.7e-17 | 1.8e-16 | +8.7e-17 (zero rows exactly 0) |
| rank 4 of 9 | 4 | 1.0e-15 | 1.3e-15 | +5.2e-16 |
| rank 5 of 9 plus zero rows | 5 | 2.4e-15 | 2.7e-15 | +8.0e-16 (zero rows exactly 0) |
| min eig 1e-10, 1e-14, 1e-18 | 8 | 1.4e-16 to 1.8e-16 | 3.1e-16 | +1.4e-16 at most |
| col 7 = 1e3 x col 1 | 6 | 1.1e-15 | 2.5e-15 | +1.1e-21 |
| rank 3, shifted by -1e-16 | 3 | 2.8e-15 | 2.8e-15 | +7.0e-16 |

All are within 1e-10, and no variance is invented beyond rounding. An
all-zero input gives a 3 x 0 factor.

### 3b. LOW (new code, latent, optional): qr.solve() tolerance in jcs_fd_root()

`Kb <- t(qr.solve(Ld, ...))` uses qr's default `tol = 1e-7`. When the drawn
parameters' standard deviations differ by 1e8 or more and they correlate,
qr calls the full-rank `Ld` singular and `fitted()` stops with "singular
matrix 'a' in solve". At sd ratios up to 1e7 it runs (`qrsolve.txt`).
`qr.solve(Ld, ..., tol = 0)` is exact to 1.9e-16 at ratios 1e8 and 1e9.
`Ld` is full column rank by construction. Not reached on a real fit here.
A theta near a boundary that is not flagged lost (sd near 1e4 against a
whitened coefficient sd near 1e-3) is the closest case. The failure is
loud, not silent.

## 4. Pre-existing MEDIUM, reproduced on rellib-r22

`Rscript dev/sparsese-cumlost.R base 2 4000` (`cumlost-base-s2.txt`):

- `lib rellib-r22`. The fit warns that `theta_1` is nonfinite (log sd
  -354.9).
- `fitted()` max Est.Error is 5.5931, and 289 of 900 entries are above
  0.5. This is the same as the lane build.
- `summary()` thresholds have NaN; the coefficients have NA.

The independent reference holds the sd at its estimate (b fixed at its
mode, about 0). It gives a max Est.Error of 0.1492 in-sample and 0.0709
for new rows, so the package is 141 and 113 times too high. The gradient
of the reference objective at the package's estimates is -1.2e-6 and
-1.3e-7 in sd units, so both use the same optimum. The finite value is
wrong, not merely large. `fitted()` gives no warning of its own.
Recommendation: raise the lane's MEDIUM to HIGH (user decision).

## 5. Tests (lane library, `NOT_CRAN=true`, `dev/release/run-tests.R`)

```
lib: frmtmb 0.83.0 C:/Users/adf44/source/r/wt-sparsese-lib
RESULT test-se-reparam-reml.R pass=193 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-fitted.R pass=129 fail=0 err=0 skip=0 warn=0
RESULT test-predict-lp-basis.R pass=30 fail=0 err=0 skip=0 warn=0
```

## Correction to the review (2026-11-17-sparsese.md, section 1)

`sparsese-review-indep.R` line 100 holds `names(of) %in% names(lost)`. But
`se_lost` names the parameter `theta_1`, while `opt$par` names it `theta`.
So the reference never held the lost sd: it kept a theta variance near
1.16e8. The b-row differences at 4500 s2 and s3 (7.45e-7, 3.38e-7 fitted
RE; 2.87e-6, 1.93e-6 in-sample) do not come from "the choice of the held
direction". They come from this mismatch. With the index stripped, those
rows agree to 8.0e-9 and 8.7e-9, and the b variances agree exactly. A held
direction choice is real only for the lost correlation (section 1, slope
seed 2).

## Pre-existing, beside the lane (file if not filed)

- LOW. `vcov(full = TRUE)` of a REML fit prints finite values for a lost
  parameter, where the message says NaN in `vcov()` and the ML fit gives
  NaN. Examples:
  - boundary sd: variance 1.157e8 (`vcovfull-base.txt`);
  - boundary correlation: variance -5.85e5, a negative variance
    (`vcovneg-base.txt`).

  `rellib-r22` and the lane give identical output. Its theta diagonal
  (0.0183, 0.0178) is also not the `cov_fixed_prop` that the summary and
  predictions use (0.242, 0.046).
- LOW, extends P2 to cumulative fits. The outer covariance `vcov(full)`
  is 7.8e-4 and 5.0e-4 (scaled) from the reference's Richardson Hessian
  (s30 seeds 1 and 3), and 2.3e-5 at 4508. A Richardson Jacobian of the
  tape's AD gradient is closer to the reference (6.4e-5, 9.0e-5) than to
  the package (8.1e-4, 4.8e-4). Effect on Est.Error: up to 1.4e-4
  relative, the same on both routes.

## Files

New in this check: `dev/reviews/sparsese-short.md`,
`dev/reviews/sparsese-short-{psd,probe,cum,f1int,slope,lostb,qrsolve,vcovfull,vcovneg}.R`,
`dev/reviews/sparsese-short-out/`.
