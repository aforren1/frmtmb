# Lane acprior: autocorrelation priors as brms's

Base 70476dc3 (frmtmb 0.75.0), worktree `frmtmb-wt-acprior`, branch
`wt-acprior`, library `C:/Users/adf44/source/r/wt-acprior-lib`. "Before"
is rellib-r13 throughout. Nothing is committed.

The two items of `dev/test-backlog.md`, "Filed at the 0.75.0 release",
"Default priors and pointwise likelihoods of residual terms", first two
entries.

## What brms 2.23.0 does (run, not recalled)

`dev/acprior-brms-defaults.R` (log `dev/acprior-out/brms-defaults.txt`)
and `dev/acprior-brms-bounds.R` (`brms-bounds.txt`), `default_prior()`
and `stancode()` only, seed 1, 4 groups of 10.

| term | class row | lb, ub | Stan declaration | default |
|---|---|---|---|---|
| `ar(p = 1)`, `ar(p = 2)` (`cov = FALSE`) | `ar` | none | `vector[Kar] ar;` | flat |
| `ma()`, `arma()` (`cov = FALSE`) | `ar`, `ma` | none | `vector[Kma] ma;` | flat |
| `ar(cov = TRUE)` | `ar` | -1, 1 | `vector<lower=-1,upper=1>[Kar] ar;` | flat |
| `ma(cov = TRUE)` | `ma` | -1, 1 | `vector<lower=-1,upper=1>[Kma] ma;` | flat |
| `arma(cov = TRUE)` | `ar`, `ma` | -1, 1 each | both bounded | flat |
| `cosy()` | `cosy` | 0, 1 | `real<lower=0,upper=1> cosy;` | flat |
| `unstr()` | `cortime` | none | `cholesky_factor_corr` | `lkj(1)` |

- The bound on `ar` is a DECLARED bound of the Stan parameter under
  `cov = TRUE`, not a stationarity constraint; under `cov = FALSE` there
  is none at all. brms refuses `cov = TRUE` above order one
  ("Covariance formulation of ARMA structures is only possible for
  effects of maximal order one").
- A user prior is a plain `lprior += normal_lpdf(ar | 0.5, 0.2)` on the
  coefficient, inside the declared bound. The bound stays under a user
  density (`cosy_beta`: still `<lower=0,upper=1>`). A user `lb`/`ub`
  replaces that side only: `lb = -0.1` on cosy gives
  `<lower=-0.1,upper=1>`, `ub = 0.8` gives `<lower=0,upper=0.8>`, and a
  bounds-only `set_prior("", class = "cosy", lb = 0.2)` gives
  `<lower=0.2,upper=1>` with no density, i.e. flat there. brms even
  takes `lb = -2` on a `cov = TRUE` ar (nonstationary); frmtmb refuses
  past the end.
- brms has no `class = "theta"`: `validate_prior(set_prior(..., class =
  "theta"), ar model)` stops "do not correspond to any model parameter".

## A. Default priors flat on the natural parameter

### The defect, measured

frmtmb.sample's `default_priors_for()` (extensions/frmtmb.sample/R/
sample.R) gave `ar`, `ma` and `cosy` no row, so the sampled density was
flat on frmtmb's INTERNAL parameter. For a first-order `cov = TRUE`
coefficient `rho = th / sqrt(1 + th^2)`, which is
`(1 - rho^2)^(-3/2)` on `rho`, with its mass at the ends; for `cosy` the
logistic onto `(-1/(d - 1), 1)`, which is `1 / ((rho + a)(1 - rho))`.
Under `cov = FALSE` the coefficients ARE the internal parameters
(identity map), so flat there was already brms's flat; nothing changes.
`unstr()` already had `lkj(1)` on `cortime` (since 0.72.0), unchanged.

`dev/acprior-agree.R`, seed 20261027, 4 chains of 4000 on both sides,
log `dev/acprior-out/agree.txt` (pasted verbatim):

    ===== ar
    ML estimate 0.9265
    frm_sample before          mean 0.9455 sd 0.0306 q05 0.8876 q50 0.9510 q95 0.9843  ESS 2539  mcse 0.0006
    frm_sample after           mean 0.9226 sd 0.0354 q05 0.8565 q50 0.9279 q95 0.9711  ESS 2875  mcse 0.0007
    brms                       mean 0.9156 sd 0.0371 q05 0.8466 q50 0.9209 q95 0.9658  ESS 6108  mcse 0.0005
    brms, weighted to frmtmb   mean 0.9214 sd 0.0355 q05 0.8554 q50 0.9265 q95 0.9691  ESS 6108  mcse 0.0005
    z(after - brms ref) 1.43; z(before - brms ref) 31.79

    ===== cosy
    ML estimate -0.0089
    frm_sample before          mean 0.0106 sd 0.0788 q05 -0.1051 q50 0.0035 q95 0.1517  ESS 6790  mcse 0.0010
    frm_sample after           mean 0.0756 sd 0.0601 q05 0.0054 q50 0.0616 q95 0.1950  ESS 3842  mcse 0.0010
    brms                       mean 0.0760 sd 0.0609 q05 0.0051 q50 0.0622 q95 0.1943  ESS 5230  mcse 0.0008
    z(after - brms ref) -0.35; z(before - brms ref) -51.39
    share of draws below 0: before 0.4844 after 0.0000 brms 0.0000

Data: ar, AR(1) with phi = 0.9 in 8 groups of 6 time points; cosy,
correlation 0.03 in 30 groups of 4. The ar comparison needs the
importance weights: brms's `cov = TRUE` ar() is on the innovation sd
`sigma_i`, frmtmb's on the marginal `sigma_m = sigma_i / sqrt(1 -
phi^2)` (documented divergence, `?frmtmb-autocor`), so the same
`student_t(3, 0, 2.5)` default on sigma is a different joint prior. The
weight `w = pt(sigma_i / s) / (pt(sigma_i) s)`, `s = sqrt(1 - phi^2)`,
is the exact density ratio after the change of variables, and the
unweighted brms row is shown beside it. This is an identity, not a fit.
cosy uses the marginal sigma on both sides; no weights.

A user prior, `dev/acprior-agree-up.R` (same data, seed, chains; log
`dev/acprior-out/agree-up.txt`, verbatim):

    ===== ar_normal
               mean     sd    q05    q50    q95    min      ess  mcse
    before   0.9081 0.0411 0.8329 0.9138 0.9648 0.6252 2920.426 8e-04
    after    0.9081 0.0411 0.8329 0.9138 0.9648 0.6252 2920.426 8e-04
    brms_ref 0.9072 0.0408 0.8314 0.9132 0.9631 0.7086 6198.511 5e-04
    z(after - brms) 1.04; z(before - brms) 1.04

    ===== ar_lb
               mean     sd    q05    q50    q95    min      ess  mcse
    before   0.9464 0.0294 0.8903 0.9516 0.9841 0.7702 3205.924 5e-04
    after    0.9219 0.0355 0.8562 0.9271 0.9701 0.6713 3491.360 6e-04
    brms_ref 0.9217 0.0364 0.8539 0.9274 0.9702 0.6757 6518.333 5e-04
    z(after - brms) 0.32; z(before - brms) 35.88

    ===== cosy_norm
               mean     sd     q05    q50    q95     min      ess  mcse
    before   0.0087 0.0610 -0.0859 0.0060 0.1152 -0.1759 7016.176 7e-04
    after    0.0540 0.0413  0.0039 0.0455 0.1336  0.0000 4362.368 6e-04
    brms_ref 0.0525 0.0414  0.0038 0.0436 0.1339  0.0000 4834.483 6e-04
    z(after - brms) 1.74; z(before - brms) -46.51

- `ar_normal`, `set_prior("normal(0.5, 0.2)", class = "ar")`: core has
  applied the user density on the NATURAL coefficient with the map's
  log Jacobian since 0.49 (`trans_dist()`, `prior_logdens()`), which is
  brms's `normal_lpdf(ar | 0.5, 0.2)` on the coefficient; before and
  after are the same draws, and both agree with brms (z 1.04). Not a
  defect.
- `ar_lb`, a bounds-only `set_prior("", class = "ar", lb = 0.5)`: brms
  is flat on [0.5, 1]. frmtmb dropped the default slot for the user's
  (`drop_superseded()`) and so fell back to flat on the internal
  parameter even after the new default: fixed by `keep_flat_defaults()`.
- `cosy_norm`, `set_prior("normal(0, 0.1)", class = "cosy")`: brms keeps
  its declared [0, 1]; frmtmb sampled the whole window. Fixed by the
  same function (the default's `lb = 0` is kept on each side the user
  leaves unset).

### ML fits are unaffected

`dev/acprior-ml.R` (seed 11, 30 groups of 5, AR(1) 0.5): `frm()` with
no prior and with a user prior on each class (ar, ar cov, ar(2) cov, ma
cov, arma cov, cosy, unstr), objective to 10 digits and every
`thetaac` to 8, before and after. `diff` of the two logs (minus the
library line): IDENTICAL (`dev/acprior-out/ml-before.txt`,
`ml-after.txt`). The fit route reads no default (`prior_table()` reads
the registry only on `route = "sample"`), and the user-prior densities
were already on the natural scale.

### What changed

- Core, `R/priors.R`: `autocor_flat_priors(fit)`, exported for the
  sampling extension on the `?frmtmb-sampling-api` page (one new line in
  the `@rawNamespace export()` list of `R/sampling-api.R`). One row per
  class and response held through a transform: the ends of a uniform
  over the whole image of the map (`(-1, 1)` at order one; cosy's
  window `(-1/(d - 1), 1)`; the registry term's `bounds`), and the
  narrower bound brms declares (`lb = 0` for cosy). For a `cov = TRUE`
  ar or ma of order `k > 1`, which brms does not have, the image is the
  stationary (invertible) region, which lies inside the box
  `|phi_j| <= choose(k, j)`, so a uniform on
  `(-choose(k, k %/% 2), choose(k, k %/% 2))` is flat on it. No row for
  `cov = FALSE` and for cortime.
- frmtmb.sample, `default_priors_for()`: one loop over
  `autocor_flat_priors()` replaces the sar-only loop over
  `frmtmb_autocor_pars()`, so core's classes and a registered term go
  through one rule. `keep_flat_defaults()` (new, beside
  `drop_superseded()`): a user spec of the same slot as a flat default
  inherits its `lb`/`ub` on each unset side, and a bounds-only spec
  inherits the flat density. The announcement prints the bound
  (`cosy uniform(-0.25, 1), lb = 0`) and `[flat on the stationary or
  invertible region]` for the order-k case.
- Why a uniform over the image and not brms's narrower interval as the
  density: a bounds-only user prior that widens cosy below 0
  (`lb = -0.1`, which brms takes) must stay flat there; a uniform(0, 1)
  density would be NaN on (-0.1, 0).

### Choice that needs the user: cosy on [0, 1]

brms 2.23.0 declares `cosy` on [0, 1]; frmtmb's `frm()` fits the whole
positive-definite window `(-1/(d - 1), 1)`. I followed brms on the
sampling route only (the project tiebreaker: brms is restrictive here,
not wrong): `frm_sample()` bounds cosy at 0 by default and under a user
cosy prior with no `lb`, and `default_prior(route = "sample")` lists
0 and 1. `frm()` is unchanged and its table lists the window. So an
ML fit with a negative cosy (the data above: -0.0089) starts its chains
at the clamped value, with the existing warning "The ML mode violates
the requested bound(s) on thetaac_1". A user who wants a negative cosy
in `frm_sample()` writes `lb`, as in brms. If the user prefers
frmtmb's window on both routes, drop the `lb` from
`autocor_flat_priors()` (the `bd` it reads from
`autocor_class_bounds(ac, "cosy", "sample")`).

## B. One prior slot per autocorrelation parameter

### The defect, measured

`dev/acprior-frm-table.R before` (`frm-table-before.txt`): every
autocorrelation model listed the class row AND `theta thetaac_k` rows
(unstr with 10 time points: `cortime` plus 45 theta rows).
`dev/acprior-double.R before` (`double-before.txt`, seed 3): the fit did
not apply both priors, as the release record said; it applied the ONE
WRITTEN LAST, silently. `ar normal(0.5, 0.2)` then `theta normal(0, 1)`
kept only the theta entry (`entry thetaac 1 normal`, thetaac 0.62568);
the reverse order kept the ar entry (`trans normal`, 0.60349, the same
as ar alone). Worst case: `unstr` with `cortime lkj(2)` plus `theta
thetaac_1 normal(0, 1)` kept only the one-position normal; the joint
LKJ on all 45 correlations was retired by `claim()`.

### What changed

- `autocor_class_coefs(frame)` (R/priors.R) lists every `thetaac`
  position that a class reaches, with the class and response:
  `ar`, `ma`, `cosy`, `cortime` through `autocor_class_idx()`, and a
  registry term's bounded parameter through `autocor_seam_bounded()`.
  It replaces `autocor_seam_bounded_coefs()` (R/autocor-seam.R,
  removed), which was the 0.75.0 special case for sar.
- `prior_table()` lists none of those internal rows, and the theta
  branch of `resolve_priorlist()` refuses them: "class = "theta", coef =
  "thetaac_1" is the internal parameter of ar, which has a class of its
  own, as in brms. Write the prior with class = "ar"; one parameter
  takes one prior", with `, resp = "y2"` in a multivariate model.
- `class = "theta"` keeps the parameters no class reaches: RE
  `theta_k`, `thetar_k`, and a registered term's parameter without an
  interval (tested with the `uscale` term in test-autocor-seam.R: the
  row is listed, and a `ub` on it holds).
- Consequence: the internal parameters of a `cov = TRUE` ar or ma of
  order >= 2 take no box at all now. The order-2 refusal used to point
  at `class = "theta", coef = "thetaac_1"`; it now says brms's cov =
  TRUE form is of order one only. The `cortime` lb/ub refusal no longer
  names thetaac_1. The `class = "theta"` with no coef refusal on an
  autocorrelation-only model names the classes.

### The prior table's lb and ub, and the round trip

Rows of `default_prior()` now carry brms's bounds:
`frm-table-after.txt` lists `ar -1 1` for `ar(cov = TRUE)`, `ma -1 1`,
both for `arma(cov = TRUE)`, nothing for `cov = FALSE` and for order 2,
and `cosy -0.1111111 1` (10 time points) on the fit route; on
`route = "sample"` cosy is 0 and 1 (test-default-priors-autocor.R).
Row ORDER equals brms's in all 9 cases (`dev/acprior-order.R`, log
`dev/acprior-out/order-vs-brms.txt`: `same=TRUE` for ar1, ar2, ma1,
arma, ar1cov, ma1cov, armacov, cosy and unstr), and the class and coef
columns of the before and after tables are identical once the trailing
theta rows are taken out (a `diff` of the two logs).

Defect found on the way, fixed: the `sar()` table did not pass back as
`prior =` at 0.75.0 (`dev/acprior-roundtrip.R before`,
`roundtrip-before.txt`): "class = "lagsar" lb = -0.999999999999999 is
outside the interval of sar(W), (-1, 1)", because the table lists the
ends of the interval and `ac_bound_theta()` refused an end. A bound AT
an end now bounds nothing (`-Inf`/`Inf`); past the end it is still
refused. After: sar, ar cov and cosy tables round-trip with the logLik
unchanged (`roundtrip-after.txt`: -36.5520620288 both, -36.2657108716
both, -36.7308045125 both). This flips one assertion of
test-autocor-seam.R that pinned "an end of the interval is not inside
it, as for ar" (ub = 2 on bscale(V)); its reason was symmetry with ar,
and ar now accepts the end too.

## Tests

Seen to fail on rellib-r13 first (one file per process, logs in
`dev/acprior-out/`):

    t-pac-before.txt:RESULT test-priors-autocor-classes.R pass=61 fail=11 err=2 skip=0 warn=0
    t-pac-after.txt:RESULT test-priors-autocor-classes.R pass=151 fail=0 err=0 skip=0 warn=0
    t-seam-before.txt:RESULT test-autocor-seam.R pass=69 fail=0 err=1 skip=0 warn=0
    t-seam-after.txt:RESULT test-autocor-seam.R pass=80 fail=0 err=0 skip=0 warn=0
    t-sar-before.txt:RESULT test-sar.R pass=80 fail=1 err=1 skip=3 warn=0
    t-sar-after.txt:RESULT test-sar.R pass=82 fail=0 err=0 skip=3 warn=0
    t-dpa-before.txt:RESULT test-default-priors-autocor.R pass=15 fail=13 err=0 skip=0 warn=0
    t-dpa-after.txt:RESULT test-default-priors-autocor.R pass=28 fail=0 err=0 skip=0 warn=0
    g-dpf-before.txt:RESULT test-default-priors-fit.R pass=21 fail=3 err=0 skip=0 warn=0
    g-dpf-after.txt:RESULT test-default-priors-fit.R pass=24 fail=0 err=0 skip=0 warn=0

The behavioral failures before (not missing symbols): the theta rows
listed and the theta prior accepted; the table's lb/ub NA; the end bound
refused; the flat-spread of the sampled density 0.677 (ar), 0.677
(ma), 0.887 (ar(2)), 0.752 (arma), 0.652 (cosy), 0.906 (bounds-only
ar), 0.421 (bounds-only registered term) against < 1e-10 after; cosy's
lb absent; and in the gated file the ar mean z 16.1, cosy draws below 0
and a cosy width ratio 5.8. The "guard absent" cases pass on both: the
class itself takes the prior, class theta bounds a parameter no class
reaches, past-the-end bounds are refused.

New or changed tests:
- tests/testthat/test-priors-autocor-classes.R: three tests rewritten
  (the cortime lb/ub message, the order-2 bound, class theta's reach),
  four new (one slot per parameter for ar, ar cov, ma, arma, cosy,
  unstr and a multivariate model; the table's bounds and round trip;
  the end bound; get_prior has no theta rows).
- tests/testthat/test-autocor-seam.R: the end bound flipped, the table
  round trip, the bounded term's internal spelling refused by name, the
  unbounded term's still listed and bounding.
- extensions/frmtmb.spatial/tests/testthat/test-sar.R: the ar assertion
  flipped (no theta row), the sar table round trip.
- extensions/frmtmb.sample/tests/testthat/test-default-priors-autocor.R
  (new, ungated): flat on the natural parameter for ar, ma, ar(2),
  arma (cov = TRUE), none under cov = FALSE, cosy's lb kept under a
  user prior, the bounds-only density, the sampling table and message,
  a registered term through the same rule.
- extensions/frmtmb.sample/tests/testthat/test-default-priors-fit.R
  (gated, FRMTMB_BRMS_FIT_TESTS): ar near 1 and cosy near 0 against
  brms's chains, 2 chains of 2000, seed 3.

### Every suite, one file per process (dev/release/run-par.sh, N = 8)

Ungated (NOT_CRAN=true), from `dev/acprior-sum.sh acprior-u`:

    frmtmb           files=279 pass=19857 fail=0 err=0 skip=204 warn=0
    frmtmb.coupling  files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam       files=32 pass=1838 fail=0 err=0 skip=5 warn=0
    frmtmb.latent    files=12 pass=681 fail=0 err=0 skip=2 warn=0
    frmtmb.learn     files=15 pass=430 fail=0 err=0 skip=13 warn=0
    frmtmb.ode       files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample    files=68 pass=3081 fail=0 err=0 skip=11 warn=0
    frmtmb.spatial   files=8 pass=454 fail=0 err=0 skip=5 warn=0
    frmtmb.survival  files=15 pass=552 fail=0 err=0 skip=0 warn=0
    logs whose frmtmb is not the lane lib: 0
    logs without a RESULT line: 0

Gated (NOT_CRAN=true, FRMTMB_BRMS_FIT_TESTS=true),
`dev/acprior-sum.sh acprior-g`:

    frmtmb           files=279 pass=22293 fail=0 err=0 skip=15 warn=0
    frmtmb.coupling  files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam       files=32 pass=1844 fail=0 err=0 skip=3 warn=0
    frmtmb.latent    files=12 pass=681 fail=0 err=0 skip=2 warn=0
    frmtmb.learn     files=15 pass=501 fail=0 err=0 skip=2 warn=0
    frmtmb.ode       files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample    files=68 pass=3363 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial   files=8 pass=529 fail=0 err=0 skip=0 warn=0
    frmtmb.survival  files=15 pass=552 fail=0 err=0 skip=0 warn=0
    logs whose frmtmb is not the lane lib: 0
    logs without a RESULT line: 0

The gated tier's 29 skips are other tiers' gates: 15 FRMTMB_SCALE_TESTS,
13 FRMTMB_DRMTMB_FIT_TESTS, 1 FRMTMB_FUZZ; none is a brms gate. The
ungated tier ran test-default-priors-autocor.R before its last three
assertions (the user density kept, the absent case) were added; the
gated tier and the one-file runs above have them.

### OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`)

The four changed ungated files, one per process, under the OpenBLAS
copy of R (`La_version()` 3.12.0, a 2000 x 2000 dgemm in 0.13 s; logs
`dev/acprior-out/ob-*.txt`): test-priors-autocor-classes.R pass=151,
test-autocor-seam.R pass=80, test-sar.R pass=82 skip=3 (brms gates),
test-default-priors-autocor.R pass=25 (before its last three
assertions); all fail=0 err=0. No assertion sits on a rounding-decided
branch (the flat-spread tests compare a ratio of order 1 before against
1e-10; the round trips compare two fits in one process), so 0.3.32 was
not run. `dev/ciharden-out/` is deleted.

### Every suite in one session, three orders (dev/release/run-sessions.sh)

    SESSION frmtmb forward files=269/279 pass=19857 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=269/279 pass=20007 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=269/279 pass=20100 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam forward files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.eam preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.latent forward files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent reverse files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.latent preload files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
    SESSION frmtmb.learn forward files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.learn reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.learn preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
    SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
    SESSION frmtmb.sample forward files=63/68 pass=3081 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample reverse files=63/68 pass=3081 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample preload files=63/68 pass=3081 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
    SESSION frmtmb.spatial forward files=8/8 pass=454 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.spatial reverse files=8/8 pass=454 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.spatial preload files=8/8 pass=454 fail=0 err=0 skip=5 warn=0
    ACPRIOR-S ran 27 of 27
    untagged LEAK lines: 0

Every LEAK line is tagged "[with load of ...]" (an extension or a
dependency loaded by a file), as at 0.74.0. frmtmb.sample's three
sessions were run again after the last three assertions of
test-default-priors-autocor.R were added (tier acprior-s2):

    SESSION frmtmb.sample forward files=63/68 pass=3084 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample reverse files=63/68 pass=3084 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample preload files=63/68 pass=3084 fail=0 err=0 skip=11 warn=0
    untagged LEAK lines: 0

The tests that register (`test-autocor-seam.R`'s terms,
`test-default-priors-autocor.R`'s `acpterm`) call
`local_frmtmb_registries()` / `frmtmb_local_registries()` first.

## Ported brms-suite rows

Read with `dev/brmsport-verdicts.tsv` and `dev/brmsport-bin2-ledger.tsv`
(not edited). The verdict rows that read these priors are `priors:92`
and `priors:94` (`class == "ar"` by resp, on the setup
`bf(y1 ~ x + (x|ID1|g) + ar()) + bf(y2 ~ 1)`; that setup and its
`default_prior()` call are the record log's setup lines 89 and 90, not
verdict rows; corrected in fix round 1, item 8). They pass before and
after: the `ar` row per response is unchanged and `cov = FALSE` gets no
bound. The gated test-brms-suite-priors.R passes on the lane library
(round 0: `dev/release/acprior-g.log`, pass=36 fail=0). None flips.
`priors:5` and `priors:123`
are divergences for frmtmb's class theta rows of RANDOM-EFFECT
parameters, which this lane leaves alone (see below); they do not
flip. The bin-2 ledger has no row that reads an autocorrelation prior;
its residual-term rows (`log_lik:130/133/150-176/230`,
`posterior_predict:93-135`) read log_lik and posterior_predict, the
loglik lane's.

## Not done, and why

- **Backlog entry: the same double slot for random-effect parameters
  and `rescor`.** Class `sd` or `cor` and `class = "theta", coef =
  "theta_k"`, and `class = "rescor"` and `theta thetar_1`, are two
  spellings of one parameter. Measured by the reviewer
  (`dev/review-acprior-re.R`, `re-before.txt` equals `re-after.txt`
  past the library line; block `(1 + x + z | g)`):
  - `frm()`: `cor lkj(2)` then `theta theta_4 normal(0, 1)` leaves only
    the theta entry on position 4, and the LKJ on theta_5 and theta_6
    is gone; in the other order the LKJ wins. The same for `rescor
    lkj(2)` and `thetar_1`.
  - `frm_sample()`: one `theta theta_4` prior retires the default
    `lkj(1)` on all three correlations, so theta_5 and theta_6 fall
    back to flat on the internal parameter (a changed sampled density).
    `sd` with `coef = "x"` plus `theta theta_2` keeps only the theta
    spec, and the sd spec is dropped with no message.
  - Cost: a filter in `prior_table()` and a refusal map in the theta
    branch of `resolve_priorlist()` from RE positions to `sd`/`coef`
    and `cor`/`group`, about 40 to 60 lines in R/priors.R, plus tests.
    The RE theta vector also holds registered covstruct parameters that
    only class theta reaches; those must keep the spelling.
  - It removes the box on a single RE correlation (class theta on one
    position), which brms does not have: a user decision.
  - It flips the ported divergences `priors:5` and `priors:123`, which
    count frmtmb's extra class theta rows.
  This is `claim()`'s documented "later wins across spellings", not a
  new defect. The decision of fix round 1 is to file it.
- The legacy named-list spelling (`prior = list(thetaac_1 =
  prior_normal())`) addresses internal parameters by design and is not
  a prior-table slot; unchanged.
- (Withdrawn in fix round 1, item 3: the bounded natural defaults
  sigma, nu and ndt now keep their bound under a user prior; see below.)
- The sampling table prints `uniform(-1, 1)` where brms prints
  `(flat)` with lb/ub (as sar's row did at 0.75.0).
- README's sentence that the prior vocabulary reaches every parameter
  is narrowed to "every one that brms reaches".

## Overlaps with the other lanes

- loglik: `R/autocor.R` (I changed one roxygen paragraph only),
  `R/autocor-seam.R` (removed `autocor_seam_bounded_coefs()`, changed
  the `frmtmb_autocor_pars()` docs), frmtmb.sample `R/sample.R`, NEWS.
  If loglik reads `autocor_seam_bounded_coefs()` it must switch to
  `autocor_class_coefs()`.
- mocode: `R/priors.R` (prior table; my edits are the autocor loop, the
  theta rows, the theta branch, the order-2 and lkj messages, and new
  functions after `autocor_class_idx()`), frmtmb.sample `R/sample.R`
  (`default_priors_for()`; my block replaced the sar loop between the
  lncor loop and the unstr loop), NEWS of core and frmtmb.sample.
- `R/sampling-api.R`: one alias and one export line, each on its own
  line.
- `R/compat.R`: the text of two `group:autocor` rows.

## R CMD check --as-cran

Built and checked in `dev/acprior-check/<pkg>` (deleted afterwards; the
logs are `dev/acprior-out/check-<pkg>.log` and `build-<pkg>.log`), with
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH and
`R_LIBS` = lane lib, rellib-r13, user library (`dev/acprior-check.sh`):

    check-frmtmb.log:Status: 1 NOTE
    check-frmtmb.sample.log:Status: OK
    check-frmtmb.spatial.log:Status: 1 NOTE

Both NOTEs are the environmental "Skipping checking math rendering:
package 'V8' unavailable" on the HTML manual. Core's tests ran 17 min
in the check, OK.

## Bump and floor

- frmtmb: minor. Behavior changes (prior table rows removed and lb/ub
  added, a refusal where a prior was taken, a new export
  `autocor_flat_priors()` on the sampling-API page).
- frmtmb.sample: minor (the default prior changes, BREAKING). Floor:
  the new frmtmb (it calls `autocor_flat_priors()`; on 0.75.0 its
  `default_priors_for()` would fail to find it).
- frmtmb.spatial: patch (tests and NEWS only, no R code). Floor: the
  new frmtmb, because its test-sar.R asserts the round trip and the
  missing theta row that only the new core gives.
- Other extensions: unchanged, and their suites pass on the new core.

## Fix round 1 (review `dev/reviews/2026-10-27-acprior.md`)

Numbered as the review. Everything is uncommitted. Round 0's sections
above stay as written; where this round corrects one of them, it says
so here.

### 1. The sigma default on brms's innovation sd (major)

Round 0 compared frm_sample with brms's draws reweighted to frmtmb's
marginal-sd parameterization, and the NEWS called that "brms's". Against
brms's own draws the round-0 build still differed: ar mean 0.9226
against 0.9156, z 9.37 (reviewer, `dev/review-acprior-unweighted.R`).
The remaining gap was the default on sigma.

- brms, read and run: `cholesky_cor_ar1()`, `cholesky_cor_ma1()` and
  `cholesky_cor_arma1()` (brms/chunks/fun_cholesky_cor_*.stan) return
  covariances with diagonal `1 / (1 - ar^2)`, `1 + ma^2` and
  `(1 + ma^2 + 2 ar ma) / (1 - ar^2)`. `dev/acprior-armacheck.R`
  (`armacheck.txt`, 3 random points per structure): frmtmb's
  correlation equals brms's covariance over its diagonal to at most
  2.2e-16, for ar, ma and arma. So brms's sigma is frmtmb's
  `sigma * c`, `c = 1 / sqrt(diagonal)`.
- The change is small (about 60 lines with comments):
  - core `R/priors.R`: `autocor_innov(ac)` and `innov_logc()`; an entry
    of a natural-scale spec with `innov_sd = TRUE` carries `innov`, and
    `neg_log_prior_fn()` hands `prior_logdens()` the log factor `lc`,
    which reads the density at `x * exp(lc)` and adds `lc`.
    `autocor_innov_sd(fit)` (new export) names the responses.
  - frmtmb.sample `default_priors_for()`: the sigma default of such a
    response carries `innov_sd`; the message says "[on brms's
    innovation sd]". Only the default: a user prior on sigma stays on
    frmtmb's marginal sd, the sigma that `frm()` reports.
  - Not covered, filed: `sigma ~ x` under a cov = TRUE term (brms's
    default is then on the log-scale intercept of the innovation sd),
    and orders above one (brms has none).
- Measured against brms's UNWEIGHTED posterior, `dev/acprior-innov.R`
  (seed 3, 4 chains of 4000 on each side; ar data of
  dev/acprior-agree.R; ma: MA(1) 0.6, 10 groups of 6; arma: 0.7 and
  0.3, 12 groups of 6; seed 20261027). Log `dev/acprior-out/innov.txt`,
  verbatim:

    
    #### frm_sample before against brms's own draws
    
    ===== ar 
    ar[1]        frm mean 0.9455 sd 0.0306 | brms mean 0.9156 sd 0.0371 | z  42.26 | sd ratio 0.825
    sigma        frm mean 0.4165 sd 0.0482 | brms mean 0.4324 sd 0.0511 | z -19.55 | sd ratio 0.942
    b_Intercept  frm mean 0.8806 sd 0.4914 | brms mean 0.8935 sd 0.3798 | z  -1.43 | sd ratio 1.294
    b_x          frm mean 0.4885 sd 0.0502 | brms mean 0.4877 sd 0.0518 | z   0.73 | sd ratio 0.968
    
    ===== ma 
    ma[1]        frm mean 1.0000 sd 0.0000 | brms mean 0.6411 sd 0.1696 | z     NA | sd ratio 0.000
    sigma        frm mean 0.7873 sd 0.0741 | brms mean 0.8692 sd 0.0911 | z -53.30 | sd ratio 0.814
    b_Intercept  frm mean 1.1784 sd 0.1935 | brms mean 1.1737 sd 0.1716 | z   1.56 | sd ratio 1.128
    b_x          frm mean 0.7095 sd 0.1016 | brms mean 0.7227 sd 0.1057 | z  -7.51 | sd ratio 0.961
    
    ===== arma 
    ar[1]        frm mean 0.6748 sd 0.0927 | brms mean 0.6811 sd 0.0926 | z  -3.59 | sd ratio 1.001
    ma[1]        frm mean 1.0000 sd 0.0000 | brms mean 0.6666 sd 0.1797 | z     NA | sd ratio 0.000
    sigma        frm mean 0.8297 sd 0.0741 | brms mean 0.9335 sd 0.0926 | z -70.14 | sd ratio 0.800
    b_Intercept  frm mean 0.9634 sd 0.4586 | brms mean 0.9758 sd 0.4393 | z  -1.50 | sd ratio 1.044
    b_x          frm mean 0.5412 sd 0.0811 | brms mean 0.5589 sd 0.0900 | z -12.09 | sd ratio 0.901
    
    #### frm_sample after against brms's own draws
    
    ===== ar 
    ar[1]        frm mean 0.9165 sd 0.0363 | brms mean 0.9156 sd 0.0371 | z   1.17 | sd ratio 0.979
    sigma        frm mean 0.4315 sd 0.0494 | brms mean 0.4324 sd 0.0511 | z  -1.12 | sd ratio 0.966
    b_Intercept  frm mean 0.8910 sd 0.3792 | brms mean 0.8935 sd 0.3798 | z  -0.34 | sd ratio 0.998
    b_x          frm mean 0.4890 sd 0.0514 | brms mean 0.4877 sd 0.0518 | z   1.31 | sd ratio 0.992
    
    ===== ma 
    ma[1]        frm mean 0.6471 sd 0.1756 | brms mean 0.6411 sd 0.1696 | z   1.19 | sd ratio 1.035
    sigma        frm mean 0.8665 sd 0.0892 | brms mean 0.8692 sd 0.0911 | z  -1.51 | sd ratio 0.979
    b_Intercept  frm mean 1.1700 sd 0.1726 | brms mean 1.1737 sd 0.1716 | z  -1.18 | sd ratio 1.006
    b_x          frm mean 0.7205 sd 0.1054 | brms mean 0.7227 sd 0.1057 | z  -1.14 | sd ratio 0.997
    
    ===== arma 
    ar[1]        frm mean 0.6801 sd 0.0941 | brms mean 0.6811 sd 0.0926 | z  -0.50 | sd ratio 1.016
    ma[1]        frm mean 0.6672 sd 0.1789 | brms mean 0.6666 sd 0.1797 | z   0.13 | sd ratio 0.996
    sigma        frm mean 0.9322 sd 0.0939 | brms mean 0.9335 sd 0.0926 | z  -0.73 | sd ratio 1.014
    b_Intercept  frm mean 0.9767 sd 0.4330 | brms mean 0.9758 sd 0.4393 | z   0.11 | sd ratio 0.986
    b_x          frm mean 0.5598 sd 0.0898 | brms mean 0.5589 sd 0.0900 | z   0.55 | sd ratio 0.997

  frmtmb's sigma draws are carried to brms's innovation sd for the
  comparison. Before (rellib-r13) the ar mean was 42.3 MC SE from
  brms's and every `ma[1]` draw sat at 1.0000 (flat on the internal
  parameter); after, every parameter of the three cases is within
  |z| 1.51 and every sd ratio within 0.966 to 1.035. brms had 0
  divergent transitions in all three; frm_sample after had none
  (before: 46, and R-hat 2.13, on ma).
- The user-prior comparisons of round 0 were weighted too. Against
  brms's own draws (`dev/acprior-agree-up.R` with `summary unweighted`,
  `agree-up-unw.txt`): `ar_normal` z 8.32 before, 1.13 after; `ar_lb`
  43.81 before, 0.48 after; `cosy_norm` -46.51 before, 1.54 after.
- frmtmb.sample's NEWS is rewritten to state only these unweighted
  comparisons.

### 2. cosy's divergent transitions (major)

Round 0's own run had 27 divergent transitions (`agree-after.txt`, 4 x
4000), and the findings did not say so. The reviewer measured 5 to 9 per
4 x 2000 run; brms had 0.

The cause, by construction (`dev/acprior-cosydiv.R`, `cosydiv.txt`;
random starts, 4 x 2000, the same objective in three arms):

    seed 3: A 5 divergent, cosy mean 0.0761 | B 0, 0.0755 | C 6, 0.0765
    seed 4: A 6 divergent, cosy mean 0.0742 | B 0, 0.0770 | C 2, 0.0760
    seed 5: A 6 divergent, cosy mean 0.0751 | B 0, 0.0759 | C 3, 0.0758

A is the round-0 construction (the window's logistic onto (-a, 1), and
Stan's `lb + exp(u)` on the internal parameter). B samples `v` with
`cosy = plogis(v)` and hands the model `qlogis((cosy + a) / (1 + a))`:
one logistic onto brms's interval, as brms's own transform. C adds a
finite upper bound to A, so Stan's logit replaces its exp. Only B has
none, so the cause is the log bound stacked on the window's logistic.

The fix: core `autocor_bound_chart(fit, lower, upper)` (new export)
builds B's chart for every bounded first-order residual-correlation
parameter whose map is not the identity (cosy, a cov = TRUE ar or ma, a
registry term's bounded parameter). frmtmb.sample samples through it
(`chart_fn()`, `chart_pars()`, `chart_draws_back()`) and drops Stan's
bound there. A `cov = FALSE` coefficient keeps Stan's bound, which is
brms's transform already. The chart's ends are snapped to the bound as
written, since the round trip through the internal bound carries a few
ulps (5.6e-17 for lb = 0).

After (`dev/acprior-cosyfix.R`, `cosyfix-frm.txt`, `cosyfix-brms.txt`;
the reviewer's seeds and starts, 4 x 2000):

    init last.par.best seed 3: 0 divergent; cosy mean 0.0750 sd 0.0600 min 1.07e-06
    init random        seed 3: 0 divergent; cosy mean 0.0772 sd 0.0618 min 2.49e-05
    init last.par.best seed 4: 0 divergent; cosy mean 0.0775 sd 0.0623 min 1.38e-04
    init random        seed 4: 0 divergent; cosy mean 0.0783 sd 0.0606 min 4.25e-05
    init last.par.best seed 5: 0 divergent; cosy mean 0.0760 sd 0.0605 min 6.81e-05
    init random        seed 5: 0 divergent; cosy mean 0.0757 sd 0.0607 min 2.24e-07
    brms seed 3: 0 divergent; cosy mean 0.0754 sd 0.0606
    brms seed 4: 0 divergent; cosy mean 0.0778 sd 0.0601
    brms seed 5: 0 divergent; cosy mean 0.0784 sd 0.0619

There are 0 divergent transitions in all six runs, and the means are
0.0750 to 0.0783 against brms's 0.0754 to 0.0784 on the same seeds.

### 3. Bounded defaults keep their bound under a user prior (major)

`keep_flat_defaults()` is now `keep_default_bounds()`: every default
with a `lb` or `ub` lends it to a user specification of the same slot,
on each side that specification leaves unset. The flat density is
still lent only by a `flat_interval` default, to a bounds-only
specification. Tests in test-default-priors-bounds.R:
- nu on an identity link under `normal(2, 2)` keeps `lb = 1`, and a user
  `lb = 1.5` replaces it;
- sigma on an identity link under `normal(1, 1)` keeps `lb = 0`;
- ndt under `normal(0.2, 0.1)` keeps `ub = log(min_Y)`;
- the absent case: nu on its logm1 link gets no bound, and the user's
  normal is the density applied;
- sampled: under `normal(2, 2)` no nu draw is below 1 (2 chains of
  1000, seed 3; the reviewer measured 3.55 percent before,
  `dev/review-acprior-nu.R`).
The four behavior assertions fail on rellib-r13.

### 4. The warning on a mode below cosy's default bound (minor)

The bound is now the chart's, so Stan's "violates the requested
bound(s)" no longer fires there; `chart_mode_warning()` does. When the
bound equals the default's (`autocor_flat_priors()`'s lb), it reads:
"The ML estimate of cosy, -0.008901, is below 0, the lower bound brms
declares for cosy, which is frm_sample()'s default there; every chain
starts just inside it. To sample the whole window frm() fits, write
set_prior("", class = "cosy", lb = -1/3)". It writes `-1/(d - 1)` as a
fraction, which parses to the window's end exactly and so bounds
nothing. A bound the user wrote gets "is outside the requested
interval".

### 5 and 6. Docs and NEWS (minor)

`R/fit.R` and core NEWS name the second kind of parameter that lost its
class-theta box: a single `unstr()` correlation. Core NEWS has a New
features bullet for the three sampling-API exports,
`autocor_flat_priors()`, `autocor_bound_chart()` and
`autocor_innov_sd()`.

### 7. Tests (minor)

- The flat-spread tests use a central-difference log determinant of
  `autocor_natural()` (`dpa_fd_logjac()`), independent of the closed
  form the prior adds, against a 1e-6 ratio. The sigma test reads its
  reference density from `stats::dt()` at `exp(eta) * c`, with c from
  the natural coefficients.
- The gated ar test compares with brms's own draws, unweighted, for
  `ar[1]` and for sigma (frmtmb's carried to the innovation sd). The
  reviewer's direct-fit Stan edit is not needed: with the sigma default
  fixed, frm_sample targets brms's own model.

### 8. Ported-row citation (nit)

Fixed in the section "Ported brms-suite rows".

### Tests of this round, seen to fail on rellib-r13 first

    t2-dpa-before: RESULT test-default-priors-autocor.R pass=20 fail=17 err=3 skip=0 warn=0
    t2-dpa-after: RESULT test-default-priors-autocor.R pass=51 fail=0 err=0 skip=0 warn=0
    t2-dpb-before: RESULT test-default-priors-bounds.R pass=17 fail=4 err=0 skip=0 warn=0
    t2-dpb-after: RESULT test-default-priors-bounds.R pass=21 fail=0 err=0 skip=0 warn=0
    g2-dpf-before: RESULT test-default-priors-fit.R pass=24 fail=4 err=0 skip=0 warn=0
    g2-dpf-after: RESULT test-default-priors-fit.R pass=28 fail=0 err=0 skip=0 warn=0

The before failures are behavioral: the sigma density off by 0.21, 0.08
and 0.15 (ratio to its parts, ar, ma, arma) where after it is below
1e-10; the flat spreads 0.96, 0.96, 0.89, 0.75, 0.69, 0.91 and 0.60;
the bounds absent; a nu draw below 1; and, gated, the ar mean z 19.7
and sigma z 10.9 against brms's own draws, cosy draws below 0 and a
width ratio of 5.8. The missing `autocor_bound_chart()` is the weak
form for three tests; their behavioral counterpart is the divergence
measurement of item 2.

### Every suite again, one file per process (tiers acprior-u2, acprior-g2)

Ungated (NOT_CRAN=true), `dev/acprior-sum.sh acprior-u2`:

    frmtmb           files=279 pass=19857 fail=0 err=0 skip=204 warn=0
    frmtmb.coupling  files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam       files=32 pass=1838 fail=0 err=0 skip=5 warn=0
    frmtmb.latent    files=12 pass=681 fail=0 err=0 skip=2 warn=0
    frmtmb.learn     files=15 pass=430 fail=0 err=0 skip=13 warn=0
    frmtmb.ode       files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample    files=68 pass=3114 fail=0 err=0 skip=11 warn=0
    frmtmb.spatial   files=8 pass=454 fail=0 err=0 skip=5 warn=0
    frmtmb.survival  files=15 pass=552 fail=0 err=0 skip=0 warn=0
    logs whose frmtmb is not the lane lib: 0
    logs without a RESULT line: 0

Gated (FRMTMB_BRMS_FIT_TESTS=true as well), `dev/acprior-sum.sh
acprior-g2`. Its 29 skips are the scale (15), drmTMB (13) and fuzz (1)
gates; none is a brms skip:

    frmtmb           files=279 pass=22293 fail=0 err=0 skip=15 warn=0
    frmtmb.coupling  files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam       files=32 pass=1844 fail=0 err=0 skip=3 warn=0
    frmtmb.latent    files=12 pass=681 fail=0 err=0 skip=2 warn=0
    frmtmb.learn     files=15 pass=501 fail=0 err=0 skip=2 warn=0
    frmtmb.ode       files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample    files=68 pass=3397 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial   files=8 pass=529 fail=0 err=0 skip=0 warn=0
    frmtmb.survival  files=15 pass=552 fail=0 err=0 skip=0 warn=0
    logs whose frmtmb is not the lane lib: 0
    logs without a RESULT line: 0

### One session, three orders: core, sample, spatial (tier acprior-s3)

    SESSION frmtmb forward files=269/279 pass=19857 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=269/279 pass=20007 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=269/279 pass=20100 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb.sample forward files=63/68 pass=3114 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample reverse files=63/68 pass=3114 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample preload files=63/68 pass=3114 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.spatial forward files=8/8 pass=454 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.spatial reverse files=8/8 pass=454 fail=0 err=0 skip=5 warn=0
    SESSION frmtmb.spatial preload files=8/8 pass=454 fail=0 err=0 skip=5 warn=0
    ACPRIOR-S3 ran 9 of 9
    untagged LEAK lines: 0

Every LEAK line is tagged "[with load of ...]".

### OpenBLAS 0.3.26 with LAPACK

`La_version()` 3.12.0; the copy is deleted. The changed and core files
(`dev/acprior-out/ob2-*.txt`): test-priors-autocor-classes.R pass=151,
test-default-priors-autocor.R pass=51, test-default-priors-bounds.R
pass=21, all fail=0 err=0. The new assertions compare a ratio of order
1 before against 1e-6 or 1e-10, or one process's draws, so no rounding
decides one.

### R CMD check --as-cran, fix round 1

    frmtmb: Status: 1 NOTE
    frmtmb.sample: Status: OK
    frmtmb.spatial: Status: 1 NOTE

Both NOTEs are the environmental V8 math-rendering NOTE. Logs
`dev/acprior-out/check2-<pkg>.log`; the check trees are deleted.

### Bump and floor, after fix round 1

Unchanged from round 0: frmtmb minor (now three new sampling-API
exports); frmtmb.sample minor, floor the new frmtmb (it calls
`autocor_flat_priors()`, `autocor_bound_chart()` and
`autocor_innov_sd()`); frmtmb.spatial patch, floor the new frmtmb.

### Overlaps added in fix round 1

- `R/priors.R`: `prior_logdens()` gained an argument `lc`, and
  `neg_log_prior_fn()` computes it for an entry with `innov`; the loglik
  lane does not touch the prior functions, but any lane that calls
  `prior_logdens()` positionally is unaffected (`lc` is last, default 0).
- frmtmb.sample `R/sample.R`: `prior_augmented_obj()` and
  `ncp_objective()` gained an argument `chart`; the draws pass through
  `chart_draws_back()` before the GEV and `mo()` back-transforms. The
  mocode lane edits `mo_*` helpers near these; my edits are the two
  signatures, the two `MakeADFun()` calls, and new functions after
  `prior_augmented_obj()`.

## Fix round 2 (review "Fix round 1 check", dev/reviews/2026-10-27-acprior.md)

### 1. A user prior on sigma is read at the innovation sd too (major)

The defect of fix round 1: the default on sigma under a first-order
`cov = TRUE` ar, ma or arma read brms's innovation sd, but a user prior
with the same text, or the sampling table passed back as `prior =`,
read frmtmb's marginal sd, up to 0.72 nats apart (reviewer,
`dev/review-acprior-r1-trip.R`, `r1-trip.txt`). So the sampling table
did not round-trip.

The fix, in core `resolve_priorlist()` (about 25 lines with the
refusal): the `innov` field is attached to every natural-scale entry on
class sigma whose response has `autocor_innov()`, whoever writes it,
instead of to specs that carry frmtmb.sample's `innov_sd` flag. So the
same reading holds on `frm()`'s MAP route and on `frm_sample()`'s.
`target_coefs()` carries the target's `dpar`. frmtmb.sample keeps the
flag for its message line only.

Tests (seen to fail on the fix-round-1 build first, logs
`dev/acprior-out/t3-*-r1build.txt`):
- core, `test-priors-autocor-classes.R`, "a sigma prior under a cov =
  TRUE term is read on brms's sigma": a typed `student_t(3, 0, 2.5)` on
  sigma has brms's density, `stats::dt()` at `exp(eta) c` with the
  Jacobian, at 5 points each for ar, ma and arma (spread ratio below
  1e-10). On the fix-round-1 build: 0.209, 0.089, 0.153. The absent
  case: cosy, order-two ar and cov = FALSE ar get no `innov`.
- frmtmb.sample, `test-default-priors-autocor.R`, "the sampling table
  passed back is the defaults' density": the table as `prior =`, and
  the default's sigma text typed, give the defaults' density, compared
  entry by entry with `expect_identical()` at 4 points each for ar, ma
  and arma. On the fix-round-1 build the summed density differed by
  0.0002 to 0.29 (24 failures). The first version of this test compared
  the summed densities; after the fix those were 1 ulp apart, from the
  order the two stacks sum their entries in, so the test compares the
  entry for each parameter exactly instead. Those per-entry densities
  are bitwise equal.
- frmtmb.sample, "frm()'s MAP under the sampling table is the sampler's
  mode": `frm(prior = default_prior(route = "sample"))` is a stationary
  point of the objective `frm_sample()` samples by default (gradient
  ratio below 1e-3; 0.0102 on the fix-round-1 build), and the two
  objectives agree there (145.1547 against 145.1505 on that build).
- On rellib-r13 these two frmtmb.sample tests pass, since there both
  readings were the marginal sd; the core sigma test fails there as on
  the fix-round-1 build.

    t3-pac-r1build: RESULT test-priors-autocor-classes.R pass=157 fail=4 err=1 skip=0 warn=0
    t3-pac-before: RESULT test-priors-autocor-classes.R pass=67 fail=15 err=3 skip=0 warn=0
    t3-pac-after: RESULT test-priors-autocor-classes.R pass=166 fail=0 err=0 skip=0 warn=0
    t3-dpa-r1build: RESULT test-default-priors-autocor.R pass=51 fail=26 err=0 skip=0 warn=0
    t3-dpa-before: RESULT test-default-priors-autocor.R pass=58 fail=17 err=3 skip=0 warn=0
    t3-dpa-after: RESULT test-default-priors-autocor.R pass=89 fail=0 err=0 skip=0 warn=0

### 2. A bound on sigma there is refused (major)

brms 2.23.0 takes it, on its innovation sd
(`dev/acprior-brms-sigmabound.R`,
`brms-sigmabound.txt`): `set_prior("", class = "sigma", lb = 0.2)` gives
`real<lower=0.2> sigma` and `ub = 3` gives `real<lower=0,upper=3>
sigma`, with `cholesky_cor_ar1(ar[1], ...)` the covariance. Matching it
would need a chart on sigma that moves with ar. frmtmb refuses a `lb`
above 0 or any `ub` on sigma under such a term, with the reason:
"class = "sigma" takes no lb here: under ar(t, g, cov = TRUE) brms's
sigma, which a prior on class "sigma" is about, is the innovation sd,
frmtmb's marginal sd times a factor of the correlation parameters, so a
bound on it is not a bound on the parameter frmtmb holds. Write the
density alone". `lb = 0` is the same box on either scale and is taken
(it is the identity-link default's own bound). Test "a bound on sigma
under a cov = TRUE term is refused", with the absent cases (cosy takes
`lb = 0.2`; identity-link sigma takes `lb = 0`).

### 3. Reporting sigma as the innovation sd: backlog, not this lane

The docs now say it plainly: core NEWS, frmtmb.sample NEWS, `?set_prior`
("Residual correlation") and `?frm_sample` ("Default priors") state
that `sigma` under `ar/ma/arma(cov = TRUE)` is the marginal sd and that
its prior is read at the innovation sd, as brms's is.

**Backlog entry: report brms's `sigma` (the innovation sd) under a
first-order `cov = TRUE` ar, ma or arma.** Measured by the reviewer
(`dev/review-acprior-r1-report.R`, `r1-report.txt`, AR(1) 0.9 data):
brms's posterior mean of `sigma` is 0.432; frmtmb reports the marginal
sd under the same name: `summary(fit)` 1.08, `VarCorr(fit)$residual__$sd`
1.078, `sigma(fit)` 1.078, `posterior_summary(ds)` and `summary(ds)`
1.158, the draws' `sigma` column 1.158; `variables(ds)` lists `sigma`.
brms's draws carried to the marginal sd give 1.149 and frmtmb's carried
to the innovation sd 0.4324, so the model is the same and only the name
differs. It touches `summary()`, `VarCorr()`, `sigma()`, the draws and
`variables()`, the predictions and `posterior_predict()` (which read
sigma), the ML standard errors (a delta-method row in the new
coordinates), `log_lik()` (the loglik lane's), `hypothesis()` names,
and the documented "Divergence from brms" in `?frmtmb-autocor`. A
larger API change; frmtmb has no users, so it is allowed, in a lane of
its own.

### Runs of this round

    Ungated, one file per process (acprior-u3):
    
    frmtmb           files=279 pass=19872 fail=0 err=0 skip=204 warn=0
    frmtmb.coupling  files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam       files=32 pass=1838 fail=0 err=0 skip=5 warn=0
    frmtmb.latent    files=12 pass=681 fail=0 err=0 skip=2 warn=0
    frmtmb.learn     files=15 pass=430 fail=0 err=0 skip=13 warn=0
    frmtmb.ode       files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample    files=68 pass=3152 fail=0 err=0 skip=11 warn=0
    frmtmb.spatial   files=8 pass=454 fail=0 err=0 skip=5 warn=0
    frmtmb.survival  files=15 pass=552 fail=0 err=0 skip=0 warn=0
    logs whose frmtmb is not the lane lib: 0
    logs without a RESULT line: 0
    
    Gated (acprior-g3; 29 skips: scale 15, drmTMB 13, fuzz 1; no brms skip):
    
    frmtmb           files=279 pass=22308 fail=0 err=0 skip=15 warn=0
    frmtmb.coupling  files=11 pass=551 fail=0 err=0 skip=5 warn=0
    frmtmb.eam       files=32 pass=1844 fail=0 err=0 skip=3 warn=0
    frmtmb.latent    files=12 pass=681 fail=0 err=0 skip=2 warn=0
    frmtmb.learn     files=15 pass=501 fail=0 err=0 skip=2 warn=0
    frmtmb.ode       files=11 pass=549 fail=0 err=0 skip=1 warn=0
    frmtmb.sample    files=68 pass=3435 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial   files=8 pass=529 fail=0 err=0 skip=0 warn=0
    frmtmb.survival  files=15 pass=552 fail=0 err=0 skip=0 warn=0
    logs whose frmtmb is not the lane lib: 0
    logs without a RESULT line: 0
    
    One session, three orders (acprior-s4):
    
    SESSION frmtmb forward files=269/279 pass=19872 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb reverse files=269/279 pass=20022 fail=0 err=0 skip=204 warn=0
    SESSION frmtmb preload files=269/279 pass=20115 fail=0 err=0 skip=206 warn=0
    SESSION frmtmb.sample forward files=63/68 pass=3152 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample reverse files=63/68 pass=3152 fail=0 err=0 skip=11 warn=0
    SESSION frmtmb.sample preload files=63/68 pass=3152 fail=0 err=0 skip=11 warn=0
    ACPRIOR-S4 ran 6 of 6
    untagged LEAK lines: 0
    
    OpenBLAS 0.3.26 LAPACK (La_version 3.12.0; copy deleted):
    RESULT test-default-priors-autocor.R pass=89 fail=0 err=0 skip=0 warn=0
    RESULT test-priors-autocor-classes.R pass=166 fail=0 err=0 skip=0 warn=0
    
    R CMD check --as-cran (V8 NOTE only):
    frmtmb: Status: 1 NOTE
    frmtmb.sample: Status: OK
    frmtmb.spatial: Status: 1 NOTE
