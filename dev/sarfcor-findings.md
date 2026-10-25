# Lane sarfcor: brms's sar() and fcor(), and the open parts of the covariance-structure seam

Base 18b536fc (frmtmb 0.74.0 plus two CI commits), worktree
`frmtmb-wt-sarfcor`, branch `wt-sarfcor`, library
`C:/Users/adf44/source/r/wt-sarfcor-lib` in front of `rellib-r12`.
Nothing is committed. Logs are under `dev/sarfcor-out/`,
`dev/sarfcor-log/` and `dev/release/sarfcor-*` (gitignored); every
number below was pasted from them, and each names its script and seed.

Part A and part B are separate patches: `dev/sarfcor-out/A.diff` is the
base commit to the tree with part A alone, and `dev/sarfcor-out/B.diff`
is that tree to the final one (`dev/sarfcor-ab-diff.sh` writes both).

## Part A. brms's sar() and fcor()

### A1. What brms 2.23.0 does (run, not recalled)

Scripts `dev/sarfcor-brms-probe.R` (log `brms-probe.txt`) and
`dev/sarfcor-brms-probe2.R` (`brms-probe2.txt`): `stancode()`,
`standata()`, `default_prior()` and `brm(empty = TRUE)`, no compile.
brms's source in `dev/brms-suite/brms` was read for the post-fit
functions the probes do not reach (`log_lik.R`, `posterior_predict.R`,
`posterior_epred.R`).

Families. `stancode()` for nine families:

    sar:  gaussian, student -> code; poisson, lognormal, skew_normal,
          bernoulli, Gamma, exgaussian, asym_laplace -> "SAR terms are
          not implemented for this family."
    fcor: gaussian, student -> code; the same seven -> "FCOR terms are
          not implemented for this family."

Likelihood. With `A = I - rho W` (`W_tilde = add_diag(-rho * W, 1)`):

    normal_lagsar_lpdf    = N/2 log(1/sigma^2) + log_det
                            - 0.5 |A y - mu|^2 / sigma^2
    normal_errorsar_lpdf  = the same with A (y - mu)
    log_det               = sum(log1m(rho * eigenW))
    student_t_lagsar_lpdf = -K/2 log(nu) + lgamma((nu+K)/2) - lgamma(nu/2)
                            + K/2 log(1/sigma^2) + log_det
                            - (nu+K)/2 log1p(|A y - mu|^2 / (sigma^2 nu))
    normal_fcor_hom_lpdf  = multi_normal_cholesky_lpdf(y | mu, sigma * L)
    normal_fcor_het_lpdf  = ... diag_pre_multiply(sigma, L)   (sigma ~ x)
    student_t_fcor_*      = multi_student_t_lpdf(y | nu, mu, Cov)

with `L = cholesky_decompose(Mfcor)`. The SAR functions leave out
`-N/2 log(2 pi)` (student: `-N/2 log(pi)`); the fcor ones do not. brms
marks all of them "normalized".

Parameters. `lagsar` (lag) and `errorsar` (error), declared
`real<lower=min_eigenMsar,upper=max_eigenMsar>`, flat prior in
`default_prior()`. `fcor()` adds no parameter: `sigma` is brms's own
`sigma`, so a covariance `M` is scaled by `sigma^2` (`standata()` keeps
`diag(Mfcor) = 2` for `M = 2 V`), and a vector `M` becomes its diagonal.

Constraints and checks on `M`:

- sar: `validate_sar_matrix()` takes a `matrix`, a `listw` (spdep
  `listw2mat()`) or an `nb` (spdep `nb2mat()`, row-standardized by
  default); anything else "must be of class 'matrix', 'listw', or
  'nb'"; no check of symmetry or row standardization (a 0/1 adjacency
  matrix passes, eigenvalues -2.83 to 2.83, and its rho range in
  brms is then (-2.83, 2.83)); `M` subset by the dropped NA rows, then
  "Dimensions of 'M' for SAR terms must be equal to the number of
  observations"; eigenvalues by `eigen(M)$values`, passed to Stan as
  they come.
- fcor: a vector becomes `diag()`; a matrix must be symmetric; the
  definiteness check is `min(eigen(M)$values <= 0)`, the minimum of a
  LOGICAL vector, so it refuses only when every eigenvalue is
  non-positive. An `M` with eigenvalue -0.228 was accepted
  (`brms-probe.txt`, "fcor indefinite M"); a negative definite one was
  refused. `M` is subset by the dropped rows as for sar.
- sar refuses a predicted `sigma` or `nu` ("SAR models are not
  implemented when predicting 'sigma' or 'nu'"); fcor takes a
  predicted `sigma` (the het form) and writes
  `student_t_fcor_hom_lpdf(Y | nu, ...)` with a vector `nu` for a
  predicted `nu`, which is not a valid call of a `real nu` function.
- Both: only on `mu` ("Explicit covariance terms can only be specified
  on 'mu'"), not with `rescor`, one per class, one covariance matrix
  of natural residuals (sar + fcor, fcor + cosy/unstr, sar +
  ar(cov = TRUE) refused), and both combine with `ar(cov = FALSE)` and
  random effects.
- Addition terms: fcor refuses `weights`, `se`, `cens` ("Invalid
  addition arguments for this model"), and writes code for `mi()`. sar
  writes code for all of them, and the code for `weights`
  (`normal_lagsar_lpdf(Y[n] | mu, ...)`), `se` (a vector `se` for
  `real sigma`), `cens` and `trunc` (`normal_lagsar_lccdf`, which
  brms never defines) does not type-check in Stan.
- New data (`brms-probe2.txt`): `standata(fit, newdata = 4 rows)`
  refuses both ("Dimensions of 'M' ... must be equal to the number of
  observations"); with `newdata2 = list(W = W[1:4, 1:4])` both give a
  4 x 4 matrix.
- `log_lik()` for both is a leave-one-out conditional density per row
  (`log_lik_gaussian_lagsar()`, whose code cites the loo vignette on
  non-factorized models), not a factor of the joint density.
  `posterior_epred()` for lag SAR is `solve(I - rho W, mu)`;
  `posterior_predict()` draws the whole response.

Defects of brms 2.23.0 found here (none posted upstream; corrected in
fix round 1, where the review showed 1 and 3 overstated):

1. The fcor definiteness check passes an indefinite `M` (above). Not a
   wrong model: brms's Stan program takes `cholesky_decompose(Mfcor)`
   in transformed data, and the fit fails there (the review's
   `fcorindef.txt`). brms fails late; frmtmb refuses at the start, with
   the same outcome. A note in `dev/upstream-bugs.md`, not an entry.
2. A row-standardized `W` of an irregular graph, the spdep default,
   comes out of `eigen()` with a rounding imaginary part (2.4e-17 on
   the 3 x 3 rook lattice of the probe), and brms passes the complex
   vector to Stan: the brms fit of `sar(Wrow)` on the 8 x 8 lattice
   failed with "variable does not exist; processing stage=data
   initialization; variable name=eigenMsar" (`brms-fit.txt`, arm
   lag-row).
3. For a `W` that is not row-standardized, brms's interval for rho,
   (min eigenvalue, max eigenvalue), reaches past 1 / max eigenvalue,
   where `log1m(rho * eigenW)` is not a number. Stan rejects such a
   proposal, so for a W with min eigenvalue <= -1 <= 1 <= max
   eigenvalue (every 0/1 adjacency) brms's support is frmtmb's
   interval: the cost is rejected proposals and initial values, not a
   different model. A note in `dev/upstream-bugs.md`, not an entry.
4. sar with weights, se, cens or trunc: code that does not compile.

Entries in `dev/upstream-bugs.md`: defect 2 is brms-36, defect 4
brms-37, and fix round 1 found brms-38 (`posterior_epred()` fails on
every error SAR fit).

### A2. What frmtmb does

**fcor() in core.** `R/autocor-seam.R` (new) holds the residual
correlation registry, the exported seam `frmtmb_register_autocor()`, and
core's own entry for `fcor()`. A registered term is one affine map,
`whiten(y, mu, sigma, theta, aux) -> list(u, logjac)` on the tape and
its inverse `unwhiten(u, mu, sigma, theta, aux)` on numbers. Core writes
from the map:

- the density (`autocor_seam_loglik()`): the gaussian of `u` plus
  `logjac`, or the multivariate t of `u` with scale `I`, with
  `lgamma_shift_diff()` for the head as the ar() family's student does;
- the mean, `unwhiten(0)`, which `fitted()`, `predict(type =
  "response")` and the response residuals read
  (`autocor_mean_dpars()`); for fcor it is `mu` bit for bit;
- the draws, `unwhiten()` at drawn innovations (one chi-square scaling
  for the multivariate t);
- the covariance over `sigma^2`, `B B'` with `B` read off `unwhiten()`
  column by column, which `autocor_matrix()` returns and the pearson
  residuals divide by (`sigma * sqrt(diag(B B'))`).

fcor's map: `u = L^-1 (y - mu) / sigma`, `logjac = -log|L| - sum(log
sigma)`, `L` the Cholesky factor of `M` over the kept rows, with
`L^-1` built once. `M`: a square numeric matrix or a vector (diagonal);
finite; `n_all` rows (every row of the data), then the kept rows; symmetric;
positive definite by its smallest eigenvalue (frmtmb refuses the
indefinite `M` brms accepts). The term takes a predicted `sigma`.

**sar() in frmtmb.spatial** (`R/sar.R`), through the seam alone. Map:
`u = (A y - mu) / sigma` (lag) or `A (y - mu) / sigma` (error),
`logjac = sum(log1p(-rho ev)) - n log(sigma)`, `A = I - rho W` with `W`
a sparse `Matrix` (a sparse matrix times an advector stays on the
tape: `dev/sarfcor-rtmb-probe.R`, value and jacobian equal to the dense
product). `rho = lo + (hi - lo) plogis(theta)` on `(1 / min ev, 1 /
max ev)`, spatialreg's default interval, starting at rho = 0. The
eigenvalues are computed once; an imaginary part above 1e-8 of the
largest modulus is refused by name ("complex"), and a matrix without a
negative and a positive eigenvalue is refused (W = I, brms's own test
matrix, is one). `listw` and `nb` are read without spdep, as
`listw2mat()` and `nb2mat()` read them (an `nb` with an empty
neighbor set is refused, as `nb2mat()` refuses it). `sigma_varies =
FALSE` refuses a predicted `sigma` with brms's words.

Why not the covariance-structure seam: an entry of
`frmtmb_register_covstruct()` is the density of a latent vector `b`
that the linear predictor ADDS, under the family's own per-row density.
The error SAR's residual IS the correlated vector, with no iid residual
beside it, and the lag SAR's mean is `A^-1 mu`, which is not `mu + b`
for any `b` that does not depend on the coefficients. So the term needs
a seam at the level of the response density; `frmtmb_register_autocor()`
is the smallest one that holds both brms terms: two functions, plus the
bookkeeping every registration has (args, parse, build, npar, start,
natural, `sigma_varies`, usage, `compat_note`, package, replace).

Choices, each against brms:

| Point | brms 2.23.0 | frmtmb | why |
|---|---|---|---|
| sar rho interval | (min ev, max ev) | (1 / min ev, 1 / max ev) | for an unstandardized W brms's effective support is the same, through rejected proposals (defect 3); for a row-standardized W brms's is narrower on the negative side, so frmtmb can estimate rho beyond it (fix round 1, F8) |
| row-standardized W of an irregular graph | fit fails (defect 2) | real part taken | the imaginary part is rounding |
| fcor with indefinite M | passes the R check, fails in Stan (defect 1) | refused | not a covariance |
| sar with weights/se/cens/trunc/mi | code that does not compile, or (mi) fits | refused | the density is joint, as for ar(cov = TRUE) |
| fcor with mi() | fits | refused | as frmtmb's other residual covariance terms |
| sar or fcor with ar(cov = FALSE) | fits | refused: one residual term per response | not built |
| new data | same rows taken; other rows refused unless `newdata2` | same rows taken (fix round 1, F3); other rows refused; frmtmb has no `newdata2`; `frm_linpred()` works | matches brms's default |
| `log_lik()` | LOO-conditional per row | refused, as for ar(cov = TRUE) | a pointwise density of the joint likelihood is not defined |
| `set_prior(class = "lagsar")` | flat on its interval, a class of its own | a class of its own, carried to the internal parameter; frmtmb.sample's default flat on its interval (fix round 1, F2) | the `bounds` slot |

### A3. Evidence

Log densities against brms's Stan program (`dev/sarfcor-lpcheck.R`,
`lpcheck.txt`; `brms_lp_check()` of helper-brms.R: brms's `log_prob()`
at frmtmb's estimates, flat priors, minus frmtmb's `logLik()`; fcor on
`fcor_data()` seed 1, n = 30; sar on a 7 x 7 lattice seed 2 with the
symmetric normalization `D^-1/2 A D^-1/2`, whose eigenvalues are real,
because of defect 2):

    fcor           lp - logLik 7.1054273576e-15  expected 0  diff 7.11e-15  max|grad| 2.52e-07
    fcor het       lp - logLik 0  expected 0  diff 0  max|grad| 2.73e-06
    fcor student   lp - logLik 8.28638683004e-07  expected 0  diff 8.29e-07  max|grad| 2.54e-07
    sar lag        lp - logLik 45.027988127  expected 45.027988127  diff -7.11e-15  max|grad| 2.46e-06
    sar error      lp - logLik 45.027988127  expected 45.027988127  diff -1.42e-14  max|grad| 0.000168
    sar student    lp - logLik 28.0458865746  expected 28.0458822033  diff 4.37e-06  max|grad| 1.17e-07
    nu at the fcor student estimate: 4.887e+09
    nu at the sar student estimate: 5.142e+09

The expected constants are `N/2 log(2 pi)` and `N/2 log(pi)` with
N = 49. The two student rows sit at nu near 5e9 (one response is ONE
multivariate-t draw, so nu is not identified and runs off), where
brms's `lgamma((nu + K)/2) - lgamma(nu/2)` loses digits to
cancellation; the tests check the student density at nu = 5 against a
dense reference instead (below). The same checks run gated in
`test-fcor.R` and frmtmb.spatial's `test-sar.R` (gated logs
`fcor-gated.log`, pass=48 fail=0; `sar-gated.log`, pass=53 fail=0).

Fit agreement with brms fits on the same data (`dev/sarfcor-brms-fit.R`,
`brms-fit.txt`; 8 x 8 lattice, 64 rows, seed 11; brms 4 chains x 2000,
seeds 101 to 103; rho 0.6, the symmetric W of above):

    lag        b_Intercept  ML   1.55067 SE  0.28248 | brms median   1.59604 sd  0.29049 | (ML - median)/SE  -0.161  sd/SE  1.028
    lag        b_x          ML   0.70862 SE  0.09378 | brms median   0.70653 sd  0.10124 | (ML - median)/SE   0.022  sd/SE  1.079
    lag        sigma        ML   0.66887 SE  0.06000 | brms median   0.69513 sd  0.06419 | (ML - median)/SE  -0.438  sd/SE  1.070
    lag        lagsar       ML   0.42060 SE  0.11130 | brms median   0.40246 sd  0.11342 | (ML - median)/SE   0.163  sd/SE  1.019
    lag        max Rhat 1.0014
    error      b_Intercept  ML   1.03251 SE  0.17858 | brms median   1.04726 sd  0.22874 | (ML - median)/SE  -0.083  sd/SE  1.281
    error      b_x          ML   0.65842 SE  0.07880 | brms median   0.65885 sd  0.08193 | (ML - median)/SE  -0.005  sd/SE  1.040
    error      sigma        ML   0.61297 SE  0.05636 | brms median   0.63327 sd  0.06218 | (ML - median)/SE  -0.360  sd/SE  1.103
    error      errorsar     ML   0.56704 SE  0.12132 | brms median   0.57703 sd  0.12835 | (ML - median)/SE  -0.082  sd/SE  1.058
    error      max Rhat 1.0018
    fcor       b_Intercept  ML   0.96433 SE  0.22331 | brms median   0.96206 sd  0.23085 | (ML - median)/SE   0.010  sd/SE  1.034
    fcor       b_x          ML   0.63040 SE  0.06253 | brms median   0.63154 sd  0.06452 | (ML - median)/SE  -0.018  sd/SE  1.032
    fcor       sigma        ML   0.77528 SE  0.06853 | brms median   0.79759 sd  0.07633 | (ML - median)/SE  -0.326  sd/SE  1.114
    fcor       max Rhat 1.0013
    lag-row brms ERROR: The model does not contain posterior draws.

The largest gap is sigma's, 0.44 SE, the ML sd's known downward bias
against a posterior median. The lag-row arm is defect 2.

Tests (one file per process, NOT_CRAN, lane library):

- core `test-fcor.R` (new): pass 42, skip 1 ungated (the brms row);
  pass 48 gated. The dense gaussian at the estimates (1e-10) and its
  own maximum (the optimized dense likelihood is not higher, 1e-7
  relative); covariance `M` scales `sigma` by exactly 1/2 for `4 V`; a
  vector `M`; the het form against `D M D`; pearson residuals; the
  multivariate t at nu = 5 against `mvtnorm::dmvt` (1e-10); dropped NA
  rows equal to the subset fit, bit for bit; a random effect against
  the dense marginal (1e-8); REML is GLS (`sigma^2` ratio `n/(n-p)` to
  1e-6); `profile = TRUE`; draws' mean and covariance against
  `sigma^2 M` at 4000 draws (6 sampling sds); every refusal with its
  absent case; a multivariate model; `frm_compat()`.
- core `test-autocor-seam.R` (new): pass 36. A term registered through
  the export alone (the help page's example, `myfix(V)`) gives
  `fcor()`'s logLik, coefficients, fitted values, draws and
  `autocor_matrix()` bit for bit; `sigma_varies = FALSE` refuses a
  predicted sigma and its absent case fits; a term with a parameter
  reports it through `natural` in `variables()` and
  `confint_varcorr()` and equals the same model written as
  `fcor() + sigma ~ h` (1e-7); eleven refused registrations leave the
  registry and the vocabulary as they were; a slot of the wrong shape
  is refused at the frame; a saved fit whose term is gone names its
  package (`frmtmb_extension_missing`); nothing leaks to the next test.
- frmtmb.spatial `test-sar.R` (new): pass 47, skip 1 ungated; pass 53
  gated. The dense gaussian for lag and error at the estimates (1e-10)
  and the lag maximum; fitted is `A^-1 mu` (lag) and `mu` (error);
  residuals and pearson residuals; `autocor_matrix()` is
  `(A'A)^-1`; the multivariate t at nu = 5, rho = 0.3 against a dense
  reference (1e-10); rho below brms's lower bound on a queen lattice
  (data drawn at rho = -0.8; min eigenvalue above -0.8; estimate below
  it); an unstandardized W keeps `rho * max ev < 1`; sparse, `nb` and
  `listw` inputs give the matrix fit's logLik (1e-12); dropped NA rows;
  draws' moments; refusals and their absent case (a row-standardized W
  whose `eigen()` values are complex by rounding fits); a random effect
  against the dense marginal of the lag model (1e-8); `frm_compat()`.
- core `test-conditions.R`: the block that pinned the two refusals now
  checks `fcor()` with an identity `M` against the independent model
  (1e-10) and the frmtmb.spatial refusal of `sar()` when it is not
  registered. `test-guardrail1.R`: one more probe, `sar(M, type =
  'error')`. `test-compat.R`: the autocor keys hold fcor, and extra
  keys must be registered terms.
- frmtmb.sample `test-autocor-seam-draws.R` (new): pass 5. A term
  registered in the test reports `lagsar` as a derived draws column,
  read off each draw's internal parameter (1e-12).

Seen to fail on rellib-r12 (0.74.0): fcor and sar were refused by name
there, so `test-fcor.R`, `test-sar.R` and `test-autocor-seam.R` fail at
their first fit (the weak form, the refusal). The behavioral form for
frmtmb.sample: before the pattern change, the draws of a sar fit had
the columns `b_Intercept b_x sigma thetaac_1 lp__` and no `lagsar`
(`dev/sarfcor-sample-probe.R`, `sample-probe-before.txt`); after,
`... thetaac_1 lagsar` (`sample-probe-after.txt`).

Post-fit probe (`dev/sarfcor-postfit-probe.R`, `postfit-probe.txt`,
seed 9, 6 x 6 lattice with `(1 | h)`): `fitted(re_formula = NA)` and
`fitted()` equal `A^-1 X b` and `A^-1 (X b + Z u)` to 6.7e-16 and
2.2e-15; `conditional_effects()` drops the term, as brms's
`incl_autocor = FALSE` does; `bf(autocor = ~ sar(W))` and
`acformula(~ fcor(V))` give the same logLik as the term in the formula
(difference 0); `update()` keeps the term; `summary()` prints
"Correlation Structures: lagsar 0.52 (0.11, 0.92)".

frmtmb.sample on the two fits (`sample-probe-after.txt`): sampling,
`posterior_predict()` (600 x 36) and `posterior_epred()` (largest
difference of the draw mean from `fitted()` 0.019 and 0.022) work, and
`log_lik()` refuses with its residual-covariance message.

### A4. Guardrail 1

`extension_special_homes` gains `sar = "frmtmb.spatial"` (on its own
line), and `extension_autocor_specials` (new) says which homes are
residual terms. `check_extension_specials()` skips a name either
registry holds. `brms_special_refusals` is empty. The refusal text that
ended "The supported residual correlation terms are ar(), ma(),
arma(), cosy() and unstr(); the car() field is in the extension
package frmtmb.spatial" is now built by `autocor_supported_sentence()`
from the terms the session resolves.

Classification of every brms name in three arms
(`dev/sarfcor-guard.sh`, which runs `dev/spatial-guard-report.R` with
test-guardrail1.R's probes and classifier; logs `guard-base.txt`,
`guard-lane.txt`, `guard-lane-spatial.txt`), the rows that moved:

    probe                          base (0.74.0)      lane, core alone      lane + frmtmb.spatial
    sar(M)                         refused            ext:frmtmb.spatial    core
    sar(M, type = 'error')         refused            ext:frmtmb.spatial    core
    fcor(M)                        refused            core                  core

`test-guardrail1.R` (one more probe): pass 35, fail 0 with frmtmb.spatial
installed (its extension arm ran); it asserts that the names refused by
name equal the declared refusals, which are now `custom` and `bhaz`.

### A5. The ported brms-suite rows

No row of `dev/brmsport-verdicts.tsv` (bin 1) uses `sar` or `fcor`: its
`standata` and `stancode` blocks that do were never ported. Bin 2
(`dev/brmsport-bin2-ledger.tsv`) has 13, all "cannot transfer / absent
/ gap" at 0.74.0 because `frm()` refused both terms. They are not
ported, so there is nothing to rerun; `dev/sarfcor-brmsrows.R`
(`brmsrows.txt`) replays what each asserts through the public API on a
fit of the same shape (10 rows, 3 draws of `frm_sample()`):

    posterior_predict:111  literal    M = diag(10) refused: sar(): 'M' needs a negative and a positive eigenvalue, ...
    posterior_predict:111  pass       path-graph W: dim(posterior_predict) = 3 x 10
    posterior_predict:113  ... the same, pass (student)
    posterior_predict:118  ... the same, pass (error)
    posterior_predict:120  ... the same, pass (error, student)
    posterior_predict:133  pass       Mfcor = diag(10): dim(posterior_predict) = 3 x 10
    posterior_predict:135  pass       Mfcor = diag(10): dim(posterior_predict) = 3 x 10
    log_lik:150, 152, 157, 159, 174, 176, 230   refused: log_lik() needs a likelihood that factors into one term per observation ...

So 6 rows can flip from gap to pass at the consolidation's regeneration
(the four sar rows with a weights matrix that has a zero diagonal: brms's
literal `Msar = diag(10)` is refused, rightly, because `I - rho I` makes
rho and the scale one parameter), and 7 move from "gap" to a refusal
by design (no pointwise density of a joint likelihood). The ledger was
not edited.

### A6. What A did not do

- Priors on `lagsar`/`errorsar` by brms's class: built in fix round 1
  (F2), through a `bounds` slot; see there.
- `log_lik()` in brms's leave-one-out form, `mi()`, combination with
  `ar(cov = FALSE)`, `newdata2`: see the table in A2.
- An interval for rho on its bounded scale: `confint_varcorr()` uses
  the delta method on the raw scale (type "raw"), so the upper end can
  pass `1 / max ev` when the estimate is near it.

## Part B. The open parts of the covariance-structure seam

### B1. frmtmb.sample asks the registry

New export `frmtmb_covstruct_info(name)` (R/covstruct-seam.R): the
structure's package, and whether it is built in, a field (reports
through `natural`), a predictor special, levels-coupled and full rank.
frmtmb.sample's `draws_cov_fast()` now reads `field` from it; its list
`c("smooth", "gp", "hsgp", "car", "spde")` is gone. For every built-in
structure the two agree (the five listed are the five built-ins and
spatial structures with `natural`).

The fast path, measured (`dev/sarfcor-b1-fastpath.R`,
`b1-fastpath.txt`, seed 21): a registered field `toyfld(A, gr = g)`
beside `(1 + x | g)` on the SAME grouping factor, as
`car(W, gr = loc) + (1 | loc)` is; 30 groups x 8 rows; 1000 draws of
`frm_sample()`; both arms in one process on the same draws, the old
function taken from 18b536fc's source:

    fast map, old: NULL (every draw through hyp_env_vals())
    fast map, new: sd_g__Intercept cor_g__Intercept__x sd_g__x
    draws: 1000  values identical: FALSE  max |diff|: 2.220446e-16
    per call, min over 5 rounds: old 0.7000 s, new 0.0527 s, control (new again) 0.0548 s; old/new 13.3, control ratio 1.04
    per draw: old 0.7 ms, new 0.0527 ms

Each block repeats the call past 1.2 s (2 calls old, 40 new), because
`proc.time()` ticks at 10 ms. The old map put the field's block on the
names `sd_g__Intercept`, which the `(1 + x | g)` block owns, failed its
own two-point check, and fell back. With the field on a grouping
factor no other block uses, both arms took the fast path (first run of
the same script, 0.04 s each), which is why the spatial lane's "slower,
not wrong" needs the shared factor to show. The values differ by one
ulp between the two paths, as the fast path's own check allows.

Test: frmtmb.sample `test-draws-sd-names.R`, new block "a registered
field over the same grouping factor keeps the fast map": pass on the
lane (16 of 16 in the file); on rellib-r12 it fails behaviorally
(`is.null(dv$fast)` TRUE) and then errors (13 pass, 1 fail, 1 error).

### B2. The weights check reads the structure's rank

The list `unbounded_full_rank_blocks` is gone. Each built-in entry
states `full_rank` (TRUE for the same 13 structures, FALSE for the
rest, so no built-in fit changes), `frmtmb_register_covstruct()` takes
`full_rank =`, and `cs_full_rank()` reads it, or, for a registered
structure that does not state it, reads it off the structure's
covariance at the start: a bar structure with uncoupled levels whose
`vcov` there is positive definite (smallest eigenvalue above
sqrt(eps) times the largest) is full rank; a field, a coupled structure
and a singular covariance are not, so the check stays open for them,
as for the built-in fields.

Constructed (`dev/sarfcor-b2-probe.R`, wu_data() seed 20261020, 72
rows, Gamma(link = "log"), weights 1.5, one level per row):

    base: y | weights(w15) ~ x + (1 | obs) -> REFUSED: The weighted likelihood of this model has no maximum. ...
    base: y | weights(w15) ~ x + myiid(1 | obs) -> FIT logLik 849.1237348 shape 5.471e+14 warnings: 2
    lane: y | weights(w15) ~ x + (1 | obs) -> REFUSED: ...
    lane: y | weights(w15) ~ x + myiid(1 | obs) -> REFUSED: ... `myiid(1 | obs)` gives every row with a positive weight ...

`myiid` is `(1 | obs)` under another name. The lane's message is the
built-in one with the term's name in it (asserted). The absent case, a
structure that is not full rank: `rankone(0 + one + one2 | obs)`, a
two-coefficient block whose covariance has rank one along (1, -1),
on columns `one = one2 = 1`, so the effects add nothing a row can use
to reach its observation. It is not refused, and its maximum is the
model without it (logLik equal to 1e-6, coefficients to 1e-4); its
own sd is flat there, which the fit's warnings say. Test:
`test-weights-unbounded.R`, "a registered structure is checked by its
own rank", pass on the lane (65 in the file); on rellib-r12 it errors
at `expect_match()` because the registered model returned a fit.

### B3. autoscale: a message, not support

autoscale rescales a random slope with its column only for `us`,
`diag`, `us_t` and `diag_t`, where the first `dim` theta entries are
per-component log standard deviations that nothing else reads, so the
rescaling is an exact change of coordinates. Nothing in a registered
entry says that of its parameters: `sd_idx` names log standard
deviations, but not that each belongs to one coefficient alone (the
seam's own toy `xcs` has one sd for every coefficient), and a wrong
guess would change the fitted model, not its coordinates. So I chose
the message: `autoscale_note_unscaled()` names each registered
structure whose component carries a planned column, once per fit (the
call that runs the pre-fit, not the refits handed its template):
"autoscale rescaled the fixed-effect columns but not the random slopes
of mydiag() on x: the structure is registered by another package, and
nothing tells frmtmb that its parameters scale with the column. If the
fit stops short, rescale the column before the fit". Frmtmb's own
unscaled structures (cs, ar1, ...) stay as documented and silent.
Test: `test-autoscale.R`, new block: the message for a registered
diagonal structure on a column of sd 1e-4; no message for `(1 + x ||
g)` and none with `autoscale = FALSE`. Fails on rellib-r12 ("did not
throw a message"; 48 pass, 1 fail).

### B4. rr(), gp()/hsgp() and the Student-t sharing: measured, not ported

`dev/sarfcor-b4-measure.R` (`b4-measure.txt`), on 18b536fc: code lines
of R/ that special-case each by name or by its own helpers (compat.R
notes counted apart), the top-level functions they sit in, and the
commits of 4f5ea39f..18b536fc (0.68.1 to 0.74.0) that changed such a
line.

| Feature | code lines (files) | functions | compat notes | commits, 5 rounds | extension R/ lines |
|---|---|---|---|---|---|
| rr() factor space | 41 (11) | 23 | 13 | 3 (9b93ae2e, ea13efc6, ed824ffb) | 1 (frmtmb.sample non-centering) |
| gp()/hsgp() specials | 65 (15) | 39 | 15 | 5 (9b93ae2e, ea13efc6, ed824ffb, 58b97a2c, 83708f40) | 4 (frmtmb.sample twice, frmtmb.latent, draws_cov_fast's list, now gone) |
| us_t/diag_t sharing | 38 (10) | 22 | 4 | 6 (9b93ae2e, df70786f, ed824ffb, 0a8d9a02, c1618670, 58b97a2c) | 3 (frmtmb.sample's ncp_reasons) |

The patterns are crude (a helper name or the quoted structure name),
so the counts are an order of magnitude, as the spatial lane's were.

What each would take:

- rr()'s factor-space expansion. `b` has `rank x levels` entries and
  the coefficients `dim x levels`, through loadings that depend on
  theta. The published `expand` is a FIXED linear map of equal
  length. Putting rr on it means an `expand` that depends on theta and
  changes length, with its jacobian in both `b` and theta, read by the
  objective (`expand_b()`), every delta method (`rr_jacobians()`),
  ranef's conditional sds, `draw_b()`, the brms names, the prior
  table, conditional_effects levels, new levels and frmtmb.sample's
  non-centering. 41 lines in 11 files, 23 functions; every fit with an
  rr term would retape through new code. Not cheap; left.
- gp() and hsgp() as specials. A field `term` is one intercept per
  level of a grouping variable `gr`; a gp is a function of continuous
  coordinates with kriging at new positions, and hsgp's design is a
  basis matrix, not level indicators. Neither fits the `term` slot
  without a new kind of term. 65 lines in 15 files, 39 functions, 5
  commits in five rounds. Left.
- The Student-t sharing of us_t/diag_t. Blocks over one grouping
  factor share one latent mixing scale per level
  (`student_share_sets()`), a CROSS-block construction; the entry
  shape is per block. It would need a slot naming a sharing key and a
  density over the set of blocks. 38 lines in 10 files, 22 functions,
  6 commits. Left.

None is cheap, and each would change the tape of every fit that uses
it, so none was ported.

## Suites, sessions, OpenBLAS and checks

Every `lib:` line of every log below names frmtmb from
`wt-sarfcor-lib` (checked: 0 of 444 and 0 of 445 per-file logs name
another), and frmtmb.sample and frmtmb.spatial from there too; the
other six extensions are rellib-r12's 0.74.0-round builds.

One file per process, NOT_CRAN, all nine suites
(`dev/release/run-par.sh sarfcor-u1`, `sarfcor-u1.log`, summed by
`dev/sarfcor-sumlog.sh`), on the build before the last compat and
message fixes:

    frmtmb files 276 pass 19013 fail 5 err 0 skip 202 warn 0
    frmtmb.coupling files 11 pass 551 fail 0 err 0 skip 5 warn 0
    frmtmb.eam files 32 pass 1838 fail 0 err 0 skip 5 warn 0
    frmtmb.latent files 12 pass 679 fail 0 err 0 skip 2 warn 0
    frmtmb.learn files 15 pass 430 fail 0 err 0 skip 13 warn 0
    frmtmb.ode files 11 pass 549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample files 65 pass 3028 fail 0 err 0 skip 9 warn 0
    frmtmb.spatial files 8 pass 375 fail 0 err 0 skip 3 warn 0
    frmtmb.survival files 14 pass 405 fail 0 err 0 skip 0 warn 0

The 5 failures were three files, all fixed and rerun (`sarfcor-u2.log`,
18 files, fail 0 err 0 warn 0): test-compat.R (the autocor keys had
fcor; and `student x fcor()` was decided by file order, a tie between
`fcor() x group:autocor_families` and `group:autocor x student`, now
the explicit row `fcor() x student`), test-compat-extensions.R (the
seam's families row and frmtmb.spatial's own tied, so one won no pair;
the seam now takes the package's note, `compat_note =`) and
test-message-uniqueness.R (the new seam repeated three message
templates of the covariance seam; both now call shared helpers,
`seam_check_package()`, `seam_check_name()`, `seam_fn()`).

One file per process, NOT_CRAN and FRMTMB_BRMS_FIT_TESTS, all nine
suites, on the final code (`sarfcor-g1.log`):

    frmtmb files 276 pass 21455 fail 0 err 0 skip 15 warn 0
    frmtmb.coupling files 11 pass 551 fail 0 err 0 skip 5 warn 0
    frmtmb.eam files 32 pass 1844 fail 0 err 0 skip 3 warn 0
    frmtmb.latent files 12 pass 679 fail 0 err 0 skip 2 warn 0
    frmtmb.learn files 15 pass 501 fail 0 err 0 skip 2 warn 0
    frmtmb.ode files 11 pass 549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample files 66 pass 3303 fail 0 err 0 skip 1 warn 0
    frmtmb.spatial files 8 pass 444 fail 0 err 0 skip 0 warn 0
    frmtmb.survival files 14 pass 405 fail 0 err 0 skip 0 warn 0

Test counts that changed, base files on rellib-r12 (`sarfcor-base1.log`,
the base commit's test files) against the lane (`sarfcor-u2.log`):
test-conditions.R 150 -> 150 (two refusal expectations became an fcor
fit and a conditional sar refusal), test-guardrail1.R 35 -> 35,
test-compat.R 763 -> 767, test-registry-helper.R 7 -> 7,
test-weights-unbounded.R 55 -> 65 (B2), test-autoscale.R 46 -> 49 (B3),
frmtmb.sample test-draws-sd-names.R 12 -> 16 (B1). New files:
test-fcor.R 42 (48 gated), test-autocor-seam.R 36, frmtmb.spatial
test-sar.R 47 (53 gated), frmtmb.sample test-autocor-seam-draws.R 5.
No expectation was dropped.

One session per suite, three orders, NOT_CRAN
(`dev/release/run-sessions.sh sarfcor-sess 8`, all nine,
`sarfcor-sess.log`; `dev/release/run-session.R` now also compares the
residual correlation registry after each file): fail 0, err 0 and warn
0 in all 27 sessions:

    SESSION frmtmb forward  files=266/276 pass=19020 fail=0 err=0 skip=202 warn=0
    SESSION frmtmb reverse  files=266/276 pass=19120 fail=0 err=0 skip=202 warn=0
    SESSION frmtmb preload  files=266/276 pass=19213 fail=0 err=0 skip=204 warn=0
    coupling 11/11 pass 551; eam 32/32 pass 1838; latent 12/12 pass 679;
    learn 15/15 pass 430; ode 11/11 pass 549; sample 61/66 pass 3033
    (the 5 brms-suite files contribute no rows in a session, as at
    0.74.0's 60/65); spatial 8/8 pass 375; survival 14/14 pass 405;
    each the same in all three orders.

The run found three untagged LEAK lines, all mine:
`test-autocor-seam.R autocor registry: 1 -> 2 entries` (forward,
reverse) and `2 -> 3` (preload). The last block of that file put the
registry back with `on.exit()`, which ran AFTER
`local_frmtmb_registries()`'s restore and so restored the state with
its own registration in it. The block now relies on the helper alone.
Rerun of frmtmb and frmtmb.spatial, three orders each
(`sarfcor-sess2.log`): the same SESSION lines, 21 LEAK lines, every one
tagged `[with load of ...]`, 0 untagged. The other 42 LEAK lines of the
first run are all tagged too.

OpenBLAS 0.3.26 with LAPACK, 4 threads (`dev/ciharden-openblas.sh
0.3.26 lapack`, then `dev/ciharden-scan.sh` with every gate on), on
the 14 files that read the new numerical code or the entries B
changed (`dev/sarfcor-jobs-ob.lst`: test-fcor, test-autocor-seam,
test-weights-unbounded, test-autoscale, test-covstruct-seam,
test-guardrail1, test-conditions, test-compat, test-autocor,
test-autocor-cond; frmtmb.spatial test-sar, test-car; frmtmb.sample
test-draws-sd-names, test-autocor-seam-draws), against the same scan
on the reference BLAS: 14 of 14 files, fail 0 in both, and the two
`.sum` files are identical line for line
(`dev/ciharden-log/sarfcor-ob.sum`, `sarfcor-ref.sum`). No new test
sits on a rounding-decided branch: the tolerances that decide (the
complex-eigenvalue cut at 1e-8 of the largest modulus, the sign of the
extreme eigenvalues at 1e-8, positive definiteness at n eps) are far
from rounding on every matrix the tests build, so 0.3.32 was not run.
The R copy was deleted after.

R CMD check --as-cran (`dev/sarfcor-check.ps1`, built and checked in
`dev/sarfcor-check/<pkg>`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`,
NOT_CRAN, pandoc and TinyTeX on PATH, R_LIBS = lane library, rellib-r12,
user library; logs `dev/sarfcor-log/check-<pkg>.log`):

    frmtmb          Status: 1 NOTE   (V8 math rendering, expected)
                    tests: [ FAIL 0 | WARN 0 | SKIP 229 | PASS 18932 ]
    frmtmb.spatial  Status: 1 NOTE   (V8 math rendering, expected)
                    tests: [ FAIL 0 | WARN 0 | SKIP 4 | PASS 329 ]
    frmtmb.sample   Status: OK
                    tests: [ FAIL 0 | WARN 0 | SKIP 17 | PASS 3019 ]

The first core check failed to build: R CMD build copies the whole
tree before it reads .Rbuildignore, and it met the frmtmb.spatial
check writing into `dev/sarfcor-check/` at the same moment
("problem copying ... 00install.out: Permission denied"). Rerun alone,
it passed as above. Run core's check alone, or outside the tree being
built.

## Files touched

`git apply --numstat` of the two patches (added, removed, file).

Part A, 49 files (A.diff):

    core R/: autocor-seam.R (new, 582), autocor.R +125 -19,
      brms-vocabulary.R +33 -8, compat.R +19 -2, covstruct-seam.R +42 -20
      (reserved names; the shared check helpers), extension-home.R +5,
      families.R +10, frame.R +14 -1, objective.R +14, parse.R +11 -7,
      predict.R +38 -12, registries.R +14 -10, simulate-newdata.R +3,
      structure.R +13 -1
    core other: NAMESPACE +2, NEWS.md +40, README.md +2 -1,
      _pkgdown.yml +1, vignettes/brms-migration.Rmd +2,
      vignettes/compatibility.Rmd +3 -2; man/autocor-terms.Rd,
      autocor_matrix.Rd, frmtmb-autocor.Rd, frmtmb-extension-api.Rd,
      frmtmb_local_registries.Rd, frmtmb_register_autocor.Rd (new)
    core tests: test-fcor.R (new, 256), test-autocor-seam.R (new, 195),
      test-compat.R +9 -4, test-conditions.R +13 -6,
      test-guardrail1.R +3 -1, test-registry-helper.R +2 -1
    frmtmb.spatial: R/sar.R (new, 306), R/zzz.R +16 -3,
      R/frmtmb.spatial-package.R +5 -3, DESCRIPTION +3 -1, NAMESPACE +2,
      NEWS.md +31, README.md +6 -3, vignettes/spatial.Rmd +45 -2,
      man/sar.Rd (new), man/frmtmb.spatial-package.Rd,
      tests/testthat/test-sar.R (new, 309), helper-brms.R +10
    frmtmb.sample: R/draws-brms.R +4 -2, NEWS.md +8,
      tests/testthat/test-autocor-seam-draws.R (new, 78)
    dev: core-boundary.md +4 -1, release/run-session.R +6

Part B, 16 files (B.diff, on top of A):

    core R/: covstruct-seam.R +65 -2, covstruct.R +44, weights-unbounded.R
      +5 -10, autoscale.R +41, fit.R +5
    core other: NAMESPACE +1, NEWS.md +24, _pkgdown.yml +1,
      man/frmtmb_covstruct_info.Rd (new), man/frmtmb_register_covstruct.Rd
    core tests: test-weights-unbounded.R +73, test-autoscale.R +49
    frmtmb.sample: R/draws-brms.R +4 -1, NEWS.md +11,
      tests/testthat/test-draws-sd-names.R +60
    dev: core-boundary.md +2

`git apply` takes A.diff on a `git archive 18b536fc` and B.diff after
it, and the result equals the worktree outside dev/ (checked with
`diff -r --strip-trailing-cr`; all 96 R files of the A tree parse). The
A tree alone was not installed and run: the lane rules allow one
private library, and it holds the final build. Lane scripts and records
(`dev/sarfcor-*`) are not in either patch.

## Bumps and floors

- frmtmb: a MINOR bump. New exports `frmtmb_register_autocor()`,
  `fcor()` (part A) and `frmtmb_covstruct_info()` (part B); a new
  argument `frmtmb_register_covstruct(full_rank =)`; `fcor()` fits
  where it was refused.
- frmtmb.spatial: a MINOR bump (new term `sar()`). Its floor must be the
  frmtmb release that exports `frmtmb_register_autocor()`: it imports
  it. DESCRIPTION still says `frmtmb (>= 0.74.0)`, because a lane does
  not bump core; set it at consolidation.
- frmtmb.sample: a PATCH bump. Part A alone (the `lagsar` pattern)
  needs no new floor; part B needs the frmtmb release that exports
  `frmtmb_covstruct_info()`, which its `draws_cov_fast()` calls.
  DESCRIPTION still says `frmtmb (>= 0.74.0)`; set it at
  consolidation. Its test `test-autocor-seam-draws.R` needs
  `frmtmb_register_autocor()` too.
- The other six extensions are unchanged and pass on the new core
  (suites above).

## What needs the user

- brms 2.23.0 defects found (A1, items 1 to 4): the fcor
  definiteness check, the complex eigenvalues of a row-standardized W,
  the rho interval of an unstandardized W, and sar code with weights,
  se, cens or trunc that does not compile. Nothing was posted upstream;
  the user decides.
- Priors on `lagsar` and `errorsar`: done in fix round 1 (F2). `ar`
  and `cosy` keep frmtmb.sample's default flat on the internal
  parameter, which is not brms's flat on the natural one; the
  coordinator put that on the backlog.
- `log_lik()` in brms's leave-one-out form for residual covariance
  terms (ar(cov = TRUE), cosy, unstr, fcor, sar): refused by design
  today; 7 bin-2 brms rows read it (A5).
- Overlaps with the other lanes: `R/brms-vocabulary.R` (bhaz changes
  `brms_aterm_refusals` beside this lane's `brms_special_refusals` and
  `extension_special_homes`; my `sar` line is on its own line),
  `R/parse.R` (brmsnames), `R/compat.R`, `R/covstruct-seam.R`,
  `R/covstruct.R`, `NEWS.md`, `_pkgdown.yml`, frmtmb.spatial,
  frmtmb.sample (`draws-brms.R`, NEWS), `tests/testthat/test-guardrail1.R`
  and `dev/core-boundary.md`. test-guardrail1.R asserts that the names
  refused by name equal the declared refusals; with bhaz moved to
  frmtmb.survival that set would be `custom` alone, which the test
  still allows.
- The A/B split puts the shared check helpers of the two seams
  (`seam_check_package()`, `seam_check_name()`, `seam_fn()`) in A,
  because part A's seam needed them for test-message-uniqueness.R.

## Fix round 1 (review `dev/reviews/2026-10-24-sarfcor.md`)

Every fix is uncommitted. A and B stay separable: A edits went to the
stage-A copy as well, part B's to the worktree only, and
`dev/sarfcor-ab-diff.sh` rebuilt `A.diff` (55 files) and `B.diff` (16
files). Both apply with `git apply` on 18b536fc in order, and the
result equals the worktree outside `dev/` (checked file by file for
`R/`, `tests/`, `man/`, `vignettes/`, `extensions/`, NEWS, NAMESPACE,
`_pkgdown.yml`, README and the three dev files). `dev/upstream-bugs.md`
is now one of the dev files the diffs carry (part A).

### Part A

- F1. With `sigma_varies = FALSE` the map gets one sigma
  (`autocor_seam_sigma()`), in the density, the mean, the draws and the
  covariance; sar()'s `whiten()` uses `n log(sigma)`. The help example
  of `frmtmb_register_autocor()` runs (`@examplesIf` withr): it
  registers `myfix()` inside `frmtmb_local_registries()` and gives the
  logLik of core's `fcor()`, -21.15892 both, with the registry put
  back. `autocor_seam_check_map()` runs once per term when the frame is
  built: `unwhiten()` of `u = sin(1:n)` must give n finite values and,
  for `cbind(u, rev(u))`, an n x 2 matrix whose first column is the
  vector result; `whiten()` must give n innovations and one finite
  log-Jacobian; `max |whiten(unwhiten(u)) - u| <= 1e-6 max |u|`; a
  `bounds` slot must match `natural` at the start. A failure is refused
  as "The '<name>' term's map is wrong: ... for <label>". Tests
  (test-autocor-seam.R): `whiten()` one innovation short ("gave 23
  innovation(s) for 24 rows"), `unwhiten()` with `t(L)` (round trip),
  an `unwhiten()` that flattens a matrix, a per-row log-Jacobian; the
  absent case: a correct toy and `fcor()` pass, and test-sar.R fits
  lag and error on a dense and a sparse W through the same check. A
  toy that stops when handed more than one sigma reproduces fcor()'s
  logLik, fitted, simulate and pearson residuals bit for bit.
- F2. A `bounds = function(aux)` slot (one row, lower and upper, per
  parameter) declares that `natural` is `lower + (upper - lower)
  plogis(theta)`. Core then resolves `set_prior(class = "<name>")`
  (the natural-scale spec is held as an Intercept with `dpar = name`)
  through `resolve_ac_class()` with a new map "bounded": the density on
  the natural value plus `log(hi - lo) + log s + log(1 - s)`, and
  `lb`/`ub` as new interval ends (an end of the interval itself is
  refused, as for ar). A class prior on a parameter of a term without
  `bounds` is refused by name with the internal spelling
  (`class = "theta", coef = "thetaac_k"`). New export
  `frmtmb_autocor_pars(fit)` lists resp, term, name, internal coef and
  interval. frmtmb.sample's default adds `uniform(lo, hi)` on each
  bounded parameter by name, flagged `flat_interval`, announced as
  `toyphi  uniform(-0.4, 0.9)  [flat on its interval, brms's (flat)]`
  (for sar on the rook lattice: `lagsar uniform(-0.999999999999998,
  0.999999999999999)`). Cost: about 55 lines in
  `R/priors.R` and 60 in `R/autocor-seam.R` (slot, check, the two
  readers, the export), 15 in frmtmb.sample; no change to any tape
  but through the existing natural-scale prior path, and it reuses the
  logistic Jacobian core already had for ar/cosy. Cheap, so supported
  rather than refused. Flat-likelihood check (a registered term whose
  density ignores its parameter, `dev/sarfcor-f2-flat-probe*.R`): 4
  chains x 5000 draws, ESS 21870, mean 0.2508 (want 0.25), sd 0.3772
  (want 0.3753), deciles within rounding, KS p = 0.58; the test
  (test-autocor-seam-draws.R, 4 x 2500) checks each decile of the
  interval holds a tenth of the draws within 4 ESS-based standard
  errors, and the mean likewise. `ar` and `cosy` keep frmtmb.sample's
  default flat on the internal parameter, which is not brms's flat on
  the natural one: not changed here (coordinator's backlog).
- F3. `newdata` with as many rows as the fit is read as the fitted rows
  (`check_autocor_seam_newdata()`): `fitted()`/`predict()` mean is the
  term's mean at newdata's `mu` (`A^-1 eta_new` for lag), `simulate()`
  and `predict()` draw through the map, other row counts are refused
  with brms's `newdata2` named. Pinned against brms's
  `posterior_epred()` with brms's draws held at frmtmb's estimates
  (Fixed_param), `dev/sarfcor-f3-epred.R`, `sarfcor-out/f3-epred.txt`:
  lag SAR 4.2e-16 relative (the counterfactual moved the mean by
  0.497), fcor 0; brms refuses 5 rows ("Dimensions of 'M' ... must be
  equal"). The same pins are gated tests (test-sar.R, test-fcor.R).
  brms 2.23.0's `posterior_epred()` and `fitted()` fail on EVERY error
  SAR fit ("'a' is 0-diml"), so the error form is pinned against the
  dense reference only; new upstream entry brms-38. Message, NEWS,
  `?frmtmb-autocor` (fcor), `?sar`, `?frmtmb_register_autocor`, the
  compat notes of both terms and both vignettes now say this.
- F4. `ac_formula_terms()` runs `check_extension_specials()` on each
  `autocor =` term, so `bf(y ~ x, autocor = ~ sar(M))` without
  frmtmb.spatial names the package. test-guardrail1.R probes it (core
  alone: `ext:frmtmb.spatial`; with the package: resolves) and
  `autocor = ~ fcor(M)` (core).
- F5. frmtmb.sample's derived names are `union(<ar/cosy pattern>,
  frmtmb_autocor_pars(fit)$name)`; the literal lagsar/errorsar pattern
  is gone. Test: a registered term whose bounded parameter is `toyphi`
  gives the column `toyphi` and a default prior by that name.
- F6. The student brms rows of test-fcor.R and test-sar.R evaluate
  brms's `log_prob` at hand-set parameters with nu = 5 (b0 0.8, b1
  0.4, sigma 0.9, rho 0.35), against frmtmb's objective at the same
  point, to 1e-10 relative; measured (`sarfcor-out/f6-lp5.txt`): lag
  and error 3.6e-15, fcor 1.4e-14. NEWS numbers replaced.
- F7. Defects 1 and 3 corrected above (A1) and in the A2 table; both
  are notes in `dev/upstream-bugs.md` ("Considered and excluded").
  Entries brms-36 (complex eigenvalues, `dev/sarfcor-upstream-eigen.R`:
  max |Im| 7.97e-17, Stan "variable does not exist ... eigenMsar", 0
  chains with draws), brms-37 (stanc rejects sar with weights, se,
  cens, trunc; plain model type-checks; `dev/sarfcor-upstream-stanc.R`)
  and brms-38 (error SAR `posterior_epred()`;
  `dev/sarfcor-upstream-epred.R`: `prep$ac$lagsar` is 200 x 0, not
  NULL). The bhaz lane added a brms-36 in its tree, so consolidation
  may renumber. Nothing posted.
- F8. NEWS and `?sar` give the queen-lattice numbers (eigenvalues
  -0.5006 to 1, data at rho -0.8: frmtmb -1.103, logLik 3.74 above
  brms's best, brms median -0.443, a fifth of its draws within 0.02 of
  its bound). Complete graph: WARN, not refuse. When the extreme
  eigenvalue repeats n - 1 times, `I - rho W` at that end keeps one
  direction; with an intercept the likelihood has no maximum there and
  rho runs to the end (K_9: false convergence, gradient 2.5e4). But
  `build` sees W, not the formula, and a model without an intercept
  (or with the intercept's direction otherwise pinned) keeps its
  maximum, so a refusal would fire on correct models. The warning names
  the eigenvalue, its count, the rho at the end and the remedy. Test:
  K_9 warns (required), a lattice does not.
- F9. `unwhiten()` takes a matrix of innovations; `autocor_seam_B()`
  makes one call with `diag(n)`, sar() does one sparse factorization
  for all columns, and pearson residuals read `rowSums(B * B)`
  (`autocor_seam_marginal()`) instead of the diagonal of `B B'`.
  Measured, n = 2500 (50 x 50 rook, row-standardized,
  `dev/sarfcor-f9-time.R`, two runs each, `sarfcor-out/f9-time.txt`;
  "before" puts the old functions back into the loaded namespaces):

      before  pearson 16.42 / 15.77 s   autocor_matrix 15.96 / 15.76 s
      after   pearson  0.47 /  0.48 s   autocor_matrix  4.32 /  4.34 s
      (frm itself 13.7 to 13.9 s in all four; results identical:
      sum(rp) 4.57336552624, sum(diag M) 3184.2867729)

  Of the 4.3 s, B takes 0.36 s and `tcrossprod(B)` 3.96 s (reference
  BLAS): the dense n x n covariance `autocor_matrix()` returns.
- F13. frmtmb.spatial test-guardrail.R saves a sar fit beside the car
  and spde fits and runs all eleven methods in a child without the
  package (unloadable and absent arms, helper-absent-package.R): 33
  refusals, "The sar() term of this model is from the extension
  package frmtmb.spatial". Three runners: see below.
- F14. `seam_fn()` accepts a function of `...` (tests for both seams).
  `frmtmb_covstruct_info()` (part B) reads the registry and errors
  "no covariance structure '<name>' is registered in this session",
  adding the extension when core lists one; its docs say that `"car"`
  and `"spde"` load the extension. Autoscale (part B): both report. The
  message now names every block component on a planned column that is
  not rescaled, core's `homdiag()`, `cs()`, `ar1()` as much as a
  registered structure, because the user's remedy (rescale the column)
  is the same and a silent skip of a core structure is the same gap.
  Complex eigenvalues: left as they are (real part when the imaginary
  part is below 1e-8 of the largest modulus, else refused); a general
  complex log-determinant is not built.

### Part B

- F10. `us_t` and `diag_t` are `full_rank = TRUE`. The review's case,
  Gamma(log), weights 1.5, `(1 | gr(obs, dist = "student"))`, is now
  refused (it fitted to logLik 4.7e7, shape 8e19 on 0.74.0 too); at
  weights 1 it fits (shape below 100). The refusal now spells the term
  as written, `(1 | gr(obs, dist = "student"))`, not `us_t(1 | obs)`
  (`block_term_as_written()`). NEWS: a bug fix.
- F11. No inference: `cs_full_rank()` is `isTRUE(entry$full_rank)`, so
  a registered structure that states nothing leaves the check open, and
  `frmtmb_covstruct_info()` reports NA. Docs say a package should state
  `full_rank = TRUE`. Tests: a stated full-rank iid block is refused
  exactly as `(1 | obs)`; the same block unstated is not; the review's
  bounded-support latent (scaled beta(2, 2) on (-0.3, 0.3)) fits
  unstated; the rank-one block fits. No subset was argued safe for
  inference: a positive definite covariance says nothing of support.

### Runs

Every `lib:` line of the 445 per-file logs names frmtmb from
`wt-sarfcor-lib`, and frmtmb.sample and frmtmb.spatial from there too.

One file per process, NOT_CRAN and FRMTMB_BRMS_FIT_TESTS, all nine
suites, 8 at a time, 16 min 29 s (`dev/release/run-par.sh sarfcor-fr1`,
`sarfcor-fr1.log`):

    frmtmb files 276 pass 21494 fail 0 err 0 skip 15 warn 0
    frmtmb.coupling files 11 pass 551 fail 0 err 0 skip 5 warn 0
    frmtmb.eam files 32 pass 1844 fail 0 err 0 skip 3 warn 0
    frmtmb.latent files 12 pass 679 fail 0 err 0 skip 2 warn 0
    frmtmb.learn files 15 pass 501 fail 0 err 0 skip 2 warn 0
    frmtmb.ode files 11 pass 549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample files 66 pass 3312 fail 0 err 0 skip 1 warn 0
    frmtmb.spatial files 8 pass 499 fail 0 err 0 skip 0 warn 0
    frmtmb.survival files 14 pass 405 fail 0 err 0 skip 0 warn 0

Per file, gated, before this round (`sarfcor-g1.log`) and after:
test-fcor.R 48 -> 53, test-autocor-seam.R 36 -> 57, test-guardrail1.R
35 -> 37, test-weights-unbounded.R 65 -> 75, test-autoscale.R 49 -> 50,
frmtmb.spatial test-sar.R 53 -> 64 and test-guardrail.R 104 -> 148,
frmtmb.sample test-autocor-seam-draws.R 5 -> 14. test-compat.R 767 and
test-conditions.R 150 unchanged. No expectation was dropped.

One session per suite, three orders, NOT_CRAN, core, frmtmb.spatial
and frmtmb.sample (`dev/release/run-sessions.sh sarfcor-fr1s 8`, 24 min
20 s, `sarfcor-fr1s.log`): fail 0, err 0, warn 0 in all nine:

    SESSION frmtmb forward  files=266/276 pass=19058 skip=204
    SESSION frmtmb reverse  files=266/276 pass=19158 skip=204
    SESSION frmtmb preload  files=266/276 pass=19251 skip=206
    SESSION frmtmb.spatial  8/8 pass=428 skip=5, all three orders
    SESSION frmtmb.sample   61/66 pass=3042 skip=9, all three orders

33 LEAK lines, every one tagged `[with load of ...]`, 0 untagged.

OpenBLAS 0.3.26 with LAPACK, 4 threads (`dev/ciharden-openblas.sh
0.3.26 lapack`, `dev/ciharden-scan.sh`, every gate on) against the
reference BLAS, on 17 files (`dev/sarfcor-jobs-ob.lst`: the earlier 14
plus frmtmb.spatial test-guardrail.R, test-brms-priors.R and
test-get-prior-route.R for the prior resolver): 17 of 17, fail 0 in
both, and `sarfcor-fr1-ob.sum` and `sarfcor-fr1-ref.sum` are identical
line for line. The new tolerances are relative and far from rounding
(the map's round trip at 1e-6 of max |u|; the uniform check at 4 ESS
standard errors). The R copy (`dev/ciharden-out`) was deleted after.

F13, the three runners for frmtmb.spatial test-guardrail.R:

- One file per process: 148 pass (both arms run; the sar fit's 11
  methods refused by name in each).
- R CMD check --as-cran: the suite FAIL 0 WARN 0 SKIP 6 PASS 360; the
  absent arm skips there as the shared helper says ("frmtmb's
  dependencies are only in a library that also holds frmtmb.spatial",
  test-guardrail.R:119), and the unloadable arm runs.
- covr: `covr::package_coverage()` could not install frmtmb.spatial on
  this machine: "lazy loading failed", with no message, in each of 13
  attempts with covr's default arguments (and in some with
  `type = "none"`), and the same for the BASE commit's frmtmb.spatial on
  rellib-r12 (so not this lane). The same R CMD INSTALL line from a
  shell installs every time; from inside an R process it failed in most
  attempts. `dev/sarfcor-covr-spatial.R` therefore does covr 3.6.5's
  two steps by hand: the install from the shell with covr's flags, then
  `tools::testInstalledPackage()` in covr's library stack (install path
  first, R_LIBS, R_LIBS_USER and R_LIBS_SITE all the whole stack,
  R_COVR = "true"): FAIL 0 WARN 0 SKIP 5 PASS 428, both arms of
  test-guardrail.R run. The R process then ended without its closing
  `proc.time()`, so testInstalledPackage() reported failure: the same
  silent exit as the lazy-loading step, after every test had passed.
  Only covr's instrumentation is missing from this runner. The CI job
  on Linux is the real covr run; this machine's silent exit at the end
  of an R child is worth its own look.

R CMD check --as-cran, one at a time (`dev/sarfcor-check.ps1`, logs
`dev/sarfcor-log/check-<pkg>.log`):

    frmtmb.sample   Status: OK       7 min 16 s
                    tests: [ FAIL 0 | WARN 0 | SKIP 17 | PASS 3028 ]
    frmtmb.spatial  Status: 1 NOTE   (V8 math rendering)  1 min 15 s
                    tests: [ FAIL 0 | WARN 0 | SKIP 6 | PASS 360 ]
    frmtmb          Status: 1 NOTE   (V8 math rendering)  30 min 19 s
                    tests: [ FAIL 0 | WARN 0 | SKIP 231 | PASS 18970 ]

## Fix round 2

The short check gave SHIP on A and B; this round made three changes,
uncommitted, all part A (synced to the stage-A copy; the shared NEWS
files edited in both). `A.diff` (55 files) and `B.diff` (16) were
rebuilt, apply on 18b536fc in order, and reproduce the worktree outside
`dev/`.

### (1) newdata after the fit drops rows

What brms 2.23.0 does, measured (`dev/sarfcor-fr2-rows.R`; 36 rows on a
6 x 6 lattice, `x` missing in rows 5 and 20, fcor(V) and lag sar(Ws),
brms's draws held at frmtmb's estimates by Fixed_param;
`sarfcor-out/fr2-rows-before.txt`, `fr2-rows-after.txt`):

    newdata                          brms 2.23.0                 frmtmb before   frmtmb now
    36 rows, same 2 x missing        36 values, NA where x is    refused         same values, same NA rows
                                     missing (sar: all 36 NA)
    the 34 kept rows                 refused ("Dimensions of     34 values       refused
                                     'M' ... must be equal")
    36 rows, every x present         36 values                   refused         sar 3.5e-16, fcor 0 (relative)
    36 rows, x missing elsewhere     36 values, NA there         refused         same values, same NA rows
    36 rows, y missing               36 values                   refused         same

So brms reads newdata against the whole matrix in `data2` and drops
no row of it: a missing predictor gives `NA` (which the lag form
spreads to every row, through the solve), and the rows the fit dropped
are predicted like any other. The fit's own row drop does not apply
to newdata. frmtmb now does the same: it takes exactly `nrow(M)` rows
(`ac$n_all`, recorded when the block is built) and refuses every other
count, the kept rows alone included, because brms refuses them. When
the fit dropped rows, `autocor_seam_for_newdata()` builds the term again
over every row of M from `fit$data2` and carries its parameters over at
their natural values, because brms's `lagsar` is rho itself while
frmtmb's internal parameter maps onto an interval that depends on the
rows: through `bounds` (`autocor_seam_theta_map()`; a value outside the
new interval is refused), or unchanged when a term without `bounds`
gives the same `natural()` on both sets of rows (otherwise refused).
`autocor_seam_mean()`, `_draw()` and `_B()` read the map. When the fit
dropped no row nothing changes.

The new message names the row count, the dropped rows, the NA rule and
brms's `newdata2`. `?frmtmb_register_autocor`, `?frmtmb-autocor`
(fcor), `?sar`, both terms' compat notes, the spatial vignette and both
NEWS files say this. Tests: test-sar.R (25 rows, 2 dropped: the mean
over all 25 rows at the fitted rho against the dense solve, predict()
and simulate() over 25, all NA with a missing predictor, 23 rows
refused with "dropped 2 of those"), test-fcor.R (10 rows, 2 dropped),
and a gated brms pin in test-sar.R (36 matched to 1e-10, 34 refused by
both).

The new tests found a gap: `predict(newdata = )` checked the rows on
the rebuilt block, but its replicates drew through the fit's own frame,
so every replicate still had the 23-row block for 25 rows and stopped
("A simulation with sar(W) needs newdata ..."). `predict_simulate()`
now puts the rebuilt block in its local copy of the fit's frame before
it makes the parameter drawer, so every replicate draws through it
(`R/predict-brms.R`). `fitted()` and `simulate()` were right from the
start.

### (2) the map check reads both columns

`autocor_seam_check_map()` compared only the first column of
`unwhiten(cbind(u, rev(u)))` with the vector call. It now compares each
column with the vector call on that column (`rev(u)` for the second),
to 1e-8 of the largest value. Test: an `unwhiten()` that is right on
one column and zero on the others is refused.

### (3) the announced interval

The default prior on a bounded parameter is announced with its ends to
6 digits, `uniform(-1, 1)` where it printed
`uniform(-0.999999999999998, 0.999999999999999)`. Only the line
changes; the prior holds the exact ends. Test (frmtmb.sample
test-autocor-seam-draws.R): the stored ends equal
`frmtmb_autocor_pars()`, the shown ones have 6 digits at most and are
within a relative 1e-5.

### Runs

One file per process (lane library):

    file                          gated (before -> now)   ungated
    test-sar.R                    64 -> 75                63, skip 3
    test-fcor.R                   53 -> 57                50, skip 3
    test-autocor-seam.R           57 -> 58                58
    test-autocor-seam-draws.R     14 -> 19                19
    test-predict-newdata.R        -                       12 (unchanged)
    test-simulate-newdata.R       -                       57 (unchanged)

All fail 0, err 0, warn 0. The brms row's "Rows containing NAs were
excluded" warnings are expected (brms drops the two rows of the fit)
and are allowed by name.

One session per suite, NOT_CRAN, three orders (reverse among them),
core and frmtmb.spatial (`run-sessions.sh sarfcor-fr2s 6`, 20 min 22
s, `sarfcor-fr2s.log`): fail 0, err 0, warn 0 in all six:

    SESSION frmtmb forward  files=266/276 pass=19063 skip=204
    SESSION frmtmb reverse  files=266/276 pass=19163 skip=204
    SESSION frmtmb preload  files=266/276 pass=19256 skip=206
    SESSION frmtmb.spatial  8/8 pass=435 skip=5, all three orders

21 LEAK lines, every one tagged `[with load of ...]`, 0 untagged.

R CMD check --as-cran, core then frmtmb.spatial, one at a time
(`dev/sarfcor-check.ps1`), because the docs changed:

    frmtmb          Status: 1 NOTE   (V8 math rendering)  22 min 37 s
                    tests: [ FAIL 0 | WARN 0 | SKIP 231 | PASS 18975 ]
    frmtmb.spatial  Status: 1 NOTE   (V8 math rendering)  1 min 4 s
                    tests: [ FAIL 0 | WARN 0 | SKIP 6 | PASS 367 ]
