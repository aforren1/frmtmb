# frmtmb specification

frmtmb is an R package for formula-based regression modeling with the brms
grammar and a frequentist backend. It fits models by maximum likelihood with
the Laplace approximation for latent effects, through RTMB. It does not use
MCMC and it does not compile code at run time.

Status: design specification; v0.19 implemented (see NEWS.md for the
consolidated changelog). The original roadmap and its extensions are
complete: v0.19 latent classes with continuous random effects (growth
mixtures, via the sum-integral swap); v0.18 gp() exact and
Hilbert-space Gaussian processes, mo()/mi() interactions, group-level
latent-class mixtures; v0.17 mixture() families, mi(sdx) measurement
error, cs() category-specific effects, equalto(), rr se.fit,
registered insight methods; v0.16 rr() reduced-rank covariance
(two-space b) and mi() one-step imputation; v0.15: mo() monotonic
terms, sratio/cratio/acat ordinal families, hetar1/homcs/homtoep and
spatial exp/gau/mat covariance structures, influence()/
cooks.distance, frm_multiple() Rubin pooling, vint()/vreal()
custom-family data, CE prediction intervals and condition sets, and a
full method surface over frm_sample() draws; v0.14: se()
meta-analysis, proportion trials(), ten families, ranef(condVar),
frm_allfit, frm_simulate, profile control; v0.13 adds the
conventional method surface:
sigma/terms/weights/model.matrix/deviance/extractAIC/ngrps accessors,
lme4-convention coef(), conditional_effects() with plotting, plot(fit)
residual diagnostics, pp_check() via bayesplot, hypothesis() with
wald/profile/boot methods over an environment that includes
natural-scale sd_/cor_/sigma names, profile(), prior_summary(),
refit(), and frm_bootstrap() - parametric bootstrap over warm-started
refits). Earlier: v0.5 implemented (v0.4.x: cumulative
ordinal with threshold extra-parameters, cens()/trunc() for
CDF-carrying families; v0.5: nl = TRUE nonlinear formulas,
custom_family() with check_custom_family(), emmeans registration,
as_tmbstan(); plus predvars freezing, rank-deficiency column dropping,
and an edge-case suite mined from lme4/glmmTMB/brms issue history -
see dev/test-backlog.md). Functional regression is complete:
scalar-on-function, function-on-scalar, and function-on-function all
work through the long-format + matrix-covariate smooth representation
(validated against mgcv). Earlier: (v0.4 core: v0.1 GLMM core;
v0.2 distributional regression, single-response family roster, newdata
prediction with delta-method SEs, simulate, confint/anova/diagnose;
v0.3 s()/t2() smooths in any dpar plus homdiag/cs/ar1 covariance
structures; v0.4 multivariate models with per-response families,
gaussian rescor, brms |ID| cross-formula random-effect correlation,
matrix-response multinomial, and zi/hurdle families with full zi/hu
formulas). The main fit function is `frm()` (cf. `brm()`/brms).
Interfaces can change until v1.0.

Multivariate design notes: every linear predictor's Z spans the full b
vector (sparse), so Z column indices are b indices and |ID| merging is
only column placement; |ID|-linked terms merge into one block whose
per-level coefficient vector concatenates the component terms. That
block is unstructured by default, and a gr_cov/gr_prec Kronecker block
of the merged dimension when every linked term names the same grouping
factor and the same relationship matrix, which makes the |ID| spelling
of a multi-trait animal model the same fit as the long-format one.
rescor
standardizes per-response residuals and evaluates one constant
correlation matrix (plus a log-sigma Jacobian), which stays vectorized
under distributional sigma. Matrix responses (multinomial) use
`primary_dpars`: families whose location predictors are mu2..muK all
receive the main formula, individually overridable as dpar formulas.

Deviations from the original plan that remain true today:
`simulate()` is a numeric R-level simulator per family instead of
`obj$simulate()` (sparse Z x simref is unsupported; see the simref
note in section 4.7). Smooth terms, mo(), mi(), cs(), and gp() calls
are extracted from the formula top level before
reformulas::splitForm runs, because splitForm silently strips ANY
term whose expression mentions a special name (even inside I() or
nested calls); plain-function uses of structure names (e.g. exp(x)
in a fixed formula) are protected by alias substitution.
fitted/residuals/simulate are univariate-only, and OBS() is applied
only to univariate non-matrix responses, because RTMB registers
observations under the deparsed argument expression, so identical
expressions in a loop over responses silently collide. Everything
else once listed here (OSA residuals, RTMBdist families, gr(cov=),
ou/toep, propto-equivalent equalto, smooth edf reporting) has since
shipped; `propto` itself is spelled `gr(g, cov = A)`; gp() now
spans up to 3 dimensions (per-dimension or iso lengthscales) and the
exact form kriges at unseen positions. Remaining deferrals:
`ar()/ma()` residual autocorrelation terms.

## 1. Thesis

glmmTMB is one fixed C++ likelihood behind a formula front end. frmtmb
inverts that design: it is a **formula compiler**. It parses a brms-style
`bf()` specification into an intermediate representation (IR), assembles
design matrices and covariance-structure blocks from data, and generates an
RTMB objective closure in R. RTMB tapes that closure once per fit.

This design unlocks four things that no current frequentist package offers
together:

1. Distributional regression with random effects and smooths in every
   distributional-parameter (dpar) formula. This is the union of gamlss2
   (all-parameter regression, no marginal-likelihood random effects) and
   glmmTMB (random effects, but only `disp` and `zi` sub-models).
2. Nonlinear formulas (`nl = TRUE`) where each nonlinear parameter has the
   full predictor grammar, including random effects and smooths.
3. Multivariate joint models with a different family per response, residual
   correlation, and `|ID|` random-effect correlation across formulas
   (glmmTMB issue #1267, open).
4. Custom families as plain R functions. brms requires Stan code; glmmTMB
   requires a C++ pull request.

## 2. Backend decision: RTMB

RTMB replaces the classic per-model C++ TMB template. The evidence:

- RTMB is maintained by the TMB author (Kasper Kristensen) and has 13 CRAN
  reverse dependencies. The glmmTMB maintainers are porting glmmTMB itself
  from TMB to RTMB (GSoC 2026, sponsored by Brooks, Bolker, and Kristensen).
- After `MakeADFun`, an RTMB object is a TMB object. `sdreport`,
  `tmbprofile`, `tmbroot`, `oneStepPredict`, `checkConsistency`, and
  `tmbstan` work unchanged. Function and gradient evaluation run on the same
  TMBad tape at the same speed as compiled TMB.
- `obj$simulate()` and one-step-ahead (OSA) residuals come free through
  `OBS()` marking. Classic TMB needs hand-written `SIMULATE` blocks.
- No compiler at run time. Users install Windows binaries without Rtools.

Constraints the implementation must respect:

- (a) **Vectorize.** Elementwise `a[i] <- x` sub-assignment inside taping
  loops is more than 1000x slower than vectorized code. Observation-length
  loops are forbidden in `R/objective.R` and `R/covstruct.R`. Use `[[i,j]]`
  only to fill small dim-sized factors.
- (b) **No branching on parameters.** Branch on data at tape time instead;
  that is free. `TapeConfig(comparison = "tape")` covers indicator needs.
- (c) **Data is baked into the tape** through the closure. A refit with new
  data re-tapes. Taping costs roughly 10-20% over classic TMB when the code
  is vectorized. Accepted.
- (d) Helper functions that sub-assign advectors need
  `ADoverload("[<-")`.
- (e) RTMB has no OpenMP parallel accumulation. This is the only real gap
  against classic TMB and it is not on the critical path. The fit object is
  backend-agnostic (it wraps a MakeADFun result), so a compiled fast path
  for hot model classes stays possible later.

## 3. Dependency posture

Imports: `RTMB (>= 1.9)`, `RTMBdist (>= 1.0.6)`, `TMB`, `reformulas
(>= 0.4.4)`, `Matrix`, `mgcv`, `methods`, `stats`. No Stan, no compiled
code of our own, no brms.

- AD log-densities on the tape come from RTMB and RTMBdist.
- Numeric d/p/q/r functions for post-processing come from gamlss.dist
  (already transitive through RTMBdist) or extraDistr (Suggests).
- brms 2.23.0 has rstan in hard Imports, so brms is never an Import. It sits
  in Suggests only, guarded by `skip_if_not_installed("brms")`, for
  parser-agreement tests. Individual GPL-2 functions from brms's
  `distributions.R` (for example exgaussian, mean-parameterized Beta) can be
  vendored with attribution; frmtmb is GPL (>= 2) to keep that option.

Reused machinery, never reimplemented:

- **reformulas**: `splitForm()` separates fixed, bar, and special terms and
  is the intended extension point for registering specials; `mkReTrms()`
  gives sparse `Zt`, `Gp`, `flist`, `cnms`; `nobars`, `expandDoubleVerts`.
- **mgcv**: `smoothCon()` and `smooth2random()` convert `s()`/`t2()` into a
  fixed part plus an iid-Gaussian wiggly block; `PredictMat()` for newdata.
- **stats**: `model.frame`/`model.matrix`/`terms`/`.getXlevels` with stored
  contrasts for correct newdata prediction.
- **glmmTMB covstruct vignette**: the documented unconstrained
  parameterizations for `us`, `toep`, `cs`, `diag`, `ar1`, `ou`,
  spatial, and `rr`, re-expressed over advectors.

## 4. Architecture

### 4.1 Pipeline

```
bf(y ~ x + (x|g), sigma ~ s(z)) + gaussian()      user grammar
        |  bf(), mvbf(), +.frmtmb_formula
        v
frmtmb_formula                                     unevaluated spec
        |  parse_spec()               [R/parse.R]
        v
frmtmb_spec                                        data-free IR
        |  assemble_frame(spec, data) [R/frame.R]
        v
frmtmb_frame                                       X, Z, blocks, par template
        |  build_objective(frame)     [R/objective.R]
        v
nll(pars) closure                                  data baked in as constants
        |  RTMB::MakeADFun(nll, template, random = "b" [+ "beta" if REML])
        v
obj -- nlminb(obj$par, obj$fn, obj$gr) -- sdreport(obj)
        v
frmtmb_fit                                         S3 fit object
```

### 4.2 IR data structures

`frmtmb_formula` (output of `bf()`, mirrors the brms shape):

```r
structure(list(
  formula = y ~ x + (x | g),        # mu formula, or nl body if nl = TRUE
  pforms  = list(sigma = ~ s(z)),   # extra dpar / nlpar formulas
  nl      = FALSE,
  family  = NULL                    # attached by `+.frmtmb_formula`
), class = "frmtmb_formula")
```

`frmtmb_spec`, the data-free IR and the contract between the parser and
everything downstream:

```r
spec := list(
  responses = list(<resp> = resp_spec, ...),
  rescor    = logical,
  re_ids    = list()                # |ID| groups shared across formulas
)
resp_spec := list(
  resp_name, family,                # frmtmb_family object
  aterms = list(weights=, trials=, cens=, trunc=, se=, offset=, rate=),
  dpars  = list(mu = dpar_spec, sigma = dpar_spec, ...),
  nlspec = NULL | list(body = quote(a * exp(-b * x)), nlpars = c("a", "b"))
)
dpar_spec := list(
  name, link,
  fixed    = ~ x,                   # bars and specials removed
  re       = list(re_term, ...),
  smooth   = list(sm_term, ...),    # s()/t2() calls, unevaluated
  offset   = NULL | language,
  constant = NULL | numeric         # dpar fixed at a value
)
re_term := list(
  term, group,
  covstruct = "us" | "diag" | "homdiag" | "cs" | "toep" | "ar1" | "ou" |
              "exp" | "gau" | "mat" | "rr" | "propto",
  by = NULL, cov = NULL,            # gr(by =, cov =)
  id = NULL                         # |ID| key
)
```

`frmtmb_frame`, spec plus data:

```r
frame := list(
  spec, data,
  y = list(<resp> = numeric | matrix),
  n_obs,
  linpreds   = list(<resp>.<dpar> = linpred),   # flat, iterable
  re_blocks  = list(re_block, ...),             # flat; |ID| blocks merged
  par_template = list(beta =, b =, theta =, ...extras),
  map        = list(),                          # MakeADFun map
  terms_info = list(<resp>.<dpar> = list(terms, xlevels, contrasts, smooths))
)
linpred := list(
  resp, dpar,
  X = dense matrix,                 # fixed cols ++ smooth fixed parts
  Z = dgCMatrix | NULL,             # RE + smooth wiggly cols for this linpred
  beta_idx, b_idx,                  # ranges into flat beta / b
  offset = numeric | NULL,
  link                              # resolved link object
)
re_block := list(
  type = "covstruct" | "smooth",
  covstruct, n_levels, dim,
  b_idx, theta_idx,
  levels, cnms,                     # for ranef()/VarCorr() labels
  aux = list()                      # times, coords, known cov, smooth trans
)
```

### 4.3 Parameter layout

Flat named-list template, glmmTMB convention:

```r
list(
  beta  = numeric(p),   # all fixed effects, every resp/dpar, smooth-fixed cols
  b     = numeric(q),   # all random-effect modes, incl. smooth wiggly coefs
  theta = numeric(t)    # all covariance / smoothing parameters
)                       # + family extras (e.g. ordinal tau), never in random
```

- `random = "b"` always. REML sets `random = c("b", "beta")` (glmmTMB
  trick; valid generally under Laplace). `logLik` is then tagged REML and
  `anova` refuses comparisons across fixed-effect changes.
- Every family auxiliary parameter is a dpar with at least an intercept-only
  formula and its own link, so `shape`, `nu`, and friends live in `beta`
  uniformly. Fixing a dpar (`bf(..., sigma = 1)`) uses `map =` with
  `factor(NA)` on that beta segment.

### 4.4 The generated objective

All family, link, and covstruct dispatch resolves in plain R before or
during taping; those are data. The closure stays vectorized and branch-free
in parameters:

```r
build_objective <- function(frame) {
  lps       <- frame$linpreds
  blocks    <- frame$re_blocks
  block_nll <- lapply(blocks, \(bk) covstruct_registry[[bk$covstruct]]$nll)
  fam_lpdf  <- lapply(frame$spec$responses, \(r) r$family$lpdf)
  y         <- frame$y

  function(pars) {
    RTMB::getAll(pars, warn = FALSE)
    "[<-" <- RTMB::ADoverload("[<-")
    nll <- 0
    for (bk in blocks)                       # loop over structure, not obs
      nll <- nll - block_nll[[...]](b[bk$b_idx], theta[bk$theta_idx], bk$aux)
    dparv <- list()
    for (lp in lps) {
      eta <- drop(lp$X %*% beta[lp$beta_idx])
      if (!is.null(lp$Z))      eta <- eta + drop(lp$Z %*% b[lp$b_idx])
      if (!is.null(lp$offset)) eta <- eta + lp$offset
      dparv[[lp$key]] <- lp$link$inv(eta)
    }
    for (r in names(y)) {
      yobs <- RTMB::OBS(y[[r]])              # free simulate() and OSA
      nll  <- nll - sum(w[[r]] * fam_lpdf[[r]](yobs, dpars_for(r, dparv),
                                               frame$spec$responses[[r]]$aterms))
    }
    nll
  }
}
```

Zero-inflation and hurdle terms branch on `y == 0`, which is data, so the
mixture density is an indicator-weighted sum, branch-free in parameters.

### 4.5 Covstruct registry

`R/covstruct.R` holds one entry per structure:

```r
covstruct_registry$us <- list(
  npar  = function(dim, aux) dim + dim * (dim - 1) / 2,
  nll   = function(U, theta, aux) ...,   # vectorized over levels
  vcov  = function(theta, aux) ...,      # numeric Sigma for VarCorr()
  start = function(dim, aux) ...
)
```

Parameterizations follow glmmTMB: log-SD plus scaled-Cholesky correlation
for `us`; `dautoreg`/`dseparable` for `ar1`; `dgmrf` reserved for later
spatial structures. `us` evaluates one row-wise `RTMB::dmvnorm` on the
`n_levels x dim` reshape of its `b` segment.

### 4.6 Smooths

`mgcv::smoothCon(absorb.cons = TRUE)` plus `smooth2random(type = 2)` split
each smooth into a fixed part (into X) and a wiggly part (into Z as a
one-variance re_block). Smoothing parameters are then variance components in
`theta`, estimated automatically by Laplace ML/REML; results must match
`mgcv::gam(method = "ML"/"REML")`. Prediction uses `PredictMat` with the
stored `trans.U`/`trans.D`. Supported smooth classes are those
`smooth2random` handles: `s` and `t2`, not `te` (same restriction as brms).

### 4.7 Family system

```r
structure(list(
  family     = "nbinom2",
  dpars      = c("mu", "shape"),
  links      = list(mu = "log", shape = "log"),
  lpdf       = function(y, dpars, aterms) ...,  # vectorized advector log-density
  lcdf       = NULL,                     # optional; enables cens()/trunc()
  valid_y    = function(y, aterms) ...,  # data check at assembly time
  extra_pars = NULL,                     # e.g. ordinal thresholds
  init_dpars = list(shape = 1),
  type       = "continuous",
  specials   = c("cens", "trunc"),
  post       = list(mean_fn =, var_fn =) # numeric helpers for fitted/residuals
), class = "frmtmb_family")
```

- `lpdf` built from RTMB/RTMBdist `d*` primitives keeps `obj$simulate()`
  and OSA residuals automatic. A hand-written lpdf sets a flag and supplies
  `simulate` explicitly. Known v0.1 limitation: simulation mode replays the
  objective with `simref` objects, and the `matrix()` reshape in the `us`
  covstruct breaks there; v0.2 needs a simref-safe formulation before
  `simulate()` ships.
- `zi_*` and `hu_*` families are generated by a wrapper that adds a
  logit-link dpar and rewrites the lpdf in indicator-mixture form.
- `custom_family()` is the same constructor with a user R lpdf.
  `check_custom_family()` tapes it on toy data and runs finite-difference
  gradient checks.
- Links live in an AD-safe registry (`R/links.R`) with `fun`, `inv`, and
  `mu_eta`. Do not call `stats::make.link` internals; their C-level clamping
  is invisible to AD.

### 4.8 Estimation

1. Start values: dpar intercepts from link-transformed response moments;
   `theta` from covstruct `start()`; `b = 0`. `start =` overrides.
2. `nlminb(obj$par, obj$fn, obj$gr)` with glmmTMB-style restarts.
   `frmtmb_control()` exposes optimizer choice and iteration limits.
3. Standard errors from `sdreport(obj)`; joint precision cached lazily for
   prediction SEs. `confint(method = c("wald", "profile", "uniroot"))` maps
   to sdreport, `TMB::tmbprofile`, and `TMB::tmbroot`.
4. `diagnose(fit)`: convergence code, max absolute gradient, Hessian
   eigenvalues, NaN SEs, near-singular fits with suggested simplifications.
5. Bayesian escape hatch: `tmbstan::tmbstan(fit$obj)` works unchanged; a
   documented `as_tmbstan(fit)` wrapper arrives in v0.5.

### 4.9 Post-processing

The standard S3 method set (`predict` with newdata through stored
terms/xlevels/contrasts, `vcov`, `coef`, `fixef`, `ranef`, `VarCorr`,
`logLik`, `model.frame`, `terms`, `simulate`, `residuals` including OSA)
makes emmeans (`recover_data`/`emm_basis`, registered in `.onLoad`) and
marginaleffects work with little extra glue. `frmtmb_adreport(fit, fn)`
re-tapes with a derived quantity under `ADREPORT` — the generic delta-method
path for `conditional_effects` ribbons.

Fits keep the tape (`fit$obj`) so profile, simulate, and OSA work without a
refit. Tapes do not serialize; `strip_tape(fit)` and `retape(fit)` handle
`saveRDS` round-trips.

## 5. Grammar compatibility

frmtmb reimplements a documented subset of the brms grammar with identical
spelling: `bf()`, `lf()`, `nlf()`, `mvbf()`, `set_rescor()`, dpar formulas,
`nl = TRUE`, aterm names (`weights`, `trials`, `cens`, `trunc`, `se`,
`rate`), RE specials (`gr`, `mm`), and `s()`/`t2()`. Unsupported brms terms
fail at parse time with a clear message naming the term. brms code with
priors removed should port mechanically. Attaching brms alongside frmtmb
masks `bf()`; ordinary R masking rules apply.

## 6. Milestones

Each milestone ends green against a reference implementation.

| Version | Scope | Acceptance |
|---|---|---|
| v0.1 | gaussian/poisson/binomial (with `trials()`), mu formula only, `(1|g)` and `(1+x|g)` with `us`/`diag`, `weights()`/`offset()`, ML + REML, core methods, `dry_run` IR printing | logLik vs glmmTMB to 1e-6; fixef 1e-5; SEs 1e-4; REML vs `lmer` (sleepstudy, cbpp) |
| v0.2 | full dpar formulas incl. REs; Gamma, nbinom1/2, beta, student, tweedie, lognormal, compois; newdata predict + se.fit; residuals/simulate/confint/diagnose/anova | `sigma ~ x` vs `glmmTMB(dispformula=)` and gamlss to 1e-5; OSA residuals uniform on simulated data |
| v0.3 | `s()`/`t2()` in any dpar; `gr(by=, cov=)`; ar1/cs/toep/homdiag/ou/propto; edf reporting | vs `mgcv::gam(method="ML"/"REML")` to 1e-4; `ar1` vs glmmTMB to 1e-6 |
| v0.4 | mvbf, per-response families, rescor (gaussian/student), `\|ID\|` cross-formula RE correlation, zi/hurdle wrappers, cumulative ordinal, cens/trunc | ordinal vs `ordinal::clmm` 1e-5; zi vs glmmTMB 1e-6; censored gaussian vs `survival::survreg` |
| v0.5 | `nl = TRUE`, `custom_family()` public, emmeans/marginaleffects polish, `conditional_effects`, `as_tmbstan`, `mm()`, pkgdown + brms-migration vignette | nl growth models vs `nlme::nlme`; custom nbinom2 matches built-in to 1e-10 |

Deferred (candidates for v0.6+): `mo()`, `gp()` (HSGP or `dgmrf`),
CAR/SAR, `me()`, `cs()`, `mixture()` (not latent-Gaussian; Laplace
inappropriate; multimodal ML). Excluded: `mi()` missing-data terms.

## 7. Testing strategy

- `helper-reference.R`: `expect_loglik_equal(fit, ref, tol = 1e-6)`,
  `expect_fixef_equal`, `expect_se_equal`; every reference comparison
  wrapped in `skip_if_not_installed`.
- One file per axis: parse (IR snapshots, no fitting), frame (X/Z/index
  invariants), one per family-milestone, predict-newdata round trips,
  method snapshots.
- Simulation-recovery tests for structures with no reference implementation
  (cross-response `|ID|`).
- Gradient sanity per family: finite-difference check on a toy fit.
- Performance canary: tape-construction time bound for an n = 1e5 GLMM,
  guarding the vectorization discipline.
- CI: GitHub Actions matrix (Windows, macOS, Linux; release and devel). The
  Windows job proves the no-Rtools claim.

## 8. Risks

- **Taping regressions** from accidental elementwise ops: canary test plus
  the style rule in section 2(a).
- **smooth2random coverage**: not all mgcv smooth classes convert; enforce
  the `s`/`t2` restriction with an assembly-time error.
- **reformulas API stability**: pinned `>= 0.4.4`; the specials registration
  path is the intended extension point and glmmTMB already depends on it.
- **REML with extra_pars**: ordinal `tau` must stay out of `random`;
  encoded in family metadata.
- **rescor beyond gaussian/student**: no clean likelihood; hard error,
  documented.
