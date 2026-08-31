# Feature gaps vs lme4 / glmmTMB / brms vignettes

Source: sweep of the installed vignettes of all three packages
(2026-08-31, post v0.13). Items already rejected or delivered are at
the bottom.

## Tier 1: cheap, real value — DONE in v0.14

All items below shipped in v0.14 (plus `frm_simulate()` from Tier 2
and the `anova()` duplicate-label fix). Notes: `profile = TRUE` is the
glmmTMB-style approximation (betas ride inside the Laplace step), and
`condVar` follows the TMB/glmmTMB convention, not lme4's
purely-conditional `postVar`.

- `ranef(fit, condVar = TRUE)`: conditional SDs of the modes from
  `sdreport`'s `diag.cov.random` (already computed for smooth edf).
  Both lme4 (`postVar` attribute, caterpillar plots) and glmmTMB
  (`condsd` column) provide this; we provide nothing.
- `as.data.frame()` methods for `ranef` (grp/term/condval/condsd) and
  `VarCorr` (grp/var1/var2/vcov/sdcor): the tidy shapes broom.mixed
  and plotting workflows expect.
- `se()` addition term: known per-observation sampling SD
  (meta-analysis). Pairs with `gr(g, cov = A)` for phylogenetic
  meta-analysis; brms treats this as core. Small: fixed additive
  component in the sigma of gaussian-type families.
- Binomial proportion response + `weights = size` (the glmer idiom).
  We currently require integer counts + `trials()`.
- Family quick wins: `bernoulli` (alias), `geometric` (negbinomial
  with shape fixed), `exponential` (Gamma shape 1), `weibull`,
  `shifted_lognormal` (ndt dpar), `hurdle_gamma`, `hurdle_lognormal`,
  `zero_inflated_binomial`, `zero_inflated_beta`, `asym_laplace`
  (quantile regression: fixed-quantile dpar).
- `frm_allfit()`: refit under every registered optimizer and tabulate
  logLik/fixef/theta agreement (lme4 `allFit`). Trivial given
  pluggable optimizers + `update()`.
- `MakeADFun(profile = "beta")`: profile fixed effects out of the
  inner problem (glmmTMB `glmmTMBControl(profile = TRUE)`, glmer
  `nAGQ = 0` relative). Speedup for many-coefficient models; expose
  as a `frmtmb_control()` flag after benchmarking.

## Tier 2: moderate effort - DONE in v0.15 except where noted

- ~~De novo simulation~~ DONE in v0.14 (`frm_simulate()`).
- ~~mo() monotonic effects~~ DONE (standalone additive terms;
  interactions and group-level mo() remain out).
- ~~sratio / cratio / acat~~ DONE. ~~cs()~~ DONE in v0.17
  (sratio/cratio/acat; refused for cumulative).
- ~~influence() + cooks.distance~~ DONE.
- ~~hetar1 / homcs / homtoep / exp / gau / mat~~ DONE. ~~equalto~~
  DONE in v0.17 (zero-parameter registry entry).
- ~~rr(..., d = k)~~ DONE in v0.16 (two-space b: factors in the
  parameter vector, expand_b() maps to the coefficient space the Z
  matrices span; exact vs glmmTMB; se.fit deferred).
- ~~conditional_effects(method = "predict") + data-frame
  conditions~~ DONE.
- ~~vint() / vreal()~~ DONE (Wiener diffusion showcase in
  test-v15.R).

## Missing data (assessed 2026-08-31)

- `frm_multiple()` (Rubin pooling over imputations, mice interop):
  DONE in v0.15.
- ~~In-model `mi()`~~ DONE in v0.16 for continuous predictors
  (gaussian/student imputation models; validated against the
  closed-form linear-gaussian marginal). Discrete predictors stay
  impossible (no discrete latents in the Laplace class - same wall
  Stan has). `me()` measurement error is a small follow-on (mi with
  a known-error observation model). mi() in dpar formulas, smooths
  of mi() variables, and mi() interactions remain out.

## Mixture models (assessed 2026-08-31)

- ~~Observation-level `mixture(fam1, fam2, ...)`~~ DONE in v0.17
  (logsumexp lpdf via RTMB::logspace_add; suffixed component dpars;
  mixing weights are multinomial-logit dpars with full linear
  predictors, so mixture-of-experts is free; quantile-spread mean
  inits; exact vs direct ML). Multimodality is inherent to mixture
  ML - use frm_allfit / bounds.
- ~~Group-level mixtures~~ DONE: v0.18 latent-class regression,
  v0.19 WITH continuous random effects (growth mixtures). The
  earlier "needs per-group per-class Laplace / crossed impossible"
  analysis was wrong: swapping sum and integral, the marginal is
  exactly the integral of a per-group class-mixture conditional on
  b, so ONE Laplace over b handles it and crossed designs are
  structurally fine. What remains approximate is the Laplace of the
  (non-Gaussian) mixture integrand: ~0.1 logLik bias in validation,
  parameters to 0.01; quadrature exact for univariate per-group
  integrands, approximate when class-specific intercepts couple.

## Tier 3: positioning decisions

- `frm_multiple`: fit over multiply-imputed datasets and pool by
  Rubin's rules (frequentist analog of `brm_multiple`; mice interop).
- Sparse X option (`sparse.model.matrix`) for many-level fixed
  factors; lme4's workaround is modular hacking, glmmTMB has
  `sparseX`.
- `autoscale`-style internal predictor scaling with back-transform
  (lme4 >= 1.1.37). We currently only advise rescaling in docs.
- Sandwich/robust SEs (`vcovHC`, `bread`/`estfun`): still skipped;
  glmmTMB does cluster-level scores. Revisit only on demand.

## Already ahead (no action)

- Custom families are plain R functions with AD verification vs
  brms's stanvars/Stan-code injection and glmmTMB's C++ enum-edit +
  recompile chapter (their "hacking" vignette is our `custom_family()`
  one-liner).
- Priors: brms-style `set_prior` with natural-sd-scale Jacobian;
  glmmTMB's priors vignette lists LKJ/elementwise as "planned".
- Post-hoc MCMC: `frm_sample`/`check_laplace` vs glmmTMB's DIY
  Metropolis vignette.
- OSA residuals, DHARMa, pp_check, quadrature matching GLMMadaptive,
  hypothesis wald/profile/boot, warm-start bootstrap.
- `||`, nesting/crossing, offsets, `numFactor`, dispformula-equivalent
  (dpar formulas), zi/hurdle with full formulas: all present.

## Explicitly out (locked or moot)

- ~~`me()` measurement error~~ DONE in v0.17 as `x | mi(sdx)` (the
  brms spelling since me() was folded into mi). Discrete-predictor
  mi(): impossible.
- ~~`gp()`~~ DONE in v0.18: exact (dense SE kernel + nugget) and
  Hilbert-space (`k =`) forms; 1-D only for now.
- ~~mo()/mi() interactions~~ DONE in v0.18 (two-way `:`/`*` with
  numeric terms; shared simplex per mo variable).
- OpenMP objective parallelism (RTMB limitation; benchmarks fine
  without it). brms threading, glmmTMB parallel vignettes are moot.
- Compilation management (precompile, cmdstanr backends): moot,
  nothing ever compiles.
- Group-level mixtures for crossed designs (see the mixture section:
  out of the Laplace class).
