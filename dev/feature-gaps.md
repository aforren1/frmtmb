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
- ~~sratio / cratio / acat~~ DONE. `cs()` category-specific effects
  still out.
- ~~influence() + cooks.distance~~ DONE.
- ~~hetar1 / homcs / homtoep / exp / gau / mat~~ DONE (`equalto`
  still out; map-trick candidate).
- `rr(..., d = k)` reduced-rank / GLVM: STILL OUT. Requires the b
  template segment (latent factors) to differ in size from the Z
  columns (loadings expand factors to coefficients inside the
  objective) - a two-sizes-of-b architecture change. Deferred until
  demanded.
- ~~conditional_effects(method = "predict") + data-frame
  conditions~~ DONE.
- ~~vint() / vreal()~~ DONE (Wiener diffusion showcase in
  test-v15.R).

## Missing data (assessed 2026-08-31)

- `frm_multiple()` (Rubin pooling over imputations, mice interop):
  DONE in v0.15.
- In-model `mi()` (brms one-step imputation): FEASIBLE for
  continuous predictors under Laplace, contrary to the original
  exclusion. Missing x-values are latent Gaussians given the
  multivariate formula's x-model - exactly the latent-Gaussian class
  Laplace integrates; mark them as random effects, joint likelihood,
  done. Discrete predictors stay impossible (no discrete latents in
  the Laplace class - same wall Stan has). Missing responses need no
  machinery (predict handles them). Implementation would ride on the
  mvbf machinery: `bf(y ~ mi(x)) + bf(x | mi() ~ z)`, with the
  missing entries appended to `b`. Substantial parse/frame work;
  roadmap item, not scheduled.

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

- `mi()` (excluded by decision); `me()` (folded into mi in brms).
- `gp()` HSGP (deferred; `gr(prec=)` covers the GMRF path).
- OpenMP objective parallelism (RTMB limitation; benchmarks fine
  without it). brms threading, glmmTMB parallel vignettes are moot.
- Compilation management (precompile, cmdstanr backends): moot,
  nothing ever compiles.
- `mixture()` (deferred: not latent-Gaussian).
