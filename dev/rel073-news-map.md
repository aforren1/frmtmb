| package | section (as the lane filed it) | lane | record | bullet |
|---|---|---|---|---|
| frmtmb | New features | loodraws | `dev/loodraws-findings.md` | The generics of brms's `loo_epred()`, `loo_linpred()`, |
| frmtmb | Breaking changes | cdfs | `dev/cdfs-findings.md` | A discrete response is censored and truncated as brms 2.23.0 |
| frmtmb | Breaking changes | cdfs | `dev/cdfs-findings.md` | `trunc()` refuses a response outside its window, and a windo |
| frmtmb | Breaking changes | cdfs | `dev/cdfs-findings.md` | The refusal of `cens()` and `trunc()` names the family and p |
| frmtmb | New features | cdfs | `dev/cdfs-findings.md` | `cens()` and `trunc()` take 23 more families, the ones brms  |
| frmtmb | New features | cdfs | `dev/cdfs-findings.md` | A family can give both tails of its CDF on the log scale, th |
| frmtmb | New features | cdfs | `dev/cdfs-findings.md` | `fitted()` reports the truncated mean of a `trunc()` respons |
| frmtmb | Bug fixes | cdfs | `dev/cdfs-findings.md` | `simulate()` of a `trunc()` fit stopped on a valid model whe |
| frmtmb | Bug fixes | cdfs | `dev/cdfs-findings.md` | `trunc(lb = )` without an upper bound divided by `1 - F(lb)` |
| frmtmb | Bug fixes | cdfs | `dev/cdfs-findings.md` | An interval-censored row lost its digits above the median an |
| frmtmb | Bug fixes | cdfs | `dev/cdfs-findings.md` | `trunc(lb = -Inf)` or `trunc(ub = Inf)` stopped on a NaN |
| frmtmb | Bug fixes | cdfs | `dev/cdfs-findings.md` | The `exgaussian()` density was a gaussian approximation from |
| frmtmb | Bug fixes | cdfs | `dev/cdfs-findings.md` | `inverse.gaussian()`'s CDF was `NaN` once `shape / mu` passe |
| frmtmb | Bug fixes | cdfs | `dev/cdfs-findings.md` | `poisson()` and `inverse.gaussian()` score a right-censored  |
| frmtmb | Breaking changes | nb2dm | `dev/nb2dm-findings.md` | The predictive draws of `dirichlet()` and `brmsfamily("diric |
| frmtmb | New features | nb2dm | `dev/nb2dm-findings.md` | `dirichlet_multinomial()`, brms's overdispersed multinomial. |
| frmtmb | New features | nb2dm | `dev/nb2dm-findings.md` | `brmsfamily("negbinomial2")`, brms's negative binomial in |
| frmtmb | New features | nb2dm | `dev/nb2dm-findings.md` | Both families are in `frm_compat()`, reach every post-fit pa |
| frmtmb | Breaking changes | fuzzfix | `dev/fuzzfix-findings.md` | `quadrature = TRUE` refuses a row that reads more than one |
| frmtmb | Breaking changes | fuzzfix | `dev/fuzzfix-findings.md` | `frm()` refuses a weighted likelihood that has no maximum. W |
| frmtmb | Breaking changes | fuzzfix | `dev/fuzzfix-findings.md` | asym_laplace() with a random effect of its own on every row  |
| frmtmb.sample | New features | loodraws | `dev/loodraws-findings.md` | brms's leave-one-out functions on draws: `loo_epred()`, |
| frmtmb.sample | New features | loodraws | `dev/loodraws-findings.md` | `add_criterion(x, "loo_R2")` stores the `loo_R2()` bootstrap |
| frmtmb.sample | New features | loodraws | `dev/loodraws-findings.md` | `psis()` honors `draw_ids`, which it ignored, and refuses an |
| frmtmb.sample | New features | loodraws | `dev/loodraws-findings.md` | `loo()` and the `loo_()` functions accept brms's default val |
| frmtmb.sample | New features | loodraws | `dev/loodraws-findings.md` | `resp` takes brms's spelling (`ya` for a column `y_a`) in al |
| frmtmb.sample | Breaking changes | loodraws | `dev/loodraws-findings.md` | On a subset of the draws (`ndraws`, `draw_ids`), `loo()`, `p |
| frmtmb.sample | Breaking changes | loodraws | `dev/loodraws-findings.md` | `posterior_average(weights = "stacking")` and `"pseudobma"`  |
| frmtmb.sample | Breaking changes | loodraws | `dev/loodraws-findings.md` | The `kfold()`, `reloo()`, `loo_moment_match()` and `loo_subs |
| frmtmb.sample | Divergences from brms | loodraws | `dev/loodraws-findings.md` | `newdata`, `re_formula` and the new-level arguments, given a |
| frmtmb.sample | Divergences from brms | loodraws | `dev/loodraws-findings.md` | `ndraws` and `draw_ids` choose the draws once, for the weigh |
| frmtmb.sample | Divergences from brms | loodraws | `dev/loodraws-findings.md` | An argument the `loo_()` functions have no use for is an err |
| frmtmb.sample | Divergences from brms | loodraws | `dev/loodraws-findings.md` | `loo_R2()` refuses an ordinal or categorical family, as `bay |
| frmtmb.sample | New features | nb2dm | `dev/nb2dm-findings.md` | `brmsfamily("negbinomial2")` and `dirichlet_multinomial()` f |
| frmtmb.sample | New features | nb2dm | `dev/nb2dm-findings.md` | A `dirichlet_multinomial()` fit's draws take the `multinomia |
| frmtmb.sample | Breaking changes | fuzzfix | `dev/fuzzfix-findings.md` | `frm_sample()` from a formula refuses a posterior that does  |
