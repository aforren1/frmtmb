| package | section (as filed) | source | record | bullet |
|---|---|---|---|---|
| frmtmb | Breaking changes | meme | `dev/meme-findings.md` | brms's prior classes `meanme`, `sdme` and `corme` of a `me() |
| frmtmb | Breaking changes | meme | `dev/meme-findings.md` | In a `default_prior(route = "sample")` table, a group row re |
| frmtmb | New features | meme | `dev/meme-findings.md` | `set_prior(class = "meanme")` is a density on the latent mea |
| frmtmb | New features | meme | `dev/meme-findings.md` | `default_prior()`, `validate_prior()` and `prior_summary()`  |
| frmtmb | New features | meme | `dev/meme-findings.md` | `default_prior()` and `validate_prior()` answer a `me(x)` wr |
| frmtmb | New features | meme | `dev/meme-findings.md` | With a proper prior on every class, the posterior density of |
| frmtmb | Sampling API | rescorpp | `dev/rescorpp-findings.md` | `predict_rescor_draw()` is exported. It is the joint predict |
| frmtmb | Sampling API | rescorpp | `dev/rescorpp-findings.md` | `rescor_row_loglik()` takes `resp`, the responses whose join |
| frmtmb.sample | Breaking changes | meme | `dev/meme-findings.md` | A `default_prior(route = "sample")` table passed as `prior = |
| frmtmb.sample | Breaking changes | meme | `dev/meme-findings.md` | A bound on one coefficient, `coef` with `lb` or `ub`, is ref |
| frmtmb.sample | Breaking changes | meme | `dev/meme-findings.md` | That spelling's bound now reaches its prior draws. Before, i |
| frmtmb.sample | Breaking changes | meme | `dev/meme-findings.md` | A `me()` model is sampled under brms 2.23.0's default priors |
| frmtmb.sample | Breaking changes | meme | `dev/meme-findings.md` | `get_prior(route = "sample")` and `prior_summary()` of draws |
| frmtmb.sample | New features | meme | `dev/meme-findings.md` | Priors on brms's classes `meanme`, `sdme` and `corme` apply  |
| frmtmb.sample | New features | meme | `dev/meme-findings.md` | `sample_prior = "only"` samples a `me()` model once `meanme` |
| frmtmb.sample | Bug fixes | meme | `dev/meme-findings.md` | The draws of a model with correlated `me(gr = )` terms carry |
| frmtmb.sample | Bug fixes | meme | `dev/meme-findings.md` | The default-prior message prints brms's `(flat)` on a bounde |
| frmtmb.sample | Breaking changes | rescorpp | `dev/rescorpp-findings.md` | `posterior_predict()` of a `set_rescor(TRUE)` model draws th |
| frmtmb.sample | Breaking changes | rescorpp | `dev/rescorpp-findings.md` | Without `resp`, `posterior_predict()`, `posterior_epred()`, |
| frmtmb.sample | Breaking changes | rescorpp | `dev/rescorpp-findings.md` | `resp` takes several responses, in the order you name them,  |
| frmtmb.sample | Breaking changes | rescorpp | `dev/rescorpp-findings.md` | Where a response predicts more than one value per row (the c |
| frmtmb.sample | Breaking changes | rescorpp | `dev/rescorpp-findings.md` | `log_lik(resp =)` of a `set_rescor(TRUE)` model is brms's: o |
| frmtmb.sample | Breaking changes | rescorpp | `dev/rescorpp-findings.md` | One response of a `student()` model with `set_rescor(TRUE)`  |
| frmtmb.sample | Breaking changes | consolidation | `dev/round-20261029.md` | Where `a ~ 1`, `set_prior(class = "Intercept", nlpar = "a")` |
| frmtmb.sample | Bug fixes | consolidation | `dev/round-20261029.md` | `bayes_R2()` and `loo_R2()` refuse a response named twice in |
