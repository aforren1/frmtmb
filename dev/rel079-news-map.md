| package | section (as filed) | source | record | bullet |
|---|---|---|---|---|
| frmtmb | Breaking changes | ncgp | `dev/ncgp-findings.md` | The exact `gp()` is non-centered: `b` of the block is z, a s |
| frmtmb | Breaking changes | ncgp | `dev/ncgp-findings.md` | The nugget is relative to `sdgp`. So the model of `y  c` is  |
| frmtmb | Breaking changes | ncgp | `dev/ncgp-findings.md` | `fit$estimates$b`, `fit$obj`, the joint precision and `start |
| frmtmb | Breaking changes | ncgp | `dev/ncgp-findings.md` | An exact `gp()` starts where the model's marginal, which int |
| frmtmb | Breaking changes | ncgp | `dev/ncgp-findings.md` | What that reaches, on rows no test had seen. On a curve with |
| frmtmb | Breaking changes | ncgp | `dev/ncgp-findings.md` | What it does not reach. Low-count Poisson read 3 of 6 short, |
| frmtmb | Breaking changes | ncgp | `dev/ncgp-findings.md` | `frm_curve_deriv()` on a fit with an exact `gp()` takes a de |
| frmtmb | Limits | ncgp | `dev/ncgp-findings.md` | `sigma` near `1e-6  sdgp`, where a nugget of `1e-12` reaches |
| frmtmb | Limits | ncgp | `dev/ncgp-findings.md` | Below about `1e-4  sdgp` nlminb also stops at false converge |
| frmtmb | Limits | ncgp | `dev/ncgp-findings.md` | On noise-free `sin(2  x)`, 60 rows, the fit is 0.003 nats sh |
| frmtmb | Costs | ncgp | `dev/ncgp-findings.md` | A gaussian fit with an exact `gp()` takes longer than the ce |
| frmtmb | Costs | ncgp | `dev/ncgp-findings.md` | On a fit with an exact `gp()`, a standard error's row under  |
| frmtmb | New features | ncgp | `dev/ncgp-findings.md` | `gp_nugget_message()`, in the sampling API: the new warning' |
| frmtmb | New features | ncgp | `dev/ncgp-findings.md` | `gp_tape_scope()` and `gp_tape_order()`, in the sampling API |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `class = "sds"` is brms's class for the standard deviation o |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `mm()` takes an `\|ID\|` key, `(x \| p \| mm(g1, g2))` or |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `gr(cov = A)` and `mm(cov = A)` refuse a matrix that is not  |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `mi()` on a response takes every family brms takes it for: ` |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | On `frm()`, a missing response value that no density reads l |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | A missing value that a density in the objective reads, a pre |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `frm()` warns, with class `frmtmb_mi_laplace_error`, when th |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `frmtmb_family(mi_chart = )` declares the chart of a family' |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `se()` on `skew_normal`, and `mi()` with `se()` on one respo |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `cens()` on a `mixture()` whose components have a CDF is the |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `gr(g, id = "p", cor = FALSE)` is accepted when every term u |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `importance` and `vcov(cluster = )` run on a fit whose missi |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | These stay refused where brms 2.23.0 builds the model, for t |
| frmtmb | New features | refuse | `dev/refuse-findings.md` | `frm_compat()` rows, each checked by running both packages:  |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | `frm(importance =)` no longer stops with "did not converge . |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | A gaussian response fitted with `importance =` now returns |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | The stop that said "ROSE" and blamed an unidentified |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | `confint(method = "profile")` on an `importance =` fit gets  |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | The warning that a correction "moved by the same amount" in |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | `frm(importance =)` starts its rounds from where the Laplace |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | `frm(importance =)` fits a response in large units. On `y x  |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | `frm(importance =)` honors a prior. It used to stop on any p |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | A constant distributional parameter, as in `bf(y ~ x + (1 \|  |
| frmtmb | Bug fixes | impsamp | `dev/impsamp-findings.md` | The numbers that `vignette("reinforcement-learning")` and |
| frmtmb | New features | latres | `dev/latres-findings.md` | BREAKING. `ar()`, `ma()`, `arma()`, `cosy()` and `unstr()` o |
| frmtmb | New features | latres | `dev/latres-findings.md` | `ar()` with `cov = FALSE` on these families is brms's moving |
| frmtmb | New features | latres | `dev/latres-findings.md` | `frm()` warns, with class `frmtmb_latres_laplace`, where the |
| frmtmb | New features | latres | `dev/latres-findings.md` | BREAKING for a call that fit in the previous development bui |
| frmtmb | New features | latres | `dev/latres-findings.md` | On a link that does not hold the mean (identity, inverse, `1 |
| frmtmb | New features | latres | `dev/latres-findings.md` | A dispersion the rows do not ask for no longer reaches an un |
| frmtmb | New features | latres | `dev/latres-findings.md` | `frm()` refuses `mi()` on the response beside a latent resid |
| frmtmb | New features | latres | `dev/latres-findings.md` | `frm()` refuses `ar()`, `ma()`, `arma()` and `cosy()` under  |
| frmtmb | New features | latres | `dev/latres-findings.md` | `frm()` refuses `asym_laplace()` with a latent residual, whe |
| frmtmb | New features | latres | `dev/latres-findings.md` | A fit at a limit of a latent residual says so: a message of  |
| frmtmb | New features | latres | `dev/latres-findings.md` | `predict()` holds `sderr` and the term's coefficients at the |
| frmtmb | New features | latres | `dev/latres-findings.md` | brms's refusals hold on these families, in brms's words: `co |
| frmtmb | New features | latres | `dev/latres-findings.md` | The methods read the latent residual as brms reads `err`: `f |
| frmtmb | New features | latres | `dev/latres-findings.md` | Above 3000 rows `fitted()`, `residuals()` and `frm_linpred(s |
| frmtmb | New features | latres | `dev/latres-findings.md` | A decision, not a defect. Under `cov = FALSE` frmtmb's `fitt |
| frmtmb | New features | latres | `dev/latres-findings.md` | `default_prior()`, `prior_summary()`, `variables()`, `summar |
| frmtmb | Breaking changes | nlfix | `dev/nlfix-findings.md` | `beta_binomial()` and `zero_inflated_beta_binomial()` have a |
| frmtmb | Breaking changes | nlfix | `dev/nlfix-findings.md` | `nl = TRUE` on `cumulative()`, `sratio()`, `cratio()`, `acat |
| frmtmb | Breaking changes | nlfix | `dev/nlfix-findings.md` | A body that is constant across the rows at a stationary poin |
| frmtmb | Breaking changes | nlfix | `dev/nlfix-findings.md` | A fit whose responses are all discrete and whose log likelih |
| frmtmb | Breaking changes | nlfix | `dev/nlfix-findings.md` | An ordinal fit with a `disc` body no longer warns that the f |
| frmtmb | Breaking changes | nlfix | `dev/nlfix-findings.md` | `frm_compat()` gives measured statuses for `nl` by family, i |
| frmtmb | New features | nlfix | `dev/nlfix-findings.md` | `bf(..., nl = TRUE, autocor = ~ term)` and `+ acformula(~ te |
| frmtmb | New features | nlfix | `dev/nlfix-findings.md` | Under `cov = FALSE`, `fitted()`, `predict()` and `logLik()`  |
| frmtmb | New features | nlfix | `dev/nlfix-findings.md` | `nl = TRUE` on `categorical()`, `multinomial()`, `dirichlet( |
| frmtmb | New features | nlfix | `dev/nlfix-findings.md` | A nonlinear body can call the functions of the Stan language |
| frmtmb | New features | nlfix | `dev/nlfix-findings.md` | `x^b` and `pow(x, b)` in a body, with a data `x` that is 0 o |
| frmtmb | New features | nlfix | `dev/nlfix-findings.md` | A nonlinear fit under a link that cannot evaluate a body at  |
| frmtmb | Bug fixes | nlfix | `dev/nlfix-findings.md` | `frm()` refuses a nonlinear body the likelihood does not rea |
| frmtmb | Bug fixes | nlfix | `dev/nlfix-findings.md` | `cs()` in a nonlinear body is refused with a sentence that s |
| frmtmb | Bug fixes | nlfix | `dev/nlfix-findings.md` | `conditional_effects()` on an ordinal body draws its Wald ba |
| frmtmb | Bug fixes | nlfix | `dev/nlfix-findings.md` | `frm_linpred()` and `fitted()` on a nonlinear body under `na |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `frm_linpred(type = "link")` with no `dpar` returns every pr |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `frm_curve()`, `frm_curve_deriv()`, `frm_curve_feature()` an |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `fitted(dpar = "mu", scale = "linear")` of a `cs()` fit retu |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `frm_linpred(dpar = "mu", type = "response")`, `fitted(dpar  |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `emmeans()` of a `cs()` fit with no `dpar` averages the thre |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `frm_curve()` under a `cov = FALSE` autocorrelation term ref |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | A new row of a `me()` term has no latent value, and brms 2.2 |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `fitted()` on such a row is the mean of the expected respons |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `fitted()`'s `Q` columns are the interval on the link's scal |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | One `me()` call is one latent value per row and draw, whiche |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `frm_curve(contrast = )` reads one draw of a `me()` value fo |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `marginaleffects`'s default `vcov` gives a warning of class |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `fitted()`'s `Q` columns of a category probability are a Wal |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `fitted()`'s `Q` columns on a `mixture_mvn()` fit are a Wald |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `conditional_effects()`'s `se__` on a lognormal or hurdle fi |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `fitted(scale = "linear")` of a distributional parameter wit |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `frm_curve_feature()` with no `resp` on a multivariate fit a |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `frm_linpred()` with no `resp` on a multivariate fit returns |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | A univariate fit does not read `resp`, as in brms: `fitted() |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `conditional_effects()` of a multivariate fit gives one pane |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `conditional_effects()` of a `me(gr = )` fit holds the group |
| frmtmb | Breaking changes | predread | `dev/predread-findings.md` | `predict(newdata = )` answers under `ar()`, `ma()`, `arma()` |
| frmtmb | New features and fixes | predread | `dev/predread-findings.md` | `pp_mixture()` takes one response of a multivariate mixture, |
| frmtmb | New features and fixes | predread | `dev/predread-findings.md` | `insight::get_predicted(predict = "link")` reads one predict |
| frmtmb | New features and fixes | predread | `dev/predread-findings.md` | `emmeans()` of a model with a `me()` term no longer leaves t |
| frmtmb | New features and fixes | predread | `dev/predread-findings.md` | `frm_curve()` no longer refuses a grid with a row that has a |
| frmtmb | New features and fixes | predread | `dev/predread-findings.md` | New sampling-API exports for frmtmb.sample: `ce_by_response( |
| frmtmb | New features and fixes | predread | `dev/predread-findings.md` | `frm_compat()`: `cs_pred()` x `frm_curve` is now "works", wh |
| frmtmb | New features | consolidation | `dev/round-20261104.md` | `?"frmtmb-gp"` is a new help topic. It states the nugget, th |
| frmtmb | New features | consolidation | `dev/round-20261104.md` | `frmtmb_control(check_mi = c("warning", "ignore", "stop"))`  |
| frmtmb | New features | consolidation | `dev/round-20261104.md` | `options(frmtmb.latres_dense_rows = )` sets the row limit ab |
| frmtmb | Bug fixes | consolidation | `dev/round-20261104.md` | `gp(x, k = )` has no nugget, but it is not a remedy by itsel |
| frmtmb | Bug fixes | consolidation | `dev/round-20261104.md` | A gap of 94 to 1127 nats between `frm()` and brms on an exac |
| frmtmb | Bug fixes | consolidation | `dev/round-20261104.md` | `fitted()` on a new row of a latent residual holds the laten |
| frmtmb | Bug fixes | consolidation | `dev/round-20261104.md` | `(1 \| p \| gr(g, cor = FALSE))` beside `(1 \| p \| g)` builds t |
| frmtmb | Bug fixes | consolidation | `dev/round-20261104.md` | `frm_curve()` at its default reads a row of `NA` as `NA`, an |
| frmtmb.sample | Breaking changes | ncgp | `dev/ncgp-findings.md` | An exact `gp()` is non-centered in core, as brms's `zgp`, at |
| frmtmb.sample | Breaking changes | ncgp | `dev/ncgp-findings.md` | `frm_sample(fit)` starts at random where the fit's noise has |
| frmtmb.sample | Breaking changes | ncgp | `dev/ncgp-findings.md` | `frm_sample(laplace = TRUE)` on an exact `gp()` takes 3.5 ti |
| frmtmb.sample | Other changes | refuse | `dev/refuse-findings.md` | brms's class `"sds"` carries a smooth's default, |
| frmtmb.sample | Other changes | refuse | `dev/refuse-findings.md` | `frm_sample()` samples `mi()` under every family core takes, |
| frmtmb.sample | Other changes | refuse | `dev/refuse-findings.md` | A missing value's draw is reported under brms's name, |
| frmtmb.sample | Other changes | refuse | `dev/refuse-findings.md` | `frm_sample(formula, laplace = TRUE)` takes the route of `fr |
| frmtmb.sample | Other changes | latres | `dev/latres-findings.md` | `frm_sample()` samples core's latent residuals: `ar()`, `ma( |
| frmtmb.sample | Other changes | latres | `dev/latres-findings.md` | `frm_sample()` fits the forms of a latent residual that `frm |
| frmtmb.sample | Other changes | latres | `dev/latres-findings.md` | `posterior_epred()`, `log_lik()` and `posterior_predict()` r |
| frmtmb.sample | Bug fixes | nlfix | `dev/nlfix-findings.md` | `frm_sample()` on `nl = TRUE` with `cumulative()`, `sratio() |
| frmtmb.sample | Bug fixes | nlfix | `dev/nlfix-findings.md` | `frm_sample()` takes `bf(..., nl = TRUE, autocor = ~ term)`, |
| frmtmb.sample | Breaking changes | predread | `dev/predread-findings.md` | A new row of a `me()` term draws a new latent value per post |
| frmtmb.sample | Breaking changes | predread | `dev/predread-findings.md` | One `me()` call is one latent value per row and draw, across |
| frmtmb.sample | Breaking changes | predread | `dev/predread-findings.md` | `posterior_linpred()` with no `dpar` on an ordinal family wi |
| frmtmb.sample | Breaking changes | predread | `dev/predread-findings.md` | `conditional_effects()` of a multivariate model gives one pa |
| frmtmb.sample | Breaking changes | predread | `dev/predread-findings.md` | A univariate model does not read `resp`, as in brms: `poster |
| frmtmb.sample | Breaking changes | predread | `dev/predread-findings.md` | `update(<draws>, prior = )` merges with the first call's pri |
| frmtmb.sample | Breaking changes | predread | `dev/predread-findings.md` | `posterior_predict(newdata = )` answers under `ar()`, `ma()` |
| frmtmb.sample | Breaking changes | predread | `dev/predread-findings.md` | A new row with no `x` or `sdx` of a `me()` term reads `NA` i |
| frmtmb.sample | Breaking changes | predread | `dev/predread-findings.md` | `conditional_effects()` of a `me(gr = )` model holds the gro |
| frmtmb.sample | New features and fixes | predread | `dev/predread-findings.md` | `bayes_R2()` and `loo_R2()` read an empty `resp` as every re |
| frmtmb.sample | New features and fixes | predread | `dev/predread-findings.md` | `pp_mixture()` of a multivariate mixture reads the one respo |
| frmtmb.sample | Bug fixes | consolidation | `dev/round-20261104.md` | The drawn route for an exact `gp()` needs long chains at low |
| frmtmb.learn | Other changes | impsamp | `dev/impsamp-findings.md` | The documented numbers for `frm(importance =)` are measured  |
| frmtmb.learn | Other changes | nlfix | `dev/nlfix-findings.md` | `frm_compat()` by `nl` is "refused" on these families: `nl = |
| frmtmb.learn | Bug fixes | consolidation | `dev/round-20261104.md` | `rlddm()` is a density, and declares `type = "continuous"`.  |
| frmtmb.eam | Other changes | nlfix | `dev/nlfix-findings.md` | `frm_compat()` by `nl`: `wiener` and `wiener_gng` are "works |
| frmtmb.survival | Other changes | nlfix | `dev/nlfix-findings.md` | `frm_compat()` by `nl`: `cox` and `royston_parmar` are "work |
