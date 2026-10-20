| package | section (as the lane filed it) | lane | record | bullet |
|---|---|---|---|---|
| frmtmb | Breaking changes | catmulti | `dev/catmulti-findings.md` | `multinomial()` takes brms's arguments, `(link, refcat, K)`, |
| frmtmb | Breaking changes | catmulti | `dev/catmulti-findings.md` | A `categorical()` dpar name is brms's `make_stan_names()`. D |
| frmtmb | Breaking changes | catmulti | `dev/catmulti-findings.md` | `fitted(scale = "linear")` on a family with several location |
| frmtmb | Breaking changes | catmulti | `dev/catmulti-findings.md` | `categorical()` takes `refcat` as its second argument, as br |
| frmtmb | New features | catmulti | `dev/catmulti-findings.md` | `categorical(refcat = )` and `multinomial(refcat = )` name t |
| frmtmb | New features | catmulti | `dev/catmulti-findings.md` | `fitted()` on `multinomial()` returns brms's `n x 4 x K` arr |
| frmtmb | New features | catmulti | `dev/catmulti-findings.md` | `logistic_normal()` fits inside `mvbf()`. Its correlations a |
| frmtmb | Bug fixes | catmulti | `dev/catmulti-findings.md` | `cumulative()` and `hurdle_cumulative()` gave a log probabil |
| frmtmb | Bug fixes | catmulti | `dev/catmulti-findings.md` | `conditional_effects()` on `multinomial()` draws one curve p |
| frmtmb | Bug fixes | catmulti | `dev/catmulti-findings.md` | `set_prior(resp = )` takes brms's Stan name of a response (` |
| frmtmb | Bug fixes | catmulti | `dev/catmulti-findings.md` | `confint(parm = )` takes the name of a correlation that is o |
| frmtmb | Bug fixes | catmulti | `dev/catmulti-findings.md` | `drop1()` on `categorical()`, `multinomial()`, the simplex f |
| frmtmb | Breaking changes | prilabel | `dev/prilabel-findings.md` | Nested and multi-membership groups carry brms's labels. |
| frmtmb | Breaking changes | prilabel | `dev/prilabel-findings.md` | The groups of each predictor are in brms's order, sorted by  |
| frmtmb | Breaking changes | prilabel | `dev/prilabel-findings.md` | `ranef()` and `coef()` name the levels of an interaction gro |
| frmtmb | Breaking changes | prilabel | `dev/prilabel-findings.md` | `variables()` lists every group-level standard deviation, th |
| frmtmb | Breaking changes | prilabel | `dev/prilabel-findings.md` | Student-t terms of one grouping factor share one latent scal |
| frmtmb | New features | prilabel | `dev/prilabel-findings.md` | A class `"sd"` prior can name one coefficient, as in brms. |
| frmtmb | New features | prilabel | `dev/prilabel-findings.md` | `gr(g, id = )` is brms's correlation id. |
| frmtmb | New features | prilabel | `dev/prilabel-findings.md` | `re_formula = ~ diag(x \| g)` names the term it wraps, as `~  |
| frmtmb | New features | prilabel | `dev/prilabel-findings.md` | A refused class `"cor"` prior names a `\|\|` or `gr(cor = FALS |
| frmtmb | New features | prilabel | `dev/prilabel-findings.md` | `validate_prior()` and `prior_summary()` list the class `"sd |
| frmtmb | Breaking changes | sampleprior | `dev/sampleprior-findings.md` | A prior density with a bounded support is `NaN` outside it o |
| frmtmb | Breaking changes | sampleprior | `dev/sampleprior-findings.md` | A bounded density on a parameter with no such bound of its o |
| frmtmb | New features | sampleprior | `dev/sampleprior-findings.md` | `set_prior()` takes brms's `uniform(lb, ub)` and `lognormal( |
| frmtmb | New features | sampleprior | `dev/sampleprior-findings.md` | A natural-scale prior on a logit-linked distributional param |
| frmtmb | New features | sampleprior | `dev/sampleprior-findings.md` | A prior on the `nu` that the Student-t responses of a |
| frmtmb | New features | sampleprior | `dev/sampleprior-findings.md` | The sampling API exports `ord_lp_blocks()`, the threshold bl |
| frmtmb | Internal | bin2port | `dev/bin2port-findings.md` | Bin 2 of the ported brms suite reaches 0.71.0's families. Th |
| frmtmb | Breaking changes | ndaterm | `dev/ndaterm-findings.md` | `predict(newdata =)` reads every addition term from newdata, |
| frmtmb | Breaking changes | ndaterm | `dev/ndaterm-findings.md` | A binomial-type family no longer reads a missing `trials()`  |
| frmtmb | Breaking changes | ndaterm | `dev/ndaterm-findings.md` | A custom family's declared `required_aterms` are checked aft |
| frmtmb | Breaking changes | ndaterm | `dev/ndaterm-findings.md` | An addition term coded by its levels is read on newdata thro |
| frmtmb | Breaking changes | ndaterm | `dev/ndaterm-findings.md` | `conditional_effects()` holds the variables of `thres()`, `v |
| frmtmb | Breaking changes | ndaterm | `dev/ndaterm-findings.md` | The `conditional_effects()` grid refuses a row without its |
| frmtmb | Bug fixes | ndaterm | `dev/ndaterm-findings.md` | A `mixture()` draws each component at its own rows' addition |
| frmtmb | Bug fixes | ndaterm | `dev/ndaterm-findings.md` | `predict(newdata =)` and `predict(re_formula =)` hand the si |
| frmtmb | Bug fixes | ndaterm | `dev/ndaterm-findings.md` | `conditional_effects(method = "predict")` codes a registered |
| frmtmb | Breaking changes | illse2 | `dev/illse2-findings.md` | REML and `profile = TRUE` fits on an ill-conditioned fixed d |
| frmtmb | Breaking changes | illse2 | `dev/illse2-findings.md` | A log standard deviation that the optimizer leaves near zero |
| frmtmb | Breaking changes | illse2 | `dev/illse2-findings.md` | Newton steps from the stop respect the bounds a prior sets.  |
| frmtmb | Breaking changes | illse2 | `dev/illse2-findings.md` | Nonlinear formulas take the coordinates of their Jacobian. A |
| frmtmb | Breaking changes | illse2 | `dev/illse2-findings.md` | Fits with up to 400 outer parameters take the design's |
| frmtmb | Breaking changes | illse2 | `dev/illse2-findings.md` | On an identity-link gaussian design the Hessian in the desig |
| frmtmb | New features | illse2 | `dev/illse2-findings.md` | `cooks.distance()` and the Wald rules of `anova()` on |
| frmtmb | Bug fixes | illse2 | `dev/illse2-findings.md` | When the standard-error check removes nothing, a log sd with |
| frmtmb | Documentation | illse2 | `dev/illse2-findings.md` | `?vcov.frmtmb_fit` says what the REML covariance of the coef |
| frmtmb | Breaking changes | adhess | `dev/adhess-findings.md` | Standard errors of a fit without random effects come from th |
| frmtmb | Breaking changes | adhess | `dev/adhess-findings.md` | With random effects, the outer Hessian's step follows the de |
| frmtmb | Breaking changes | adhess | `dev/adhess-findings.md` | A variance that is zero, negative or not finite is shown as  |
| frmtmb | Breaking changes | drop1raw | `dev/drop1raw-findings.md` | `update()` takes what it does not replace from the fit, as b |
| frmtmb | Breaking changes | drop1raw | `dev/drop1raw-findings.md` | A variable the fit does not hold is read only from data that |
| frmtmb | Breaking changes | drop1raw | `dev/drop1raw-findings.md` | `drop1()`, and so `step()`, and `anova()` refuse a fit whose |
| frmtmb | Breaking changes | drop1raw | `dev/drop1raw-findings.md` | `drop1()` keeps the fit's bases. A submodel keeps the knots  |
| frmtmb | Breaking changes | drop1raw | `dev/drop1raw-findings.md` | `influence()` keeps the fit's bases in every deletion refit, |
| frmtmb | Bug fixes | drop1raw | `dev/drop1raw-findings.md` | `drop1()` refits a term that reads a raw variable. `z ~ x +  |
| frmtmb | Bug fixes | drop1raw | `dev/drop1raw-findings.md` | `drop1()` refits a model with an interval `cens(cc, y2)`. Th |
| frmtmb | Bug fixes | drop1raw | `dev/drop1raw-findings.md` | `drop1()` takes the family and the control from the fit. It |
| frmtmb | Bug fixes | drop1raw | `dev/drop1raw-findings.md` | `step()` removes a term, on the fit's data. It passes its up |
| frmtmb | Bug fixes | drop1raw | `dev/drop1raw-findings.md` | `influence()` reads the raw variables by default. Its defaul |
| frmtmb | Bug fixes | drop1raw | `dev/drop1raw-findings.md` | `fitted(fit, re_formula = ~ (1 + log(x2) \| g))`, the fit's o |
| frmtmb | Bug fixes | drop1raw | `dev/drop1raw-findings.md` | A fit on a tibble with a transformed term no longer warns "S |
| frmtmb.sample | Breaking changes | catmulti | `dev/catmulti-findings.md` | `posterior_linpred()` with no `dpar` on a family with severa |
| frmtmb.sample | New features | catmulti | `dev/catmulti-findings.md` | Draws of a `logistic_normal()` response inside `mvbf()` repo |
| frmtmb.sample | Breaking changes | prilabel | `dev/prilabel-findings.md` | A default steps aside only for a prior on its own slot, as i |
| frmtmb.sample | Breaking changes | prilabel | `dev/prilabel-findings.md` | A class-wide prior does not replace a default written per |
| frmtmb.sample | Breaking changes | sampleprior | `dev/sampleprior-findings.md` | `frm_sample()` puts brms 2.23.0's default prior on every |
| frmtmb.sample | Breaking changes | sampleprior | `dev/sampleprior-findings.md` | The ordinal thresholds take brms's `student_t(3, 0, 2.5)`, c |
| frmtmb.sample | Breaking changes | sampleprior | `dev/sampleprior-findings.md` | `sigma ~ 1` takes brms's default on sigma's intercept. It go |
| frmtmb.sample | Breaking changes | sampleprior | `dev/sampleprior-findings.md` | An `unstr()` residual correlation takes brms's `lkj(1)`, cla |
| frmtmb.sample | Breaking changes | sampleprior | `dev/sampleprior-findings.md` | `frm_sample(seed = )` fixes every chain. The init jitter of |
| frmtmb.sample | Breaking changes | sampleprior | `dev/sampleprior-findings.md` | A natural-scale default carries brms's bound on its class wh |
| frmtmb.sample | New features | sampleprior | `dev/sampleprior-findings.md` | The draws carry brms's `sd_<group>__<coef>` and |
| frmtmb.sample | Breaking changes | ndaterm | `dev/ndaterm-findings.md` | `posterior_predict(newdata =)` reads every addition term fro |
| frmtmb.sample | Bug fixes | ndaterm | `dev/ndaterm-findings.md` | `posterior_predict(newdata =)` and `posterior_predict(re_for |
| frmtmb.eam | Breaking changes | ndaterm | `dev/ndaterm-findings.md` | `wiener()` refuses a boundary indicator that is not 0 or 1 a |
| frmtmb.eam | Breaking changes | ndaterm | `dev/ndaterm-findings.md` | `dec()` on newdata is coded as the fit coded it, with frmtmb |
| frmtmb.eam | Bug fixes | ndaterm | `dev/ndaterm-findings.md` | `predict(newdata =, summary = FALSE)` and frmtmb.sample's |
| frmtmb.eam | Bug fixes | ndaterm | `dev/ndaterm-findings.md` | `conditional_effects()` reads a `dec()` factor with this pac |
| frmtmb.latent | Bug fixes | ndaterm | `dev/ndaterm-findings.md` | `hmm()` draws each state at its own rows' addition terms. Th |
| frmtmb.latent | Bug fixes | drop1raw | `dev/drop1raw-findings.md` | `hmm_starts()` refits with `update()`, which takes the data, |
| frmtmb.coupling | Breaking changes | adhess | `dev/adhess-findings.md` | `frm_coherence()` refuses, by name, an interval that comes w |
