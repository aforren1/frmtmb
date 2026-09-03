# Feature compatibility

This package accepts a large grammar: about thirty families, a dozen
addition terms, twenty covariance structures, six predictor specials,
eight estimation modes, and several model structures. Most of the
combinations are meaningful. Some are refused. A few are accepted and
give a wrong answer.

The worst defects in this package’s history were all of the third kind.
`rescor = TRUE` once accepted `cens()` and then dropped the censoring
from the likelihood. `frmtmb_control(profile = TRUE)` once accepted
parameter bounds and applied them to the wrong parameters. Neither
raised an error, so neither looked like a defect.

The lesson is that a guard which does not exist looks exactly like a
guard that passed. The absence of an error is not evidence of support.
This page therefore gives every pair of features one of five states, and
it treats “nobody has checked” as a state of its own.

## The five states

| Symbol | Status | Meaning |
|:---|:---|:---|
| \+ | works | Supported and exercised by the test suite. |
| ~ | conditional | Supported, but the combination must meet a stated condition. |
| x | refused | frm() or the post-fit method stops with an error. The refusal is deliberate. |
| ! | broken | Accepted, but the result is wrong or the failure does not explain itself. Do not use the pair. |
| ? | untested | Nothing checks this pair. It may work. Treat a silent success as unverified. |

`untested` is the state that matters most. It does not mean the
combination fails. It means that if it fails, it will fail quietly,
because nothing in the test suite would notice.

## Reading the registry from code

The tables on this page are generated from the registry, so you can ask
it the same questions directly.

``` r

frm_compat("rescor", "cens()")
#>   feature_a    kind_a feature_b kind_b  status
#> 1    rescor structure    cens()  aterm refused
#>                                                                        note
#> 1 Refused. This pair was once accepted with the censoring silently dropped.

frm_compat("quadrature", "trunc()")$status
#> [1] "refused"
```

Give one feature to get every pair it takes part in, and `status` to
filter:

``` r

head(frm_compat("se()", status = "works"), 3)
#>   feature_a kind_a feature_b    kind_b status
#> 1      se()  aterm  gaussian    family  works
#> 2      se()  aterm   student    family  works
#> 3      se()  aterm        us covstruct  works
#>                                                                                      note
#> 1 Supported. se(x, sigma = TRUE) keeps the estimated residual SD alongside the known one.
#> 2 Supported. se(x, sigma = TRUE) keeps the estimated residual SD alongside the known one.
#> 3       Addition terms change the likelihood, covariance structures change the predictor.
```

[`frm_compat_features()`](https://aforren1.github.io/frmtmb/reference/frm_compat_features.md)
lists the vocabulary.
[`frm_compat_rules()`](https://aforren1.github.io/frmtmb/reference/frm_compat_rules.md)
shows the underlying rules, which are written per kind or per feature
group rather than per pair.

## Combinations to avoid

These pairs are accepted and should not be. Each one is a defect, not a
design decision.

## Response distributions and addition terms

Addition terms are the most family-sensitive part of the grammar.
`cens()` and [`trunc()`](https://rdrr.io/r/base/Round.html) need a
family that supplies a cumulative distribution function, and `cens()`
additionally refuses discrete responses. `se()` and `mi()` need a
gaussian or student model.

|  | weights() | trials() | cens() | trunc() | se() | mi() | vint() | vreal() |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| gaussian | \+ | ? | \+ | \+ | \+ | \+ | ? | ? |
| student | \+ | ? | x | x | \+ | \+ | ? | ? |
| lognormal | \+ | ? | \+ | \+ | x | x | ? | ? |
| shifted_lognormal | \+ | ? | x | x | x | x | ? | ? |
| skew_normal | \+ | ? | x | x | x | x | ? | ? |
| exgaussian | \+ | ? | x | x | x | x | ? | ? |
| asym_laplace | \+ | ? | x | x | x | x | ? | ? |
| Gamma | \+ | ? | x | x | x | x | ? | ? |
| weibull | \+ | ? | \+ | \+ | x | x | ? | ? |
| exponential | \+ | ? | \+ | \+ | x | x | ? | ? |
| inverse.gaussian | \+ | ? | \+ | \+ | x | x | ? | ? |
| beta | \+ | ? | x | x | x | x | ? | ? |
| tweedie | \+ | ? | x | x | x | x | ? | ? |
| poisson | \+ | ? | x | ~ | x | x | ? | ? |
| negbinomial | \+ | ? | x | x | x | x | ? | ? |
| nbinom1 | \+ | ? | x | x | x | x | ? | ? |
| geometric | \+ | ? | x | x | x | x | ? | ? |
| compois | \+ | ? | x | x | x | x | ? | ? |
| binomial | \+ | \+ | x | x | x | x | ? | ? |
| bernoulli | \+ | ? | x | x | x | x | ? | ? |
| beta_binomial | \+ | \+ | x | x | x | x | ? | ? |
| multinomial | \+ | ~ | x | x | x | x | ? | ? |
| zero_inflated_poisson | \+ | ? | x | x | x | x | ? | ? |
| zero_inflated_negbinomial | \+ | ? | x | x | x | x | ? | ? |
| zero_inflated_binomial | \+ | \+ | x | x | x | x | ? | ? |
| zero_inflated_beta | \+ | ? | x | x | x | x | ? | ? |
| hurdle_poisson | \+ | ? | x | x | x | x | ? | ? |
| hurdle_gamma | \+ | ? | x | x | x | x | ? | ? |
| hurdle_lognormal | \+ | ? | x | x | x | x | ? | ? |
| cumulative | \+ | ? | x | x | x | x | ? | ? |
| sratio | \+ | ? | x | x | x | x | ? | ? |
| cratio | \+ | ? | x | x | x | x | ? | ? |
| acat | \+ | ? | x | x | x | x | ? | ? |
| categorical | \+ | ? | x | x | x | x | ? | ? |
| von_mises | \+ | ? | x | x | x | x | ? | ? |
| cox | \+ | ? | \+ | ~ | x | x | ? | ? |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | poisson + trunc() | Discrete truncation needs a lower bound of at least 1; trunc(lb = 0) is not truncation and is refused. |
| ~ | multinomial + trials() | Required. The row sums of the response matrix must equal the trials. |
| ~ | cox + trunc() | Runs through the same log-CDF the censoring uses, so trunc(lb = ) is delayed entry. The truncation bound is evaluated against the SAME spline basis the response is, which means a bound outside the boundary knots is clamped to them rather than extrapolated. Untested against an external left-truncated reference. |
| x | student + cens(); shifted_lognormal + cens(); skew_normal + cens(); exgaussian + cens(); asym_laplace + cens(); and 8 more | Refused: cens() needs a family with an AD log-CDF. |
| x | student + trunc(); shifted_lognormal + trunc(); skew_normal + trunc(); exgaussian + trunc(); asym_laplace + trunc(); and 19 more | Refused: trunc() needs a family with an AD log-CDF. |
| x | lognormal + se(); shifted_lognormal + se(); skew_normal + se(); exgaussian + se(); asym_laplace + se(); and 29 more | Refused: known standard errors are added to the residual variance, which only the gaussian and student families have. |
| x | lognormal + mi(); shifted_lognormal + mi(); skew_normal + mi(); exgaussian + mi(); asym_laplace + mi(); and 29 more | Refused: an imputation model must be gaussian or student. |
| x | poisson + cens() | Refused: censoring is not supported for discrete families yet, even though poisson carries a CDF. |
| x | negbinomial + cens(); nbinom1 + cens(); geometric + cens(); compois + cens(); binomial + cens(); and 6 more | Refused: censoring is not supported for discrete families yet. |
| x | cumulative + cens(); cumulative + trunc(); sratio + cens(); sratio + trunc(); cratio + cens(); and 3 more | Refused: ordinal families carry no AD log-CDF over the response scale. |
| x | categorical + cens() | Refused: a nominal response carries no order, so it has no CDF for a censored row to contribute. |
| x | categorical + trunc() | Refused for the same reason: no order, no CDF, no truncation window. |

## Response distributions and post-fit methods

What you can do with a fit depends on what the family supplies. A family
without a simulator cannot be simulated from. A family without a mean
function has no [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)
value.

|  | fitted | predict | simulate | residuals | residuals_osa | emmeans | frm_sample | confint_profile | hypothesis_profile |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| gaussian | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| student | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| lognormal | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| shifted_lognormal | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| skew_normal | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| exgaussian | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| asym_laplace | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| Gamma | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| weibull | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| exponential | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| inverse.gaussian | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| beta | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| tweedie | ~ | ~ | x | ~ | ~ | ~ | ~ | ~ | ~ |
| poisson | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| negbinomial | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| nbinom1 | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| geometric | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| compois | ~ | ~ | x | ~ | ~ | ~ | ~ | ~ | ~ |
| binomial | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| bernoulli | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| beta_binomial | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| multinomial | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| zero_inflated_poisson | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| zero_inflated_negbinomial | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| zero_inflated_binomial | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| zero_inflated_beta | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| hurdle_poisson | ~ | ~ | x | ~ | ~ | ~ | ~ | ~ | ~ |
| hurdle_gamma | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| hurdle_lognormal | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| cumulative | ~ | ~ | \+ | ~ | \+ | ~ | ~ | ~ | ~ |
| sratio | ~ | ~ | \+ | ~ | \+ | ~ | ~ | ~ | ~ |
| cratio | ~ | ~ | \+ | ~ | \+ | ~ | ~ | ~ | ~ |
| acat | ~ | ~ | \+ | ~ | \+ | ~ | ~ | ~ | ~ |
| categorical | ~ | ~ | \+ | x | x | ~ | ~ | ~ | ~ |
| von_mises | ~ | ~ | \+ | ~ | x | ~ | ~ | ~ | ~ |
| cox | x | ~ | x | ~ | ~ | ~ | ~ | ~ | ~ |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | gaussian + fitted; student + fitted; lognormal + fitted; shifted_lognormal + fitted; skew_normal + fitted; and 24 more | Needs a family with a mean function. |
| ~ | gaussian + predict; student + predict; lognormal + predict; shifted_lognormal + predict; skew_normal + predict; and 25 more | Rank-deficient designs drop aliased columns at fit time. New data that is not estimable from the retained columns predicts NA and warns. |
| ~ | gaussian + residuals; gaussian + confint_profile; gaussian + hypothesis_profile; student + residuals; student + confint_profile; and 98 more | Depends on which post-fit ingredients the family supplies (CDF, simulator, variance function). |
| ~ | gaussian + residuals_osa; student + residuals_osa; lognormal + residuals_osa; shifted_lognormal + residuals_osa; skew_normal + residuals_osa; and 25 more | One-step-ahead residuals need the family to register its observation through OBS(). |
| ~ | gaussian + emmeans; student + emmeans; lognormal + emmeans; shifted_lognormal + emmeans; skew_normal + emmeans; and 27 more | Univariate fits only, and the mu predictor must be linear. |
| ~ | gaussian + frm_sample; student + frm_sample; lognormal + frm_sample; shifted_lognormal + frm_sample; skew_normal + frm_sample; and 31 more | Chains start jittered around the fitted mode. Use init_jitter to widen the spread, or init = “random” for a multimodal posterior. |
| ~ | cumulative + fitted; sratio + fitted; cratio + fitted; acat + fitted | Returns the same n x K matrix of category probabilities predict(type = “response”) returns, not a vector: an ordinal response has no mean, so the modelled response is the category distribution. The predict(type = “response”) == fitted() identity holds. The latent linear predictor is predict(fit, type = “link”), which is also what emmeans and insight see. |
| ~ | cumulative + predict; sratio + predict; cratio + predict; acat + predict | type = “response” returns an n x K matrix of category probabilities (rows summing to 1, columns named by the response’s own levels), not a vector: an ordinal response has no mean. It equals fitted(). cs() terms are honored and re-evaluated on newdata. type = “link” gives the latent predictor, which is where se.fit is available; se.fit is refused on the response scale. |
| ~ | cumulative + residuals; sratio + residuals; cratio + residuals; acat + residuals | “response” and “pearson” score the categories by the same codes 1..K the likelihood uses: y - sum_k k \* P(y = k), standardized by that distribution’s own sd. That is a residual on a SCORE, not on the ordinal scale; “osa” and dharma_residuals() use only the order. “deviance” is refused, as for every family without a standard unit deviance. |
| ~ | cumulative + emmeans; sratio + emmeans; cratio + emmeans; acat + emmeans | Works on the LATENT linear predictor, emmeans’s mode = “latent” convention for clm-like models: the intercept is dropped there (the K-1 thresholds take its place), so contrasts are on the latent scale and absolute means carry no threshold offset. For category probabilities use predict(fit, type = “response”) or conditional_effects(), which are on a different scale from these means. |
| ~ | categorical + fitted | Returns the n x K matrix of category probabilities, not a vector: a nominal response has no mean, so the modelled response is the category distribution. The predict(type = “response”) == fitted() identity holds. |
| ~ | categorical + predict | type = “response” returns an n x K matrix of category probabilities, columns named by the response’s own levels and rows summing to 1, exactly as for the ordinal families; it equals fitted(). se.fit is refused there. Each category’s latent predictor is predict(type = “link”, dpar = “mu”), which is where se.fit works. |
| ~ | von_mises + fitted | Returns the mean DIRECTION in radians on (-pi, pi\], which is what brms’s posterior_epred() reports for this family; a circular response has no arithmetic mean. |
| ~ | cox + predict | type = “response” and fitted() are refused: a survival time has no mean the censored rows identify. type = “link” gives the log hazard ratio, and cox_baseline() the fitted baseline weights. |
| x | tweedie + simulate; compois + simulate; hurdle_poisson + simulate | Refused: this family has no simulator yet. |
| x | categorical + residuals | Refused: the categories carry no order, so no residual has a scale to live on. Compare fitted(fit), the n x K category probabilities, against the observed categories instead. |
| x | categorical + residuals_osa | Refused with residuals() as a whole: a one-step-ahead residual is a CDF value, and a nominal response has no CDF. |
| x | von_mises + residuals_osa | Refused upstream: RTMBdist::dvm() rejects the osa observation object, because a wrapped support has no one-step CDF on the line. |
| x | cox + fitted | Refused: a survival time has no mean on the response scale here. Use predict(type = “link”) for the log hazard ratio. |
| x | cox + simulate | Refused: drawing a survival time means inverting the cumulative baseline hazard, which this family does not carry a quantile function for. simulate(), posterior_predict() and frm_simulate() each say so in their own words and then repeat the family’s reason. |

## Estimation modes

The modes constrain each other more than they constrain anything else.
`REML`, `quadrature`, and `profile` all move the fixed effects out of
the outer optimization problem, and no two of them can do that at once.
Bounds name outer parameters, so bounds on fixed effects are refused
under `REML` and under `profile`.

`autoscale`, `sparse_x`, and `verbose` are different in kind. They
change the arithmetic or the output, not the model.

|            | REML | quadrature | profile | autoscale | sparse_x | prior | bounds | verbose |
|:-----------|:----:|:----------:|:-------:|:---------:|:--------:|:-----:|:------:|:-------:|
| REML       |      |     x      |    x    |    \+     |    \+    |   ~   |   x    |   \+    |
| quadrature |  x   |            |    x    |    \+     |    \+    |  \+   |   \+   |   \+    |
| profile    |  x   |     x      |         |    \+     |    \+    |   ~   |   x    |   \+    |
| autoscale  |  \+  |     \+     |   \+    |           |    \+    |  \+   |   \+   |   \+    |
| sparse_x   |  \+  |     \+     |   \+    |    \+     |          |  \+   |   \+   |   \+    |
| prior      |  ~   |     \+     |    ~    |    \+     |    \+    |       |   \+   |   \+    |
| bounds     |  x   |     \+     |    x    |    \+     |    \+    |  \+   |        |   \+    |
| verbose    |  \+  |     \+     |   \+    |    \+     |    \+    |  \+   |   \+   |         |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | REML + prior | Priors on fixed effects are accepted under REML while bounds on the same parameters are refused. The two surfaces are inconsistent. Priors on variance parameters behave as expected. |
| ~ | profile + prior | Priors on the profiled fixed effects are accepted, while bounds on the same parameters are refused. Treat priors under profile with care. |
| x | REML + quadrature | Refused: quadrature already marginalizes the random effects, so there is no inner problem left for REML to integrate. |
| x | REML + profile | Refused: profile = TRUE and REML both remove the fixed effects from the outer problem. |
| x | REML + bounds | Refused: under REML the fixed effects leave the outer parameter vector, so bounds naming them are rejected as unknown parameters. |
| x | quadrature + profile | Refused: profile = TRUE cannot be combined with quadrature = TRUE. |
| x | profile + bounds | Refused: with the fixed effects profiled out they leave the outer parameter vector, so bounds naming them are rejected. This pair was once accepted and the bounds were then applied to the wrong parameters. |

### Modes and addition terms

|           | REML | quadrature | profile | autoscale | sparse_x | prior | bounds | verbose |
|:----------|:----:|:----------:|:-------:|:---------:|:--------:|:-----:|:------:|:-------:|
| weights() |  ?   |     \+     |   \+    |    \+     |    \+    |   ?   |   ?    |   \+    |
| trials()  |  ?   |     ?      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| cens()    |  \+  |     \+     |   \+    |    \+     |    \+    |   ?   |   ?    |   \+    |
| trunc()   |  \+  |     x      |   \+    |    \+     |    \+    |   ?   |   ?    |   \+    |
| se()      |  \+  |     \+     |   \+    |    \+     |    \+    |   ?   |   ?    |   \+    |
| mi()      |  \+  |     x      |   \+    |    \+     |    \+    |   ?   |   ?    |   \+    |
| vint()    |  ?   |     ?      |    ?    |    \+     |    \+    |   ?   |   ?    |    ?    |
| vreal()   |  ?   |     ?      |    ?    |    \+     |    \+    |   ?   |   ?    |    ?    |

| Status | Pairs | Note |
|:---|:---|:---|
| x | trunc() + quadrature | Refused: the truncation normalizer is log(F(ub) - F(lb)) over plain CDFs, and the Gauss-Kronrod nodes reach random-effect values where that difference underflows to exactly zero while the density is still representable. The integrand is then +Inf and the marginalized objective is -Inf, at the Laplace optimum as well as at the starting values, so the fit used to report logLik = +Inf as converged. Laplace stays near the mode and is unaffected: use quadrature = FALSE, REML, or profile for truncated responses. |
| x | mi() + quadrature | Refused: the imputed values are themselves random effects that quadrature cannot marginalize. |

### Modes and covariance structures

Only `quadrature` reads the covariance structures. It marginalizes one
scalar random intercept at a time, so every block in the model must be a
one-dimensional `us`, `diag`, or `homdiag` term.

|         | REML | quadrature | profile | autoscale | sparse_x | prior | bounds | verbose |
|:--------|:----:|:----------:|:-------:|:---------:|:--------:|:-----:|:------:|:-------:|
| us      |  \+  |     ~      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| diag    |  \+  |     ~      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| homdiag |  \+  |     ~      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| cs      |  \+  |     x      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| ar1     |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| hetar1  |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| ou      |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| toep    |  \+  |     x      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| homtoep |  \+  |     x      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| homcs   |  \+  |     x      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| exp     |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| gau     |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| mat     |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| rr      |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| equalto |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| gr_cov  |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| gr_prec |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| smooth  |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| gp      |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| hsgp    |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| car     |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| spde    |  \+  |     x      |    ?    |    \+     |    \+    |   ~   |   ~    |    ~    |
| us_t    |  \+  |     ~      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| diag_t  |  \+  |     ~      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | us + quadrature; diag + quadrature; homdiag + quadrature | Allowed only when the block is one-dimensional, that is a scalar random intercept. Correlated slopes are refused. Several such blocks are fine, nested ones included: (1 \| ga/gb) becomes an iterated one-dimensional integral. |
| ~ | ar1 + prior; ar1 + bounds; ar1 + verbose | ar1() reads the level order. Over gapped integer levels it keeps positional semantics and warns. Use ou() over num_factor() for irregular spacing. |
| ~ | hetar1 + prior; hetar1 + bounds; hetar1 + verbose | hetar1() reads the level order. Over gapped integer levels it keeps positional semantics and warns. |
| ~ | ou + prior; ou + bounds; ou + verbose | ou() is the irregular-spacing structure. Build the levels with num_factor() so the metric distance is recoverable. |
| ~ | exp + prior; exp + bounds; exp + verbose; gau + prior; gau + bounds; and 4 more | Spatial structures need coordinates built with num_factor(x, y). |
| ~ | rr + prior; rr + bounds; rr + verbose | rr() gives a reduced-rank block; the rank d must not exceed the term dimension. |
| ~ | equalto + prior; equalto + bounds; equalto + verbose | equalto(x + 0 \| g, V) fixes the term covariance to V, which must be square and match the term dimension, and belongs in data2 = list(V = V). |
| ~ | gr_cov + prior; gr_cov + bounds; gr_cov + verbose | gr(cov = A) accepts correlated slopes; the block covariance is the Kronecker product of A and the term covariance. A needs dimnames covering every grouping level, and belongs in data2 = list(A = A). Terms sharing an \|ID\| key over the same factor and the same A merge into one such block, which is the same model as writing the traits long with a single gr() term. |
| ~ | gr_prec + prior; gr_prec + bounds; gr_prec + verbose | gr(prec = Q) takes correlated slopes; the block precision is the Kronecker product of Q and the inverse term covariance, so it stays as sparse as Q. Q needs dimnames covering every grouping level, and belongs in data2 = list(Q = Q). Terms sharing an \|ID\| key over the same factor and the same Q merge into one such block. |
| ~ | smooth + prior; smooth + bounds; smooth + verbose | smooth is the internal structure behind s() and t2(); it is not written directly in a formula. |
| ~ | gp + prior; gp + bounds; gp + verbose; hsgp + prior; hsgp + bounds; hsgp + verbose | gp() and hsgp() are predictor specials, not bar terms. Write gp(x), not (gp(x) \| g). |
| ~ | car + prior; car + bounds; car + verbose | car(M, gr = g, type = ) is a predictor special, not a bar term. M is a symmetric adjacency matrix with dimnames (rownames, colnames, or both, which then have to agree) covering every location; entries must be present and non-negative, and non-zero weights are binarized. type = “escar” is the proper CAR, “icar”/“esicar” the intrinsic one under a soft sum-to-zero constraint (con_sd), “bym2” the scaled mixture; escar needs every location to have a neighbor. M belongs in data2 = list(M = M). |
| ~ | spde + prior; spde + bounds; spde + verbose | spde(fem, gr = node) is a predictor special taking a mesh’s finite-element matrices (M0/M1/M2 or c0/g1/g2) as fixed data; gr maps observations onto mesh nodes BY ROW NUMBER (whole numbers in 1..nrow(M0), as integers or as a factor/character spelling of them), because the matrices carry no dimnames to match labels against. Unobserved nodes keep their column; a general projector matrix is not supported yet. The matrices belong in data2 = list(fem = fem). |
| ~ | us_t + quadrature; diag_t + quadrature | Allowed for one-dimensional blocks, and RECOMMENDED there: the Gauss-Kronrod rule marginalizes a scalar t latent EXACTLY, where the Laplace default is approximate. Verified against adaptive Gauss-Hermite quadrature to 1e-6 in the log-likelihood and in every estimate. Correlated slopes are refused, as they are for a gaussian block. |
| x | cs + quadrature; ar1 + quadrature; hetar1 + quadrature; ou + quadrature; toep + quadrature; and 14 more | Refused: quadrature marginalizes one scalar random intercept at a time. Every block must be a dimension-1 us, diag, or homdiag term. |

## Model structures

|             | weights() | trials() | cens() | trunc() | se() | mi() | vint() | vreal() |
|:------------|:---------:|:--------:|:------:|:-------:|:----:|:----:|:------:|:-------:|
| mvbf        |    \+     |    \+    |   \+   |   \+    |  \+  |  \+  |   ?    |    ?    |
| rescor      |     x     |    ?     |   x    |    x    |  x   |  x   |   ?    |    ?    |
| \|ID\|      |     ?     |    ?     |   ?    |    ?    |  ?   |  ?   |   ?    |    ?    |
| nl          |     ?     |    ?     |   \+   |    ?    |  ?   |  ?   |   ?    |    ?    |
| mixture     |    \+     |    ?     |   x    |    x    |  x   |  x   |   ?    |    ?    |
| mixture_mvn |     ?     |    ?     |   x    |    x    |  ?   |  ?   |   ?    |    ?    |
| hmm         |     x     |    \+    |   x    |    x    |  x   |  x   |   ?    |    ?    |
| lca         |     x     |    x     |   x    |    x    |  x   |  x   |   x    |    x    |

| Status | Pairs | Note |
|:---|:---|:---|
| x | rescor + weights(); rescor + se() | Refused. |
| x | hmm + weights() | Refused: weights() scales a per-row log-density, and an HMM’s contribution is per sequence. |
| x | lca + weights(); lca + trials(); lca + cens(); lca + trunc(); lca + se(); and 3 more | Refused: an lca() response is a matrix of item codes with no per-row weight, censoring window, trial count or known standard error to attach an addition term to. One message covers the whole set. |
| x | rescor + cens() | Refused. This pair was once accepted with the censoring silently dropped. |
| x | mixture + cens() | Refused: mixture() has no CDF, so the CDF guard rejects cens(). |
| x | mixture_mvn + cens(); mixture_mvn + trunc() | Refused: mixture_mvn() has no CDF. |
| x | hmm + cens() | Refused: censoring replaces a row’s density with a CDF difference, which the forward recursion has no room for. |
| x | rescor + trunc() | Refused. This pair was once accepted with the truncation silently dropped. |
| x | mixture + trunc() | Refused: mixture() has no CDF. |
| x | hmm + trunc() | Refused for the same reason as cens(). |
| x | mixture + se() | Refused: se() is supported for gaussian and student families only. |
| x | hmm + se() | Refused: a known measurement SD is a per-row modification of the emission density; combining it with the state sum is not implemented. |
| x | rescor + mi() | Refused: mi() cannot be combined with rescor = TRUE. |
| x | mixture + mi() | Refused: mi() on the mixture response is not supported. |
| x | hmm + mi() | Refused: an NA response is handled by hmm() itself - the row is kept and its emission masked, so the chain keeps its length - and needs no latent parameter. |

### Structures and post-fit methods

Several post-fit methods are univariate-only. A multivariate fit
predicts one response at a time, through `predict(fit, resp = )`, and
refuses [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html). The inference
surface is not univariate-only:
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
and
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
read the outer parameter vector, which a multivariate fit has like any
other.

|  | fitted | predict | simulate | residuals | residuals_osa | emmeans | frm_sample | confint_profile | hypothesis_profile |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| mvbf | x | \+ | x | x | x | x | \+ | \+ | \+ |
| rescor | x | \+ | x | x | x | x | \+ | \+ | \+ |
| \|ID\| | ? | ? | ? | ? | ? | ? | ~ | ? | ? |
| nl | \+ | ~ | ? | ? | ? | x | ? | ? | ? |
| mixture | ? | ? | ~ | ? | ? | ? | ~ | ? | ? |
| mixture_mvn | ? | ? | \+ | ? | ? | ? | ~ | ? | ? |
| hmm | \+ | ~ | ~ | ~ | x | ? | \+ | ? | ? |
| lca | x | ~ | \+ | x | x | ? | \+ | \+ | ? |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | \|ID\| + frm_sample; mixture_mvn + frm_sample | Chains start jittered around the fitted mode. Use init_jitter to widen the spread, or init = “random” for a multimodal posterior. |
| ~ | nl + predict | Point predictions work. se.fit is not supported for the nonlinear predictor; request a nonlinear parameter with dpar instead. |
| ~ | mixture + simulate | Works only when every component family has a simulator. |
| ~ | mixture + frm_sample | Mixture posteriors are multimodal. Sample with init = “random” rather than the mode-anchored default. |
| ~ | hmm + predict | type = “link” and dpar = work normally, including the transition logits. type = “response” equals fitted() in sample; it is refused for newdata (state occupancy conditions on the observed responses of a whole sequence) and se.fit is refused on the response scale. |
| ~ | lca + predict | predict() returns the gating linear predictor (theta1 by default, any theta with dpar =), including on newdata. type = “response” is refused with the fitted() message. |
| ~ | hmm + simulate | A draw walks the chain forward per sequence and then emits, so it needs the emission family to have a simulator. re.form and censored = TRUE are refused. Since v0.36 the chain walk is the family’s structured simulator (fam\$sim_ctx), so posterior_predict() and frm_simulate() reach it too; posterior_predict(newdata =) is refused, because the sequence structure indexes the fitted rows. |
| ~ | hmm + residuals | type = “response” and “pearson” are computed against the occupancy-weighted mean, with the pearson scale the law-of-total-variance mixture variance. type = “deviance” is refused: there is no per-row likelihood to saturate. |
| x | mvbf + fitted; rescor + fitted | Refused: fitted() calls single_response() and stops with ‘fitted() is not supported yet for multivariate fits’. Predict one response at a time instead: predict(fit, resp = ). |
| x | mvbf + simulate; mvbf + residuals; mvbf + emmeans | Refused: the post-fit methods below are univariate-only for now. |
| x | mvbf + residuals_osa | Refused: residuals() is not supported for multivariate fits yet, one-step-ahead residuals included. |
| x | rescor + simulate | Refused: simulate() is not supported for multivariate fits yet. |
| x | rescor + residuals; rescor + residuals_osa | Refused: residuals() is not supported for multivariate fits yet. |
| x | rescor + emmeans | Refused: emmeans support is univariate-only for now. |
| x | nl + emmeans | Refused: emmeans support needs a linear mu predictor. |
| x | lca + fitted | Refused: the response is a matrix of nominal item codes, so there is no mean to fit. lca_probs() and lca_profiles() are the post-fit surface. |
| x | lca + residuals | Refused for the same reason as fitted(): no fitted mean, so no residual. |
| x | hmm + residuals_osa | Refused: one-step prediction needs the taped density of one observation given the earlier ones, and the tape holds a forward recursion over each whole sequence with no registered observation vector. |
| x | lca + residuals_osa | Refused: one observation is a subject’s whole item response pattern, not a value with a univariate conditional CDF to step through. |

## Within-group residual correlation

[`ar()`](https://rdrr.io/r/stats/ar.html), `ma()`, `arma()`, `cosy()`
and `unstr()` replace the response’s per-row density with one joint
density per group. That is what decides almost every pair below: a
family needs a real residual to correlate, and an addition term that
reshapes a per-row contribution has nothing left to reshape.

|  | gaussian | student | lognormal | shifted_lognormal | skew_normal | exgaussian | asym_laplace | Gamma | weibull | exponential | inverse.gaussian | beta | tweedie | poisson | negbinomial | nbinom1 | geometric | compois | binomial | bernoulli | beta_binomial | multinomial | zero_inflated_poisson | zero_inflated_negbinomial | zero_inflated_binomial | zero_inflated_beta | hurdle_poisson | hurdle_gamma | hurdle_lognormal | cumulative | sratio | cratio | acat | categorical | von_mises | cox |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| ar() | \+ | ~ | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x |
| ma() | \+ | ~ | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x |
| arma() | \+ | ~ | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x |
| cosy() | \+ | ~ | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x |
| unstr() | \+ | ~ | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x | x |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | ar() + student; ma() + student; arma() + student; cosy() + student; unstr() + student | The multivariate-t has one shape parameter per group, so nu must be constant; a predicted nu ~ … is refused. The density is brms’s multi_student_t with scale matrix D R D, verified against mvtnorm::dmvt exactly. |
| x | ar() + lognormal; ma() + lognormal; arma() + lognormal; cosy() + lognormal; unstr() + lognormal; and 165 more | Refused: a residual correlation needs a family with a real residual. brms accepts the same spelling for other families but fits a different model there - a latent gaussian AR process added to the linear predictor - which is spelled here as a random effect over the time factor: + ar1(factor(week) + 0 \| subj), or toep()/us() for a freer lag structure. |

|         | weights() | trials() | cens() | trunc() | se() | mi() | vint() | vreal() |
|:--------|:---------:|:--------:|:------:|:-------:|:----:|:----:|:------:|:-------:|
| ar()    |     x     |    x     |   x    |    x    |  x   |  x   |   ?    |    ?    |
| ma()    |     x     |    x     |   x    |    x    |  x   |  x   |   ?    |    ?    |
| arma()  |     x     |    x     |   x    |    x    |  x   |  x   |   ?    |    ?    |
| cosy()  |     x     |    x     |   x    |    x    |  x   |  x   |   ?    |    ?    |
| unstr() |     x     |    x     |   x    |    x    |  x   |  x   |   ?    |    ?    |

| Status | Pairs | Note |
|:---|:---|:---|
| x | ar() + weights(); ma() + weights(); arma() + weights(); cosy() + weights(); unstr() + weights(); and 25 more | Refused: the group’s density is joint, so there is no per-row contribution for a frequency weight to repeat, a censoring indicator to replace with a tail probability, a truncation bound to renormalize, or a known standard error to add to. brms refuses weights(), cens() and trunc() here with ‘Invalid addition arguments for this model’. |

|         | REML | quadrature | profile | autoscale | sparse_x | prior | bounds | verbose |
|:--------|:----:|:----------:|:-------:|:---------:|:--------:|:-----:|:------:|:-------:|
| ar()    |  \+  |     x      |   \+    |    \+     |    \+    |   ~   |   \+   |   \+    |
| ma()    |  \+  |     x      |   \+    |    \+     |    \+    |   ~   |   \+   |   \+    |
| arma()  |  \+  |     x      |   \+    |    \+     |    \+    |   ~   |   \+   |   \+    |
| cosy()  |  \+  |     x      |   \+    |    \+     |    \+    |   ~   |   \+   |   \+    |
| unstr() |  \+  |     x      |   \+    |    \+     |    \+    |   ~   |   \+   |   \+    |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | ar() + prior; ma() + prior; arma() + prior; cosy() + prior; unstr() + prior | Priors on the fixed effects and on random-effect covariance parameters work as usual. set_prior() cannot target the residual-correlation parameters themselves yet; bounds on thetaac\_\* are the available lever. |
| x | ar() + quadrature; ma() + quadrature; arma() + quadrature; cosy() + quadrature; unstr() + quadrature | Refused: the Gauss-Kronrod rule integrates a random effect against per-observation densities, and this residual is a joint density over each group. |

|         | mvbf | rescor | \|ID\| | nl  | mixture | mixture_mvn | hmm | lca |
|:--------|:----:|:------:|:------:|:---:|:-------:|:-----------:|:---:|:---:|
| ar()    |  \+  |   x    |   ?    |  x  |    x    |      x      |  x  |  ?  |
| ma()    |  \+  |   x    |   ?    |  x  |    x    |      x      |  x  |  ?  |
| arma()  |  \+  |   x    |   ?    |  x  |    x    |      x      |  x  |  ?  |
| cosy()  |  \+  |   x    |   ?    |  x  |    x    |      x      |  x  |  ?  |
| unstr() |  \+  |   x    |   ?    |  x  |    x    |      x      |  x  |  ?  |

| Status | Pairs | Note |
|:---|:---|:---|
| x | ar() + rescor; ma() + rescor; arma() + rescor; cosy() + rescor; unstr() + rescor | Refused: both describe the residual covariance - one across time, one across responses - and the joint structure is their Kronecker product, which is not implemented. brms refuses the same pair. |
| x | ar() + nl; ma() + nl; arma() + nl; cosy() + nl; unstr() + nl | Refused: a nonlinear mu is arbitrary R code, so the term would be evaluated rather than read. brms reaches this model through acformula(), which has no analog here. |
| x | ar() + mixture; ma() + mixture; arma() + mixture; cosy() + mixture; unstr() + mixture | Refused: a mixture likelihood has no single residual to correlate. The term is rejected as sitting on mu1 rather than mu, which is also how brms rejects it. |
| x | ar() + mixture_mvn; ma() + mixture_mvn; arma() + mixture_mvn; cosy() + mixture_mvn; unstr() + mixture_mvn | Refused for the same reason as mixture(). |
| x | ar() + hmm; ma() + hmm; arma() + hmm; cosy() + hmm; unstr() + hmm | Refused: a residual correlation term is rejected as sitting on mu1 rather than mu, exactly as it is for mixture(), and for the same reason - there is no single residual to correlate. |

## Predictor specials

|           | REML | quadrature | profile | autoscale | sparse_x | prior | bounds | verbose |
|:----------|:----:|:----------:|:-------:|:---------:|:--------:|:-----:|:------:|:-------:|
| s()       |  ?   |     x      |    ?    |     ?     |    \+    |   ?   |   ?    |    ?    |
| t2()      |  ?   |     x      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| mo()      |  \+  |     ?      |   \+    |     ?     |    ?     |   ?   |   ?    |    ?    |
| mi_pred() |  ?   |     ?      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |
| gp_pred() |  \+  |     x      |    ?    |     ?     |    ?     |   ?   |   ?    |    ?    |
| cs_pred() |  ?   |     ?      |    ?    |    \+     |    \+    |   ?   |   ?    |   \+    |

| Status | Pairs | Note |
|:---|:---|:---|
| x | s() + quadrature | Refused in practice: a smooth is a wide random-effect block, so the scalar-intercept guard rejects it. |
| x | t2() + quadrature | Refused in practice: a smooth is a wide random-effect block. |
| x | gp_pred() + quadrature | Refused in practice: gp() builds a wide block, so the scalar-intercept guard rejects it. |

The specials `mi()`, `gp()`, and `cs()` share a name with a covariance
structure. The registry writes the predictor forms as `mi_pred()`,
`gp_pred()`, and `cs_pred()` to keep them apart.

## Formula grammar

Five spellings have restrictions of their own. `mm()` and `mmc()` are
the multi-membership pair: `(1 | mm(g1, g2))` averages the member
levels’ effects into one design row, and `mmc()` gives that row a
member-specific covariate. Their full rule set is longer than one table
row; see `?frmtmb-multimembership`.

| Spelling | Status | Note |
|:---|:---|:---|
| x \* (1 \| g) | x | Refused: a bar term crossed with \* or : (as in x \* (1 \| g)) is not a random-effect specification (lme4#196). Write the crossing inside the bar: (x \| g). This spelling was once accepted with the crossing silently dropped. |
| (1 \| factor(x)) | \+ | Call-valued grouping factors are supported: (1 \| factor(x)) and (1 \| interaction(a, b)) both build the grouping factor from the model frame. |
| (x \|\| g) | \+ | (x \|\| g) gives uncorrelated terms. With a factor on the left, (f \|\| g) routes to diag, that is one independent effect per factor level. |
| (1 \| mm(g1, g2)) | \+ | The membership design is built before the family sees it, exactly as an ordinary grouping factor’s is. |
| (mmc(x1, x2) \| mm(g1, g2)) | ~ | mmc() only means something on the left of a multi-membership bar, where it supplies one covariate value per member. Anywhere else it is refused, including over a single-membership grouping factor. Inside an mm() term it composes like any other random-slope column. |

## Coverage

| Status      | Pairs | Share |
|:------------|------:|:------|
| works       |  1664 | 32%   |
| conditional |  1861 | 35%   |
| refused     |   764 | 15%   |
| broken      |     0 | 0%    |
| untested    |   967 | 18%   |

The untested share is the honest measure of what this registry does not
yet know. It shrinks as pairs are tested, not as the code is trusted. To
close one, add a test and change the rule in `R/compat.R`.
