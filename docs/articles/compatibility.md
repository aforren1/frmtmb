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
#>   feature_a kind_a feature_b    kind_b  status
#> 1    cens()  aterm    rescor structure refused
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

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | poisson + trunc() | Discrete truncation needs a lower bound of at least 1; trunc(lb = 0) is not truncation and is refused. |
| ~ | multinomial + trials() | Required. The row sums of the response matrix must equal the trials. |
| x | student + cens(); shifted_lognormal + cens(); skew_normal + cens(); exgaussian + cens(); asym_laplace + cens(); and 7 more | Refused: cens() needs a family with an AD log-CDF. |
| x | student + trunc(); shifted_lognormal + trunc(); skew_normal + trunc(); exgaussian + trunc(); asym_laplace + trunc(); and 18 more | Refused: trunc() needs a family with an AD log-CDF. |
| x | lognormal + se(); shifted_lognormal + se(); skew_normal + se(); exgaussian + se(); asym_laplace + se(); and 26 more | Refused: known standard errors are added to the residual variance, which only the gaussian and student families have. |
| x | lognormal + mi(); shifted_lognormal + mi(); skew_normal + mi(); exgaussian + mi(); asym_laplace + mi(); and 26 more | Refused: an imputation model must be gaussian or student. |
| x | poisson + cens() | Refused: censoring is not supported for discrete families yet, even though poisson carries a CDF. |
| x | negbinomial + cens(); nbinom1 + cens(); geometric + cens(); compois + cens(); binomial + cens(); and 6 more | Refused: censoring is not supported for discrete families yet. |
| x | cumulative + cens(); cumulative + trunc(); sratio + cens(); sratio + trunc(); cratio + cens(); and 3 more | Refused: ordinal families carry no AD log-CDF over the response scale. |

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
| multinomial | ~ | ~ | x | ~ | ~ | ~ | ~ | ~ | ~ |
| zero_inflated_poisson | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| zero_inflated_negbinomial | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| zero_inflated_binomial | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| zero_inflated_beta | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| hurdle_poisson | ~ | ~ | x | ~ | ~ | ~ | ~ | ~ | ~ |
| hurdle_gamma | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| hurdle_lognormal | ~ | ~ | \+ | ~ | ~ | ~ | ~ | ~ | ~ |
| cumulative | ~ | ~ | x | ~ | \+ | ~ | ~ | ~ | ~ |
| sratio | ~ | ~ | x | ~ | \+ | ~ | ~ | ~ | ~ |
| cratio | ~ | ~ | x | ~ | \+ | ~ | ~ | ~ | ~ |
| acat | ~ | ~ | x | ~ | \+ | ~ | ~ | ~ | ~ |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | gaussian + fitted; student + fitted; lognormal + fitted; shifted_lognormal + fitted; skew_normal + fitted; and 28 more | Needs a family with a mean function. |
| ~ | gaussian + predict; student + predict; lognormal + predict; shifted_lognormal + predict; skew_normal + predict; and 28 more | Rank-deficient designs drop aliased columns at fit time. New data that is not estimable from the retained columns predicts NA and warns. |
| ~ | gaussian + residuals; gaussian + confint_profile; gaussian + hypothesis_profile; student + residuals; student + confint_profile; and 94 more | Depends on which post-fit ingredients the family supplies (CDF, simulator, variance function). |
| ~ | gaussian + residuals_osa; student + residuals_osa; lognormal + residuals_osa; shifted_lognormal + residuals_osa; skew_normal + residuals_osa; and 24 more | One-step-ahead residuals need the family to register its observation through OBS(). |
| ~ | gaussian + emmeans; student + emmeans; lognormal + emmeans; shifted_lognormal + emmeans; skew_normal + emmeans; and 28 more | Univariate fits only, and the mu predictor must be linear. |
| ~ | gaussian + frm_sample; student + frm_sample; lognormal + frm_sample; shifted_lognormal + frm_sample; skew_normal + frm_sample; and 28 more | Chains start jittered around the fitted mode. Use init_jitter to widen the spread, or init = “random” for a multimodal posterior. |
| x | tweedie + simulate; compois + simulate; multinomial + simulate; hurdle_poisson + simulate; cumulative + simulate; and 3 more | Refused: this family has no simulator yet. |

## Estimation modes

The modes constrain each other more than they constrain anything else.
`REML`, `quadrature`, and `profile` all move the fixed effects out of
the outer optimization problem, and no two of them can do that at once.
Bounds name outer parameters, so bounds on fixed effects are refused
under `REML` and under `profile`.

`autoscale`, `sparse_x`, and `verbose` are different in kind. They
change the arithmetic or the output, not the model.

|            | REML | quadrature | profile | autoscale | sparse_x | priors | bounds | verbose |
|:-----------|:----:|:----------:|:-------:|:---------:|:--------:|:------:|:------:|:-------:|
| REML       |      |     x      |    x    |    \+     |    \+    |   ~    |   x    |   \+    |
| quadrature |  x   |            |    x    |    \+     |    \+    |   \+   |   \+   |   \+    |
| profile    |  x   |     x      |         |    \+     |    \+    |   ~    |   x    |   \+    |
| autoscale  |  \+  |     \+     |   \+    |           |    \+    |   \+   |   \+   |   \+    |
| sparse_x   |  \+  |     \+     |   \+    |    \+     |          |   \+   |   \+   |   \+    |
| priors     |  ~   |     \+     |    ~    |    \+     |    \+    |        |   \+   |   \+    |
| bounds     |  x   |     \+     |    x    |    \+     |    \+    |   \+   |        |   \+    |
| verbose    |  \+  |     \+     |   \+    |    \+     |    \+    |   \+   |   \+   |         |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | REML + priors | Priors on fixed effects are accepted under REML while bounds on the same parameters are refused. The two surfaces are inconsistent. Priors on variance parameters behave as expected. |
| ~ | profile + priors | Priors on the profiled fixed effects are accepted, while bounds on the same parameters are refused. Treat priors under profile with care. |
| x | REML + quadrature | Refused: quadrature already marginalizes the random effects, so there is no inner problem left for REML to integrate. |
| x | REML + profile | Refused: profile = TRUE and REML both remove the fixed effects from the outer problem. |
| x | REML + bounds | Refused: under REML the fixed effects leave the outer parameter vector, so bounds naming them are rejected as unknown parameters. |
| x | quadrature + profile | Refused: profile = TRUE cannot be combined with quadrature = TRUE. |
| x | profile + bounds | Refused: with the fixed effects profiled out they leave the outer parameter vector, so bounds naming them are rejected. This pair was once accepted and the bounds were then applied to the wrong parameters. |

### Modes and addition terms

|           | REML | quadrature | profile | autoscale | sparse_x | priors | bounds | verbose |
|:----------|:----:|:----------:|:-------:|:---------:|:--------:|:------:|:------:|:-------:|
| weights() |  ?   |     \+     |   \+    |    \+     |    \+    |   ?    |   ?    |   \+    |
| trials()  |  ?   |     ?      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| cens()    |  \+  |     \+     |   \+    |    \+     |    \+    |   ?    |   ?    |   \+    |
| trunc()   |  \+  |     x      |   \+    |    \+     |    \+    |   ?    |   ?    |   \+    |
| se()      |  \+  |     \+     |   \+    |    \+     |    \+    |   ?    |   ?    |   \+    |
| mi()      |  \+  |     x      |   \+    |    \+     |    \+    |   ?    |   ?    |   \+    |
| vint()    |  ?   |     ?      |    ?    |    \+     |    \+    |   ?    |   ?    |    ?    |
| vreal()   |  ?   |     ?      |    ?    |    \+     |    \+    |   ?    |   ?    |    ?    |

| Status | Pairs | Note |
|:---|:---|:---|
| x | trunc() + quadrature | Refused: the truncation normalizer is log(F(ub) - F(lb)) over plain CDFs, and the Gauss-Kronrod nodes reach random-effect values where that difference underflows to exactly zero while the density is still representable. The integrand is then +Inf and the marginalized objective is -Inf, at the Laplace optimum as well as at the starting values, so the fit used to report logLik = +Inf as converged. Laplace stays near the mode and is unaffected: use quadrature = FALSE, REML, or profile for truncated responses. |
| x | mi() + quadrature | Refused: the imputed values are themselves random effects that quadrature cannot marginalize. |

### Modes and covariance structures

Only `quadrature` reads the covariance structures. It marginalizes one
scalar random intercept at a time, so every block in the model must be a
one-dimensional `us`, `diag`, or `homdiag` term.

|         | REML | quadrature | profile | autoscale | sparse_x | priors | bounds | verbose |
|:--------|:----:|:----------:|:-------:|:---------:|:--------:|:------:|:------:|:-------:|
| us      |  \+  |     ~      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| diag    |  \+  |     ~      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| homdiag |  \+  |     ~      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| cs      |  \+  |     x      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| ar1     |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |    ~    |
| hetar1  |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |    ~    |
| ou      |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |    ~    |
| toep    |  \+  |     x      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| homtoep |  \+  |     x      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| homcs   |  \+  |     x      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| exp     |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |   \+    |
| gau     |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |   \+    |
| mat     |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |   \+    |
| rr      |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |    ~    |
| equalto |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |    ~    |
| gr_cov  |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |    ~    |
| gr_prec |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |    ~    |
| smooth  |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |    ~    |
| gp      |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |   \+    |
| hsgp    |  \+  |     x      |    ?    |    \+     |    \+    |   ~    |   ~    |   \+    |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | us + quadrature; diag + quadrature; homdiag + quadrature | Allowed only when the block is one-dimensional, that is a scalar random intercept. Correlated slopes are refused. Several such blocks are fine, nested ones included: (1 \| ga/gb) becomes an iterated one-dimensional integral. |
| ~ | ar1 + priors; ar1 + bounds; ar1 + verbose | ar1() reads the level order. Over gapped integer levels it keeps positional semantics and warns. Use ou() over num_factor() for irregular spacing. |
| ~ | hetar1 + priors; hetar1 + bounds; hetar1 + verbose | hetar1() reads the level order. Over gapped integer levels it keeps positional semantics and warns. |
| ~ | ou + priors; ou + bounds; ou + verbose | ou() is the irregular-spacing structure. Build the levels with num_factor() so the metric distance is recoverable. |
| ~ | exp + priors; exp + bounds; gau + priors; gau + bounds; mat + priors; mat + bounds | Spatial structures need coordinates built with num_factor(x, y). |
| ~ | rr + priors; rr + bounds; rr + verbose | rr() gives a reduced-rank block; the rank d must not exceed the term dimension. |
| ~ | equalto + priors; equalto + bounds; equalto + verbose | equalto(x + 0 \| g, V) fixes the term covariance to V, which must be square and match the term dimension. |
| ~ | gr_cov + priors; gr_cov + bounds; gr_cov + verbose | gr(cov = A) accepts correlated slopes; the block covariance is the Kronecker product of A and the term covariance. A needs dimnames covering every grouping level. |
| ~ | gr_prec + priors; gr_prec + bounds; gr_prec + verbose | gr(prec = Q) supports intercept-only terms: (1 \| gr(g, prec = Q)). |
| ~ | smooth + priors; smooth + bounds; smooth + verbose | smooth is the internal structure behind s() and t2(); it is not written directly in a formula. |
| ~ | gp + priors; gp + bounds; hsgp + priors; hsgp + bounds | gp() and hsgp() are predictor specials, not bar terms. Write gp(x), not (gp(x) \| g). |
| x | cs + quadrature; ar1 + quadrature; hetar1 + quadrature; ou + quadrature; toep + quadrature; and 12 more | Refused: quadrature marginalizes one scalar random intercept at a time. Every block must be a dimension-1 us, diag, or homdiag term. |

## Model structures

|             | weights() | trials() | cens() | trunc() | se() | mi() | vint() | vreal() |
|:------------|:---------:|:--------:|:------:|:-------:|:----:|:----:|:------:|:-------:|
| mvbf        |    \+     |    \+    |   \+   |   \+    |  \+  |  \+  |   ?    |    ?    |
| rescor      |     x     |    ?     |   x    |    x    |  x   |  x   |   ?    |    ?    |
| \|ID\|      |     ?     |    ?     |   ?    |    ?    |  ?   |  ?   |   ?    |    ?    |
| nl          |     ?     |    ?     |   \+   |    ?    |  ?   |  ?   |   ?    |    ?    |
| mixture     |    \+     |    ?     |   x    |    x    |  x   |  x   |   ?    |    ?    |
| mixture_mvn |     ?     |    ?     |   x    |    x    |  ?   |  ?   |   ?    |    ?    |

| Status | Pairs | Note |
|:---|:---|:---|
| x | rescor + weights(); rescor + se() | Refused. |
| x | rescor + cens() | Refused. This pair was once accepted with the censoring silently dropped. |
| x | mixture + cens() | Refused: mixture() has no CDF, so the CDF guard rejects cens(). |
| x | mixture_mvn + cens(); mixture_mvn + trunc() | Refused: mixture_mvn() has no CDF. |
| x | rescor + trunc() | Refused. This pair was once accepted with the truncation silently dropped. |
| x | mixture + trunc() | Refused: mixture() has no CDF. |
| x | mixture + se() | Refused: se() is supported for gaussian and student families only. |
| x | rescor + mi() | Refused: mi() cannot be combined with rescor = TRUE. |
| x | mixture + mi() | Refused: mi() on the mixture response is not supported. |

### Structures and post-fit methods

Several post-fit methods are univariate-only. A multivariate fit
supports [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) and refuses the
rest.

|  | fitted | predict | simulate | residuals | residuals_osa | emmeans | frm_sample | confint_profile | hypothesis_profile |
|:---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| mvbf | \+ | \+ | x | x | ? | x | x | x | x |
| rescor | \+ | \+ | x | ? | x | x | ? | ? | ? |
| \|ID\| | ? | ? | ? | ? | ? | ? | ~ | ? | ? |
| nl | \+ | ~ | ? | ? | ? | x | ? | ? | ? |
| mixture | ? | ? | ~ | ? | ? | ? | ~ | ? | ? |
| mixture_mvn | ? | ? | x | ? | ? | ? | ~ | ? | ? |

| Status | Pairs | Note |
|:---|:---|:---|
| ~ | \|ID\| + frm_sample; mixture_mvn + frm_sample | Chains start jittered around the fitted mode. Use init_jitter to widen the spread, or init = “random” for a multimodal posterior. |
| ~ | nl + predict | Point predictions work. se.fit is not supported for the nonlinear predictor; request a nonlinear parameter with dpar instead. |
| ~ | mixture + simulate | Works only when every component family has a simulator. |
| ~ | mixture + frm_sample | Mixture posteriors are multimodal. Sample with init = “random” rather than the mode-anchored default. |
| x | mvbf + simulate; mvbf + residuals; mvbf + emmeans; mvbf + frm_sample; mvbf + confint_profile; mvbf + hypothesis_profile | Refused: the post-fit methods below are univariate-only for now. |
| x | rescor + simulate | Refused: simulate() is not supported for multivariate fits yet. |
| x | rescor + residuals_osa | Refused: residuals() is not supported for multivariate fits yet. |
| x | rescor + emmeans | Refused: emmeans support is univariate-only for now. |
| x | nl + emmeans | Refused: emmeans support needs a linear mu predictor. |
| x | mixture_mvn + simulate | Refused: mixture_mvn() has no simulator yet. |

## Predictor specials

|           | REML | quadrature | profile | autoscale | sparse_x | priors | bounds | verbose |
|:----------|:----:|:----------:|:-------:|:---------:|:--------:|:------:|:------:|:-------:|
| s()       |  ?   |     x      |    ?    |     ?     |    \+    |   ?    |   ?    |    ?    |
| t2()      |  ?   |     x      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| mo()      |  \+  |     ?      |   \+    |     ?     |    ?     |   ?    |   ?    |    ?    |
| mi_pred() |  ?   |     ?      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |
| gp_pred() |  \+  |     x      |    ?    |     ?     |    ?     |   ?    |   ?    |    ?    |
| cs_pred() |  ?   |     ?      |    ?    |    \+     |    \+    |   ?    |   ?    |   \+    |

| Status | Pairs | Note |
|:---|:---|:---|
| x | s() + quadrature | Refused in practice: a smooth is a wide random-effect block, so the scalar-intercept guard rejects it. |
| x | t2() + quadrature | Refused in practice: a smooth is a wide random-effect block. |
| x | gp_pred() + quadrature | Refused in practice: gp() builds a wide block, so the scalar-intercept guard rejects it. |

The specials `mi()`, `gp()`, and `cs()` share a name with a covariance
structure. The registry writes the predictor forms as `mi_pred()`,
`gp_pred()`, and `cs_pred()` to keep them apart.

## Formula grammar

Three spellings have restrictions of their own.

| Spelling | Status | Note |
|:---|:---|:---|
| x \* (1 \| g) | x | Refused: a bar term crossed with \* or : (as in x \* (1 \| g)) is not a random-effect specification (lme4#196). Write the crossing inside the bar: (x \| g). This spelling was once accepted with the crossing silently dropped. |
| (1 \| factor(x)) | \+ | Call-valued grouping factors are supported: (1 \| factor(x)) and (1 \| interaction(a, b)) both build the grouping factor from the model frame. |
| (x \|\| g) | \+ | (x \|\| g) gives uncorrelated terms. With a factor on the left, (f \|\| g) routes to diag, that is one independent effect per factor level. |

## Coverage

| Status      | Pairs | Share |
|:------------|------:|:------|
| works       |  1708 | 46%   |
| conditional |   891 | 24%   |
| refused     |   391 | 10%   |
| broken      |     0 | 0%    |
| untested    |   760 | 20%   |

The untested share is the honest measure of what this registry does not
yet know. It shrinks as pairs are tested, not as the code is trusted. To
close one, add a test and change the rule in `R/compat.R`.
