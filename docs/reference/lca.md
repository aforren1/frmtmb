# Latent class analysis

`lca(K)` fits the classic latent class measurement model, the one
'poLCA' fits (Linzer and Lewis 2011): `J` polytomous indicator items per
subject, one latent class `c` in `1..K` per subject, and the items
conditionally independent given the class,

## Usage

``` r
lca(K, ncat = NULL, na.rm = TRUE)
```

## Arguments

- K:

  Number of latent classes (at least 2).

- ncat:

  Number of categories per item: a vector of length `J`, a single value
  for equally-sized items, or `NULL` (the default) to infer each item's
  count as its largest observed code, which is what poLCA does.

- na.rm:

  If `TRUE` (default) subjects with any missing item are dropped by
  `na.action`; if `FALSE` they are kept and each missing item's term is
  masked out of that subject's likelihood.

## Value

A `frmtmb_family`.

## Details

\$\$P(y_i) = \sum_c w_c \prod_j \pi\_{j,c,y\_{ij}}.\$\$

The response is a MATRIX, one row per subject and one column per item,
holding integer category codes `1..C_j`. Write it as
`cbind(item1, item2, ...)` on the left of the formula, or attach a
matrix column to the data (`dd$Y <- data.matrix(dd[items])`) and name
it. Items may have different numbers of categories.

The class-membership weights are the `theta1 ... theta{K-1}` dpars,
multinomial logit against class `K`, and the main model formula applies
to every one of them. So `bf(cbind(a, b, c) ~ 1)` is the plain
measurement model and `bf(cbind(a, b, c) ~ age + educ)` is poLCA's
latent class regression: covariates on class membership come from the
ordinary linear-predictor machinery, with the usual
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
on top. Individual gating predictors are overridable as
`bf(Y ~ x, theta2 ~ 1)` (all but `theta1`, which the main formula owns).

The item profiles are family extra parameters, held as
reference-category logits, one vector `pi<j>` per item (item `j`'s
`K * (C_j - 1)` free logits, class-major), and reported on the
probability scale by
[`lca_profiles()`](https://aforren1.github.io/frmtmb/reference/lca_profiles.md).
They take no linear predictor - a covariate acts on class membership,
never on an item's conditional response probability, which is what makes
the classes interpretable as measurement.

## Labeling and starting values

A latent class likelihood is invariant to relabeling the classes and is
genuinely multimodal, so the answer depends on where the optimizer
starts. poLCA uses random restarts (`nrep`). The starting values here
are deterministic instead: subjects are scored by the mean of their item
codes rescaled to `[0, 1]`, cut into `K` equal-count slices by that
score, and each slice's smoothed empirical category proportions become
one class's starting profile. Class 1 is therefore the low-score end and
class `K` the high-score end, and re-running the same data gives the
same labeling.

That fixes reproducibility, not multimodality. Do what
`poLCA(..., nrep = 10)` does and compare starts before reading a
solution: perturb the item parameters and refit,

    p0 <- fit$frame$par_template[fit$frame$extra_names]
    refits <- replicate(10, simplify = FALSE,
      frm(bf(Y ~ 1), family = lca(K = 3), data = dd,
          start = lapply(p0, function(v) v + rnorm(length(v)))))
    sapply(refits, logLik)

and keep the best.
[`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md)
is a different check: it re-runs the four optimizers from the SAME
start, so it tests the optimizer, not the surface. Ordering the gating
intercepts through `lower` and `upper` is the way to pin the labeling
itself.

Asking for more classes than the data hold drives item probabilities to
0 and 1. The optimizer then reports singular convergence and the
standard errors come back `NaN`, which is the boundary showing through
rather than a fault; poLCA lands on the same solutions. Compare
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) across `K` and read
[`lca_probs()`](https://aforren1.github.io/frmtmb/reference/lca_probs.md)'s
entropy before committing to a `K`.

## Missing item responses

With `na.rm = TRUE` (the default, poLCA's default) a subject with any
missing item is dropped by the usual `na.action`, and the message naming
the dropped rows is the ordinary one. With `na.rm = FALSE` the subject
is kept and only the missing item's factor leaves that subject's
likelihood, which is poLCA's `na.rm = FALSE` behavior; the mask is data,
so the tape never sees a branch. A subject missing EVERY item
contributes a constant and is better dropped.

## What is and is not supported

Post-fit,
[`lca_probs()`](https://aforren1.github.io/frmtmb/reference/lca_probs.md)
gives posterior class-membership probabilities per subject (with the
relative-entropy classification diagnostic attached),
[`lca_profiles()`](https://aforren1.github.io/frmtmb/reference/lca_profiles.md)
gives the item profile table, and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws a class per
subject and then its items.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
`predict(type = "response")` and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) are refused: the
response is a matrix of nominal codes, so an "expected item code" would
be an average of arbitrary labels. Read
[`lca_probs()`](https://aforren1.github.io/frmtmb/reference/lca_probs.md)
and
[`lca_profiles()`](https://aforren1.github.io/frmtmb/reference/lca_profiles.md)
instead. [`predict()`](https://rdrr.io/r/stats/predict.html) itself
returns the gating linear predictor (`theta1` by default, any `theta`
with `dpar =`), on `newdata` as well.

The gating coefficients are ordinary fixed effects, so
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md),
[`confint()`](https://rdrr.io/r/stats/confint.html) (Wald, profile and
uniroot),
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md),
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md),
`lower`/`upper` bounds and
[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.html)
all work on them; [`anova()`](https://rdrr.io/r/stats/anova.html)
compares nested gating formulas at one `K`.

Refused in this version: random effects and smooths anywhere in the
model (latent classes plus continuous random effects is the
growth-mixture shape - use `mixture(..., groups = ~g)` for that), `REML`
and `frmtmb_control(profile = TRUE)` (the classes are exchangeable, so
there is no single inner mode), `quadrature`, every addition term
([`weights()`](https://rdrr.io/r/stats/weights.html), `cens()`,
[`trunc()`](https://rdrr.io/r/base/Round.html), `se()`, `mi()`,
`trials()`),
[`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md), and
`residuals(type = "osa")`.

## References

Linzer, D. A. and Lewis, J. B. (2011). poLCA: An R Package for
Polytomous Variable Latent Class Analysis. *Journal of Statistical
Software*, 42(10), 1-29.

## See also

[`lca_probs()`](https://aforren1.github.io/frmtmb/reference/lca_probs.md),
[`lca_profiles()`](https://aforren1.github.io/frmtmb/reference/lca_profiles.md),
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md)
for a mixture of continuous responses,
[`mixture_mvn()`](https://aforren1.github.io/frmtmb/reference/mixture_mvn.md)
for model-based clustering of a numeric matrix.

## Examples

``` r
set.seed(1)
# four binary items measuring two well-separated classes
n <- 300
cl <- rbinom(n, 1, 0.4) + 1
pr <- rbind(c(0.85, 0.80, 0.75, 0.90), c(0.15, 0.20, 0.25, 0.10))
Y <- matrix(0L, n, 4)
for (j in 1:4) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
dd <- data.frame(x = rnorm(n))
dd$Y <- Y

fit <- frm(bf(Y ~ 1), family = lca(K = 2), data = dd)
lca_profiles(fit)
#> <lca profiles> 2 classes, 4 items
#> 
#> Estimated class sizes (mean prior probability):
#> class1 class2 
#> 0.2863 0.7137 
#> 
#> item1:
#>          cat1   cat2
#> class1 0.9661 0.0339
#> class2 0.1869 0.8131
#> 
#> item2:
#>          cat1   cat2
#> class1 0.8423 0.1577
#> class2 0.2413 0.7587
#> 
#> item3:
#>          cat1   cat2
#> class1 0.7982 0.2018
#> class2 0.2777 0.7223
#> 
#> item4:
#>          cat1   cat2
#> class1 0.9520 0.0480
#> class2 0.1459 0.8541
head(lca_probs(fit))
#>           class1      class2
#> [1,] 0.009336046 0.990663954
#> [2,] 0.009336046 0.990663954
#> [3,] 0.006291359 0.993708641
#> [4,] 0.992685029 0.007314971
#> [5,] 0.000915022 0.999084978
#> [6,] 0.061162256 0.938837744

# latent class regression: covariates gate class membership
frm(bf(Y ~ x), family = lca(K = 2), data = dd)
#> frmtmb fit: Y ~ x 
#> Family: lca(K = 2)   Method: ML 
#> logLik: -697.815  AIC: 1415.63  nobs: 300 
#> 
#> Fixed effects:
#>  theta1:
#> (Intercept)           x 
#>    -0.91700     0.04477 
```
