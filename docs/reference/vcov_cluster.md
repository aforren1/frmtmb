# Cluster-robust (sandwich) covariance

The misspecification-robust covariance of the estimates when
observations are correlated within clusters in a way the model does not
describe. With `B` the inverse observed information of the marginal
likelihood at the optimum (that is, `vcov(object, full = TRUE)`) and
`s_g` the score of cluster `g`
([`cluster_scores()`](https://aforren1.github.io/frmtmb/reference/cluster_scores.md)),

## Usage

``` r
vcov_cluster(
  object,
  cluster,
  type = c("CR0", "CR1", "CR1p", "CR1S"),
  full = FALSE
)
```

## Arguments

- object:

  A `frmtmb_fit`.

- cluster:

  The clustering factor, as in
  [`cluster_scores()`](https://aforren1.github.io/frmtmb/reference/cluster_scores.md).

- type:

  Small-sample correction: `"CR0"` (none), `"CR1"`, `"CR1p"` or
  `"CR1S"`, spelled and defined as in `clubSandwich`.

- full:

  If `TRUE`, return the whole outer parameter block (covariance
  parameters included), named as in `vcov(full = TRUE)`; otherwise the
  fixed-effect block, as
  [`vcov.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/vcov.frmtmb_fit.md)
  returns.

## Value

A covariance matrix, with attributes `type`, `nclusters` and `df` (the
`G - 1` reference degrees of freedom).

## Details

\$\$V = a\_{G} \\ B \left(\sum_g s_g s_g'\right) B\$\$

the Liang-Zeger / White cluster sandwich for an M-estimator. The
small-sample factor \\a_G\\ is the clubSandwich one, with `G` clusters,
`N` rows and `p` estimated parameters: `"CR0"` is 1, `"CR1"` is
`G / (G - 1)`, `"CR1p"` is `G / (G - p)`, and `"CR1S"` is
`G * (N - 1) / ((G - 1) * (N - p))`, the Stata `vce(cluster)` factor.

## What it is a covariance for

The sandwich is taken over the WHOLE outer parameter vector, so the
returned fixed-effect block already carries the cost of estimating the
covariance parameters.
[`clubSandwich::vcovCR()`](http://jepusto.github.io/clubSandwich/reference/vcovCR.md)
on an `lmerMod` instead conditions on the variance parameters - it
sandwiches only the mean-model bread - so the two agree exactly only
when the fixed-effect / covariance-parameter block of the observed
information vanishes. To reproduce that conditional form, build it from
[`cluster_scores()`](https://aforren1.github.io/frmtmb/reference/cluster_scores.md)
and the fixed-effect block of `solve(vcov( object, full = TRUE))`.

`"CR2"` (Bell-McCaffrey) and `"CR3"` are refused rather than
approximated. Both are defined through the hat matrix of a linear (or
GLS) model; a Laplace-marginal likelihood with a nonlinear link has no
such representation, and nothing in the literature defines the
adjustment for one. Use `clubSandwich` on a matched linear model if
`"CR2"` is what the analysis needs.

## When it is refused

The marginal likelihood factors over clusters only when every random
effect is nested in (or equal to) the clustering factor, and when every
likelihood term is a product over rows. Fits that fail either condition
are refused with the reason: a random effect whose level spans two
clusters (including crossed effects, `mm()` pooled levels, a global
smooth, `gp()`, `car()` and the SPDE), a group-level mixture whose
groups span clusters, an `autocor()` residual, a family whose structure
does not declare `cluster_robust`, `rescor = TRUE`, `mi()`/`me()`,
`REML = TRUE`, `profile = TRUE`, `quadrature = TRUE`, and any fit made
with priors.
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
is the fallback in every one of those cases.

Few clusters make the estimator badly biased whatever the correction;
the usual advice is at least 30-50, and `"CR1"` with a `t(G - 1)`
reference distribution below that.
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`hypothesis()`](https://aforren1.github.io/frmtmb/reference/hypothesis.md)
use the `t(G - 1)` reference automatically when they are handed a matrix
from this function.

## References

Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis using
generalized linear models. *Biometrika* 73, 13-22.

White, H. (1982). Maximum likelihood estimation of misspecified models.
*Econometrica* 50, 1-25.

Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
cluster-robust inference. *Journal of Human Resources* 50, 317-372.

Pustejovsky, J. E. and Tipton, E. (2018). Small-sample methods for
cluster-robust variance estimation and hypothesis testing in fixed
effects models. *Journal of Business & Economic Statistics* 36, 672-683.

## See also

[`cluster_scores()`](https://aforren1.github.io/frmtmb/reference/cluster_scores.md)
for the pieces,
[`vcov.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/vcov.frmtmb_fit.md)
which accepts `cluster` and forwards here, and
[`frm_bootstrap()`](https://aforren1.github.io/frmtmb/reference/frm_bootstrap.md)
for the cases this refuses.

## Examples

``` r
set.seed(1)
G <- 40
dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                 x = rnorm(G * 6))
# cluster-specific error scale the model does not describe
dd$y <- 1 + 0.5 * dd$x +
  rnorm(G * 6, 0, rep(runif(G, 0.3, 2.5), each = 6))
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd,
           REML = FALSE)

sqrt(diag(vcov(fit)))                          # model-based
#>       (Intercept)                 x sigma_(Intercept) 
#>        0.10787700        0.11225600        0.04564353 
sqrt(diag(vcov_cluster(fit, ~ g, "CR1")))      # cluster-robust
#>       (Intercept)                 x sigma_(Intercept) 
#>        0.09237769        0.12179568        0.08048348 

# feeds straight into the inference methods, which need the whole
# outer parameter vector and pick up the t(G - 1) reference from it
confint(fit, parm = "x",
        vcov = vcov_cluster(fit, ~ g, "CR1", full = TRUE))
#>         lwr       upr       est
#> x 0.4379456 0.9306556 0.6843006
```
