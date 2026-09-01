# Multivariate gaussian mixture family

`mixture_mvn(K, D)` does model-based clustering of an n x D matrix
response (mclust-style): K classes, each with its own D-dimensional mean
and unstructured D x D covariance. Every class mean is a full linear
predictor - the main model formula applies to all of them - so cluster
means may depend on covariates, which mclust cannot do. The location
dpars are named `mu<k>d<j>` (class k, response column j) and are
individually overridable, e.g. `bf(Y ~ x, mu2d1 ~ 1)` (all except the
first, `mu1d1`). Mixing weights are `theta1 ... theta{K-1}`, multinomial
logit against class K, each with its own linear predictor - so gating on
covariates works like [`mixture()`](mixture.md).

## Usage

``` r
mixture_mvn(K, D)
```

## Arguments

- K:

  Number of mixture classes (at least 2).

- D:

  Number of response columns (at least 2; for `D = 1` use
  `mixture(gaussian(), ...)`).

## Value

A `frmtmb_family`.

## Details

Class covariances are family-level extra parameters in the `us`
parameterization (D log-SDs, then scaled-Cholesky correlation entries),
stored as `sigmaraw<k>` in `fit$estimates`; they are covariate-free.
Log-SDs start at the per-column response SDs and correlations at zero;
class means start on spread-out per-column response quantiles to break
the label symmetry. The usual finite-mixture ML caveats apply: the
likelihood is invariant to relabeling and can be multimodal (compare
starts via [`frm_allfit()`](frm_allfit.md)).
[`mixture_probs()`](mixture_probs.md) returns posterior class
probabilities per row;
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the n x
D mixture-mean matrix.
`cens()`/[`trunc()`](https://rdrr.io/r/base/Round.html),
[`mvbf()`](mvbf.md), and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) are not supported.

## Examples

``` r
set.seed(1)
Y <- rbind(matrix(rnorm(60, 0), ncol = 2),
           matrix(rnorm(60, 4), ncol = 2))
dd <- data.frame(row = seq_len(nrow(Y)))
dd$Y <- Y
fit <- frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd)
fixef(fit)
#> $mu1d1
#> (Intercept) 
#>  0.08342536 
#> 
#> $mu1d2
#> (Intercept) 
#>   0.1339542 
#> 
#> $mu2d1
#> (Intercept) 
#>    4.111133 
#> 
#> $mu2d2
#> (Intercept) 
#>     4.11396 
#> 
#> $theta1
#>  (Intercept) 
#> 0.0008996756 
#> 
head(mixture_probs(fit))
#>         class1       class2
#> [1,] 0.9999999 9.151779e-08
#> [2,] 1.0000000 1.152650e-09
#> [3,] 1.0000000 1.619378e-10
#> [4,] 0.9999990 1.000145e-06
#> [5,] 1.0000000 3.772682e-12
#> [6,] 1.0000000 2.427070e-12
```
