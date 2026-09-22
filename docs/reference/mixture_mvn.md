# Multivariate gaussian mixture family

`mixture_mvn(K, D)` does model-based clustering of an n x D matrix
response (mclust-style): K classes, each with its own D-dimensional mean
and a D x D covariance from mclust's model taxonomy. Every class mean is
a full linear predictor - the main model formula applies to all of
them - so cluster means may depend on covariates, which mclust cannot
do. The location dpars are named `mu<k>d<j>` (class k, response column
j) and are individually overridable, e.g. `bf(Y ~ x, mu2d1 ~ 1)` (all
except the first, `mu1d1`). Mixing weights are `theta1 ... theta{K-1}`,
multinomial logit against class K, each with its own linear predictor -
so gating on covariates works like
[`mixture()`](https://aforren1.github.io/frmtmb/reference/mixture.md).

## Usage

``` r
mixture_mvn(K, D, model = "VVV")
```

## Arguments

- K:

  Number of mixture classes (at least 2).

- D:

  Number of response columns (at least 2; for `D = 1` use
  `mixture(gaussian(), ...)`).

- model:

  Covariance model name from mclust's vocabulary: one of `"EII"`,
  `"VII"`, `"EEI"`, `"VEI"`, `"EVI"`, `"VVI"`, `"EEE"`, `"VVV"` (the
  default, a free covariance per class).

## Value

A `frmtmb_family`.

## Details

Class covariances are family-level extra parameters, covariate-free, and
their structure follows `model`, mclust's volume-shape-orientation
taxonomy for `Sigma_k = lambda_k * D_k * A_k * D_k'`:

|  |  |  |
|----|----|----|
| `EII` | spherical, equal volume | `sigmaraw`, one log-SD |
| `VII` | spherical, varying volume | `sigmaraw<k>`, one log-SD each |
| `EEI` | diagonal, equal volume and shape | `sigmaraw`, D log-SDs |
| `VEI` | diagonal, varying volume, equal shape | `sigmavol<k>` plus `sigmashape` (D - 1) |
| `EVI` | diagonal, equal volume, varying shape | `sigmavol` plus `sigmashape<k>` (D - 1 each) |
| `VVI` | diagonal, free | `sigmaraw<k>`, D log-SDs each |
| `EEE` | one shared full covariance | `sigmaraw`, a `us` block |
| `VVV` | free full covariance per class (default) | `sigmaraw<k>`, one `us` block each |

A `us` block is D log-SDs then the scaled-Cholesky correlation entries,
as in [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)'s
`us()` covariance structure. The `sigmashape` vectors hold the first D -
1 log-shape entries; the last is minus their sum, which is what fixes
`det(A_k) = 1`. Log-SDs start at the per-column response SDs and
correlations at zero; class means start on spread-out per-column
response quantiles to break the label symmetry. The usual finite-mixture
ML caveats apply: the likelihood is invariant to relabeling and can be
multimodal (compare starts via
[`frm_allfit()`](https://aforren1.github.io/frmtmb/reference/frm_allfit.md)).
[`mixture_probs()`](https://aforren1.github.io/frmtmb/reference/mixture_probs.md)
returns posterior class probabilities per row;
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the n x
D mixture-mean matrix. Covariances take no linear predictor (no
covariance regression), and the models with a class-varying eigenvector
basis (`EEV`, `VEV`, `EVE`, `VEE`, `VVE`, `EVV`) are not available.
`cens()`/[`trunc()`](https://rdrr.io/r/base/Round.html),
[`mvbf()`](https://aforren1.github.io/frmtmb/reference/mvbf.md), and
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
#>                   Estimate Est.Error       Q2.5     Q97.5
#> mu1d1_Intercept 0.08342536 0.1662316 -0.2423826 0.4092333
#> mu1d2_Intercept 0.13395425 0.1434356 -0.1471744 0.4150829
#> mu2d1_Intercept 4.11113311 0.1724837  3.7730712 4.4491950
#> mu2d2_Intercept 4.11395998 0.1584580  3.8033881 4.4245319
head(mixture_probs(fit))
#>         class1       class2
#> [1,] 0.9999999 9.151779e-08
#> [2,] 1.0000000 1.152650e-09
#> [3,] 1.0000000 1.619378e-10
#> [4,] 0.9999990 1.000145e-06
#> [5,] 1.0000000 3.772682e-12
#> [6,] 1.0000000 2.427070e-12
# a shared spherical covariance (mclust's EII, k-means-like)
frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2, model = "EII"), data = dd)
#>  Family: mixture_mvn(K = 2, D = 2, model = "EII") 
#>  Links: mu1d1 = identity; mu1d2 = identity; mu2d1 = identity; mu2d2 = identity;
#>         theta1 = identity
#> 
#> Formula: Y ~ 1 
#>    Data: dd (Number of observations: 60) 
#>  Method: ML   logLik: -196.172   AIC: 404.344   BIC: 416.91 
#> 
#> Regression Coefficients:
#>                 Estimate Est.Error l-95% CI u-95% CI z value Pr(>|z|)
#> mu1d1_Intercept     0.08      0.16    -0.23     0.40    0.52      0.6
#> mu1d2_Intercept     0.14      0.16    -0.18     0.45    0.84      0.4
#> mu2d1_Intercept     4.11      0.16     3.80     4.43   25.57   <2e-16
#> mu2d2_Intercept     4.11      0.16     3.80     4.43   25.62   <2e-16
#> 
#> Further Distributional Parameters:
#>        Estimate Est.Error l-95% CI u-95% CI
#> theta1      0.5        NA       NA       NA
#> theta2      0.5        NA       NA       NA
```
