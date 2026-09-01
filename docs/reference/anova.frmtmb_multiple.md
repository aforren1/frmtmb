# Pooled model comparison across imputations (D1, D2, D3)

The multiply-imputed counterpart of
[`anova.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/anova.frmtmb_fit.md).
A likelihood-ratio test per imputation is not a test: the m statistics
have to be combined, and the three standard combining rules
([`mice::D1()`](https://amices.org/mice/reference/D1.html),
[`mice::D2()`](https://amices.org/mice/reference/D2.html),
[`mice::D3()`](https://amices.org/mice/reference/D3.html)) each
reference an F distribution whose denominator degrees of freedom absorb
the extra between-imputation uncertainty.

## Usage

``` r
# S3 method for class 'frmtmb_multiple'
anova(
  object,
  ...,
  method = c("D3", "D1", "D2"),
  use = c("likelihood", "wald"),
  constraint = NULL,
  dfcom = NULL
)
```

## Arguments

- object:

  A `frmtmb_multiple` from
  [`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md).

- ...:

  A second `frmtmb_multiple`, nested with `object` and fit to the same
  imputed datasets.

- method:

  Combining rule: `"D3"`, `"D1"` or `"D2"`.

- use:

  For `method = "D2"`, whether the per-imputation statistics are
  likelihood-ratio (`"likelihood"`) or Wald (`"wald"`) chi-squares.

- constraint:

  Instead of a second model, the names of the coefficients to test
  jointly against zero (`"D1"` and `method = "D2", use = "wald"` only).

- dfcom:

  Complete-data residual degrees of freedom for the `"D1"` reference
  distribution. Defaults to
  [`df.residual()`](https://rdrr.io/r/stats/df.residual.html) of the
  larger model's first fit; `Inf` selects the unadjusted form.

## Value

A one-row `frmtmb_pooled_anova` table: `statistic`, `df1`, `df2`, `p`,
`riv` (the relative increase in variance).

## Details

- `"D1"`:

  Multivariate Wald (Li, Raghunathan and Rubin 1991). Let `Qbar` be the
  pooled estimate of the `k` coefficients that the larger model adds,
  `Ubar` their average within-imputation covariance and `B` the
  between-imputation covariance. With the average relative increase in
  variance r, which is `(1 + 1/m) tr(B Ubar^-1) / k`, the statistic is
  `Qbar' ((1 + r) Ubar)^-1 Qbar / k` on `k` and `v` degrees of freedom.
  `v` uses the Reiter (2007) small-sample form when `dfcom` is finite
  and the Li et al form otherwise.

- `"D2"`:

  Pooled test statistics (Li, Meng, Raghunathan and Rubin 1991). From
  the `m` per-imputation statistics `d` - either likelihood-ratio
  chi-squares (`use = "likelihood"`) or Wald chi-squares
  (`use = "wald"`) - and r set to `(1 + 1/m) var(sqrt(d))`, the
  statistic is `(mean(d)/k - (m + 1)/(m - 1) r) / (1 + r)` on `k` and
  `k^(-3/m) (m - 1) (1 + 1/r)^2` degrees of freedom. D2 needs no
  covariance matrix, so it is the cheapest rule and the least efficient.

- `"D3"`:

  Pooled likelihood ratio (Meng and Rubin 1992). The average deviance
  difference is computed twice: at the imputation-specific estimates
  (`dbar`) and at the pooled estimates (`dtilde`), the second by
  evaluating each imputation's own likelihood at one common parameter
  vector. With r set to `(m + 1) / (k (m - 1)) (dbar - dtilde)`, the
  statistic is `dtilde / (k (1 + r))`. D3 is the default: it needs no
  covariance matrix either, and unlike D2 it stays valid when the
  coefficients are far from normal.

Pooled parameter values are the plain across-imputation means of the
optimizer's own parameter vector, so variance components are pooled on
the internal (log / Cholesky) scale, matching `$pooled`. No model is
refit: each imputation's objective is re-evaluated at the pooled vector,
which costs one Laplace solve.

[`mice::D3()`](https://amices.org/mice/reference/D3.html) instead fixes
the pooled coefficients as an offset and re-estimates every remaining
parameter, so it can differ by a few percent;
`mitml::testModels(method = "D3")` uses the plug-in form implemented
here. Note also that
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) counts the
dispersion parameter of a gaussian fit, which
[`lm()`](https://rdrr.io/r/stats/lm.html) does not, so `dfcom` is one
smaller than mice's default; pass `dfcom` explicitly to reproduce mice
exactly.

## References

Li, K. H., Raghunathan, T. E. and Rubin, D. B. (1991). Large-sample
significance levels from multiply imputed data using moment-based
statistics and an F reference distribution. *JASA* 86, 1065-1073.

Li, K. H., Meng, X.-L., Raghunathan, T. E. and Rubin, D. B. (1991).
Significance levels from repeated p-values with multiply-imputed data.
*Statistica Sinica* 1, 65-92.

Meng, X.-L. and Rubin, D. B. (1992). Performing likelihood ratio tests
with multiply-imputed data sets. *Biometrika* 79, 103-111.

## See also

[`frm_multiple()`](https://aforren1.github.io/frmtmb/reference/frm_multiple.md),
[`anova.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/anova.frmtmb_fit.md).

## Examples

``` r
set.seed(4)
n <- 60
imps <- lapply(1:4, function(i) {
  x <- rnorm(n)
  data.frame(y = rnorm(n, 1 + 0.5 * x), x = x, z = rnorm(n))
})
m1 <- frm_multiple(bf(y ~ x + z) + gaussian(), data = imps)
m0 <- frm_multiple(bf(y ~ x) + gaussian(), data = imps)
anova(m1, m0)
#> Pooled model comparison over 4 imputations (D3)
#>   Model 1: y ~ x + z
#>   Model 2: y ~ x
#> 
#>  statistic df1   df2     p   riv
#>     0.2167   1 7.978 0.654 1.585
anova(m1, m0, method = "D1")
#> Pooled model comparison over 4 imputations (D1)
#>   Model 1: y ~ x + z
#>   Model 2: y ~ x
#> 
#>  statistic df1   df2      p  riv
#>     0.3354   1 10.16 0.5751 1.19
```
