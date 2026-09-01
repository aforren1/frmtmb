# Fit a model across multiply imputed datasets

The frequentist counterpart of brms's `brm_multiple()`: fits the model
on every imputed dataset and pools by Rubin's rules with Barnard-Rubin
adjusted degrees of freedom. Accepts a plain list of data frames or a
[`mice::mids`](https://amices.org/mice/reference/mids.html) object.

## Usage

``` r
frm_multiple(formula, data, level = 0.95, ...)
```

## Arguments

- formula, ...:

  As in [`frm()`](frm.md).

- data:

  A list of completed data frames, or a
  [`mice::mids`](https://amices.org/mice/reference/mids.html).

- level:

  Confidence level for the `$pooled_varcorr` interval.

## Value

A `frmtmb_multiple` object: `pooled` (the Rubin table for the fixed
effects), `pooled_varcorr` (grp/term/type/estimate/ lwr/upr/df/fmi for
the random-effect SDs and correlations; `NULL` without random effects),
and `fits` (the per-imputation fits).

## Details

Fixed effects are pooled on the link scale, so distributional
coefficients like `sigma` appear in `$pooled` on their link (log) scale,
exactly as in [`vcov()`](https://rdrr.io/r/stats/vcov.html).
Random-effect SDs and correlations are pooled on the scales where a Wald
argument is defensible (log for SDs and GP ranges, Fisher z for
correlations, the [`confint_varcorr()`](confint_varcorr.md) convention)
and back-transformed, giving `$pooled_varcorr`. For missing-data
mechanisms beyond imputation, see in-model `mi()`.

## Examples

``` r
set.seed(8)
n <- 80
x <- rnorm(n)
y <- rnorm(n, 1 + 0.5 * x, 1)
x[sample(n, 15)] <- NA
imps <- lapply(1:3, function(i) {
  xi <- x
  xi[is.na(xi)] <- sample(x[!is.na(x)], sum(is.na(xi)), TRUE)
  data.frame(y = y, x = xi)
})
frm_multiple(bf(y ~ x) + gaussian(), data = imps)
#> Pooled over 3 imputations (Rubin's rules):
#> 
#>                   estimate     se    df     t         p     fmi
#> (Intercept)         0.9240 0.1330 73.93 6.946 1.244e-09 0.01093
#> x                   0.3156 0.1496 21.49 2.110 4.674e-02 0.24080
#> sigma_(Intercept)   0.1603 0.0803 70.37 1.996 4.980e-02 0.03063
```
