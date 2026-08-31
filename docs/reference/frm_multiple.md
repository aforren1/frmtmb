# Fit a model across multiply imputed datasets

The frequentist counterpart of brms's `brm_multiple()`: fits the model
on every imputed dataset and pools the fixed-effect estimates by Rubin's
rules with Barnard-Rubin adjusted degrees of freedom. Accepts a plain
list of data frames or a `mice::mids` object.

## Usage

``` r
frm_multiple(formula, data, ...)
```

## Arguments

- formula, ...:

  As in [`frm()`](frm.md).

- data:

  A list of completed data frames, or a `mice::mids`.

## Value

A `frmtmb_multiple` object: `pooled` (the Rubin table) and `fits` (the
per-imputation fits).

## Details

Only the fixed effects are pooled (covariance parameters are reported
per fit through `$fits`). For missing-data mechanisms beyond imputation,
see the roadmap note on latent-variable `mi()` in dev/feature-gaps.md.

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
