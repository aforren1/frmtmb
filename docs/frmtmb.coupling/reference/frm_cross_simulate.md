# Draw whole cross-spectral matrices from a fitted coupling model

[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) refuses on
a
[`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
fit, because a draw is a Hermitian matrix and the response slot carries
one of its four numbers. This returns all four, so a simulated data
frame can be refitted or compared against the observed one.

## Usage

``` r
frm_cross_simulate(fit, nsim = 1L, seed = NULL, newdata = NULL, re.form = NULL)
```

## Arguments

- fit:

  A
  [`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
  fit.

- nsim:

  Number of replicate data frames.

- seed:

  Optional seed, set with
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before drawing.

- newdata:

  Optional data frame of predictor values.

- re.form:

  Passed to [`stats::predict()`](https://rdrr.io/r/stats/predict.html);
  `NULL` keeps the random effects.

## Value

A list of `nsim` data frames, each with the columns
[`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md)
produces: `w11`, `w22`, `w12r`, `w12i`, `n`.

## How a draw is made

Each row's fitted spectral matrix is factored as `L L^H` and `n`
independent complex normal vectors are drawn through `L`; the returned
matrix is their outer-product sum. That is the definition of the complex
Wishart rather than a transformation of it, so a simulated row and the
density are one statement of the model. Checked against the density's
own first moment: over 20000 draws the mean of the simulated matrix
divided by `n` matches the fitted matrix to 0.0011 at `n = 16`.

## Examples

``` r
set.seed(6)
src <- rnorm(2048)
xs <- frm_cross_spectrum(src + rnorm(2048), 0.9 * src + rnorm(2048),
                         segments = 16)
xs <- xs[xs$freq < 0.1, ]
fit <- frmtmb::frm(
  frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
             pow2 ~ 1, coh ~ 1, phase ~ 1),
  family = cross_wishart(), data = xs)
str(frm_cross_simulate(fit, nsim = 2, seed = 1)[[1]])
#> 'data.frame':    12 obs. of  5 variables:
#>  $ w11 : num  27.7 28.1 37.8 37.6 38 ...
#>  $ w22 : num  24.3 19.9 28.2 33.5 33.2 ...
#>  $ w12r: num  13.3 12.1 12.9 18.5 23.5 ...
#>  $ w12i: num  -5.8 2.4 -1.17 -13.35 3.03 ...
#>  $ n   : num  16 16 16 16 16 16 16 16 16 16 ...
```
