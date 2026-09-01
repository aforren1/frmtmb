# DHARMa residual diagnostics

Builds a DHARMa object from the fit's family simulator: scaled quantile
residuals that are uniform under a correctly specified model. Plot the
result with [`plot()`](https://rdrr.io/r/graphics/plot.default.html), or
run
[`DHARMa::testUniformity()`](https://rdrr.io/pkg/DHARMa/man/testUniformity.html),
[`DHARMa::testDispersion()`](https://rdrr.io/pkg/DHARMa/man/testDispersion.html),
[`DHARMa::testZeroInflation()`](https://rdrr.io/pkg/DHARMa/man/testZeroInflation.html)
and friends on it.

## Usage

``` r
dharma_residuals(fit, nsim = 250, re.form = NULL, seed = NULL, ...)
```

## Arguments

- fit:

  A `frmtmb_fit` (univariate; the family needs a simulator).

- nsim:

  Number of simulated response vectors.

- re.form:

  Passed to
  [`simulate.frmtmb_fit()`](https://aforren1.github.io/frmtmb/reference/simulate.frmtmb_fit.md):
  `NULL` (default) conditions on the estimated random effects; `NA`
  redraws them.

- seed:

  Optional RNG seed for the simulations.

- ...:

  Passed to
  [`DHARMa::createDHARMa()`](https://rdrr.io/pkg/DHARMa/man/createDHARMa.html).

## Value

A `DHARMa` object.

## Examples

``` r
if (requireNamespace("DHARMa", quietly = TRUE)) {
  set.seed(1)
  dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
  dd$y <- rpois(100, exp(0.3 + 0.4 * dd$x + rnorm(10, 0, 0.6)[dd$g]))
  fit <- frm(bf(y ~ x + (1 | g)) + poisson(), data = dd)

  # scaled quantile residuals: uniform under a correct model, which
  # is what Pearson residuals cannot give for a discrete family
  res <- dharma_residuals(fit, nsim = 100, seed = 1)
  plot(res)
  DHARMa::testUniformity(res, plot = FALSE)
  DHARMa::testDispersion(res, plot = FALSE)
}

#> 
#>  DHARMa nonparametric dispersion test via sd of residuals fitted vs.
#>  simulated
#> 
#> data:  simulationOutput
#> dispersion = 0.99641, p-value = 0.88
#> alternative hypothesis: two.sided
#> 
```
