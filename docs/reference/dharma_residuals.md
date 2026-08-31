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

  Passed to [`simulate.frmtmb_fit()`](simulate.frmtmb_fit.md): `NULL`
  (default) conditions on the estimated random effects; `NA` redraws
  them.

- seed:

  Optional RNG seed for the simulations.

- ...:

  Passed to
  [`DHARMa::createDHARMa()`](https://rdrr.io/pkg/DHARMa/man/createDHARMa.html).

## Value

A `DHARMa` object.
