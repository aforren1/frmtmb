# Simulate responses from a formula and parameters

Builds the design from `formula` and `data` exactly as [`frm()`](frm.md)
would, sets the parameters from `newparams`, and simulates responses.
Random effects are redrawn from their covariance for every simulation
unless `newparams$b` supplies them.

## Usage

``` r
frm_simulate(formula, data, family = NULL, newparams, nsim = 1, seed = NULL)
```

## Arguments

- formula:

  A [`bf()`](bf.md) formula (with a family attached) or a plain formula
  plus `family`.

- data:

  Model data, including a dummy response column.

- family:

  Family, when `formula` does not carry one.

- newparams:

  Named list of parameter vectors matching the template: `beta` (and
  `betad`, `theta`, `thetar` as the model requires); optionally `b` to
  fix the random effects across simulations.

- nsim, seed:

  As in [`simulate()`](https://rdrr.io/r/stats/simulate.html).

## Value

A data frame with `nsim` columns of simulated responses.

## Details

`data` must contain a response column with values that are valid for the
family (any dummy values do; they only anchor the design). Inspect the
required parameter layout with
`frm(formula, data, dry_run = "frame")$par_template`.

## Examples

``` r
# power analysis: simulate from a design with chosen parameters
dd <- data.frame(x = rnorm(60), g = factor(rep(1:6, 10)), y = 0)
frm(bf(y ~ x + (1 | g)) + gaussian(), dd,
    dry_run = "frame")$par_template   # the required layout
#> $beta
#> (Intercept)           x 
#>           0           0 
#> 
#> $betad
#> sigma_(Intercept) 
#>                 0 
#> 
#> $b
#> [1] 0 0 0 0 0 0
#> 
#> $theta
#> [1] 0
#> 
sims <- frm_simulate(bf(y ~ x + (1 | g)) + gaussian(), dd,
                     newparams = list(beta = c(1, 0.5),
                                      betad = log(0.7),
                                      theta = log(0.5)),
                     nsim = 3, seed = 1)
head(sims)
#>      sim_1      sim_2      sim_3
#> 1 1.007314 0.50436944  0.7886221
#> 2 1.731551 1.20181120  0.3238773
#> 3 0.817993 0.03184368 -0.3540123
#> 4 1.327337 2.03378663  0.7988518
#> 5 1.936121 0.64057031 -0.6028489
#> 6 1.028016 0.81116086  0.5767491
```
