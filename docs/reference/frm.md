# Fit a model

Fits a model specified with [`bf()`](bf.md) by maximum likelihood, using
the Laplace approximation for random effects through RTMB.

## Usage

``` r
frm(
  formula,
  data,
  family = NULL,
  REML = FALSE,
  start = NULL,
  control = frmtmb_control(),
  se = FALSE,
  na.action = stats::na.omit,
  dry_run = NULL
)
```

## Arguments

- formula:

  A `frmtmb_formula` from [`bf()`](bf.md) (with a family attached via
  `+`), or a plain formula combined with the `family` argument.

- data:

  A data frame.

- family:

  A family: a `frmtmb_family`, a
  [stats::family](https://rdrr.io/r/stats/family.html) object or
  constructor (for example `gaussian`, `poisson`, `binomial`), or a
  family name as a string. Overrides a family already attached to
  `formula`.

- REML:

  If `TRUE`, integrate the `mu` fixed effects out of the likelihood
  along with the random effects (restricted maximum likelihood).

- start:

  Optional named list of starting values; components must match the
  parameter template (`beta`, `betad`, `theta`).

- control:

  A list from [`frmtmb_control()`](frmtmb_control.md).

- se:

  If `TRUE`, run `RTMB::sdreport()` at fit time. The default (`FALSE`)
  defers it until standard errors are first needed (`summary`, `vcov`,
  `confint`, `predict(se.fit = TRUE)`), which cuts roughly a quarter off
  fit time in fit-and-predict or bootstrap loops. The deferred report is
  cached, so nothing is computed twice.

- na.action:

  How to handle missing values, as in
  [`stats::lm()`](https://rdrr.io/r/stats/lm.html) (default
  [stats::na.omit](https://rdrr.io/r/stats/na.fail.html)).

- dry_run:

  `"spec"` returns the parsed intermediate representation without
  touching `data`; `"frame"` returns the assembled design matrices and
  parameter template without fitting.

## Value

An object of class `frmtmb_fit`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(sleepstudy, package = "lme4")
fit <- frm(bf(Reaction ~ Days + (Days | Subject)) + gaussian(),
              data = sleepstudy)
summary(fit)
} # }
```
