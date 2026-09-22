# Refit a model to a new response

Reuses the assembled design (no formula parsing, no frame assembly) and
warm-starts the optimizer at the previous estimates, so a refit costs
one re-tape plus the optimization. This is the engine for parametric
bootstrap: simulate responses with
[`simulate()`](https://rdrr.io/r/stats/simulate.html), refit to each.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
refit(object, newresp, start = NULL, ...)
```

## Arguments

- object:

  A `frmtmb_fit` for a univariate model.

- newresp:

  Replacement response: a vector of the original length, or a matrix of
  the original dimensions for matrix responses.

- start:

  Optional named start list (as in
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)); when
  given it replaces the warm start.

- ...:

  Refused: an argument the method does not have is an error naming it,
  rather than silently changing nothing.

## Value

A new `frmtmb_fit`.

## Examples

``` r
set.seed(2)
dd <- data.frame(x = rnorm(80), g = factor(rep(1:8, 10)))
dd$y <- rnorm(80, 1 + 0.5 * dd$x + rnorm(8, 0, 0.5)[dd$g], 1)
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = dd)
# refit to a simulated response (the parametric-bootstrap step)
ysim <- simulate(fit, nsim = 1, re_formula = NA)[[1]]
rf <- refit(fit, ysim)
fixef(rf)
#>            Estimate Est.Error      Q2.5     Q97.5
#> Intercept 0.9330103 0.2342002 0.4739864 1.3920342
#> x         0.5717605 0.1101704 0.3558305 0.7876906
```
