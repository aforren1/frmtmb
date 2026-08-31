# Refit a model to a new response

Reuses the assembled design (no formula parsing, no frame assembly) and
warm-starts the optimizer at the previous estimates, so a refit costs
one re-tape plus the optimization. This is the engine for parametric
bootstrap: simulate responses with
[`simulate()`](https://rdrr.io/r/stats/simulate.html), refit to each.

## Usage

``` r
refit(object, newresp, ...)

# S3 method for class 'frmtmb_fit'
refit(object, newresp, start = NULL, ...)
```

## Arguments

- object:

  A `frmtmb_fit` for a univariate model.

- newresp:

  Replacement response: a vector of the original length, or a matrix of
  the original dimensions for matrix responses.

- ...:

  Unused.

- start:

  Optional named start list (as in [`frm()`](frm.md)); when given it
  replaces the warm start.

## Value

A new `frmtmb_fit`.
