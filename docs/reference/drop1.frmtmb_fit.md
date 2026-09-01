# Single-term deletions

Drops each fixed-effect term of the primary (`mu`) formula in turn,
refits, and tabulates AIC (and likelihood-ratio tests with
`test = "Chisq"`), following
[`stats::drop1()`](https://rdrr.io/r/stats/add1.html) and lme4's
`drop1.merMod`. Random-effect, smooth, and `mo()`/`mi()` terms are not
part of the deletion scope.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
drop1(object, scope, test = c("none", "Chisq"), k = 2, ...)
```

## Arguments

- object:

  A `frmtmb_fit` from an ML fit (`REML = FALSE`) of a univariate model.

- scope:

  Terms to drop: a character vector or a right-hand-side formula.
  Defaults to all fixed-effect terms that marginality allows
  ([`stats::drop.scope()`](https://rdrr.io/r/stats/factor.scope.html)).

- test:

  `"Chisq"` adds likelihood-ratio tests.

- k:

  AIC penalty per parameter.

- ...:

  Unused.

## Value

An `anova` table with one row per dropped term.
