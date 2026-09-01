# Per-group coefficients (fixed effects plus conditional modes)

Follows the lme4/glmmTMB/brms convention: for each random-effect
grouping factor, the fixed effects of its linear predictor broadcast
over the group levels, with the conditional modes added to the matching
columns. Use
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) for
the fixed effects alone.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
coef(object, ...)
```

## Arguments

- object:

  A `frmtmb_fit`.

- ...:

  Unused.

## Details

The result is a list of data frames keyed by grouping factor. When
random effects appear in more than one dpar (or response), an outer
layer keyed like
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) is
added. Smooth terms are excluded. A fit without random effects returns
[`fixef()`](https://aforren1.github.io/frmtmb/reference/fixef.md) (the
single coefficient vector when there is one linear predictor).
