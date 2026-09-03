# Update and refit a model

Re-evaluates the stored
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) call with
the given arguments replaced. Any
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) argument
can be updated by name.

## Usage

``` r
# S3 method for class 'frmtmb_fit'
update(object, formula., ..., evaluate = TRUE)
```

## Arguments

- object:

  A `frmtmb_fit`.

- formula.:

  A complete formula or
  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md), or a
  delta such as `~ . + z` or `. ~ . + z`.

- ...:

  Arguments of
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) to
  replace, e.g. `data`, `family`, `REML`. `newdata` is accepted for
  `data`.

- evaluate:

  If `FALSE`, return the updated call instead of the refitted model.

## Value

A `frmtmb_fit`, or the updated call when `evaluate = FALSE`.

## Details

The formula argument is `formula.`, as in
[`stats::update()`](https://rdrr.io/r/stats/update.html) and in brms;
`formula = ` reaches it by partial matching. A formula carrying a `.` is
a delta applied to the stored `mu` formula with
[stats::update.formula](https://rdrr.io/r/stats/update.formula.html)
semantics - one-sided `~ . + z`, dotted `. ~ . + z`, or a changed
response `z ~ . + x` - and keeps the dpar formulas, the fixed dpar
values and the family. A formula with no `.` replaces the stored one.
brms's `newdata` is accepted as a synonym for `data`.

## Examples

``` r
set.seed(3)
dd <- data.frame(x = rnorm(60), z = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
fit <- frm(bf(y ~ x), family = gaussian(), data = dd)

# a delta on the stored formula, in either spelling
fit2 <- update(fit, ~ . + z)
formula(fit2)
#> y ~ x + z
#> <environment: 0x0000021e154607e8>
formula(update(fit, . ~ . + z))
#> y ~ x + z
#> <environment: 0x0000021e154607e8>
fit3 <- update(fit, formula. = ~ . - x, newdata = dd[1:40, ])
nobs(fit3)
#> [1] 40
```
