# The non-decision time, in the response's own units

Without `ndt_group()` this is
`predict(fit, dpar = "ndt", type = "response")`, which already reports a
time. With `ndt_group()` the `ndt` link is a plain logit on a FRACTION
of the row's own bound, so
[`predict()`](https://rdrr.io/r/stats/predict.html) reports that
fraction and this multiplies it back out. Written so that a caller need
not know which of the two a fit is.

## Usage

``` r
ndt_time(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted
  [`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md),
  [`lba()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/lba.md),
  [`rdm()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/rdm.md)
  or
  [`wiener_gng()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener_gng.md)
  model.

- newdata:

  Optional data frame. It must carry the model's `ndt_group()` column
  when the model has one.

- ...:

  Passed to [`stats::predict()`](https://rdrr.io/r/stats/predict.html);
  `re.form` and `allow_new_levels` are the useful ones.

## Value

A numeric vector, one non-decision time per row.

## Details

The bound is the fastest response of the row's `ndt_group()`, or of the
whole data set when the model has no `ndt_group()`, or the `max_ndt` the
family was given. Whichever it is, it is a property of the data the
model was FITTED to, so a prediction on new rows is scaled by the same
bound the fit used.

## See also

[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.eam/reference/wiener.md)
for what the two parameterizations are and why there are two.

## Examples

``` r
set.seed(1)
d <- ddm_simulate(300, mu = 1.2, bs = 1.5, ndt = 0.25)
fit <- frm(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5),
           family = wiener(), data = d)
# with no ndt_group() the two agree: `ndt` is already a time
head(predict(fit, dpar = "ndt", type = "response"), 3)
#>         1         2         3 
#> 0.2569868 0.2569868 0.2569868 
head(ndt_time(fit), 3)
#> [1] 0.2569868 0.2569868 0.2569868
```
