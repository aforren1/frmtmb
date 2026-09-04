# Feature metadata for the compatibility registry

The vocabulary the compatibility registry talks about: every family,
addition term, covariance structure, predictor special, estimation mode,
model structure, and post-fit method that has a declared compatibility
status.

## Usage

``` r
frm_compat_features()
```

## Value

A data frame with columns `name`, `key`, and `kind`.

## Details

`name` is how the feature is written in a formula or a call. `key` is
the identifier the package uses internally, which is what lets the tests
check the registry against
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md)'s real
vocabulary. Three specials share a name with a covariance structure, so
they carry the display names `mi_pred()`, `gp_pred()`, and `cs_pred()`.

## See also

[`frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.md),
[`frm_compat_rules()`](https://aforren1.github.io/frmtmb/reference/frm_compat_rules.md)

## Examples

``` r
head(frm_compat_features())
#>                name               key   kind
#> 1          gaussian          gaussian family
#> 2           student           student family
#> 3         lognormal         lognormal family
#> 4 shifted_lognormal shifted_lognormal family
#> 5       skew_normal       skew_normal family
#> 6        exgaussian        exgaussian family
table(frm_compat_features()$kind)
#> 
#>     aterm   autocor covstruct    family   grammar    method      mode   special 
#>         8         5        24        36         5         9         8         6 
#> structure 
#>         6 
```
