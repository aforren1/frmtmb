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

`key` is the PARSER NAME: the identifier
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) itself
matches when it reads a formula, which is what lets the tests check the
registry against the real vocabulary. A caller deciding whether a
formula could be written with a feature should read `key`, not `name`.

`name` is the DISPLAY name, which is what the printed table shows and
which carries `"()"` for a callable feature. It is usually also what a
formula writes, and for three rows it is not: `mi()`, `gp()` and `cs()`
each name both a predictor special and a covariance structure, so the
specials carry the display names `mi_pred()`, `gp_pred()` and
`cs_pred()` while keeping the keys `mi`, `gp` and `cs` that a formula
actually writes. Those three are the only rows where the two columns
differ, and `frmtmb.sample`'s pre-flight refusal asserts exactly that
set, so a fourth would fail its tests rather than being read as a name
no formula writes.

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
#>         8         5        24        36         5        10         9         7 
#> structure 
#>         6 
```
