# Family object field access

A family object is a list, and `$` on a list PARTIAL-matches: before
frmtmb answered brms's field names, `student()$link` returned the
`links` list rather than the link name brms code expects.

## Usage

``` r
# S3 method for class 'frmtmb_family'
x$name
```

## Arguments

- x:

  A `frmtmb_family`.

- name:

  The field name.

## Value

The field.

## Details

`$` on a `frmtmb_family` answers brms's link fields, `link`, `linkfun`,
`linkinv` and `link_<dpar>` for every parameter but the mean, by reading
`links` at the moment of the call, so they always describe the links the
fit uses, however the object was edited and whenever it was saved. Every
other name matches exactly. A name that is only the prefix of one field
is an error naming that field, and a name that matches nothing is
`NULL`, as on any list.

The link fields are not elements of the list, so
[`names()`](https://rdrr.io/r/base/names.html), `[[` and
[`str()`](https://rdrr.io/r/utils/str.html) do not show them. An element
stored under one of those names is refused when it is read, because `$`
would otherwise have to choose between it and the links the fit actually
applies.

## Examples

``` r
fam <- beta_binomial()
fam$link
#> [1] "logit"
fam$link_phi
#> [1] "log"
try(fam$link_p)   # a prefix of link_phi, refused
#> Error : A family object has no field `link_p`. `$` would have partial-matched `link_phi`, and it does not do that on a family: write `$link_phi` if that is the field you mean
```
