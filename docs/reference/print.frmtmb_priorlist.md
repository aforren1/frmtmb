# Print a prior specification

brms's layout: one specification prints as the parameter it addresses
followed by its density, `b_x ~ normal(0, 1)`, with any bounds in front
as `<lower=0>`; several print as a table with one row each, the table
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns,
where a row with no density of its own shows its class row's density as
`"(vectorized)"`. An empty prior prints nothing, as in brms.

## Usage

``` r
# S3 method for class 'frmtmb_priorlist'
print(x, show_df = NULL, ...)
```

## Arguments

- x:

  A `frmtmb_priorlist`.

- show_df:

  Print as a table (`TRUE`) or as one line per specification (`FALSE`).
  `NULL`, the default, prints a table when there is more than one
  specification, as brms does.

- ...:

  Must be empty.

## Value

`x`, invisibly.

## Examples

``` r
set_prior("normal(0,1)", coef = "x")
#> b_x ~ normal(0,1)
set_prior("cauchy(0,1)", class = "sd", group = "g")
#> sd_g ~ cauchy(0,1)
set_prior("normal(0, 2)", class = c("b", "sd"))
#>         prior class coef group resp dpar nlpar   lb   ub source
#>  normal(0, 2)     b                            <NA> <NA>   user
#>  normal(0, 2)    sd                            <NA> <NA>   user
```
