# Column access on a prior specification

A `frmtmb_priorlist` holds parsed specifications, and brms's prior
object is a data frame of strings. `$` reads the columns brms's object
has, from the table
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) builds,
so `set_prior("normal(0, 2)", class = c("b", "sd"))$class` is
`c("b", "sd")` in both packages. The name must match a column exactly;
anything else is `NULL`, as it is on a data frame. Assigning a column
rebuilds every specification from the edited table with
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md),
so an edit is checked as a new call would be.

## Usage

``` r
# S3 method for class 'frmtmb_priorlist'
x$name

# S3 method for class 'frmtmb_priorlist'
x$name <- value
```

## Arguments

- x:

  A `frmtmb_priorlist`.

- name:

  A column: `prior`, `class`, `coef`, `group`, `resp`, `dpar`, `nlpar`,
  `lb`, `ub` or `source`.

- value:

  The new column, recycled as a data frame column is.

## Value

For `$`, a character vector with one element per specification. For
`$<-`, the rebuilt `frmtmb_priorlist`.

## Examples

``` r
pr <- set_prior("normal(0, 2)", class = c("b", "sd"))
pr$class
#> [1] "b"  "sd"
pr$prior[2] <- "exponential(1)"
pr
#>           prior class coef group resp dpar nlpar   lb   ub source
#>    normal(0, 2)     b                            <NA> <NA>   user
#>  exponential(1)    sd                            <NA> <NA>   user
```
