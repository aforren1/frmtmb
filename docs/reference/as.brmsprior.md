# Transform an object into a prior specification

brms's `as.brmsprior()`. A data frame (or anything
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) accepts)
with a `prior` column becomes the object
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
returns, one specification per row. Missing columns take
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)'s
defaults (`class = "b"`, empty `coef`, `group`, `resp`, `dpar` and
`nlpar`, no bounds) and columns
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
does not take are dropped.

## Usage

``` r
as.brmsprior(x)
```

## Arguments

- x:

  A data frame with a `prior` column, a `brmsprior`, a
  [`default_prior()`](https://aforren1.github.io/frmtmb/reference/default_prior.md)
  table or a `frmtmb_priorlist`.

## Value

A `frmtmb_priorlist`.

## Details

The name is brms's, and so is the job: turn a table into the prior
object this package fits with. That object is a `frmtmb_priorlist`
rather than a `brmsprior`, because this package parses a density when it
is written rather than when the model is compiled. A `brmsprior` built
by brms is translated row by row, as
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) translates
one, and a table from
[`default_prior()`](https://aforren1.github.io/frmtmb/reference/default_prior.md)
or
[`validate_prior()`](https://aforren1.github.io/frmtmb/reference/validate_prior.md)
is read as written.

A row whose `prior` is empty or `"(flat)"` and that carries no bound
applies nothing, and is left out: it is the flat default a
[`default_prior()`](https://aforren1.github.io/frmtmb/reference/default_prior.md)
table lists for a slot nobody has set. So is a row whose `source` is
`"(vectorized)"`, which repeats the density of the class row above it.

## Examples

``` r
as.brmsprior(data.frame(prior = "normal(0,1)", coef = c("a", "b")))
#>        prior class coef group resp dpar nlpar   lb   ub source
#>  normal(0,1)     b    a                       <NA> <NA>   user
#>  normal(0,1)     b    b                       <NA> <NA>   user
```
