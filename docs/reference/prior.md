# Set up priors with brms's quoting spelling

`prior()` is
[`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)
with the distribution given UNQUOTED, as brms's `prior()` takes it:
`prior(normal(5000, 1000), nlpar = "ult")` is
`set_prior("normal(5000, 1000)", nlpar = "ult")`. Every argument is
deparsed rather than evaluated, so `class = b` and `class = "b"` mean
the same thing, and a variable holding a distribution is deparsed to its
NAME rather than its value: use `prior_string()` to build a prior from
strings computed at run time.

## Usage

``` r
prior(prior, ...)

prior_(prior, ...)

prior_string(prior, ...)
```

## Arguments

- prior:

  Distribution string, e.g. `"normal(0, 5)"`, or a
  [`prior_normal()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)/[`prior_t()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)/[`prior_lkj()`](https://aforren1.github.io/frmtmb/reference/frmtmb-priors.md)
  object, or `""` for bounds only.

- ...:

  Any of
  [`set_prior()`](https://aforren1.github.io/frmtmb/reference/set_prior.md)'s
  remaining arguments: `class`, `coef`, `group`, `resp`, `dpar`,
  `nlpar`, `lb`, `ub`.

## Value

A `frmtmb_priorlist`.

## Details

`prior_()` takes one-sided formulas, calls, names or constants
(`prior_(~normal(0, 10), class = ~b)`) and `prior_string()` takes plain
strings; both exist so that priors can be built programmatically, and
both are brms's.

A frmtmb prior and a brms prior are different objects, and with brms
attached after frmtmb its `prior()` masks this one. Nothing breaks:
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) and
[`frm_sample()`](https://aforren1.github.io/frmtmb/reference/frm_sample.md)
accept a `brmsprior` object and translate its rows, so
`c(prior(...), prior(...))` copied out of a brms script works whichever
`prior()` was in scope.

## Examples

``` r
# the brms nonlinear vignette's spelling
prior(normal(5000, 1000), nlpar = "ult")
#> normal(5000, 1000) class=b nlpar=ult

# combine with c() or `+`, as with set_prior()
c(prior(normal(1, 2), nlpar = "omega"),
  prior(normal(45, 10), nlpar = "theta"))
#> normal(1, 2) class=b nlpar=omega
#> normal(45, 10) class=b nlpar=theta

# the programmatic spellings
prior_(~normal(0, 10), class = ~b)
#> normal(0, 10) class=b
prior_string(paste0("normal(0, ", 2 * 5, ")"), class = "b")
#> normal(0, 10) class=b
```
