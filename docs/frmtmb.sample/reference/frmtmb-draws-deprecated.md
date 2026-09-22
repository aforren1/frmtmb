# Deprecated brms draws accessors

`posterior_samples()`, `nsamples()` and `parnames()` are brms spellings
that brms itself deprecated and still answers. Each warns and then does
what brms's does, so a ported script runs.

- `posterior_samples(x)` is `as.data.frame(x)`; `pars` is a regular
  expression unless `fixed = TRUE`, as in brms.

- `nsamples(x)` is
  [`ndraws()`](https://mc-stan.org/posterior/reference/draws-index.html).

- `parnames(x)` is
  [`variables()`](https://mc-stan.org/posterior/reference/variables.html).

Use `as_draws(x)` for a posterior draws object, `as.data.frame(x)`,
[`ndraws()`](https://mc-stan.org/posterior/reference/draws-index.html)
and
[`variables()`](https://mc-stan.org/posterior/reference/variables.html)
in new code.

## Usage

``` r
posterior_samples(x, pars = NA, ...)

# S3 method for class 'frmtmb_draws'
posterior_samples(
  x,
  pars = NA,
  fixed = FALSE,
  add_chain = FALSE,
  subset = NULL,
  as.matrix = FALSE,
  as.array = FALSE,
  ...
)

nsamples(object, ...)

# S3 method for class 'frmtmb_draws'
nsamples(object, subset = NULL, incl_warmup = FALSE, ...)

parnames(x, ...)

# S3 method for class 'frmtmb_draws'
parnames(x, ...)
```

## Arguments

- x, object:

  A `frmtmb_draws`.

- pars:

  Variable names. A regular expression unless `fixed = TRUE`; `NA` (the
  default) takes all of them.

- ...:

  Refused: an argument the method does not have is an error naming it.

- fixed:

  If `TRUE`, `pars` names variables exactly.

- add_chain:

  If `TRUE`, add the `chain` and `iter` columns brms adds.

- subset:

  Draw indices to keep.

- as.matrix, as.array:

  Return a matrix or a draws-by-chains-by-variables array instead of a
  data frame.

- incl_warmup:

  Refused:
  [`frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.md)
  discards the warmup, so there is no warmup draw to count.

## Value

A data frame (or matrix, or array) of draws for `posterior_samples()`,
one integer for `nsamples()`, and a character vector for `parnames()`.
