# Add an addition term from another package

Registers an addition term (a brms "aterm": the `trials(n)` in
`y | trials(n) ~ x`) that frmtmb will accept on the left-hand side of a
formula. A family that lives outside frmtmb uses this to give its
per-row data the spelling its literature uses, instead of asking users
for `vint()`. Register from the contributing package's `.onLoad()`.

## Usage

``` r
frmtmb_register_aterm(name, arity = 1L, coerce = as.numeric)
```

## Arguments

- name:

  The term's name, as it is written in a formula, without parentheses:
  `"dec"` accepts `y | dec(response) ~ x`.

- arity:

  How many arguments the term takes. With `arity = 1` the value arrives
  as `aterms[[name]]`; above one it arrives as
  `aterms[[paste0(name, i)]]` for each argument, which is the `vint()`
  convention.

- coerce:

  Function applied to each evaluated argument, once, at frame assembly.
  It must return a numeric vector, because the value is baked into the
  tape as data. The default is
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html); a term whose
  natural spelling is a factor gives a function that maps the factor to
  the numbers the density wants.

## Value

`NULL`, invisibly. Called for the registration.

## Details

The registered term carries DATA and nothing else. Its evaluated,
coerced value reaches the family's `lpdf`, `lcdf`, `post$mean_fn` and
`sim` in the `aterms` list under `name` (or under `name1`, `name2`, ...
when `arity` is above one, following `vint()`), it is required on
`newdata`, and it is what `frmtmb_family(required_aterms =)` names. It
cannot reshape the likelihood: censoring, truncation and case weights
are core terms because the core acts on them, and a contributed term is
not a route to that.

Registering a name twice replaces the earlier entry, so that reloading
the contributing package is not an error. The eight core terms cannot be
replaced.

## See also

[`frmtmb_family()`](https://aforren1.github.io/frmtmb/reference/frmtmb_family.md)
for `required_aterms`, which is how a family says the term is not
optional,
[`frmtmb_register_frame_check()`](https://aforren1.github.io/frmtmb/reference/frmtmb_register_frame_check.md)
for the seam that refuses a data problem, and
[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)
for the accessors an extension may use

## Examples

``` r
# a decision indicator, given as a factor and delivered as 0/1
if (FALSE) { # \dontrun{
frmtmb_register_aterm("dec", arity = 1, coerce = function(x) {
  as.integer(factor(x)) - 1L
})
} # }
```
