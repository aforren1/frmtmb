# Combine formulas into a multivariate model

Each response keeps its own formula, family, dpar formulas, and addition
terms. Residual correlation between gaussian responses is requested with
`rescor = TRUE` or `set_rescor()`. Random-effect correlation across
responses uses the brms `|ID|` syntax, e.g. `(1 | p | g)` in several
formulas correlates their `g` effects.

## Usage

``` r
mvbf(..., rescor = FALSE)

set_rescor(rescor = arg_unset(), rescor_value = arg_unset())
```

## Arguments

- ...:

  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) formulas,
  each with a family attached (or supply one `family` to
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) for all
  of them).

- rescor:

  Model residual correlation between the responses (gaussian only).
  `mvbf()` defaults to `FALSE`, `set_rescor()` to `TRUE`. It is brms's
  spelling on both, and brms's default on `set_rescor()`; brms's
  `mvbf()` defaults to `NULL` and decides later, which frmtmb settles at
  formula-assembly time instead.

- rescor_value:

  The spelling `set_rescor()` shipped with, still accepted as an alias
  of `rescor`. It existed only because this Rd page documents two
  functions and could not carry two `rescor` entries; brms spells the
  argument `rescor`, and so does this function now. Give one spelling or
  the other, not both.

## Value

An object of class `frmtmb_mvformula`.

## Details

The linked terms merge into one covariance block, so they must all name
the same grouping specification. When they all write `gr(g, cov = A)`
(or all `gr(g, prec = Q)`) with the same matrix, the merged block keeps
it: its covariance is `A (x) Sigma`, with `Sigma` unstructured across
the merged coefficients. A two-trait animal model is therefore the same
fit whether written across two responses with
`(1 | q | gr(id, cov = A))` or in long format as a single
`(0 + trait | gr(id, cov = A))`. Mixing structures under one key - a
plain `g` in one formula and `gr(g, cov = A)` in another, or `cov =`
against `prec =` - is refused, because a merged block has room for one
structure.

## Examples

``` r
set.seed(2)
n <- 160
dd <- data.frame(x = rnorm(n), g = factor(rep(1:16, 10)))
u <- cbind(rnorm(16, 0, 0.8), rnorm(16, 0, 0.8))
e <- rnorm(n)                      # a disturbance both responses see
dd$y1 <- 1 + 0.5 * dd$x + u[dd$g, 1] + e + rnorm(n, 0, 0.5)
dd$y2 <- 2 - 0.3 * dd$x + u[dd$g, 2] + e + rnorm(n, 0, 0.5)

# each response keeps its own formula and family
fit <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x)) + gaussian(), data = dd)
fixef(fit)
#>                Estimate  Est.Error       Q2.5      Q97.5
#> y1_Intercept  1.0775542 0.10391095  0.8738925  1.2812159
#> y2_Intercept  2.0235294 0.09885933  1.8297687  2.2172901
#> y1_x          0.5517536 0.09447715  0.3665817  0.7369254
#> y2_x         -0.2793594 0.08988416 -0.4555292 -0.1031897

# rescor estimates the correlation of the residuals
fit_rc <- frm(mvbf(bf(y1 ~ x), bf(y2 ~ x), rescor = TRUE) + gaussian(),
              data = dd)
rescor_matrix(fit_rc)
#>           y1        y2
#> y1 1.0000000 0.6720766
#> y2 0.6720766 1.0000000

# set_rescor() turns it on after the fact, and `+` also combines bf()s
mvbf(bf(y1 ~ x), bf(y2 ~ x)) + set_rescor(TRUE)
#> y1 ~ x
#> y2 ~ x
#> rescor: TRUE 
bf(y1 ~ x) + bf(y2 ~ x)
#> y1 ~ x
#> y2 ~ x
#> rescor: FALSE 

# |ID| correlates the random effects of the two responses
fit_id <- frm(mvbf(bf(y1 ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
                gaussian(), data = dd)
VarCorr(fit_id)
#> $g
#> $g$sd
#>               Estimate Est.Error      Q2.5    Q97.5
#> y1_Intercept 0.8606044 0.1732261 0.5210875 1.200121
#> y2_Intercept 0.6111245 0.1438907 0.3291040 0.893145
#> 
#> $g$cor
#> , , y1_Intercept
#> 
#>               Estimate Est.Error      Q2.5    Q97.5
#> y1_Intercept 1.0000000 0.0000000 1.0000000 1.000000
#> y2_Intercept 0.5987236 0.2196369 0.1682431 1.029204
#> 
#> , , y2_Intercept
#> 
#>               Estimate Est.Error      Q2.5    Q97.5
#> y1_Intercept 0.5987236 0.2196369 0.1682431 1.029204
#> y2_Intercept 1.0000000 0.0000000 1.0000000 1.000000
#> 
#> 
#> $g$cov
#> , , y1_Intercept
#> 
#>               Estimate Est.Error        Q2.5     Q97.5
#> y1_Intercept 0.7406399 0.2981582  0.15626049 1.3250193
#> y2_Intercept 0.3148905 0.1797248 -0.03736366 0.6671447
#> 
#> , , y2_Intercept
#> 
#>               Estimate Est.Error        Q2.5     Q97.5
#> y1_Intercept 0.3148905 0.1797248 -0.03736366 0.6671447
#> y2_Intercept 0.3734732 0.1758702  0.02877383 0.7181725
#> 
#> 
#> 
#> $residual__
#> $residual__$sd
#>     Estimate  Est.Error      Q2.5    Q97.5
#> y1 0.9987735 0.05886056 0.8834089 1.114138
#> y2 1.0927242 0.06440943 0.9664841 1.218964
#> 
#> 
```
