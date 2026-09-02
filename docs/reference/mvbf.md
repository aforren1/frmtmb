# Combine formulas into a multivariate model

Each response keeps its own formula, family, dpar formulas, and addition
terms. Residual correlation between gaussian responses is requested with
`rescor = TRUE` or `set_rescor()`. Random-effect correlation across
responses uses the brms `|ID|` syntax, e.g. `(1 | p | g)` in several
formulas correlates their `g` effects.

## Usage

``` r
mvbf(..., rescor = FALSE)

set_rescor(rescor_value = TRUE)
```

## Arguments

- ...:

  [`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) formulas,
  each with a family attached (or supply one `family` to
  [`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) for all
  of them).

- rescor:

  Model residual correlation between the responses (gaussian only).

- rescor_value:

  For `set_rescor()`: turn residual correlation on or off.

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
#> $y1_mu
#> (Intercept)           x 
#>   1.0775542   0.5517536 
#> 
#> $y1_sigma
#> (Intercept) 
#>   0.2721586 
#> 
#> $y2_mu
#> (Intercept)           x 
#>   2.0235294  -0.2793594 
#> 
#> $y2_sigma
#> (Intercept) 
#>   0.2223223 
#> 

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
#>   y1 1 | g + y2 1 | g [ID] 
#>               Name Std.Dev. y1.mu:(Intercept)
#>  y1.mu:(Intercept)  0.86060                  
#>  y2.mu:(Intercept)  0.61112             0.599
```
