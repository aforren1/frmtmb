# The design of a linear predictor over the coefficient vector

`predict(se.fit = TRUE)` builds a matrix `A` with one row per prediction
and one column per contributing coefficient, forms `A V A'` and keeps
only its diagonal. Every delta-method quantity over a fitted curve needs
the whole thing: a contrast between two grids, an average marginal
effect with a correct standard error, a simultaneous band, a derivative,
the time of a peak. This returns the pieces so that an extension does
not have to rebuild `A` by perturbation, one
[`predict()`](https://rdrr.io/r/stats/predict.html) call per
coefficient.

## Usage

``` r
frm_lp_basis(
  object,
  newdata = NULL,
  dpar = NULL,
  resp = NULL,
  re.form = NULL,
  allow_new_levels = FALSE
)
```

## Arguments

- object:

  A fitted `frmtmb_fit`.

- newdata:

  Data to build the design at, or `NULL` for the training data.

- dpar, resp:

  The distributional parameter and response to take the linear predictor
  of. Both default the way
  [`predict()`](https://rdrr.io/r/stats/predict.html) defaults them.

- re.form:

  `NULL` keeps every random effect, `NA` drops them all, a one-sided
  formula keeps the ones it names.

- allow_new_levels:

  Whether a grouping level the fit never saw is allowed.

## Value

A list with

- `eta`:

  the linear predictor, exactly `predict(type = "link")`.

- `A`:

  `n x p`; `d eta / d coef`.

- `coef_pos`:

  length `p`; the rows of `V` the columns of `A` belong to, in `V`'s own
  order.

- `V`:

  `p x p`;
  [`frm_joint_cov()`](https://aforren1.github.io/frmtmb/reference/frm_joint_cov.md)
  subset to `coef_pos`.

- `coef_names`:

  length `p`; the labels of those rows.

- `extra_var`:

  length `n`; variance that is NOT coefficient uncertainty, kept
  separate rather than folded into `A V A'` because an exact `gp()`'s
  kriging variance and a new grouping level's marginal variance are not.

- `nonest`:

  length `n`; rows that load on a direction the rank-deficient design
  could not identify.

`var(eta)` is `rowSums((A %*% V) * A) + extra_var`, and
`predict(se.fit = TRUE)` is written that way.

## Details

The name follows
[`emmeans::emm_basis()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.html),
which is the same idea for the fixed block alone.

## A nonlinear body

For a nonlinear predictor `A` is a JACOBIAN rather than a design, and it
is computed by taping the body against the coefficients it reaches
through. That includes a
[`ps()`](https://aforren1.github.io/frmtmb/reference/ps.md) block, whose
coefficients enter the body through a spline evaluated at an argument
the parameters move. `predict(se.fit = TRUE)` stays refused for a
nonlinear predictor; this is the route.

The Jacobian is exact, and the delta method built on it is still a
first-order approximation, which for a warped curve is a stronger
assumption than it is for a linear one. `allow_new_levels = TRUE` is
refused there, and so is a contributing exact `gp()`, because neither
variance has a chain rule through the body that has been measured.

## A reduced-rank block

An `rr()` block's loadings live in `theta`, so a design over `(beta, b)`
alone is INCOMPLETE. `A` carries the loading columns too, through
`rr_jacobians()`, and `coef_pos` names their `theta` rows; a caller
therefore gets the whole delta method rather than discovering a missing
piece.

## See also

[`frm_joint_cov()`](https://aforren1.github.io/frmtmb/reference/frm_joint_cov.md)
for the covariance alone,
[frmtmb-extension-api](https://aforren1.github.io/frmtmb/reference/frmtmb-extension-api.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(120), g = factor(rep(1:12, each = 10)))
dd$y <- rnorm(120, 1 + 2 * dd$x + rnorm(12, 0, 0.5)[dd$g], 0.4)
fit <- frm(bf(y ~ x + (1 | g)), data = dd)
nd <- data.frame(x = c(-1, 0, 1), g = factor(1, levels = levels(dd$g)))
lb <- frm_lp_basis(fit, newdata = nd, re.form = NA)
str(lb$A)
#>  num [1:3, 1:2] 1 1 1 -1 0 1
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:3] "1" "2" "3"
#>   ..$ : chr [1:2] "(Intercept)" "x"
lb$coef_names
#> [1] "beta.(Intercept)" "beta.x"          

# the covariance of the WHOLE grid, which predict() reduces to its
# diagonal
Sigma <- lb$A %*% lb$V %*% t(lb$A)
all.equal(sqrt(diag(Sigma)),
          predict(fit, newdata = nd, re.form = NA, se.fit = TRUE)$se.fit)
#> [1] TRUE
```
