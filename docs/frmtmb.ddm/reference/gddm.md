# The generalized drift-diffusion family

The generalized drift-diffusion model of Shinn, Lam and Murray (2020).
Evidence accumulates between two decision boundaries and the response
time is the first time it reaches one, as in
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md);
unlike
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md),
the drift may depend on the state of the accumulator and on a covariate,
and the boundaries may move. That generality costs the closed form:
there is no first-passage density to evaluate, so each likelihood
evaluation solves the Fokker-Planck equation forward in time and reads
the probability flux through each boundary.

## Usage

``` r
gddm(
  drift = gddm_drift_constant(),
  bound = gddm_bound_constant(),
  start = gddm_start_point(),
  lapse = c("none", "uniform"),
  control = gddm_control()
)
```

## Arguments

- drift:

  A drift component, or a list of them to be summed. The first must be a
  base term. See
  [gddm-drift](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md).

- bound:

  A boundary component. See
  [gddm-bound](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-bound.md).

- start:

  A starting-point component. See
  [gddm-start](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-start.md).

- lapse:

  `"none"`, or `"uniform"` to mix the first-passage density with a lapse
  distribution uniform over the modeled window and split evenly between
  the two responses. `"uniform"` adds the free parameter `lapse` on a
  logit link; the published models fix it rather than fitting it, which
  is spelled `bf(..., lapse = 0.05)`.

- control:

  Numerical controls, from
  [`gddm_control()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_control.md).

## Value

A `frmtmb_family`.

## The model

The accumulator `x` follows

\$\$dx = a(x, t) \\ dt + dW\$\$

from a starting distribution `X_0`, absorbed at moving boundaries \\\pm
B(t)\\. The response time is the absorption time plus a non-decision
time `ndt`. The drift `a`, the boundary `B` and the start `X_0` are
chosen by argument; each brings its own free parameters, and each of
those takes a formula like any other distributional parameter. See
[gddm-drift](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-drift.md),
[gddm-bound](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-bound.md)
and
[gddm-start](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm-start.md)
for the catalogue and for how to add to it.

The likelihood is solved after the substitution \\y = x / B(t)\\, which
pins the boundaries at \\\pm 1\\ and leaves the grid fixed while the
boundary collapses. That is what makes the model differentiable: the
usual treatment sandwiches a moving bound between two integer grid
indices, which makes the objective a function of where the bound falls
between nodes and rules out gradient-based fitting. With the walls
stationary the scheme is Crank-Nicolson.

## What the data must carry

Which boundary a trial ended at is data, and so is the condition a trial
belongs to. Two spellings carry the boundary, as they do for
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md):

    frm(bf(rt | dec(response) + vint(cond) ~ 1), family = gddm(), data = dat)
    frm(bf(rt | vint(upper, cond) ~ 1),          family = gddm(), data = dat)

`dec()` takes what brms takes, a factor whose second level is the upper
boundary or a 0/1 column. `upper` in the `vint()` spelling is the same
0/1. Two is all there is: one accumulator in one dimension between two
absorbing boundaries can end a trial at the upper wall or the lower wall
and nowhere else, so multi-alternative choice is outside this family,
and a third level is refused and pointed at
[`lba()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/lba.md)
rather than folded into one of the two.

`cond` is an integer labelling the distinct parameter settings in the
design: one solve serves every trial that shares one. Note where it
sits. `vint()` numbers its values positionally, so the condition is the
first `vint()` value when `dec()` carries the boundary and the second
when `vint()` carries it. The family reads whichever it is; you write
the pair in the order the two lines above show. A drift term that reads
a covariate needs `vreal()` as well, for example
`bf(rt | dec(response) + vint(cond) + vreal(coh) ~ 1)`.

**Every row sharing a condition must share every parameter value.** The
family cannot check this, because checking would mean comparing values
on the tape, so it is your side of the contract: build the index from
every variable that appears on the right-hand side of any of the
family's formulas, which is what
[`gddm_conditions()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_conditions.md)
does. What the family can check, and does, is that the `vreal()`
covariates are constant within a condition.

## Cost, honestly

This is much slower than
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)
and you should use
[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)
whenever it applies, which is whenever the drift is constant in the
state and in time and the boundaries do not move. There the
first-passage density is a known pair of series and costs arithmetic;
here every evaluation solves a PDE per condition. What you buy is the
class of models the analytic density cannot express at all: collapsing
boundaries, leaky or unstable integration, a drift that varies within a
trial. Cost scales with the number of conditions, so a design with many
distinct parameter settings is where this becomes painful. See
[`gddm_control()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_control.md).

## References

Shinn, M., Lam, N. H. and Murray, J. D. (2020). A flexible framework for
simulating and fitting generalized drift-diffusion models. *eLife*, 9,
e56938.

## See also

[`wiener()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/wiener.md)
for the analytic special case,
[`gddm_control()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_control.md)
for the grid,
[`gddm_conditions()`](https://aforren1.github.io/frmtmb/frmtmb.ddm/reference/gddm_conditions.md)
for building the condition index, and
[`vignette("ddm")`](https://aforren1.github.io/frmtmb/frmtmb.ddm/articles/ddm.md).

## Examples

``` r
# \donttest{
set.seed(1)
# a coarse grid, so that the example is quick; see gddm_control()
ctl <- gddm_control(t_max = 2, dt = 0.02, ny = 101)
dat <- gddm_simulate(400, mu = 2.5, bs = 3, ndt = 0.25, tau = 1,
                     bound = gddm_bound_exponential(), control = ctl)
dat$cond <- 1L
fit <- frm(bf(rt | vint(upper, cond) ~ 1, bias = 0.5),
           family = gddm(bound = gddm_bound_exponential(), control = ctl),
           data = dat)
fixef(fit)
#> $mu
#> (Intercept) 
#>    2.339251 
#> 
#> $bs
#> (Intercept) 
#>    1.273025 
#> 
#> $tau
#> (Intercept) 
#>  -0.2133438 
#> 
#> $bias
#> (Intercept) 
#>           0 
#> 
#> $ndt
#> (Intercept) 
#>   0.4489101 
#> 
# }
```
