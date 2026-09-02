# Solve an ODE once per group inside a nonlinear predictor

Evaluates a system of ordinary differential equations separately for
each group of rows and returns the solution aligned with the rows of the
data, so that a compartment model can be written directly in the body of
a `bf(..., nl = TRUE)` formula. The dynamics parameters and the initial
states are ordinary nonlinear parameters, which means they take fixed
effects, random effects and covariates like any other linear predictor,
and the Laplace approximation is exact through the solver's adjoint.

## Usage

``` r
frm_ode(
  dynamics,
  init,
  times,
  parms,
  group = NULL,
  output = NULL,
  states = NULL,
  t0 = 0,
  method = "lsoda",
  atol = 1e-08,
  rtol = 1e-08,
  on_error = c("penalize", "error"),
  penalty = 1e+06,
  ...
)
```

## Arguments

- dynamics:

  A function `function(t, y, parms)` giving the derivatives, following
  the deSolve convention: `t` is the scalar time, `y` the state vector,
  `parms` the parameter vector, and the return value is `list(dydt)`. A
  bare derivative vector is also accepted. Index `y` and `parms` by
  position; use `"c" <- RTMB::ADoverload("c")` inside the function so
  that [`c()`](https://rdrr.io/r/base/c.html) keeps the
  automatic-differentiation class.

- init:

  Initial states, one column per state: a list, a matrix, or a single
  vector for a one-state system. Each column is either one value per
  observation (constant within group) or one value shared by every
  group. The number of columns sets the number of states.

- times:

  Observation times, one per row of the data.

- parms:

  Dynamics parameters, one column per parameter, in the order `dynamics`
  expects them. Same shape rules as `init`.

- group:

  Grouping vector naming the unit that owns one system, one value per
  row. `NULL` treats the whole data as one group.

- output:

  Which states to return: `NULL` (the default) returns every state, an
  integer or character vector selects some. A single selected state is
  returned as a vector, otherwise a matrix with one column per selected
  state. Character selection requires `states`.

- states:

  Optional state names, one per column of `init`.

- t0:

  Initial time, a scalar or one value per row (constant within group).
  Every observation time must be at or after it.

- method:

  Integrator, passed to
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html). Must be
  adaptive; fixed-step integrators such as `"rk4"` and `"euler"` return
  a different likelihood and are refused.

- atol, rtol:

  Absolute and relative solver tolerances.

- on_error:

  What to do about a solve that fails: `"penalize"` (the default) fills
  that group's rows with `penalty` and warns, naming the group;
  `"error"` stops instead, also naming it. Read "Failed solves" below
  first: on the automatic-differentiation tape most failures cannot be
  detected at all, so neither setting has the reach it appears to have
  during a fit.

- penalty:

  The filler value for `on_error = "penalize"`. It is on the scale of
  the nonlinear body's result, before the response link, so lower it for
  a model whose link exponentiates.

- ...:

  Further arguments for
  [`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html).

## Value

A numeric vector of length `nrow(data)` when one state is selected,
otherwise a matrix with `nrow(data)` rows. On the
automatic-differentiation tape both carry the `advector` class.

## What the group is

`group` names the unit that owns one ODE system: a subject in a
population pharmacokinetic model, a reactor, a patient. Every row of one
group is one observation of that group's trajectory, at the time given
by `times`. One solve is performed per group, over that group's own
sorted times, and the results are scattered back into the input row
order. Ragged designs, unsorted rows, repeated times and an observation
at `t0` itself are all fine.

Groups are never stacked into one large system. That is a hard
constraint, not an implementation detail: the second-order derivative
path the Laplace approximation needs returns `NaN` above roughly eight
states in a single system, so a stacked solve gives silently wrong or
missing gradients. `frm_ode()` warns when one system alone exceeds that
many states.

## What is constant within a group

`init` and `parms` are read off each group's **first row**. A dynamics
parameter must therefore not vary inside a group: a covariate that
changes between an early and a late observation of the same subject
describes a model this helper cannot solve. Such a covariate is refused,
by name, rather than silently ignored. Covariates that are constant
within a group (a subject's weight, dose or treatment arm) are the
intended case and are unrestricted.

## Failed solves

`on_error` and `penalty` reach less than they look like they do, and the
difference matters:

- **While fitting**, the body is evaluated on the
  automatic-differentiation tape.
  [`RTMBode::ode()`](https://rdrr.io/pkg/RTMBode/man/ode.html) returns
  an `advector` and does not raise an R error when the trajectory goes
  bad, and a non-finite AD value cannot be tested for (RTMB refuses
  comparison on AD types). So a diverging region of the parameter space
  is **not** caught here. It surfaces instead as the optimizer's own
  `NA/NaN function evaluation` or `NA/NaN gradient evaluation` warning,
  and `on_error = "error"` will not name the group. The two failures
  that are still caught on the tape are an integrator that gives up and
  returns fewer time points than were asked for, and any R error raised
  by `dynamics` itself.

- **Everywhere else** the evaluation is ordinary numeric arithmetic, and
  every check applies:
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html), a direct
  call, and a fit whose `init` and `parms` contain no estimated
  parameter (which is evaluated once, numerically, as the tape is
  built). Here `on_error = "penalize"` writes `penalty` into the failed
  group's rows, and `on_error = "error"` stops and names the group.

A penalty is never written silently: `frm_ode()` warns, naming the
groups, and
[`frm_ode_failures()`](https://aforren1.github.io/frmtmb/reference/frm_ode_failures.md)
reports them afterwards. Treat those rows as missing, not as
predictions.

If a fit reports `NA/NaN gradient evaluation`, run
[`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
on it, then call `frm_ode()` directly at the suspect parameter values
with `on_error = "error"`: numerically it will name the group that
cannot be solved.

Solver warnings from deSolve ("corrector convergence failed repeatedly",
"exceeded maxsteps") during a fit come from the optimizer's probing
steps and are usually not fatal. Judge the fit by the gradient at the
optimum, not by whether the solver complained.

## Boundaries

Dosing event tables (`evid`/`amt` records, repeated doses, infusions)
are out of scope. Only models driven by their initial conditions are
supported. For event-driven population pharmacokinetics use a dedicated
tool such as `nlmixr2`.

`predict(se.fit = TRUE)` is not available for a nonlinear predictor,
including one containing `frm_ode()`; request a nonlinear parameter with
`predict(dpar = )` instead.

## Installation

`frm_ode()` needs RTMBode, which is not on CRAN:

    install.packages("RTMBode", repos = c(
      "https://kaskr.r-universe.dev",
      "https://cloud.r-project.org"))

## See also

[`frm_ode_failures()`](https://aforren1.github.io/frmtmb/reference/frm_ode_failures.md)
for the groups a penalty was written into,
[`bf()`](https://aforren1.github.io/frmtmb/reference/bf.md) for the
nonlinear formula grammar, and
[`vignette("ode")`](https://aforren1.github.io/frmtmb/articles/ode.md)
for a worked population pharmacokinetic model.

## Examples

``` r
# One-compartment oral pharmacokinetics with between-subject
# variability on the absorption and elimination rates.
#   dA/dt = -ka A            A(0) = dose
#   dC/dt =  ka A / V - ke C C(0) = 0
pk_dyn <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}

set.seed(2026)
tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
n_id <- 6
dd <- data.frame(id = factor(rep(seq_len(n_id), each = length(tt))),
                 time = rep(tt, n_id), dose = 100)
ka <- exp(rnorm(n_id, 0, 0.3))[as.integer(dd$id)]
ke <- exp(rnorm(n_id, log(0.2), 0.25))[as.integer(dd$id)]
dd$conc <- 100 * ka / (10 * (ka - ke)) *
  (exp(-ke * dd$time) - exp(-ka * dd$time)) + rnorm(nrow(dd), 0, 0.3)

if (requireNamespace("RTMBode", quietly = TRUE)) {
  # \donttest{
  fit <- frm(
    bf(conc ~ frm_ode(pk_dyn,
                      init   = list(dose, 0),
                      times  = time,
                      parms  = list(exp(lka), exp(lke), exp(lV)),
                      group  = id,
                      states = c("depot", "central"),
                      output = "central"),
       lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
      gaussian(),
    data = dd, start = list(beta = c(0, log(0.25), log(8))))
  fixef(fit)
  # }
}
#> $lka
#> (Intercept) 
#>  -0.2014577 
#> 
#> $lke
#> (Intercept) 
#>   -1.743724 
#> 
#> $lV
#> (Intercept) 
#>    2.298439 
#> 
#> $mu
#> numeric(0)
#> 
#> $sigma
#> (Intercept) 
#>   -1.293758 
#> 
```
