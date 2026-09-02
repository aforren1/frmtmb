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
  events = NULL,
  event_scale = 1,
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

- events:

  Optional dosing table: a data.frame with columns `time`, `value` and
  `state`, and optional `group`, `method` and `duration`, or a function
  of no arguments returning one. See "Dosing events" below. `NULL` (the
  default) is a model driven only by its initial conditions.

- event_scale:

  A multiplier on every `events$value`, one value per observation
  (constant within group) or a single value shared by every group. This
  is the one estimated quantity that can reach a dose, and it is how a
  bioavailability is written. Only for a table whose rows are all
  `method = "add"`.

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

## Dosing events

`events` is a data.frame of doses, one row per dose, with columns:

- `group` (optional): which group the row applies to, matching the
  values of `group`. Leave the column out and the schedule applies to
  every group.

- `time`: when the dose happens. At or after `t0`.

- `state`: the state it goes into, by name (requires `states`) or by
  position, resolved exactly as `output` is. Optional for a one-state
  system.

- `value`: how much.

- `method` (optional, default `"add"`): `"add"` puts `value` into the
  state, `"replace"` sets the state to `value`, `"multiply"` scales it.
  These are the deSolve event methods.

- `duration` (optional, default `0`): a positive value makes the row an
  infusion, delivering `value` at the constant rate `value / duration`
  over `[time, time + duration]`. Infusions must use `"add"`.

Inside a `bf(nl = TRUE)` body, name the table:

    doses <- data.frame(time = seq(12, 48, by = 12), state = "depot",
                        value = 100)
    conc ~ frm_ode(pk_dyn, ..., events = doses)

A name in a nonlinear body is normally a request for a column of `data`,
and a column of that name still wins. A data.frame is not something a
column could hold, so `doses` is read from the formula environment
instead. Writing the table **inline**, or holding it in a function of no
arguments, works the same way and is what a schedule read at fit time
wants:

    conc ~ frm_ode(pk_dyn, ..., events = data.frame(
             time = seq(12, 48, by = 12), state = "depot", value = 100))

    schedule <- function() read.csv("doses.csv")
    conc ~ frm_ode(pk_dyn, ..., events = schedule)

In NONMEM terms an `"add"` row is a dosing record (`evid = 1`) with
`amt = value` into `cmt = state`; a row with `duration` is the same
record with `rate = amt / duration`. `frm_ode()` does not read NONMEM
column names, and there is no `evid` column: observation rows are the
rows of `data`, and dose rows are the rows of `events`, which is a
separate table. A NONMEM-shaped dataset has to be split into the two.

An observation at exactly a dose time reads the state **before** the
dose, which is the trough, matching both the deSolve convention and the
usual reading of a pre-dose sample. That includes an observation at `t0`
with a dose at `t0`: it reads `init`.

The doses are not handed to deSolve as events. `frm_ode()` splits the
integration at the event times and chains one solve per interval,
carrying the state across the break itself. That is a correctness
requirement, not a style choice: RTMBode solves an augmented system
carrying the derivatives of the states with respect to the parameters,
and a deSolve event jumps the state without jumping those derivatives,
which gives a wrong gradient for `"replace"` and `"multiply"` (measured
at 42% and 59% relative error). Splitting the solve is exact for every
method. It costs one solve per dosing interval per group, so an
intensively dosed design is proportionally slower.

## Doses that depend on a parameter

`events` is data: `value` is a numeric column, so it cannot hold an
estimated quantity. `event_scale` is the way in. It is one value per
observation (constant within group, like `init` and `parms`), and it
multiplies the `value` of every event in that group, so a
bioavailability written as a nonlinear parameter estimates a dose scale
that carries covariates and random effects like any other:

    frm_ode(pk_dyn, init = list(0, 0), times = time, parms = ...,
            group = id, events = doses, event_scale = plogis(logitF))

Because scaling only makes sense for a dose, `event_scale` is refused on
a table containing `"replace"` or `"multiply"` rows.

## Boundaries

Time-varying input other than dosing is out of scope, and the reason is
worth knowing. RTMBode tapes `dynamics` once, so `t` is an
automatic-differentiation value inside it. A branch on time
(`if (t < t_end) rate else 0`) raises "Comparison is generally unsafe
for AD types", and an
[`approxfun()`](https://rdrr.io/r/stats/approxfun.html) forcing table
silently returns the value at the taping point instead of failing.
Smooth arithmetic in `t` is fine. A piecewise-constant input belongs in
`events` as an infusion, where it is carried as a parameter over each
interval and differentiated exactly; deSolve's own `forcings` argument
is not reachable, because RTMBode's compiled derivative shim has no
forcing hook.

Estimated event times, lag times and inter-dose intervals are not
supported: the event times decide where the solve is split, which is
settled before the tape is built.

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

  # Repeated dosing: 100 into the depot every 12 hours. The dose at
  # time 0 is the initial condition, the rest are events.
  doses <- data.frame(time = c(12, 24, 36), state = "depot",
                      value = 100)
  frm_ode(pk_dyn, init = list(100, 0), times = c(6, 18, 30, 42),
          parms = list(1, 0.2, 10), states = c("depot", "central"),
          output = "central", events = doses)
}
#> [1] 3.733943 4.075490 4.106474 4.109285
```
