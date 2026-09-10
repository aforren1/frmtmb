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
  parms = list(),
  group = NULL,
  output = NULL,
  states = NULL,
  t0 = 0,
  events = NULL,
  event_scale = 1,
  tv = NULL,
  tv_break = NULL,
  n_ss = 20L,
  ss_tol = 1e-06,
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
  position. RTMB's tape-safe [`c()`](https://rdrr.io/r/base/c.html),
  `[<-` and `diag<-` are put in scope automatically, so the function
  needs no `"c" <- RTMB::ADoverload("c")` boilerplate; only a HELPER it
  calls, defined elsewhere, still needs its own (lexical scope does not
  travel into other functions).

- init:

  Initial states, one column per state: a list, a matrix, or a single
  vector for a one-state system. Each column is either one value per
  observation (constant within group) or one value shared by every
  group. The number of columns sets the number of states.

- times:

  Observation times, one per row of the data.

- parms:

  Dynamics parameters, one column per parameter, in the order `dynamics`
  expects them. Same shape rules as `init`. May be empty when every
  input is given in `tv`.

- group:

  Grouping vector naming the unit that owns one system, one value per
  row. `NULL` treats the whole data as one group.

- output:

  Which states to return: `NULL` (the default) returns every state, an
  integer or character vector selects some. A single selected state is
  returned as a vector, otherwise a matrix with one column per selected
  state. Character selection requires `states`. One value per row -
  which needs more rows than there are states - instead says which state
  each row reads, over one shared solve; see "One solve, a different
  state per row".

- states:

  Optional state names, one per column of `init`.

- t0:

  Initial time, a scalar or one value per row (constant within group).
  Every observation time must be at or after it.

- events:

  Optional dosing table: a data.frame with columns `time`, `value` and
  `state`, and optional `group`, `method`, `duration`, `ii`, `addl` and
  `ss`, or a function of no arguments returning one. See "Dosing events"
  below. `NULL` (the default) is a model driven only by its initial
  conditions.

- event_scale:

  A multiplier on every `events$value`, one value per observation
  (constant within group) or a single value shared by every group. This
  is the one estimated quantity that can reach a dose, and it is how a
  bioavailability is written. Only for a table whose rows are all
  `method = "add"` or `"reset"`; a reset is not scaled.

- tv:

  Dynamics inputs that vary with time, one column per input, appended to
  `parms` in the vector `dynamics` receives. Unlike `parms` these may
  vary within a group, as a step function of time: the solve is split at
  each change point. See "Time-varying inputs".

- tv_break:

  A data column, one value per row, whose changes within a group are the
  change points of every `tv` column. Required when a `tv` column
  carries estimated values, because those cannot be compared. Optional
  otherwise, when the change points are read off the `tv` columns
  themselves.

- n_ss:

  How many dosing cycles a steady-state run-in simulates before an
  `events` row marked `ss = TRUE`. The steady state is approached, not
  solved for; see "Steady-state dosing" for the size of the
  approximation.

- ss_tol:

  Relative change between the last two run-in cycles that `frm_ode()`
  will accept without warning. Checked on the numeric path only.

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

A numeric vector of length `nrow(data)` when one state is selected, or
when `output` selects one state per row, otherwise a matrix with
`nrow(data)` rows. On the automatic-differentiation tape both carry the
`advector` class.

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

`tv` is the exception, and the next section is about it.

## Time-varying inputs

`tv` carries dynamics inputs that DO change inside a group, as long as
they change in steps: a creatinine clearance measured again at each
visit, a dose-dependent rate that switches at a protocol amendment, a
temperature held at one level and then another. The solve is split at
each change point, exactly as it is split at a dose time, and each
segment's dynamics see that segment's value.

The values follow `parms` in the vector `dynamics` is handed, so with
`parms = list(ka, V)` and `tv = list(ke)` the derivative function reads
`p[1]`, `p[2]` and `p[3]`. Nothing in `dynamics` learns that `p[3]` ever
changed: within one segment it is an ordinary constant, and because it
is a parameter it is a tape input, so an estimated time-varying value is
differentiated exactly.

The step function is **last observation carried forward**, which is
rxode2's convention for covariates. A row's value is in force from that
row's time until the next change, and the first row's value reaches back
to `t0`. A consequence worth knowing: the state is continuous across a
change point, because a covariate moves the derivative and not the
state, so an observation exactly at a change point reads what the
PRE-change dynamics produced. The difference shows up only afterwards.
(This is unlike a dose, which jumps the state, and where the same
observation reads the trough.)

The values may be estimated; the change points may not. They decide
where the solve is split, which is settled before the tape is built.
When a `tv` column carries estimated values there is nothing to
compare - RTMB refuses comparison on AD types - so `tv_break` has to
name the data column whose changes within a group mark them:

    bf(conc ~ frm_ode(pk_dyn, init = list(dose, 0), times = time,
                      parms = list(exp(lka), exp(lV)),
                      tv = list(exp(lke)), tv_break = visit,
                      group = id, output = 2L),
       lka ~ 1, lV ~ 1, lke ~ 1 + crcl, nl = TRUE)

Here `lke` is a linear predictor that varies within a subject, which
`parms` would refuse. `tv_break = visit` says the value is constant
inside a visit, and a plain-numeric `tv` column that disagrees with that
is refused. When every `tv` column is plain data, `tv_break` can be left
out and the change points are read off the columns.

A time-varying input other than a step function is not available. See
"Boundaries".

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
[`frm_ode_failures()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode_failures.md)
reports them afterwards. Treat those rows as missing, not as
predictions.

If a fit reports `NA/NaN gradient evaluation`, run `diagnose()` on it,
then call `frm_ode()` directly at the suspect parameter values with
`on_error = "error"`: numerically it will name the group that cannot be
solved.

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
  system, and ignored on a `"reset"` row, which names no compartment
  because it sets every one of them.

- `value`: how much.

- `method` (optional, default `"add"`): `"add"` puts `value` into the
  state, `"replace"` sets the state to `value`, `"multiply"` scales it -
  the deSolve event methods - and `"reset"` sets **every** state to
  `value`.

- `duration` (optional, default `0`): a positive value makes the row an
  infusion, delivering `value` at the constant rate `value / duration`
  over `[time, time + duration]`. Infusions must use `"add"`.

- `ii` (optional, default `0`): the interdose interval, for `addl` and
  `ss`.

- `addl` (optional, default `0`): how many further doses follow this
  one, at `time + ii`, `time + 2 * ii`, and so on. The rows are written
  out internally, so the table means exactly what the hand-expanded one
  means.

- `ss` (optional, default `FALSE`): the row's cycle has already reached
  steady state. See "Steady-state dosing".

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
record with `rate = amt / duration`; `ii` and `addl` are spelled the
same way there and in rxode2. A `"reset"` row is `evid = 3`, which both
NONMEM and rxode2 read as "set every compartment to zero and carry on" -
write `value = 0` for that reading. `evid = 4`, a reset followed by a
dose, is a `"reset"` row and an `"add"` row at the same time; the reset
is always applied first, so the order of the two rows in the table does
not matter.

`frm_ode()` does not read NONMEM column names, and there is no `evid`
column: observation rows are the rows of `data`, and dose rows are the
rows of `events`, which is a separate table. A NONMEM-shaped dataset has
to be split into the two.

An observation at exactly a dose time reads the state **before** the
dose, which is the trough, matching both the deSolve convention and the
usual reading of a pre-dose sample. That includes an observation at `t0`
with a dose at `t0`: it reads `init`. The one exception is an `ss` row,
which is not a jump from anything.

## Steady-state dosing

`ss = TRUE` on an `events` row says the system has already been given
this dose every `ii` for long enough to settle. It is reached by
simulation, not by a closed form, so it works for any dynamics: every
compartment is set to zero (NONMEM's and rxode2's reading of a
steady-state record, and the reason an `ss` row overrides `init`), and
the cycle is then repeated `n_ss` times before the record's time. The
record's own dose is applied afterwards, so an observation at the
record's time reads the steady-state **trough**.

    # 100 into the depot every 12 hours, already at steady state at t = 0
    data.frame(time = 0, state = "depot", value = 100, ii = 12, ss = TRUE)

`n_ss` is the honest part of this. A real steady state is the limit of
repeating until the cycle-start state stops moving, and that test
branches on a value, which the tape cannot do: the number of solves has
to be fixed before the tape is built. So `n_ss` cycles is an
**approximation**, and its error is geometric - for linear kinetics the
shortfall after `n` cycles is the accumulation factor to the power `n`,
which is `exp(-n * k * ii)` for a one-compartment system. The rate `k`
in that expression is the SLOWEST disposition eigenvalue, `lambda_z`, so
the shortfall is `exp(-n_ss * lambda_z * ii)` and it is set by the
terminal half-life measured in dosing intervals. **Choose `n_ss` so that
`n_ss * lambda_z * ii` is at least 20**, which puts the shortfall at
2e-09. The default of 20 does that only when the terminal half-life is
under about three dosing intervals. Measured, on a two-compartment oral
model at `n_ss = 20`: a 23 hour half-life dosed every 8 hours is 4e-03
short, a 107 hour half-life dosed daily is 1.2e-02 short, and a 265 hour
half-life dosed daily is 9e-02 short. A half-life of ten dosing
intervals is 25 percent short. It moves estimates: on one 30-subject
dataset simulated from the exact steady state, fitting at `n_ss = 20`
rather than at the limit moved `k21` by a factor of 3.2 and `ke` by 11
percent.

Off the tape - a direct call,
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html), a body holding no
estimated parameter - the last two cycles are compared and `frm_ode()`
warns when they still differ by more than `ss_tol`. During a fit that
check cannot run at all. **Do not choose `n_ss` from that warning**: it
reports the CYCLE-TO-CYCLE movement, which understates the distance to
the limit by about `1 / (lambda_z * ii)` and so understates it most
exactly where the error is largest (measured at 3.1x, 6.6x and 10.8x as
the half-life grows). Compare `n_ss` against `2 * n_ss` numerically
instead, or, for a linear compartment model, use
[`frm_lincmt()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_lincmt.md),
whose default sums the series and has nothing to truncate.

The cost is `n_ss` extra solves per group (two per cycle for an
infusion), so a steady-state population fit is several times a plain
one. One `ss` row per group is allowed; write later doses out with `ii`
and `addl`.

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
a table containing `"replace"` or `"multiply"` rows. A `"reset"` row is
allowed beside scaled doses and is simply not scaled: it names a level
for the states, not an amount.

## One solve, a different state per row

`output` normally selects states. Given one value per row instead, it
says which state THAT row reads, and the whole group still takes one
solve. This is the spelling for a parent and its metabolite reported in
one assay column, or for two species of a predator-prey series stacked
long:

    bf(y ~ frm_ode(pm_dyn, init = list(dose, 0, 0), times = time,
                   parms = list(exp(lka), exp(lkm), exp(lke)),
                   group = id, output = cmt),          # cmt is 2 or 3
       lka ~ 1, lkm ~ 1, lke ~ 1, nl = TRUE)

The two readings are told apart by length: a per-row selection needs one
value per row AND more rows than the system has states. A selection of
states is shorter than that. The per-row values are resolved by name or
position exactly as a scalar `output` is, and they have to be data - a
state index cannot be differentiated.

## Boundaries

Time-varying input has to be piecewise constant, and the reason is worth
knowing. RTMBode tapes `dynamics` once, so `t` is an
automatic-differentiation value inside it. A branch on time
(`if (t < t_end) rate else 0`) raises "Comparison is generally unsafe
for AD types", and an
[`approxfun()`](https://rdrr.io/r/stats/approxfun.html) forcing table
silently returns the value at the taping point instead of failing.
Smooth arithmetic in `t` is fine. A step function belongs in `tv`, or in
`events` as an infusion, where it is carried as a parameter over each
interval and differentiated exactly; deSolve's own `forcings` argument
is not reachable, because RTMBode's compiled derivative shim has no
forcing hook. Linear interpolation between measured covariate values,
which rxode2 offers, is not available: it would need the value inside
the segment to depend on `t`.

Estimated event times, lag times, inter-dose intervals and `tv` change
points are not supported: they decide where the solve is split, which is
settled before the tape is built.

`predict(se.fit = TRUE)` is not available for a nonlinear predictor,
including one containing `frm_ode()`; request a nonlinear parameter with
`predict(dpar = )` instead.

## Installation

`frm_ode()` needs RTMBode, which is not on CRAN:

    install.packages("RTMBode", repos = c(
      "https://kaskr.r-universe.dev",
      "https://cloud.r-project.org"))

## Sampling an ODE fit

[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
and
[`frmtmb.sample::as_tmbstan()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/as_tmbstan.html)
both refuse a model that contains `frm_ode()`, and the refusal is a
registered row of the compatibility registry
([`frmtmb::frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.html))
rather than a guess. Both doors are guarded, and the second one was the
worse of the two: it returned an empty `stanfit` with no error at all,
and the abort behind it left rstan unusable for the rest of the R
session.

The defect is upstream and is reproduced in
`dev/upstream/rtmbode-issues.md` without frmtmb, in bare RTMB, RTMBode
and deSolve.
[`RTMBode::ode()`](https://rdrr.io/pkg/RTMBode/man/ode.html) calls
[`deSolve::ode()`](https://rdrr.io/pkg/deSolve/man/ode.html) with no
guard. deSolve answers an extreme parameter either by raising
`illegal input detected before taking any integration steps` or by
returning fewer rows than were asked for, which the fixed-length adjoint
node reports as `Wrong output length`. An optimizer shortens its step
and a sampler rejects its proposal only when a failure arrives as `NaN`;
as an error it is fatal. Stan reaches those parameters by construction,
because its first warmup step is of size 1 on the unconstrained scale,
so the chain aborts at iteration 1 even when it starts at the fitted
optimum. The abort also crosses Stan's C++ boundary and leaves rstan's
nested autodiff arena unbalanced for the rest of the R session.

Two ways forward. Fit by maximum likelihood and read the Wald intervals,
which is what the rest of this page describes; or apply the three-patch
series in `dev/upstream/patches/` to RTMBode, after which the same
models sample. The refusal disappears on its own once a fixed RTMBode is
released and this package drops the row.

## See also

[`frm_ode_failures()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode_failures.md)
for the groups a penalty was written into,
[`frmtmb::bf()`](https://aforren1.github.io/frmtmb/reference/bf.html)
for the nonlinear formula grammar, and
[`vignette("ode")`](https://aforren1.github.io/frmtmb/frmtmb.ode/articles/ode.md)
for a worked population pharmacokinetic model.

## Examples

``` r
# One-compartment oral pharmacokinetics with between-subject
# variability on the absorption and elimination rates.
#   dA/dt = -ka A            A(0) = dose
#   dC/dt =  ka A / V - ke C C(0) = 0
pk_dyn <- function(t, y, p) {
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

  # The same schedule written compactly, and already at steady state
  frm_ode(pk_dyn, init = list(0, 0), times = c(6, 18, 30, 42),
          parms = list(1, 0.2, 10), states = c("depot", "central"),
          output = "central",
          events = data.frame(time = 0, state = "depot", value = 100,
                              ii = 12, ss = TRUE))

  # An elimination rate that doubles after hour 6, carried forward
  # from the row it appears on. `tv` values follow `parms`, so this
  # dynamics reads ka at p[1], V at p[2] and the time-varying ke at
  # p[3].
  pk_tv <- function(t, y, p) {
    list(c(-p[1] * y[1], p[1] * y[1] / p[2] - p[3] * y[2]))
  }
  tt <- c(1, 3, 6, 9, 12)
  frm_ode(pk_tv, init = list(100, 0), times = tt, parms = list(1, 10),
          tv = list(ifelse(tt < 6, 0.2, 0.4)),
          states = c("depot", "central"), output = "central")
}
#> [1] 5.6356414 6.2378071 3.7339432 1.1350284 0.3423811
```
