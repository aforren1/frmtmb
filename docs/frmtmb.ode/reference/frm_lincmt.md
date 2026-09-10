# Analytic one- to three-compartment pharmacokinetics

Evaluates a linear compartment model in closed form, as an alternative
to
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
for the systems where one exists. The dosing grammar, the grouping and
the place in a `bf(nl = TRUE)` formula are
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)'s,
so a linear model written for the solver becomes a `frm_lincmt()` call
with the same `events` table. There is no integrator: the states are a
superposition of the impulse responses of the doses that precede each
observation, and a steady-state record is the limit of that
superposition rather than a run-in.

## Usage

``` r
frm_lincmt(
  parms,
  times,
  group = NULL,
  ncmt = 1L,
  depot = FALSE,
  output = c("conc", "central", "depot"),
  t0 = 0,
  init = NULL,
  events = NULL,
  event_scale = 1,
  n_ss = Inf,
  tv = NULL,
  tv_break = NULL
)
```

## Arguments

- parms:

  A named list of parameters. See "Parameters". Inside a `bf(nl = TRUE)`
  body a bare NAME is a request for a column of the data. The arguments
  that are not columns (`ncmt`, `depot`, `output`, `n_ss`) therefore
  have to be written as literals there rather than held in a variable.

- times:

  Observation times, one per row of the data.

- group:

  Grouping vector naming the unit that owns one system, one value per
  row. `NULL` treats the whole data as one group.

- ncmt:

  Number of compartments in the disposition model: 1, 2 or 3. The depot,
  if any, is not counted.

- depot:

  Whether a depot compartment feeds the central one by first-order
  absorption at rate `ka`.

- output:

  What to return: `"conc"` (the default) is the central amount divided
  by `V`, `"central"` and `"depot"` are amounts. One column only;
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)'s
  per-row `output` is refused by name.

- t0:

  Initial time, a scalar or one value per row (constant within group).
  Every observation time must be at or after it.

- init:

  Amounts present at `t0`, as a named list with elements `depot` and
  `central`, each one value per observation (constant within group) or a
  single value. `NULL` starts empty.

- events:

  Optional dosing table, exactly
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)'s.
  See "Dosing".

- event_scale:

  A multiplier on every `events$value`, one value per observation
  (constant within group) or a single shared value.

- n_ss:

  How many dosing cycles a steady-state record carries. `Inf`, the
  default, is the exact geometric limit and costs one term. A whole
  number writes out that many cycles, which is what
  [`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
  simulates.

- tv, tv_break:

  Refused, and present so that the refusal names them. See "What is
  refused, and why".

## Value

A numeric vector of length `nrow(data)`. On the automatic
differentiation tape it carries the `advector` class.

## The model

One central compartment, up to two peripheral compartments exchanging
with it, and an optional depot feeding it by first-order absorption:

\$\$ \frac{dA_d}{dt} = -k_a A_d, \qquad \frac{dA_1}{dt} = k_a A_d -
(k_e + k\_{12} + k\_{13}) A_1 + k\_{21} A_2 + k\_{31} A_3 \$\$ \$\$
\frac{dA_j}{dt} = k\_{1j} A_1 - k\_{j1} A_j \$\$

The compartments are named `"depot"`, `"central"`, `"peripheral1"` and
`"peripheral2"`, in that order, and `events$state` may use those names
or the positions they stand at. States are AMOUNTS; the default
`output = "conc"` divides the central amount by `V`.

## Parameters

`parms` is a NAMED list, because which parameters a model needs depends
on `ncmt` and `depot`. Two spellings are read, and mixing them is
refused:

- rate constants: `ke`, `V`, and `ka`, `k12`, `k21`, `k13`, `k31` as the
  shape requires;

- clearances: `CL`, `V`, and `ka`, `Q2`, `V2`, `Q3`, `V3`, with \\k_e =
  CL/V\\, \\k\_{12} = Q_2/V\\, \\k\_{21} = Q_2/V_2\\.

Each element is one value per observation, constant within a group, or a
single value shared by every group, exactly as
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)'s
`parms` is. They are ordinary nonlinear parameters, so they take fixed
effects, random effects and covariates.

## Dosing

`events` is
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)'s
table and goes through the same validation, so `time`, `value`, `state`,
`method`, `duration`, `ii`, `addl`, `ss` and `group` mean what
[`?frm_ode`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
says they mean, including that an observation at a dose time reads the
trough. `event_scale` is the same estimated multiplier on every dose,
which is how a bioavailability is written.

A steady-state row is exact here.
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
approaches a steady state by simulating `n_ss` cycles; `frm_lincmt()`
sums the geometric series, so the default `n_ss = Inf` is the limit
itself and costs one term. A finite `n_ss` writes out that many cycles
as ordinary doses, which is what reproduces a
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
fit exactly and is otherwise not worth paying for. The limit is what it
says: a drug that is not eliminated over its dosing interval accumulates
without bound, and the series says so, where a run-in of twenty cycles
reports twenty cycles worth. Nothing overflows, so an optimizer that
wanders into that region gets a large finite objective and backs out of
it.

## What is refused, and why

Superposition adds impulse responses. It never forms the state vector,
so anything that reads or sets that vector is refused by name rather
than approximated:

- `method = "replace"` and `method = "multiply"`. A `"reset"` to zero IS
  supported, because "forget every dose before this time" needs no
  state; a reset to a non-zero level is refused.

- a dose into a peripheral compartment, and `output` naming one.

- an infusion into the depot, which is zero-order absorption.

- `tv` and `tv_break`. A rate constant that changes with time makes the
  system time-varying, and superposition over the whole history is then
  wrong rather than approximate.

- an infusion that is still running at a `ss` or `"reset"` time.

[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
fits every one of them.

## Accuracy

The textbook closed form divides by differences of rate constants and
loses all its precision where two of them meet: \\k_a = k_e\\ is the
ordinary flip-flop case of an oral model, and it is also the starting
value of any fit that starts two log rates at the same number. **The
ABSORPTION cancellation is removed completely**: every response is built
from `(exp(-p) - exp(-q)) / (q - p)` evaluated as
`exp(-max(p, q)) (1 - exp(-|p - q|)) / |p - q|`, which is exact at
`p == q` and never overflows, where the textbook difference of
exponentials returns `NaN` there and has lost most of its digits by a
separation of 1e-12. The disposition coefficients of a mammillary model
are non-negative and sum to one, so the sum of exponentials is a
positive combination.

One cancellation is NOT removed, and it is worth knowing which. The
smaller of the two disposition coefficients of a two-compartment model
is spelled `1 - c1`, and when the peripheral compartment is nearly
decoupled that subtraction loses about as many digits as `k12 k21 / g^2`
has decades, much as the textbook spelling does. It is bounded where it
matters and unbounded where it does not, and the two statements are the
same statement: the ABSOLUTE error is held at a few parts in 1e10 of the
trajectory's own maximum, so the POINTWISE relative error at a point
`10^-d` below `Cmax` is at most about that times `10^d`. An observation
within four decades of `Cmax` therefore carries at most about 1e-11 and
the 1e-8 the feature is held to is crossed only beyond seven decades
below `Cmax`, which is far below any assay's limit of quantification.

Measured: over 118 schedules the disagreement with
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
run at `atol = rtol = 1e-12` is at most 2.6e-12 of the trajectory's own
scale, and tightening the solver drives it down rather than leaving it,
so it is the solver's tolerance rather than the closed form's error.
Against a 240-bit reference that never forms an eigenvalue, the worst
relative error through every coalescence down to exact equality is
7.2e-15, and over a seven-decade box of rate constants the worst error
relative to the trajectory's maximum found by random search is 3.6e-10.
The repository's `dev/lincmt-findings.md` has the tables.

## Sampling

[`frmtmb.sample::frm_sample()`](https://aforren1.github.io/frmtmb/frmtmb.sample/reference/frm_sample.html)
refuses a model containing
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md),
because RTMBode calls deSolve unguarded and the chain aborts at warmup
iteration 1. `frm_lincmt()` calls neither, carries no refusal in
[`frmtmb::frm_compat()`](https://aforren1.github.io/frmtmb/reference/frm_compat.html),
and samples.

## Cost

One term per (observation, dose) pair, plus one term per group for a
steady-state record, all evaluated in one vectorized pass over the whole
data. There is no per-group loop on the tape and no solver, so the cost
is linear in the number of doses each observation follows. A schedule
with hundreds of `addl` doses before each observation is the case where
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)'s
per-segment cost can win.

## The accurate domain

**Up to a six-decade spread of rate constants, nothing degrades.** Over
60 random three-compartment draws per cell, the worst disagreement with
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
at `atol = rtol = 1e-12` is 9.9e-11 at a two-decade spread, 1.3e-09 at
four and 1.4e-07 at six, and in every one of those cells it is BELOW
what the solver's own tolerance costs on the same draws (1.3e-07,
1.6e-06, 1.9e-05). The direction an optimizer actually travels is
cleaner still: driving the four peripheral rate constants from one to
twelve decades below the absorption rate, which is what a fit does to a
compartment the data do not support, leaves the worst error at 6.9e-11
and does not degrade at all.

**Beyond an eight-decade spread that includes the FAST rates, both the
value and the gradient go wrong.** There the trajectory can be wrong by
up to 1.3e-04 of its own maximum against a 300-bit reference, with
nothing warning, and the gradient can be `NaN`. That region needs a
terminal half-life of seconds beside an absorption rate of hundreds per
hour, which is not a drug; and
[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
is not a usable reference there either, because DLSODA reaches its step
limit and disagrees with itself by 1e+05 and more.

## Boundary

A three-compartment model returns `NaN` in the GRADIENT, though not in
the value, when the cubic's two smaller roots have lost their precision
against each other. `lincmt_disp()` solves the cubic through
[`acos()`](https://rdrr.io/r/base/Trig.html), and a double root puts its
argument at exactly one, where
[`acos()`](https://rdrr.io/r/base/Trig.html) has an infinite derivative.

**The derivative itself is finite there.** What `frm_lincmt()` returns
is a symmetric function of the three roots, hence a function of the
characteristic polynomial's coefficients, which are polynomials in the
rate constants, so the residue formula's poles at a double root are
removable and the trajectory is analytic in the rate constants AT the
tangency. Measured at an exact double root (`k13` at 1e-300, `k31` on
the slow root of the reduced quadratic): `d/d(log ke)` is -9.2779361,
stable to eight significant digits at step sizes 1e-4 through 1e-7 and
equal to what the tape itself returns one part in 1e8 away. The `NaN` is
the implementation routing through the eigenvalues, not a derivative
that does not exist.

**How often it happens.** Not a knife edge. Over 4000 random draws per
box: 0 on a pharmacokinetic box of 1e-2 to 10 per hour, 0 on 1e-4 to 20,
38 on 1e-8 to 50, and 749 on 1e-20 to 50. By spread, 3000 draws per
cell: 0 at six decades or fewer, 0.07 percent at eight decades and 1.33
percent at ten, which is one draw in 75.

**And the `NaN` is the loud end of a soft region rather than a boundary
guarding a right answer.** On ten-decade draws where the gradient comes
back FINITE, so nothing fires at all, the value is wrong by up to
1.3e-04 of the trajectory's own maximum. A silent value error outranks a
loud `NaN`, and it is the reason this is stated as a domain above rather
than patched.

**Why it is not clamped.** Holding the
[`acos()`](https://rdrr.io/r/base/Trig.html) argument below `1 - 1e-14`
makes the gradient finite and right at a genuine double root and costs
nothing there (value error 5.177e-13 against the shipped 5.181e-13, and
no change at all over 400 ordinary draws). But on 40 draws from the wide
region it changes the VALUE by a factor of 2.8e+02 to 2.5e+06. It is the
right fix for a real collision and destructive where the roots were
already garbage, and telling the two apart needs a test on whether the
cubic's coefficients still carry relative precision, which is a
comparison on an automatic-differentiation value and is the constraint
this whole function is written around. So it is left, and said.

**What a user sees.** An optimizer cannot walk INTO the region: fitting
three compartments to two-compartment data from six starting points,
including `log(k13)` at -20, -30, -35 and -700, every fit returned
normally and `nlminb` left `log(k13)` exactly where it began, because
the gradient in that direction is zero once `k13` stops mattering. What
a user can do is START inside it, and then the fit stops with `nlminb`'s
`NA/NaN gradient evaluation`, which is loud and names neither this
function nor the cause.

A three-compartment model whose five disposition rate constants have all
underflowed to exactly zero returns `NaN`, in the value and in the
gradient, because the cubic is then a triple root at zero and its
trigonometric solution reads `0 / 0`. It fails loudly rather than
quietly, and reaching it needs five log rates below about -745 at once.
It is not offset, because offsetting it by this file's own 1e-150
returns a finite answer at eight-ninths of the true height, a silent
factor of 1.125.

## See also

[`frm_ode()`](https://aforren1.github.io/frmtmb/frmtmb.ode/reference/frm_ode.md)
for a system with no closed form, and for every feature named in "What
is refused, and why".

## Examples

``` r
# One-compartment oral pharmacokinetics, 100 into the depot every 12
# hours, and the same schedule already at steady state.
doses <- data.frame(time = 0, state = "depot", value = 100,
                    ii = 12, addl = 3L)
frm_lincmt(parms = list(ka = 1, ke = 0.2, V = 10),
           times = c(6, 18, 30, 42), ncmt = 1, depot = TRUE,
           events = doses)
#> [1] 3.733943 4.075490 4.106474 4.109285

frm_lincmt(parms = list(ka = 1, ke = 0.2, V = 10),
           times = c(6, 18, 30, 42), ncmt = 1, depot = TRUE,
           events = data.frame(time = 0, state = "depot", value = 100,
                               ii = 12, ss = TRUE))
#> [1] 4.109565265 0.375622018 0.034075678 0.003091276

# The flip-flop case, where the textbook form returns NaN
frm_lincmt(parms = list(ka = 0.2, ke = 0.2, V = 10), times = 1:4,
           ncmt = 1, depot = TRUE, init = list(depot = 100))
#> [1] 1.637462 2.681280 3.292870 3.594632

# In a formula, with between-subject variability on both rates
set.seed(2026)
tt <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
dd <- data.frame(id = factor(rep(1:6, each = length(tt))),
                 time = rep(tt, 6), dose = 100)
ka <- exp(rnorm(6, 0, 0.3))[as.integer(dd$id)]
ke <- exp(rnorm(6, log(0.2), 0.25))[as.integer(dd$id)]
dd$conc <- 100 * ka / (10 * (ka - ke)) *
  (exp(-ke * dd$time) - exp(-ka * dd$time)) + rnorm(nrow(dd), 0, 0.3)
# \donttest{
fit <- frm(
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, init = list(depot = dose)),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
    gaussian(),
  data = dd, start = list(beta = c(0, log(0.25), log(8))))
fixef(fit)
#> $lka
#> (Intercept) 
#>  -0.2014589 
#> 
#> $lke
#> (Intercept) 
#>   -1.743725 
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
#>   -1.293759 
#> 
# }
```
