# Differential equation models

Some models say what the *rate of change* of a quantity is, not what the
quantity is. A drug moves from a depot into the blood and is cleared
from it; a population grows and is eaten.
[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
puts such a system inside a nonlinear formula, so the constants of the
system are ordinary nonlinear parameters with fixed effects, random
effects and covariates, and the model is fitted by maximum likelihood
like any other
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) model.

## Setup

[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
needs the **RTMBode** package, which supplies a differentiable interface
to the **deSolve** integrators. It is not on CRAN:

``` r

install.packages("RTMBode", repos = c(
  "https://kaskr.r-universe.dev",
  "https://cloud.r-project.org"))
```

Everything below is skipped if it is not installed.

## A population pharmacokinetic model

The standard one-compartment model with first-order oral absorption has
two states: the amount `A` still in the depot, and the concentration `C`
in the central compartment.

``` math
\frac{dA}{dt} = -k_a A, \qquad
  \frac{dC}{dt} = \frac{k_a A}{V} - k_e C
```

A single dose `D` enters the depot at time zero, so $`A(0) = D`$ and
$`C(0) = 0`$. Three constants describe a subject: the absorption rate
$`k_a`$, the elimination rate $`k_e`$, and the apparent volume of
distribution $`V`$. All three must be positive, so they are estimated on
the log scale.

The system itself is a plain R function in the **deSolve** convention:
it takes the time, the state vector and the parameter vector, and it
returns the derivatives in a list. Index `y` and `p` by position.

``` r

pk_dyn <- function(t, y, p) {
  # base c() drops the automatic-differentiation class; this restores it
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1],
         p[1] * y[1] / p[3] - p[2] * y[2]))
}
```

[`datasets::Theoph`](https://rdrr.io/r/datasets/Theoph.html) is the
classic example: 12 subjects, one oral dose of theophylline each, 11
serum concentrations over 25 hours.

``` r

d <- datasets::Theoph
d$Subject <- factor(as.character(d$Subject))
head(d, 3)
#>   Subject   Wt Dose Time conc
#> 1       1 79.6 4.02 0.00 0.74
#> 2       1 79.6 4.02 0.25 2.84
#> 3       1 79.6 4.02 0.57 6.57
```

The model is one `bf(nl = TRUE)` formula. The body says how a
concentration is produced; the three nonlinear parameter formulas say
how each constant varies between subjects.

``` r

form <- bf(
  conc ~ frm_ode(pk_dyn,
                 init   = list(Dose, 0),
                 times  = Time,
                 parms  = list(exp(lka), exp(lke), exp(lV)),
                 group  = Subject,
                 states = c("depot", "central"),
                 output = "central"),
  lka ~ 1 + (1 | Subject),
  lke ~ 1 + (1 | Subject),
  lV  ~ 1,
  nl = TRUE
)

fit <- frm(form + gaussian(), data = d, se = TRUE,
           start = list(beta = c(0.5, log(0.08), log(0.5))))
summary(fit)
#> Family: gaussian 
#> Formula: conc ~ frm_ode(pk_dyn, init = list(Dose, 0), times = Time, parms = list(exp(lka), exp(lke), exp(lV)), group = Subject, states = c("depot", "central"), output = "central") 
#> Method: ML   nobs: 132 
#> Groups: Subject, 12 
#> logLik: -191.192  AIC: 394.384  BIC: 411.68 
#> 
#> Random effects:
#>   lka: 1 | Subject 
#>         Name Std.Dev.
#>  (Intercept)  0.68534
#>   lke: 1 | Subject 
#>         Name Std.Dev.
#>  (Intercept)  0.35835
#> 
#> Coefficients (lka):
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  0.37847    0.20822  1.8177  0.06911
#> 
#> Coefficients (lke):
#>             Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -2.37288    0.11884 -19.968 < 2.2e-16
#> 
#> Coefficients (lV):
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.821161   0.026642 -30.822 < 2.2e-16
#> 
#> Coefficients (mu):
#>      Estimate Std. Error z value Pr(>|z|)
#> 
#> Coefficients (sigma):
#>              Estimate Std. Error z value  Pr(>|z|)
#> (Intercept) -0.229063   0.068604 -3.3389 0.0008411
```

The typical subject’s constants come back on the natural scale, and the
derived quantities a pharmacokineticist reads off them follow:

``` r

ka <- exp(fixef(fit)$lka); ke <- exp(fixef(fit)$lke)
V <- exp(fixef(fit)$lV)
c(ka = ka, ke = ke, V = V,
  half_life = log(2) / ke,       # hours
  clearance = ke * V)            # L/h per kg of body weight
#>        ka.(Intercept)        ke.(Intercept)         V.(Intercept) 
#>            1.46004473            0.09321157            0.43992080 
#> half_life.(Intercept) clearance.(Intercept) 
#>            7.43627852            0.04100571
```

Everything downstream works as usual:
[`predict()`](https://rdrr.io/r/stats/predict.html) on new times,
[`ranef()`](https://aforren1.github.io/frmtmb/reference/ranef.md) for
the subject deviations,
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html), `REML = TRUE`.

``` r

nd <- data.frame(Subject = factor("1", levels = levels(d$Subject)),
                 Time = seq(0, 25, length.out = 5), Dose = 4.02)
predict(fit, newdata = nd)
#> [1] 0.000000 7.257109 5.663743 4.420213 3.449711
```

## How the arguments work

[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
runs **one solve per group**, over that group’s own observation times,
then puts the results back in the row order of the data.

| Argument | What it is |
|----|----|
| first argument | The derivative function, `function(t, y, parms)`. |
| `init` | Initial states, one list element per state. |
| `times` | The observation time of every row. |
| `parms` | The constants of the system, one list element per constant, in the order the derivative function indexes them. |
| `group` | The unit that owns one system: a subject, a reactor, a batch. |
| `states`, `output` | State names, and which of them the body returns. |
| `t0` | The initial time. Defaults to 0. |
| `method`, `atol`, `rtol` | Integrator and tolerances. |

`init` and `parms` are given as lists of **columns**, not as one row of
values per state. Each element is either one value per observation or a
single value shared by every group. `init = list(Dose, 0)` is therefore
“the depot starts at this row’s `Dose`, the central compartment starts
at zero”. A bare vector such as `c(100, 0)` is refused, because it
cannot be told apart from one column of two observations.

With no `output` the result is a matrix with one column per state, so
`frm_ode(...)[, 2]` also works. Naming the states and selecting one by
name is clearer and is what the example above does.

## Three things that will bite

### A covariate must not vary inside a group

[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
reads each constant off the group’s **first row**, because one solve has
one set of constants. A covariate that changes between an early and a
late observation of the same subject therefore describes a model this
helper cannot solve, and it is refused by name at frame assembly:

``` r

d$phase <- factor(ifelse(d$Time <= 2, "early", "late"))
frm(bf(conc ~ frm_ode(pk_dyn, init = list(Dose, 0), times = Time,
                      parms = list(exp(lka), exp(lke), exp(lV)),
                      group = Subject, output = 2L),
       lka ~ 1 + (1 | Subject),
       lke ~ 1 + phase + (1 | Subject),      # varies within Subject
       lV ~ 1, nl = TRUE) + gaussian(), data = d)
#> Error: Nonlinear parameter 'lke' is a dynamics input of frm_ode() but
#> is not constant within 'Subject': phaselate. ...
```

A covariate that is constant within a subject, such as weight or
treatment arm, is the intended case: `lV ~ 1 + Wt` is fine.

Without the check the fit is quietly wrong rather than loud. The
coefficient stays at its starting value, the log-likelihood does not
move, and the only symptom is an indefinite Hessian.

### Only adaptive integrators

A fixed-step integrator run at the default number of steps does not
solve the system to the tolerance the likelihood is defined at, so it
returns a *different* likelihood, not a noisier one: on the probe data
`rk4` reported -63.93 where `lsoda`, `lsode`, `adams` and `ode45` all
reported -60.46 to seven digits.
[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
refuses `method = "rk4"` and `method = "euler"`. The default `"lsoda"`
is both the fastest and the quietest of the adaptive methods.

### Keep the system small

Each group’s system must be small. Above roughly eight states in one
system, the second-order derivative path that the Laplace approximation
needs starts returning `NaN` gradients, and with some integrators the
process crashes.
[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
warns past that point. This is why groups are solved one at a time and
never stacked into one large system, even though stacking would look
like the obvious optimization.

Compartment models live far below the ceiling, so in practice this only
means: do not try to fold the subjects into the state vector yourself.

## Cost

The tape is built once, but every gradient evaluation replays one
adjoint solve per group, so the cost is linear in the number of groups
and does not shrink with tape reuse. For the two-state model above,
eight timepoints per subject, a whole
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.md) call
including `sdreport()`:

| subjects | rows | with [`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md) | with the closed form |
|----|----|----|----|
| 6 | 48 | 1.6 s | 0.05 s |
| 12 | 96 | 3.5 s | 0.06 s |
| 25 | 200 | 8.3 s | 0.08 s |
| 50 | 400 | 19.9 s | 0.11 s |

About 0.4 s per subject. A 200-subject population model is a couple of
minutes. If your system has a closed-form solution, write that instead:
it gives the same answer roughly a hundred times faster.

## What is out of scope

Dosing **event tables** are not supported: repeated doses, infusions,
`evid`/`amt` records.
[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
covers models driven by their initial conditions, which is the
single-dose case. For event-driven population pharmacokinetics use a
dedicated tool such as **nlmixr2**.

`predict(se.fit = TRUE)` is not available for any nonlinear predictor,
[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
included. Ask for a nonlinear parameter instead, with
`predict(dpar = "lka")`.

## Failed solves

Mid-optimization the optimizer probes extreme rate constants, and a
solve can fail there. `on_error` and `penalty` say what to do about
that, but they reach less far than they look like they do, and the
difference is worth knowing before you rely on them.

**While fitting, most failures cannot be caught.** The body runs on the
automatic-differentiation tape. There
[`RTMBode::ode()`](https://rdrr.io/pkg/RTMBode/man/ode.html) returns an
AD object and raises no R error when a trajectory goes bad, and RTMB
refuses comparison on AD types, so nothing can test the result for
`NaN`. A diverging region of the parameter space therefore reaches you
as the optimizer’s own warning,

    Warning: NA/NaN gradient evaluation

and `on_error = "error"` will not name the group. Two failures are still
caught on the tape: an integrator that gives up and returns fewer time
points than it was asked for, and an R error raised by your derivative
function.

**Everywhere else every check applies.**
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), a direct call
to
[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md),
and a fit whose `init` and `parms` contain no estimated parameter are
ordinary numeric arithmetic. There a failed solve is detected, and by
default its rows are filled with `penalty`.

A penalty is never written silently.
[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
warns and names the groups, and
[`frm_ode_failures()`](https://aforren1.github.io/frmtmb/reference/frm_ode_failures.md)
reads the record back afterwards:

``` r

p <- predict(fit, newdata = nd)
#> Warning: frm_ode(): the solve failed for 1 of 3 groups (2). Their rows
#> hold penalty = 1e+06, not a solution. ...

frm_ode_failures()
#> $groups
#> [1] "2"
```

Those rows are missing values dressed as numbers. Treat them as missing.

To localize an `NA/NaN gradient` from a fit, call
[`frm_ode()`](https://aforren1.github.io/frmtmb/reference/frm_ode.md)
yourself at the suspect parameter values with `on_error = "error"`.
Numerically it will name the group that cannot be solved.

Solver warnings from **deSolve** (“corrector convergence failed
repeatedly”, “exceeded maxsteps”) during a fit are normal noise from
those probing steps. Judge the fit by
[`diagnose()`](https://aforren1.github.io/frmtmb/reference/diagnose.md)
and the gradient at the optimum, not by whether the solver complained on
the way.
