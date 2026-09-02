# ODE models via RTMBode: feasibility probe

Date: 2026-09-01. Branch `wt-odeprobe`, frmtmb 0.29.0, RTMBode 1.0,
deSolve 1.42, R 4.6.1, Windows 11.

> **Status: implemented.** `frm_ode()` (`R/ode.R`) is the "one exported
> helper" of section 8, with all seven guards. `tests/testthat/test-ode.R`
> reproduces the probe A2 numbers exactly (logLik -60.46293086 against
> -60.462931, AIC 132.9259 against 132.926), and `vignettes/ode.Rmd` is
> the tutorial. Two departures from the sketch and one addition:
>
> - `init` and `parms` are lists of **columns**, not `cbind()`. `cbind()`
>   on advectors does work, but a bare `c(100, 0)` cannot be told apart
>   from one column of two observations, so lists remove the ambiguity
>   and read better next to the state names.
> - `output = "central"` (with `states =`) replaces `frm_ode(...)[, 2]`.
>   Character indexing of an advector matrix fails - RTMB drops
>   `colnames` - so the name is resolved to a position by the helper.
>   The matrix form still works when `output` is omitted.
> - Guard 4 (within-group constancy) had to move. It cannot be done on
>   the tape for a linear predictor: RTMB refuses comparison on AD types,
>   and a value-based check would miss the probe G trap anyway, because
>   the offending coefficient starts at exactly 0. It is done instead at
>   frame assembly (`check_ode_constancy()`), against the design matrices
>   and the grouping column, which is exact and start-value independent.
>   Plain data columns are still checked inside `frm_ode()`.
>
> Guard 7 (`tryCatch` a failed solve into a penalty) is real but reaches
> much less far than section 6 assumed, and the help page and the
> vignette now say so. On the tape `RTMBode::ode()` returns an advector
> and raises nothing when a trajectory diverges, and RTMB refuses
> comparison on AD types, so a non-finite solution cannot be tested for:
> a bad region of the parameter space still surfaces only as the
> optimizer's `NA/NaN gradient evaluation`, and `on_error = "error"`
> cannot name the group. What the guard does catch on the tape is a
> short solution matrix and an R error raised by the dynamics function.
> Off the tape - `predict()`, `simulate()`, `residuals()`, and a body
> whose `init` and `parms` hold no estimated parameter, which is
> evaluated numerically as the tape is built - every check applies, and
> there the penalty warns and names the groups; `frm_ode_failures()`
> reads the record back. A third numeric-path check had to be added:
> deSolve reports "integration was not successful" as a *warning* and
> still returns a full-length matrix of finite garbage (1e85 in the
> test), so neither the row count nor `is.finite()` sees that one.
>
> One hook was needed in `R/frame.R`: `all.vars()` on the nl body puts
> the *dynamics function* in `datavars`, so `model.frame()` was asked for
> a column named `pk_dyn`. `drop_nl_lexical_datavars()` leaves a body
> name that resolves to a function and is not a column of `data` to
> resolve lexically, which is what the body's own environment does
> anyway. A column always wins over a same-named function. The cost is
> that a misspelled column sharing a name with a base function (`t`,
> `c`, `df`) is no longer caught by `model.frame()`, so the objective
> re-raises the late coercion error with those names attached
> (`nl_body_error()`).

**No package code was changed.** Every result below comes from user
code in the `nl = TRUE` body or in a `custom_family()`, both of which
are already arbitrary R evaluated on the AD tape. Probe scripts are in
`dev/ode/`. They run in order from the worktree root with
`Rscript dev/ode/<script>.R`; `probeC` reads the fit `probeA2` saves.
The probes load frmtmb 0.29.0 from a scratch library, so set
`.libPaths()` to wherever the worktree build is installed.

## Verdict

Population ODE models already work in frmtmb today, on the nl route,
with correct Laplace likelihoods, correct standard errors, and a
working `predict()`/`simulate()`/`confint()` surface. The one blocking
condition is architectural, not fixable by us: **solve one small system
per group, never one big stacked system** (see "The stacked-solve
ceiling"). The remaining gap is ergonomics, not capability.

Recommendation: ship **one exported helper** that hides the per-group
solve loop and its RTMB gotchas, plus a vignette. Do **not** build an
`odefun =` grammar on `bf()`. Sizing in "Recommended design".

## 1. Backend works standalone

`dev/ode/probe0-lotka.R` reproduces RTMBode's own `?ode`
Lotka-Volterra example verbatim: `MakeADFun` over a likelihood calling
`ode()`, `nlminb` (0.17 s, convergence 0, objective 557.2158),
`sdreport` returning finite SEs for all 8 parameters (5 dynamics
parameters, 2 initial states, 1 observation SD). The backend
differentiates through both `parms` and `y` as documented.

## 2. Probe A: the nl route (the target use case)

Model: one-compartment first-order oral absorption, single bolus dose
D into the depot,

    dA/dt = -ka A                 A(0) = D
    dC/dt =  ka A / V - ke C      C(0) = 0

with between-subject variability on log(ka) and log(ke). 12 subjects x
8 timepoints, additive gaussian error. The closed form
`C(t) = D ka / (V (ka - ke)) (exp(-ke t) - exp(-ka t))` is both the
simulation truth and the reference the solver is checked against, so
an ODE failure is distinguishable from a specification failure
(`dev/ode/pk-common.R`).

### What the nl body actually sees

Read from `R/parse.R:563-611` and `R/objective.R:69-74`:

- `parse_one_response()` stores the RHS expression as
  `dpars$mu$nl_body`, its non-nlpar variables as `datavars`, and the
  formula environment as `nl_env`.
- `build_frame()` (`R/frame.R:617-637`) pulls each `datavar` out of
  the combined model frame into `data_list`, verbatim - a factor stays
  a factor, a numeric stays numeric.
- The objective evaluates `eval(lp$nl_body, c(dparv[[resp]], lp$data_list),
  lp$nl_env)`. So the body's variables are exactly: one ROW-WISE
  advector per nlpar (length `n_obs`, already carrying the fixed
  effects, random effects and offsets), plus the raw data columns.
  Everything else - including `RTMBode::ode` and any user helper -
  resolves lexically through `nl_env`, which is the formula's
  environment (the global environment in a script).

Confirmed empirically by `dev/ode/probeA2-frmtmb-nl.R`:
`data_list` names `time, id, dose` with classes `numeric, factor,
numeric`, and `nl_env` identical to `globalenv()`.

### The spelling that works

A per-subject solve is expressible because a nlpar that is constant
within a group can be read off that group's first row:

```r
pk_ode <- function(ka, ke, V, time, id, dose) {
  "[<-" <- RTMB::ADoverload("[<-")
  "c"   <- RTMB::ADoverload("c")
  ids <- as.integer(id)
  out <- numeric(length(ids))
  for (s in unique(ids)) {
    idx <- which(ids == s)
    idx <- idx[order(time[idx])]
    sol <- RTMBode::ode(
      y     = c(dose[idx[1]], 0),
      times = c(0, time[idx]),          # deSolve wants times[1] == t0
      func  = pk_dyn,
      parms = c(ka[idx[1]], ke[idx[1]], V[idx[1]]),
      method = "lsoda", atol = 1e-8, rtol = 1e-8)
    out[idx] <- sol[-1, 3]
  }
  out
}

frm(bf(conc ~ pk_ode(exp(lka), exp(lke), exp(lV), time, id, dose),
       lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
      gaussian(),
    data = d, start = list(beta = c(0, log(0.25), log(8))), se = TRUE)
```

This fits with no errors and no warnings. Verbose breakdown, 96 rows /
12 subjects: parse 0.01 s, frame 0.00 s, tape 0.22 s (6 outer, 24
inner parameters), optimize 3.30 s, sdreport 0.93 s, total 4.50 s,
`max|grad| 0.00026`.

### Correctness

Four independent anchors agree:

1. Numeric mode: `max |ode - closed form| = 2.1e-08`.
2. The hand-rolled `MakeADFun` reference (`dev/ode/probeA1-reference.R`)
   and the frmtmb fit reach the **same objective, 60.462931**, and the
   same estimates and SEs to 6 digits.
3. Replacing the ODE solve with the closed form in the same nl body
   gives **identical logLik at every size tested** (n_id 6/12/25/50:
   -21.490 / -51.056 / -134.547 / -257.099 both ways). This is the
   decisive check that the Laplace approximation is correct *through*
   the adjoint ODE node, not just the mode.
4. nlmixr2 FOCEi on the same data (`dev/ode/probeF-nlmixr2.R`):

   | quantity | truth (pop) | truth (sample) | frmtmb | nlmixr2 FOCEi |
   |---|---|---|---|---|
   | lka      | 0.000  | -0.1735 | -0.2292 | -0.2269 |
   | lke      | -1.609 | -1.6375 | -1.5932 | -1.5966 |
   | lV       | 2.303  | 2.3026  | 2.2704  | 2.2726  |
   | sd(lka)  | 0.300  | 0.2349  | 0.2533  | 0.2539  |
   | sd(lke)  | 0.250  | 0.2823  | 0.2624  | 0.2628  |
   | sigma    | 0.300  | 0.3066  | 0.3040  | 0.3041  |
   | logLik   |        |         | -60.4629 | -60.4804 |
   | AIC      |        |         | 132.926 | 132.961 |

   frmtmb finds a marginally *better* optimum; nlmixr2 itself reports
   "last objective function was not at minimum". The "sample truth"
   column (empirical mean/sd of the simulated subject deviations) is
   the honest anchor - with 12 subjects the Monte Carlo offset in
   `lka` is larger than the estimation error. `cor(ranef, simulated
   deviations)` is 0.981 (lka) and 0.987 (lke).

   Cost note: nlmixr2 needed 19.7 s of setup (it compiles C for the
   model, so it needs Rtools) against 0.6 s of optimization. The
   frmtmb route needs no compiler at all.

## 3. Probe B: the custom_family fallback

Not needed, but it works and is worth documenting as an alternative
(`dev/ode/probeB-custom-family.R`, `probeB-family-def.R`). A
`custom_family()` with `dpars = c("lka","lke","lV","lsigma")`,
`primary_dpars = "lka"`, and an `lpdf` that calls the same `pk_ode()`
with times/ids/doses carried as `vreal(time, dose) + vint(idn)`
reaches the **identical objective 60.462931** and identical estimates.

Differences from route A:

- Slower: optimize 5.36 s vs 3.30 s, sdreport 1.64 s vs 0.93 s (four
  real linear predictors instead of three plus a free-form body).
- `predict()` defaults to `type = "link"`, which on this family is the
  *lka* linear predictor, not the concentration. The response mean
  needs `predict(type = "response")` or `fitted()`; both then match
  route A exactly, including on `newdata`. This is documented
  behavior, not a bug, but it is a sharp edge for a family whose
  `mean_fn` is the whole model.
- It needs `mean_fn`, `var_fn` and `sim` written out by hand.

Route A is better for a mu-shaped ODE model. Route B is the right
answer only when the ODE enters something other than the mean.

## 4. Probe C: downstream survey (route A)

`dev/ode/probeC-downstream.R`. Everything works except one
pre-existing nl limitation:

| feature | result |
|---|---|
| `summary()`, `vcov()`, `logLik()`, `AIC()` | ok |
| `confint(method = "wald")` | ok |
| `confint(method = "profile")` | ok |
| `predict()` in-sample, `fitted()`, `residuals()` | ok |
| `predict(dpar = "lka")` | ok (12 distinct values) |
| `predict(newdata =)`, dense grid, existing subjects | **ok** |
| `predict(newdata =)`, reproduces in-sample on training rows | ok |
| `predict(newdata =)`, new subject + `allow_new_levels` | ok |
| `predict(newdata =)`, new subject + `re.form = NA` | ok |
| `predict(newdata =)`, one row per subject | ok |
| `simulate(nsim = 2)` and refit on the simulated response | ok |
| `ranef()`, `VarCorr()`, `coef()` | ok |
| `REML = TRUE` | ok |
| `predict(se.fit = TRUE)` | **fails** (see below) |

`newdata` was expected to fail and does not. The reason is that
`predict.frmtmb_fit()`'s nl branch (`R/predict.R:690-715`) predicts
each nlpar on `newdata` and then re-evaluates the same body against
`newdata`'s own columns, so the helper re-derives its groups from the
`newdata` id column and solves that grid. A `newdata` slice with a
different row set, different times, or a different number of rows per
subject all work, and a subset of subjects reproduces the in-sample
values to 1e-6.

`predict(se.fit = TRUE)` errors with the existing message "se.fit is
not supported for the nonlinear predictor yet; request the nonlinear
parameters (dpar = ...) instead" (`R/predict.R:688`). This is a
general nl gap, not ODE-specific, and the delta method through an
adjoint ODE node is the awkward part of closing it.

## 5. Probe D: edge cases

`dev/ode/probeD-edges-and-scaling.R`, all on route A.

| case | result |
|---|---|
| an observation at `t = 0` (grid becomes `c(0, 0, ...)`) | ok - deSolve tolerates the duplicate and returns `C(0) = 0` |
| duplicate times within a subject (replicate assays) | ok |
| ragged design, unsorted times, shuffled rows | ok (the helper sorts and scatters) |
| `NA` responses (`na.omit` drops rows first) | ok, `nobs` 93, prediction length 93 |
| `REML = TRUE` | ok |
| a covariate that varies WITHIN a solve group | **wrong model, but detectable** |

The last row is the real trap and belongs in any documentation. The
helper reads each group's parameters off its first row, so a
within-subject covariate on a dynamics parameter cannot affect the
likelihood. Adding `lke ~ 1 + phase + (1 | id)` with a time-varying
`phase` leaves the coefficient at its starting value (exactly 0.000),
the logLik unchanged to 6 digits, and produces "Hessian is not
positive definite; standard errors are unreliable" with `NaN` SEs. So
it is not silent - but the signal is indirect and a user could miss
it. A shipped helper should assert within-group constancy and error
by name.

## 6. Cost

Per-subject loop, full `frm()` call including sdreport, 8 timepoints
per subject:

| subjects | rows | ODE route | closed form | ratio |
|---|---|---|---|---|
| 6  | 48  | 1.60 s  | 0.05 s | 32x  |
| 12 | 96  | 3.47 s  | 0.06 s | 58x  |
| 25 | 200 | 8.26 s  | 0.08 s | 103x |
| 50 | 400 | 19.91 s | 0.11 s | 181x |

Linear in the number of subjects, about 0.4 s per subject for this
2-state system. The tape build itself is cheap (0.22 s at 12
subjects); the cost is replaying N adjoint ODE nodes on every gradient
evaluation. A 200-subject population PK fit is therefore a couple of
minutes, not seconds - usable, and roughly the ballpark nlmixr2 lands
in once its compile time is counted, but it is not free.

### Integrator choice matters

Same 2-state per-subject model, regular and ragged designs
(`dev/ode/probeG-integrator-and-trap.R`):

| method | regular | ragged | logLik |
|---|---|---|---|
| lsoda  | 3.53 s  | 3.22 s  | -60.462931 / -42.812831 |
| lsode  | 12.53 s | 21.36 s | -60.462928 / -42.812831 |
| adams  | 6.63 s  | 7.18 s  | -60.462930 / -42.812832 |
| ode45  | 8.07 s  | 8.29 s  | -60.462931 / -42.812831 |
| rk4    | 3.16 s  | 3.36 s  | **-63.930599 / -93.855616** |

Every adaptive integrator agrees to 7 digits. The fixed-step ones
(`rk4`, `euler`) give a *wrong likelihood* at the default step count
and must never be the default. `lsoda` is both fastest and quietest.

Solver warnings ("corrector convergence failed repeatedly",
"exceeded maxsteps") do appear mid-optimization when the optimizer
probes extreme rate constants, and one such step produced a single
`NA/NaN function evaluation` warning from nlminb. Every fit still
converged to the same optimum, but the noise is real and a shipped
helper should consider `tryCatch`ing a failed solve into a large
finite penalty rather than letting `NaN` reach the optimizer.

## 7. The stacked-solve ceiling (the one hard constraint)

The obvious optimization - stack all subjects into one system of
`2 x n_subject` states and make **one** `ode()` call, collapsing N
adjoint nodes into one - **does not work under Laplace**, and fails in
ways that range from silent to fatal.

`dev/ode/probeE5-random.R`: the stacked PK gives the identical
objective and gradient to the loop when `random = NULL`, and `NaN` for
both when `random = c("u_ka","u_ke")`. So the failure is specifically
in the second-order path the Laplace inner problem needs.

`dev/ode/probeE8-laplace-limit.R` is a minimal, frmtmb-free
reproduction: `n` independent decays `dy_i/dt = -exp(mu + u_i) y_i`
with `u ~ N(0,1)` integrated out. First-order AD is fine at every size
tested; the Laplace gradient degrades with the state count, and the
threshold depends on the integrator:

| method | 8 states | 16 | 24 | 32 | 48 | 64 |
|---|---|---|---|---|---|---|
| lsoda  | ok  | NaN | NaN | NaN | - | - |
| lsode  | NaN | NaN | NaN | NaN | - | - |
| adams  | ok  | ok  | ok  | ok  | ok | hangs (>10 min) |
| ode45  | ok  | ok  | ok  | **crash (exit 127)** | - | - |
| rk4    | ok  | ok  | ok  | **crash (exit 127)** | - | - |
| euler  | ok  | ok  | ok  | **crash (exit 127)** | - | - |

Stacked PK under Laplace breaks at 4 subjects (8 states) with `lsoda`,
which is why the whole stacked variant died.

Two consequences:

- **Design constraint**: solve per group. Small systems (a handful of
  states) are exactly the regime the adjoint second-order path
  handles, and pharmacometric compartment models live there. Any
  future sugar must not "optimize" by stacking.
- **Upstream report**: this is worth filing against RTMBode.
  `probeE8-laplace-limit.R` is a 25-line reproduction. Silent `NaN`
  gradients (lsoda/lsode) and a process crash (ode45/rk4/euler) are
  both bad failure modes for a package a user might reach for with a
  20-state epidemic model.

Related session hazard found while bisecting: **a `NaN` ODE tape
poisons later `MakeADFun` objects in the same R session.** After the
stacked failure, an unrelated 1-state `ode()` model in the same
session threw `'*this' must be 'advector' (lost class attribute?)`
from `$gr()`; the identical code in a fresh session returns correct
`fn`, `gr` and `he` (`probeE5` vs `probeE6` vs a fresh run). Probe
scripts that sweep sizes must fork one process per size, and any test
we ship must not run a deliberately-failing ODE model before a real
one.

## 8. Recommended design

### What to build

**One exported helper, no parser changes.** The nl body already
carries the entire model; the only thing frmtmb can usefully add is
the loop and the gotcha-proofing. Something like:

```r
frm_ode(dynamics, init, times, parms, group,
        method = "lsoda", atol = 1e-8, rtol = 1e-8, on_error = "penalize")
```

called from the nl body:

```r
bf(conc ~ frm_ode(pk_dyn, init = cbind(dose, 0),
                  times = time, parms = cbind(exp(lka), exp(lke), exp(lV)),
                  group = id)[, 2],
   lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
```

It must own, because every one of these is a way for a user to get a
wrong answer quietly:

1. the `ADoverload("[<-")` and `ADoverload("c")` locals (gotchas 3 and
   8 - lexically scoped, so they must be inside the helper);
2. prepending `t0` to each group's time vector and dropping the
   corresponding solution row;
3. sorting within group and scattering back to row order;
4. **asserting that every `parms` and `init` column is constant within
   `group`**, erroring by name instead of leaving a non-PD Hessian;
5. refusing fixed-step integrators, or at least warning that they
   change the likelihood;
6. a state-count guard that warns when a single system is large enough
   to hit the Laplace ceiling of section 7;
7. optional `tryCatch` around each solve so a failed step becomes a
   finite penalty rather than `NaN` in the optimizer.

### What NOT to build

An `odefun =` argument or an `ode()` term on `bf()` with its own
grammar (state names, dosing records, `evid`/`amt` event tables). It
would require new machinery in `parse.R`, `frame.R`, `objective.R`,
`predict.R` and `simulate-new.R` for capability the nl body already
has, and it would end up a worse nlmixr2. The union-coverage argument
(brms has no ODE grammar) is satisfied by the helper plus a vignette.

The one thing the helper does not reach is **dosing events** -
multiple doses, infusions, `evid`/`amt` records. That needs deSolve's
`events =` plumbed through `RTMBode::ode()` and was not probed. It is
the honest boundary of this feature: single-dose and
initial-condition-driven models now, event tables never (point users
at nlmixr2).

### Packaging

As already planned in `feature-gaps.md`, and unchanged by the probe:

- `Suggests: RTMBode`, `Additional_repositories: https://kaskr.r-universe.dev`
  (the brms/bayesplot/bridgesampling posture for cmdstanr).
- `deSolve` also to `Suggests` (RTMBode depends on it, but the helper
  should fail informatively if either is missing).
- `requireNamespace("RTMBode", quietly = TRUE)` guard in the helper
  with an install message naming the r-universe repo.
- `skip_if_not_installed("RTMBode")` on every test; the vignette
  chunks `eval = requireNamespace(...)`.
- Conditional compat-registry row.
- CRAN's check farm will not install it, so examples must degrade to
  nothing.

### Effort

| piece | size |
|---|---|
| `frm_ode()` helper + the seven guards above | 1-1.5 days |
| tests (skip-guarded; per-process ODE tests, see the session hazard) | 0.5 day |
| packaging (Suggests, Additional_repositories, guards, registry) | 0.5 day |
| population-PK vignette (diataxis: one tutorial, one how-to) | 0.5-1 day |
| upstream RTMBode issue with `probeE8` reproduction | 1 hour |
| **total** | **3-4 days** |

For comparison, a dedicated `bf()` ODE grammar with event tables is
2+ weeks and is not recommended.

## Probe scripts

| file | what |
|---|---|
| `dev/ode/pk-common.R` | simulation, closed form, dynamics, per-subject solver |
| `dev/ode/probe0-lotka.R` | RTMBode `?ode` example, standalone |
| `dev/ode/probeA1-reference.R` | hand-rolled RTMB population PK with ODE |
| `dev/ode/probeA2-frmtmb-nl.R` | the frmtmb nl spelling |
| `dev/ode/probeB-custom-family.R`, `probeB-family-def.R` | custom_family route |
| `dev/ode/probeB2-predict-check.R` | `predict()` link-vs-response on route B |
| `dev/ode/probeC-downstream.R` | post-processing survey, truth recovery |
| `dev/ode/probeD-edges-and-scaling.R` | edge cases, cost scaling, stacked attempt |
| `dev/ode/probeE2-state-limit.R` | first-order AD state sweep (all fine) |
| `dev/ode/probeE3-bisect.R`, `probeE4-bisect2.R` | ruling out structural causes |
| `dev/ode/probeE5-random.R` | Laplace vs no-Laplace, stacked vs loop |
| `dev/ode/probeE6-drop.R` | ruling out solution-matrix subsetting |
| `dev/ode/probeE7-stacked-laplace.R` | stacked Laplace by subject count (one per process) |
| `dev/ode/probeE8-laplace-limit.R` | minimal upstream reproduction |
| `dev/ode/probeE9-integrator-laplace.R` | the same, swept over integrators |
| `dev/ode/probeF-nlmixr2.R` | nlmixr2 FOCEi cross-check |
| `dev/ode/probeG-integrator-and-trap.R` | integrator sweep, within-group trap |
