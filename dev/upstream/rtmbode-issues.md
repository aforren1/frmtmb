# RTMBode defects: root causes, patches, verification

Date: 2026-09-02. Written for filing against `kaskr/RTMB`, subdirectory
`RTMBode/`. Everything here is reproduced without frmtmb, with bare
`RTMB` + `RTMBode` + `deSolve`.

Environment: R 4.6.1, Windows 11, gcc 14.3.0 (Rtools 4.5), RTMB 1.9,
RTMBode 1.0 (r-universe build of `5242257`, identical to the source in
this clone), deSolve 1.42, tmbstan 1.2.0, rstan 2.32.7.

The patches are in `patches/`, as a three-commit series on top of
`5242257`:

| patch | subject |
|---|---|
| `01-solver-failure-nan.patch` | a failed solve must give NaN, not an R error |
| `02-events-ad-path.patch` | name the ODE states, and refuse events that are not `add` |
| `03-augmented-workspace-guard.patch` | warn when the augmented system outgrows the integrator |

The reproduction scripts are in `repro/`. Each takes a library path as an
argument, so the same script runs against the stock build (`""`) and
against a patched one.

---

## 1. A failed solve escapes as an R error

### 1.1 What was reported

Sampling an RTMB objective that contains `RTMBode::ode()` nodes with
`tmbstan` aborts at warmup iteration 1, deterministically, even when the
chain starts at the fitted optimum where `obj$fn` and `obj$gr` both
succeed. Two signatures appear: DLSODA reports illegal input at
`TCUR = 0`, and RTMB's `ADjoint` raises `Wrong output length`. The abort
also corrupts rstan's nested autodiff arena for the rest of the R
session, so every later call fails with `empty_nested() must be true
before calling recover_memory()`.

Three causes were suspected: static pointer clobbering between two ODE
tapes in `src/odesolve.c`, tmbstan evaluating through a copied or
retaped TMB object, and garbage collection invalidating the raw
addresses the C statics hold.

### 1.2 Root cause

None of the three. The sampler is not involved at all.

`ODEadjoint()$updateSolution()` calls `deSolve::ode()` with no guard.
deSolve fails in two ways when the parameters are extreme, and both used
to escape into the caller:

1. It signals an R error, `illegal input detected before taking any
   integration steps`. lsoda raises this when the derivative is not
   finite at `t0`, which happens as soon as a rate parameter reaches
   `Inf` through `exp()`.
2. It returns early. The solution matrix then has fewer rows than
   `times`, so `sol[, solcols]` is shorter than the length `ADjoint`
   recorded at tape time, and RTMB reports `Wrong output length`.

Neither is a sampler interaction. Both are ordinary consequences of an
extreme parameter value, and Stan reaches those values by construction:
its initial step size is 1 on the unconstrained scale, and the gradient
of this objective at the mode is of order 1e2 to 1e5, so the first
leapfrog step of warmup lands many orders of magnitude away from the
mode.

### 1.3 The discriminating experiment

`repro/03-infparm.R` produces both signatures from a plain `obj$fn()`
call, with no sampler, no random effects and one ODE node:

| parameter | `exp(par)` | stock `obj$fn` |
|---|---|---|
| `logpar[2] = 50` | 5.2e+21 | finite |
| `logpar[2] = 100` | 2.7e+43 | non-finite |
| `logpar[2] = 500` | 1.4e+217 | `ERROR: illegal input detected before taking any integration steps` |
| `logpar[2] = 710` | `Inf` | same error |
| `mode + 100` (all coordinates) | | `ERROR: Wrong output length` |
| `mode + 1000` | | `ERROR: Wrong output length` |

`repro/04-hypotheses.R` shows the mechanism behind the second signature
directly. It traces `deSolve::ode` and reports the shape returned:

```
length(times) = 31
nrow(sol) returned by deSolve: 2
objective result: ERROR: Wrong output length
```

deSolve warned `Returning early. Results are accurate, as far as they
go` and handed back 2 rows where 31 were asked for.

### 1.4 The three suspected causes, eliminated

- **Static pointer clobbering between two tapes.** `repro/01-sampler.R`
  takes the number of ODE nodes as an argument. One node fails under
  tmbstan exactly as two do, so no interleaving is needed to trigger it.
  `repro/04-hypotheses.R` also runs 20 pairs of repeated two-tape
  gradients at perturbed points: every pair agrees exactly. The C
  statics are re-set by `setTape()` immediately before each solve, and
  the early return on the cache hit skips the solve as well, so the
  sequencing is in fact safe.
- **tmbstan evaluating through a different object.** Ruled out by 1.3:
  the direct path reproduces both signatures.
- **Garbage collection moving or collecting the tape.** Under
  `gctorture2(step = 1)` the direct gradient is bit-identical to the
  untortured one, and identical again after an explicit `gc()`.

### 1.5 The fix

Patch 01 wraps the `deSolve::ode()` call in `tryCatch`, detects a short
return, and substitutes a solution of the right shape filled with `NaN`.
The objective then evaluates to `NaN` instead of raising, which an
optimizer treats as a step to shorten and Stan treats as a proposal to
reject. Nothing crosses Stan's C++ boundary, so the arena stays
balanced.

Two details matter:

- The `NaN` fill uses `array()`, not `matrix()`. RTMB masks `matrix()`
  and its `num.` method returns an advector inside an active AD context,
  which reaches `EvalOp` as a complex vector and fails with
  `EvalOp: Function must return 'real' or 'integer'`.
- A solution with **more** rows than `times` is not a failure. deSolve
  inserts event times that are not already in `times`. Patch 01 maps
  that case back onto the requested times. Section 2.4 covers what this
  fixes.

### 1.6 Verification

`repro/01-sampler.R`, Lotka-Volterra with lognormal observations, one or
two `ode()` calls, `tmbstan(obj, chains = 1, iter = 30, init =
list(mode))`:

| ODE nodes | stock | patched |
|---|---|---|
| 1 | no draws, chain aborted | 15 draws |
| 2 | no draws, chain aborted | 15 draws |

`repro/05-sampler-long.R`, two ODE nodes, 400 iterations, patched build:

| init | outcome | divergences | seconds |
|---|---|---|---|
| at the mode | 200 draws | 0 | 23.6 |
| over-dispersed, `sd = 1.5` | 200 draws | 0 | 25.8 |

Both chains reach the same posterior mean to three digits, and the run
log shows deSolve `Returning early` warnings during warmup, so the
sampler is passing over failed solves rather than avoiding them.

`repro/03-infparm.R` on the patched build: every one of the nine
parameter settings above returns non-finite rather than raising.

Regression, `repro/12-regression.R`, the example in `?ode`, stock
against patched:

| quantity | stock | patched |
|---|---|---|
| optimum | 557.215804 | 557.215804 |
| `nlminb` iterations | 59 | 59 |
| estimates (`pars`, `yini`, `sdobs`) | identical | identical |
| numeric path vs `deSolve::ode` | max abs difference 0 | 0 |

---

## 2. Events on the autodiff path

### 2.1 The state vector has no names

`func2tape()` names its input `c("t", names(y), names(parms))`, but
`MakeTape(...)$par()` returns that vector unnamed. `addInfo()` therefore
stored an unnamed `info$state`, and `updateSolution()` passed an unnamed
vector to `deSolve::ode()`. deSolve reads state names off that vector
(`vars <- attr(y, "names")`), so `checkevents()` bounded the event index
by `length(vars)`, which was 0:

```
Error: too many state variables in 'event'; should be < 0
```

and, when the caller's `y` carried names:

```
Error: unknown state variable in 'event': depot
```

Only the order-0 solve was affected. The order-1 augmented state built
by `augment()` is named (`y1 y2 dy1 ...`), because it is assembled with
`c(y = ..., dy = ...)`.

Patch 02 carries `names(y)` through `addInfo()`, and keeps those names on
the leading block of the augmented state, so an event addressing a state
by name resolves at every differentiation order rather than only at
order 0.

### 2.2 With the names fixed, two of the three methods are silently wrong

RTMBode differentiates by integrating an augmented system that carries
the states together with their derivatives. A deSolve event row jumps the
state block and leaves the sensitivity block alone.

| method | state jump | correct sensitivity jump | what deSolve does |
|---|---|---|---|
| `add` | `y := y + a` | `dy := dy` | correct |
| `replace` | `y := v` | `dy := 0` | leaves `dy` alone |
| `multiply` | `y := f y` | `dy := f dy` | leaves `dy` alone |

Measured against central differences of the numeric solve, by emulating a
fixed event table with an event function that touches only the first
`nstate` entries. These four rows and the segmented-solve numbers below
come from the frmtmb probe `dev/ode/probeH2-events-adjoint.R`, which
predates this work; everything else in this document was measured here:

| method | AD `d/d(lka)` | finite differences | max relative error |
|---|---|---|---|
| `add` | 1.9931038 | 1.993096 | 3.9e-06 |
| `replace` | 0.85369576 | 1.460373 | **0.415** |
| `multiply` | -0.27782693 | -0.67921062 | **0.591** |

The solution values agree to ten digits in all three cases. Only the
gradient is wrong. An optimum found with those gradients is not the
maximum likelihood estimate, and nothing warns. `replace` is also
deSolve's fallback when no method is given, so the wrong case is the one
a user reaches by omission.

Patch 02 refuses `replace` and `multiply` by name, with an error that
names the working alternative. `events$func` is arbitrary R code and
cannot be inspected, so it warns rather than stops: it is correct where
it adds a constant to a state and wrong otherwise.

The workaround the message names is exact for every method: split the
integration at the event times, chain one `ode()` call per interval, and
apply the jump in ordinary RTMB arithmetic between the calls. RTMBode is
already differentiable through the initial state, so the adjoint is
correct by construction and the jump amount may itself be an estimated
quantity. Measured max relative gradient errors for that route, against
finite differences, over repeated bolus, bioavailability, `replace`,
`multiply` and infusion-rate cases: 9.5e-11 to 1.4e-10.

### 2.3 Verification

`repro/06-events.R`, one-compartment oral PK, three bolus doses,
`atol = rtol = 1e-10`:

| case | stock | patched |
|---|---|---|
| `add` | `unknown state variable in 'event': depot` | works, max relative error vs FD **6.97e-07** |
| `replace` | same error | refused by name |
| `multiply` | same error | refused by name |
| no `method` column | deSolve structure error | refused by name |
| `events$func` | accepted silently | accepted with a warning |
| no events at all | max relative error 4.97e-08 | 4.97e-08, unchanged |
| returned column names | `time depot central` | unchanged |

### 2.4 A third defect, found while fixing the first two

Event times that are not already in `times` make deSolve return **more**
rows than `times`, because it inserts them. The adjoint node has a fixed
output length, so the node fails. On the stock build this is masked by
2.1; with the names fixed and no row mapping it gives

```
Error: EvalOp: Function must return 'real' or 'integer'
```

The type in that message, rather than a length complaint, is the
`matrix()` masking noted in 1.5.

Patch 01 maps the returned solution back onto the requested times, taking
the first row where a time appears twice, which is the pre-event value
and matches the convention the plain numeric solve follows.

`repro/11-event-times.R`, `times = seq(0, 24, by = 2)`:

| dose times | `nrow` returned by deSolve | patched, max relative error vs FD |
|---|---|---|
| 6, 12, 18 (on the grid) | 13 | 6.04e-07 |
| 5, 11, 17 (off the grid) | 16 | 1.08e-05 |

---

## 3. The augmented-state ceiling

### 3.1 What was reported

Above about 8 states in one system, `lsoda` gives `NaN` gradients under
Laplace with no message, and some other integrators crash the process.

### 3.2 Root cause

The augmented system at differentiation order `k` has

```
neq = sum(nstate * (nstate + nparms)^(0:k))
```

states. The order needed is 1 for a gradient, **2 for a Laplace
objective** and **3 for the gradient of a Laplace objective**.
`repro/08-ceiling-mechanism.R` traces the sizes actually integrated and
confirms the split:

```
n = 10, method = lsoda
augmented sizes by order: 0:10  1:210  2:4210  3:84210
systems integrated during fn(): 10, 210, 4210
systems integrated during gr(): 84210
```

That is why `fn` stays correct while `gr` goes `NaN`: the objective needs
order 2, the gradient needs order 3, and only order 3 is too big.

The Livermore integrators that carry a full Jacobian size their real work
array as `lrw = 22 + neq * max(16, neq + 9)`, and deSolve computes that
in R **integer** arithmetic. `repro/09-lrw-overflow.R`:

```
integer.max = 2147483647
predicted first overflowing neq: 46337

neq=40000  lrw=1.6e+09    (11.9 GB)  lsoda: ok
neq=46000  lrw=2.116e+09  (15.8 GB)  lsoda: ok
neq=46340  lrw=2.148e+09  (16.0 GB)  lsoda: ERROR: cannot allocate memory block of size 134217728 Tb
neq=84210  lrw=7.092e+09  (52.8 GB)  lsoda: ERROR: cannot allocate memory block of size 134217728 Tb
```

with `adams` succeeding at every size, because its work array is linear
in `neq`. The R-level warning that accompanies the failure is `NAs
introduced by coercion to integer range`.

So the ceiling has two parts. The overflow is a **deSolve** defect and
should be filed there. The reason RTMBode reaches those sizes at all is
the cubic growth of the augmented system, which is by design and only
needs to be said out loud.

### 3.3 Characterization

`repro/07-state-ceiling.R`, `n` independent decays
`dy_i/dt = -exp(mu + u_i) y_i` with `u ~ N(0, 1)` integrated out, one
process per size:

| n | `neq` order 2 | `neq` order 3 | lsoda `fn` | lsoda `gr` | adams `gr` | ode45 `gr` |
|---|---|---|---|---|---|---|
| 2 | 42 | 170 | 0.50318823 | 0.213053 | 0.213043 | 0.213055 |
| 4 | 292 | 2340 | 0.98023318 | 0.573291 | 0.573268 | 0.573293 |
| 6 | 942 | 11310 | 1.533038 | 0.758475 | 0.758434 | 0.758479 |
| 8 | 2184 | 34952 | 2.0308163 | 1.01396 | 1.01393 | 1.01397 |
| 10 | 4210 | 84210 | 2.5397296 | **NaN** | 1.17951 | 1.17961 |
| 12 | 7212 | 173100 | 3.0472225 | **NaN** | 1.36591 | 1.36605 |
| 16 | 16912 | 541200 | 3.9875801 | **NaN** | 1.90124 | 1.90155 |

The break is between n = 8 and n = 10, which matches `neq` crossing
46337 between 34952 and 84210. `lsoda` at n = 8 does work, but only by
allocating 9.1 GB.

### 3.4 The fix

Patch 03 adds `checkWorkspace()`, called once per `ODEadjoint()`
construction, that is once per differentiation order per node. It warns
above 1 GB of requested workspace and warns more sharply above the
integer limit, naming the growth law and the two ways out.

It warns rather than stops. An error raised at that point is swallowed
by the Laplace machinery and reappears as the same unexplained `NaN`; a
warning reaches the user in every path. Verified in
`repro/10-guard.R`:

```
n = 10, method = lsoda
fn: 2.5397296
gr: NaN
WARNING: the order-3 augmented system has 84210 states, so 'lsoda' asks for a
real work array of 7.09e+09 doubles (52.8 GB), which overflows R's integer
range: deSolve cannot allocate it and the solve will fail, so this derivative
will be NaN. The augmented system grows as nstate*(nstate+nparms)^order, and
the order needed is 1 for a gradient, 2 for a Laplace objective and 3 for the
gradient of one. Either reduce the system (solve per group rather than stacking
groups into one system) or choose an integrator whose workspace is linear in
the number of states, such as method = "adams".
```

At n = 6 no warning is issued; at n = 8 the softer 9.1 GB warning is
issued and the gradient is still correct.

---

## 4. What is not explained

- **The `ode45`, `rk4` and `euler` crashes at 32 states** reported
  earlier (process exit 127) were not reproduced here. The sizes above
  cover n up to 16 for `lsoda`, `adams` and `ode45`, and all three
  either succeed or fail cleanly. The crash may be a separate
  allocation failure in the fixed-step methods, which build a dense
  output matrix. It is worth a separate probe before filing.
- **`ODEadjoint()` is reconstructed on every gradient**, but `Df` is
  cached per node, so `checkWorkspace()` warns once per node per order
  and not once per evaluation. This was observed rather than proved from
  the code; if a caller builds many nodes the warning could become
  noisy, and `warning(..., call. = FALSE)` plus a `once` guard may be
  preferable.
- **deSolve prints DLSODA diagnostics from Fortran to the console**
  during a sampler run. Patch 01 stops the failures from aborting the
  chain but not from printing. There is no deSolve argument that
  suppresses them.
- **The `?ode` example's gradient disagrees with central differences by
  10%** at the starting value, on both the stock and the patched build.
  The Lotka-Volterra likelihood over 200 output times is very
  ill-conditioned there, and the step size used was not tuned. This is
  not caused by any patch here, but it is worth a look.
- **`forcings` is unusable** through RTMBode (`initforc should be loaded
  if there are forcing functions`), because `desolve_derivs` has no
  forcing hook, and `approxfun()` in the dynamics is **silently frozen
  at t = 0** because `func2tape()` tapes the derivative function once at
  an all-zero point. That is a fourth report, not patched here.

---

## 5. Issue and pull request text

### 5.1 Issue: a failed ODE solve raises instead of returning NaN

> **Title**: `RTMBode`: a failed `deSolve` call escapes the adjoint node as an R error
>
> `ODEadjoint()$updateSolution()` calls `deSolve::ode()` with no guard.
> When the parameters are extreme the call either raises
> `illegal input detected before taking any integration steps` or returns
> early with fewer rows than `times`, which `RTMB::ADjoint` then reports
> as `Wrong output length`. Both are reachable from a plain `obj$fn()`
> with no sampler and one ODE node:
>
> ```r
> # repro/03-infparm.R
> p <- mode; p[2] <- 500          # a growth rate of exp(500)
> obj$fn(p)   # Error: illegal input detected before taking any integration steps
> obj$fn(mode + 100)              # Error: Wrong output length
> ```
>
> An optimizer can shorten its step and a sampler can reject its proposal
> only if the failure arrives as `NaN`. As an error it is fatal, and when
> the caller is Stan through `tmbstan` it also leaves Stan's nested
> autodiff arena unbalanced for the rest of the session
> (`empty_nested() must be true before calling recover_memory()`), so the
> whole session has to be restarted.
>
> Consequence: **`tmbstan` cannot sample any objective containing an
> `RTMBode::ode()` node.** The chain aborts at warmup iteration 1 even
> when it starts at the fitted optimum, because Stan's initial step size
> is 1 on the unconstrained scale and this objective's gradient is of
> order 1e2 to 1e5 at the mode.
>
> Patch attached. It converts both failure modes to a `NaN` solution of
> the right shape. After it, a 400-iteration chain over a two-node
> Lotka-Volterra model completes with 0 divergences from both a
> mode-centered and an over-dispersed init, and the `?ode` example
> optimizes to a bit-identical answer.

### 5.2 Issue: events on the autodiff path

> **Title**: `RTMBode`: `events` are unusable on the AD path, and two of
> the three methods would be silently wrong
>
> Two defects, the second reachable only once the first is fixed.
>
> 1. `addInfo()` stores an unnamed state vector, because
>    `MakeTape(...)$par()` drops the names `func2tape()` attached.
>    `deSolve` reads state names off that vector, so `checkevents()`
>    bounds the event index by 0 and any event table fails with
>    `too many state variables in 'event'; should be < 0`.
>
> 2. With names attached, `replace` and `multiply` events return wrong
>    gradients and correct values. The augmented system's sensitivity
>    block is not jumped. Measured against central differences: `add`
>    agrees to 3.9e-06, `replace` is off by 42%, `multiply` by 59%.
>    `replace` is also the fallback when no method is given.
>
> A third, found while fixing these: event times that are not already in
> `times` make `deSolve` return more rows than `times`, which the adjoint
> node cannot accept.
>
> Patch attached. It carries `names(y)` through `addInfo()` and onto the
> leading block of the augmented state, refuses `replace` and `multiply`
> by name, warns on `events$func` (which cannot be inspected), and maps
> an expanded solution back onto the requested times. After it, an `add`
> event agrees with finite differences to 6.97e-07 on and off the output
> grid.

### 5.3 Issue: the augmented system outgrows the integrator without saying so

> **Title**: `RTMBode`: silent `NaN` gradients above ~9 states under
> Laplace; `lsoda`'s work array overflows R's integer range
>
> The augmented system at order `k` has
> `neq = sum(nstate * (nstate + nparms)^(0:k))` states, and a Laplace
> *gradient* needs order 3. A 10-state, 10-parameter model therefore
> integrates 84210 equations. `lsoda`'s real work array is
> `22 + neq * max(16, neq + 9)`, computed by `deSolve` in R integer
> arithmetic, which overflows at `neq = 46337`. Past that the solver dies
> with `cannot allocate memory block of size 134217728 Tb` and RTMBode
> reports only a `NaN` gradient.
>
> ```
> n states  neq(order 3)  lsoda gr  adams gr
>    6          11310     0.758475  0.758434
>    8          34952     1.01396   1.01393
>   10          84210     NaN       1.17951
>   16         541200     NaN       1.90124
> ```
>
> Below the overflow it succeeds only by allocating `8*neq^2` bytes:
> 9.1 GB for an 8-state model under Laplace.
>
> Patch attached. It warns above 1 GB and more sharply above the integer
> limit, naming the growth law and the two ways out (solve per group, or
> use an integrator whose workspace is linear in `neq`, such as
> `method = "adams"`). It warns rather than stops because an error at
> that point is swallowed by the Laplace machinery and reappears as the
> same unexplained `NaN`.

### 5.4 Issue to file against deSolve, not RTMBode

> **Title**: `lsoda`: `lrw` is computed in integer arithmetic and
> overflows above `neq = 46337`
>
> ```r
> deSolve::ode(rep(1, 46000), c(0, 1), function(t, y, p) list(-y), NULL)  # ok
> deSolve::ode(rep(1, 46340), c(0, 1), function(t, y, p) list(-y), NULL)
> # Warning: NAs introduced by coercion to integer range
> # Error: cannot allocate memory block of size 134217728 Tb
> ```
>
> `lrw = 22 + neq * max(16, neq + 9)` exceeds `.Machine$integer.max` at
> `neq = 46337`. Computing it as a double, and refusing with a clear
> message when it cannot be allocated, would turn a nonsense allocator
> message into an actionable one. `adams` is unaffected because its work
> array is linear in `neq`.

---

## 6. Reproduction index

Every script takes the library path holding the build under test as its
first or second argument. Pass `""` for the installed build.

| script | what it establishes |
|---|---|
| `repro/lv-model.R` | shared Lotka-Volterra model, one `ode()` call per group |
| `repro/01-sampler.R <ngroup> <lib> repro` | the tmbstan abort, and that one node fails as two do |
| `repro/02-direct-error.R <lib> repro` | random offsets from the mode give non-finite, not errors |
| `repro/03-infparm.R <lib> repro` | both error signatures from a plain `obj$fn()` |
| `repro/04-hypotheses.R <lib> repro` | eliminates static pointers, retaping and gc; shows the short return |
| `repro/05-sampler-long.R <lib> repro` | 400-iteration chains, two nodes, two inits |
| `repro/06-events.R <lib> repro` | events by method, with finite-difference gradients |
| `repro/07-state-ceiling.R <n> <method> <lib>` | the ceiling table, one process per size |
| `repro/08-ceiling-mechanism.R <n> <method> <lib>` | which orders `fn` and `gr` each integrate |
| `repro/09-lrw-overflow.R` | the deSolve integer overflow, in bare deSolve |
| `repro/10-guard.R <n> <method> <lib>` | the warning text the guard produces |
| `repro/11-event-times.R <lib>` | event times off the output grid |
| `repro/12-regression.R <lib>` | the `?ode` example, stock against patched |

To build the patched package:

```sh
cd dev/upstream
gh repo clone kaskr/RTMB
cd RTMB && git am ../patches/*.patch && cd ..
mkdir -p lib
R CMD INSTALL --no-multiarch --library=lib RTMB/RTMBode
```

Then pass `dev/upstream/lib` as the library argument to any script.

`R CMD check` on the patched package is clean: examples pass, and the two
warnings and three notes it reports are all present on the stock package
(build artifacts in `src/`, a period in the `Title` field, `:::` calls
into RTMB, and the `local()` bindings in `ODEadjoint`).
