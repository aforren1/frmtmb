# Phase 0: what a realistic model costs

Date: 2026-09-08. frmtmb 0.55.0; frmtmb.eam 0.5.1, frmtmb.sample 0.3.2,
frmtmb.latent 0.2.1, frmtmb.learn 0.2.1, frmtmb.ode 0.1.1,
frmtmb.spline 0.4.0, frmtmb.coupling 0.2.0. R 4.6.1 on Windows 11, 16
logical cores, 32 GB.

This is Phase 0 of `dev/extension-gaps-plan.md`: one fit per package at
the design that plan's "Realistic scale, per field" table names, timed
in a fresh process, with the estimate against the truth the simulator
used. How the run was made, and what to distrust about it, is in
`dev/scale-lane-notes.md`.

The tier is gated on `FRMTMB_SCALE_TESTS=true` and lives in each
package's `tests/testthat/test-scale.R`, beside a `helper-scale.R` that
is byte-identical in all seven. Each package's `dev/scale-row.md`
points at its own row.

## The short version

Four things came out of this that the survey could not have guessed.

1. **`frmtmb.ode` cannot do a population fit.** 65.8 minutes for one
   `frm()` call on 800 rows. Item 5.4's closed form is a prerequisite,
   not an option, which is the question the plan deferred to Phase 0.
2. **A random effect on a non-decision time is broken, in two
   packages, and one of the ways it breaks is silent.** `wiener()`
   bounds `ndt` by the GLOBAL fastest response, and at the plan's own
   eam design two thirds of the subjects have a true `ndt` above that
   bound. The plain arm does not converge. The `sv` arm CONVERGES, and
   its population `ndt` comes back PINNED at the bound, 0.42 of a
   standard error below it, with a standard error of 7.2e-06 on a
   number that is wrong by ten percent, and `diagnose()` reports
   nothing. `rlddm()` inherits it and takes 12.3 minutes to fail the
   same way.
3. **`frmtmb.sample`'s per-draw loops are not a problem.**
   `posterior_epred()` over 4000 draws and 2000 rows is 1.01 s against
   120 s of sampling. Item 4.6 is not a prerequisite for anything.
4. **And the diagnostic built to catch the second of these has a
   directional blind spot.** `diagnose()`'s `unbounded_dpar` looks
   straight at the saturated `ndt` intercept, passes its magnitude
   threshold, and then rejects it because the standard error did not
   EXPLODE. On a bounded link it collapses instead. That is filed as a
   core seam.

| row | rows | outer par | tape build | 1 gradient (start) | 1 gradient (opt) | frm() with sdreport() | peak R heap | peak process | control | passes |
|---|---|---|---|---|---|---|---|---|---|---|
| `eam` | 12000 | 7 | 2.32 s | 1.59 s | 291.8 ms | 269.2 s | 318 MB | 4780 MB | 2.90 | 4, 269 / 309 / 383 / 422 |
| `eam-sv` | 12000 | 8 | 2.28 s | 2.44 s | 302.4 ms | 117.0 s | 245 MB | 5027 MB | 1.39 | 4, 117 / 117 / 184 / 193 |
| `eam-unbounded` | 12000 | 7 | 2.33 s | 1.85 s | 300.6 ms | 128.5 s | 234 MB | 4819 MB | 1.18 | 3, 128 / 136 / 223 |
| `learn` | 20000 | 5 | 205.3 ms | 227.3 ms | 37.0 ms | 5.72 s | 284 MB | 565 MB | 1.24 | 4, 5.72 / 6 / 7.01 / 7.11 |
| `learn-rlddm` | 20000 | 14 | 4.91 s | 6.47 s | 536.3 ms | 738 s (12.3 min) | 365 MB | 7995 MB | 1.67 | 3, 738 / 1073 / 1175 |
| `hmm` | 25000 | 13 | 839.1 ms | 939.4 ms | 259.7 ms | 117.8 s | 333 MB | 1207 MB | 1.02 | 5, 118 / 121 / 121 / 142 / 192 |
| `lca` | 2000 | 49 | 32.2 ms | 1.4 ms | 1.5 ms | 537.2 ms | 232 MB | 310 MB | 1.19 | 5, 0.537 / 0.964 / 1.25 / 1.98 / 2.07 |
| `ode` | 800 | 6 | 40.1 s | 68.7 s | 14.2 s | 3948 s (65.8 min) | 641 MB | 1094 MB | 1.31 | 1 |
| `spline` | 8000 | 7 | 41.5 ms | 7.8 ms | 7.6 ms | 3.83 s | 346 MB | 1322 MB | 1.57 | 6, 3.83 / 6.02 / 6.43 / 6.8 / 8.31 / 10.7 |
| `coupling-coh-intercept` | 4800 | 4 | 38.1 ms | 545 us | 495 us | 110.3 ms | 224 MB | 310 MB | 1.11 | 4, 0.11 / 0.162 / 0.2 / 0.245 |
| `coupling-coh-cond` | 4800 | 5 | 43.3 ms | 488 us | 476 us | 149.5 ms | 211 MB | 292 MB | 1.11 | 4, 0.15 / 0.224 / 0.231 / 0.354 |
| `coupling-coh-smooth` | 4800 | 9 | 89.9 ms | 8.8 ms | 2.6 ms | 1.39 s | 248 MB | 397 MB | 1.15 | 4, 1.39 / 1.54 / 1.63 / 3.21 |
| `coupling-coh-id` | 4800 | 10 | 121.5 ms | 157.9 ms | 28.5 ms | 8.20 s | 250 MB | 663 MB | 1.14 | 2, 8.2 / 11.3 |
| `coupling-coh-full` | 4800 | 11 | 119.9 ms | 167.6 ms | 31.5 ms | 11.5 s | 239 MB | 704 MB | 1.72 | 4, 11.5 / 14.6 / 14.7 / 30.6 |
| `sample` | 2000 | n/a | n/a | n/a | n/a | 119.3 s | 1129 MB | 1338 MB | n/a | 4, 119 / 120 / 131 / 213 |

## Post-fit calls, timed separately

| row | call | seconds |
|---|---|---|
| `hmm` | `hmm_probs()` | 102.1 ms |
| `hmm` | `hmm_viterbi()` | 55.1 ms |
| `lca` | `lca_probs()` | 11.4 ms |
| `lca` | `lca_profiles()` | 5.0 ms |
| `spline` | `the memoized solve, paid once per fit by the first curve call` | 477.3 ms |
| `spline` | `frm_curve(simultaneous = TRUE), after the solve` | 71.2 ms |
| `spline` | `frm_curve(simultaneous = FALSE), after the solve` | 3.1 ms |
| `spline` | `frm_curve_feature(), after the solve` | 9.6 ms |
| `sample` | `frm() maximum likelihood first` | 428.2 ms |
| `sample` | `posterior_epred(), all draws` | 1.01 s |
| `sample` | `posterior_epred(), a tenth of the draws` | 94.9 ms |
| `sample` | `posterior_predict()` | 1.33 s |
| `sample` | `ranef()` | 280.0 ms |
| `sample` | `log_lik()` | 1.87 s |
| `sample` | `loo()` | 7.27 s |

## How to read the table

**Tape build** is `frm(dry_run = "objective")` less
`frm(dry_run = "frame")`, two interleaved rounds with the minimum of
each taken, because frame assembly is not what a tape build costs and
because the first call through `frm()` in a process pays one-time costs
that belong to neither arm.

**One gradient** is `obj$gr(par)` over a batch that doubles until the
block passes 1.2 s or reaches 256 calls, reported as the best of three
such blocks. It is given at the STARTING values and again at the
optimum, and the two are not the same operation: for a model with
random effects the inner Newton solve starts cold at one and warm at
the other. The cold number is also the less reproducible of the two.
Between the two passes of this table the cold number moved by 3.3x on
hmm (4.35 s against 1.33 s) and 2.5x on eam (4.71 s against 1.88 s),
where the warm number moved by 1.16x and 2.20x. Read the warm gradient
as the measurement and the cold one as an indication.

**`frm()` with `sdreport()`** is one `frm(se = TRUE)` call, which is
what a user pays. For the sample row it is the `frm_sample()` call, and
the maximum-likelihood fit it needs first is in the post-fit table.

**Peak R heap** is `gc()`'s peak since a reset at the top of the row.
It does not see RTMB's tape, which is C++ memory, so it is a FLOOR.
**Peak process** is the process working set, polled twice a second from
outside; it sees everything, R's own startup included. Where the two
differ by a factor of fifteen, as they do on eam, the tape is what the
difference is.

**Control** is two blocks of the same gradient work through the same
path, reported as `max/min`, and the column carries the WORST of the
passes. It must be near 1. Nine of the fourteen rows are at or below
1.31. The four above 1.5 are eam at 2.90, coupling's top rung at 1.72,
rlddm at 1.67 and spline at 1.57, every one of them measured while an
R CMD check or a sibling lane was running; their gradient numbers are
upper bounds and nothing in this document rests on them.

**Passes** is how many times the row was run and, when more than once,
the whole-call times from each. Every number in the table is the
MINIMUM over passes, because contention can only make a measurement
longer. Where two passes disagree by more than a factor of two the
difference is the machine and not the model: the estimates from every
repeated row are identical to every digit recorded, including the two
that did not converge.

The post-fit columns are not single timings at all. `hmm`, `lca`,
`spline` and `sample` run their post-fit calls INTERLEAVED, one of each
per round for three rounds, and report the minimum per call with that
call's round-to-round spread. That machinery exists because the first
version of this document reported a ratio taken from one unreplicated
pair and the ratio was wrong; `dev/scale-lane-notes.md` records it.

## Seeds

ONE seed per row, `20260908`. This is a cost measurement, not a
recovery study: the plan puts recovery at Phase 2, at 60 replicates.
Every estimate below is one draw, and the honest reading of an estimate
near its truth here is "nothing is obviously broken", not "it
recovers". The one place a single draw is decisive is where a fit fails
STRUCTURALLY, which is what the eam and rlddm rows do, and there the
finding was reproduced in a second independent process to every digit.

## The ten-minute rule under contention

Sibling lanes ran on this machine through most of the queue, and the
runner records how many `Rterm` processes were alive at each row's
start and end. Contention can only make a wall clock LONGER, so the
rule is applied one-sidedly:

- a row measured UNDER ten minutes while the machine was loaded is
  under ten minutes;
- a row measured OVER ten minutes is checked against its own
  reproduction before the rule is applied to it.

## What each row decided

### eam: it is minutes, not hours, and the non-decision time is broken

The whole `frm(se = TRUE)` call is 269 s, the minimum of four passes
spanning 269 to 422 s, on 12,000 rows with 90 random effects and 7 outer parameters.
Peak process working set 4.78 GB, against an R heap peak of 318 MB:
almost all of it is the tape. That answers the plan's first question.
Hierarchical DDM at the realistic design is MINUTES.

The second question, "whether the bounded ndt link tolerates a random
effect at all", is answered no, and not subtly. The fit stops with
`nlminb` code 1, a maximum absolute gradient of 1.25e11, a Hessian that
is not positive definite and NaN for all seven standard errors, so the
Wald interval on the condition effect does not exist and the tier's own
assertion on it FAILS. That failure is deliberate and is left in place:
this is the row's finding, not an accident of the harness. The
`(1 | s)` standard deviation on `ndt` is 3.115 on its link, the edge of
a scaled logit, and the between-subject spread it implies on the
natural scale is 0.0058 against the simulator's 0.0263.

The mean structure is unharmed: drift intercept 0.434 against 0.4,
condition effect 0.891 against 0.9, drift standard deviation 0.334
against 0.35, boundary standard deviation 0.173 against 0.20. It is the
non-decision-time component alone that runs away. Two independent
processes produced the same log-likelihood, -7148.81, and the same
estimates to every digit recorded.

### eam with `sv`: it converges, and reports 7.2e-06 on a wrong number

The `variability = "sv"` arm is FASTER than the plain one, 117 s
against 269 s, and it converges: code 0, maximum gradient 7.5e-05,
positive definite Hessian, no bad standard errors. On its own it looks
like the healthy row of the pair.

It is not. The `ndt` random-effect standard deviation is 7.804 on its
link, further out than the 3.115 the non-convergent arm reached, and
the natural-scale spread it implies is 0.0077 against 0.0263. The
population non-decision time comes back as 0.2236 against a truth of
0.25, with a standard error of 7.2e-06, a z of 3644 against the truth.

**It is not merely ten percent low. It is pinned at the bound.** The
global fastest response on this data is 0.22359909 and the fitted
population non-decision time is 0.2235960355, which is 3.05e-06 below
it, that is 0.42 of its own standard error. The estimate has been
pushed against the wall and stopped there.

That standard error is the finding. The scaled logit's derivative
vanishes as the linear predictor leaves, so a delta-method standard
error at the boundary is not small because the estimate is precise; it
is small because the map has flattened: what is reported is the width
of the link's flattened image, not sampling uncertainty. A user sees a
converged model with a tight interval on a parameter that is wrong by
ten percent. By the plan's own Rule 3 that is the worst class of
defect: a silent wrong answer.

### And the diagnostic that should catch it has a directional blind spot

Core 0.55.0's `unbounded_dpar` check does not fire here, and the reason
matters more than the fact, because it is not that the check fails to
look. It looks straight at the right number and then rejects it.

`diagnose_unbounded_dpar()` (`R/confint.R:812`) tests a PAIR on a
dpar's fixed-effect columns, on the link scale: `abs(est) > 10` AND
`(!is.finite(se) | se > abs(est))`. Measured on the two arms, taking
the estimate from `fixef()` and the standard error from the width of
`confint()`'s Wald interval, so nothing here reaches into core:

| arm | `ndt_(Intercept)` on the link | its link-scale se | `plogis(est)`, the fraction of the bound | `abs(est) > 10` | `se > abs(est)` | fires |
|---|---|---|---|---|---|---|
| `sv` | 11.20098624 | 2.372441145 | 0.9999863 | TRUE | FALSE | no |
| plain | 7.166026912 | NaN | 0.9992282 | FALSE | TRUE | no |

Each arm fails a different half, and neither failure is about variance
components.

In the `sv` arm the fixed `ndt` intercept is not fine. It stands at
11.2 on a scaled logit whose upper bound is `min(y)`, that is 0.99999
of the bound. The check's FIRST condition passes. It is the second that
fails: the standard error is 2.37, which is not larger than 11.2.

The premise the check is built on, stated in its own comment, is "the
pair, not either half: a dpar legitimately far out on its link keeps a
small se". That is true of an UNBOUNDED link, where running away makes
the estimate large and the standard error explode. On a BOUNDED link
the same failure does the opposite: the scaled logit saturates at a
modest linear predictor and its derivative vanishes, so an
ILLEGITIMATELY saturated dpar also keeps a small standard error. The
check is calibrated for a parameter running away, and this one is
pinned, so its evidence points the wrong way.

The plain arm shows the other half of the same miscalibration. There
the standard error is NaN, so the second condition passes; it is the
MAGNITUDE threshold that fails, at 7.17 against a threshold of 10, even
though 7.17 on this link is already 0.9992 of the bound. A scaled logit
saturates long before its linear predictor reaches 10, so a threshold
in linear-predictor units cannot mean the same thing on a bounded link
as on an unbounded one.

`extreme_theta` DOES look at variance components, so "nothing looks at
this" would also have been wrong. Its threshold is 8 on a log standard
deviation, and the `ndt` block's log sd is 2.0546284 on the `sv` arm
and 1.136288 on the plain arm, that is 7.80 and 3.12 on the link. The
threshold is link-agnostic; on a scaled logit a standard deviation of
7.8 already pins most subjects at one end or the other.

So the plan's note on item 2.1, that the check "reports a dpar running
to its link's edge, which is the shape this failure takes", is half
right in the most misleading way: the check is built for exactly this
shape and inverts on a bounded link. Four bounded `ndt` links ship
today, in `wiener()`, `lba()`, `rdm()` and `gddm()`. This is filed as a
core seam.

### eam: why, measured

`wiener()` bounds `ndt` by `min(y)`, the global fastest response
(`frmtmb.eam/R/ddm-shared.R:57`, `ddm_ndt_finalize()`), through a
scaled logit. On this row's data, at the seed the row uses:

| quantity | value |
|---|---|
| global fastest response, and so the bound | 0.2262 |
| true subject non-decision times | 0.1996 to 0.2947 |
| subjects whose TRUE ndt is above the bound | 20 of 30 |
| subjects whose true ndt is above their OWN fastest response | 0 of 30 |
| largest true ndt as a fraction of the bound | 1.303 |

Two thirds of the subjects have a true non-decision time the
parameterization cannot express at any value of the random effect, and
not one of them is inconsistent with its own data: every subject's true
`ndt` is below that subject's own fastest response. The bound is global
where the constraint is per subject, which is what item 2.1 predicted
and what the per-subject bound now filed as item 1.0a is for.

The spread the simulator used, 0.0263, is not extreme. It is a
between-subject non-decision-time standard deviation of 26 ms on a mean
of 250 ms, which is at the low end of what the response-time literature
reports, and a quarter of what a uniform below the same bound could
hold.

### eam with the bound lifted: the same data, and everything works

`wiener(max_ndt = 0.45, allow_unreachable = TRUE)` on the SAME 12,000
rows. The two arms differ in the bound and in nothing else.

| | `eam`, bound = min(rt) = 0.226 | `eam-unbounded`, bound = 0.45 |
|---|---|---|
| convergence code | 1 | 0 |
| maximum absolute gradient | 1.25e11 | 0.0029 |
| positive definite Hessian | no | yes |
| standard errors that are NaN | 7 of 7 | 0 of 7 |
| log-likelihood | -7148.81 | -7027.37 |
| `sd(ndt)` on the link | 3.115 | 0.209 |
| `sd(ndt)` natural, truth 0.0263 | 0.0058 | 0.0225 |
| population `ndt`, truth 0.25 | 0.2261, se NaN | 0.2437, se 0.0045 |
| condition effect, truth 0.9 | 0.891, no interval | 0.909 (0.857, 0.962) |
| whole `frm()` call | 269 s | 128 s |

121.4 log-likelihood units better, on the same data, with the SAME
number of parameters, and faster. The bounded arm was not near the
optimum; it was against a wall.

That is the causal evidence: the bound is the cause, not the
optimizer, not the starting values, not the design.

**This arm is the diagnosis and NOT the fix, and nobody should copy
it.** `allow_unreachable = TRUE` does not install a per-subject bound;
it replaces one global bound with a laxer global bound and lifts the
refusal that guards it. Under that flag, with no variability
parameter, the density evaluates a row the model cannot reach at a
decision time of `delta = 1e-9 * min(y)` instead of returning `-Inf`,
so a single floored row would make the reported number not a log
density at all. Checked on this arm's own fit at its own optimum: the
smallest `y - ndt` over all 12,000 rows is 0.0219 s, zero rows are at
or below the floor, and all 30 fitted subject non-decision times are
below their own fastest response. The comparison is valid because the
optimum happened to land inside the per-subject-valid region, which is
a property of this fit and not a guarantee. `?wiener` says the flag
exists for mixtures, where another component covers the unreachable
rows; here there is no other component. Item 1.0a's per-subject bound
is what would make it a guarantee.

### learn: the recursion scales

`bandit2arm_delta()` with `(1 | p | id)` on both parameters, 100
subjects x 200 trials, 20,000 rows: the whole call is 5.72 s, and the
three passes agree to 5.72, 7.01 and 7.11. The tape build is 0.21 s and
one gradient is 0.037 s warm. Peak process 565 MB.

Recovery on one seed is clean. The learning-rate intercept is -0.631
with an interval of -0.761 to -0.501 covering the truth of -0.619; the
two standard deviations are 0.533 and 0.242 against 0.5 and 0.3, each
on the link the family estimates it on; the largest absolute
correlation in the block is 0.148 against a truth of zero.

That answers the plan's question for the delta learner: the per-trial
recursion is not what stops anyone running the population designs in
the literature.

### learn, rlddm: 20 minutes, and the same non-decision-time failure

`rlddm()` with `(1 | p | id)` on `alpha`, `drift`, `bs` and `ndt` at
the same 100 x 200 design, the whole call is 738 s, 12.3 minutes,
taken as the minimum of two passes that read 738 s and 1175 s. The tape
build is 5.45 s and one gradient is 0.54 s warm. Peak process 7.99 GB,
the largest anywhere in this table, and identical to the megabyte
across passes.

It does not converge: code 1, maximum absolute gradient 6.36e9, Hessian
not positive definite, four standard errors NaN.

The correlated block reads `alpha 0.520, drift 1.138, bs 0.212, ndt
2.236` against link-scale truths of 0.50, 1.00 and 0.20, and no stated
truth for `ndt`, whose link is a scaled logit bounded by the data's own
fastest response. The first three recover. The fourth is the eam row's
failure again: `rlddm()` takes its diffusion parameterization from
frmtmb.eam and inherits the same globally bounded non-decision time.
The largest absolute correlation in the block is 0.396 against a truth
of zero, which is what a component pinned at a boundary does to the
rest of the block.

The contrast with the same package's delta learner is the useful part.
`bandit2arm_delta()` on the same subject and trial counts is 7.0 s and
clean. The recursion is not the cost. The diffusion with a globally
bounded non-decision time is.

### hmm: the post-fit passes are not the problem

The plan sent this row to measure `hmm_probs()` and `hmm_viterbi()`,
"which are R loops". At 50 sequences of 500 steps, over three
interleaved rounds, they cost 0.102 s and 0.055 s against a 118 s fit,
so together they are a seventh of one percent of the call. Nothing in
the post-fit surface moves for cost.

The tape build is 0.84 s at 25,000 rows and K = 3, which sits where
`dev/hmm-feasibility.md` put it: that probe measured 1.86 s at 20,000
rows and K = 2 and called T = 20,000 the point where the build becomes
noticeable in an interactive workflow. Still true, and still not the
binding constraint next to a 118 s fit.

Recovery on one seed: the three state means come back as -2.001, 0.004
and 3.008 against -2, 0 and 3; the `tr12` intercept is -1.448 with an
interval of -1.613 to -1.283 covering -1.5; the random-effect standard
deviation on `tr12` is 0.541 against 0.6. `init = "uniform"` was used,
because the default stationary initial distribution is refused when a
transition carries a predictor, and because the simulator draws the
first state uniformly.

### lca: cheap, as the plan said

0.54 s for the whole call at n = 2000, 10 items, K = 4 and two
covariates on membership, the minimum of four passes spanning 0.54 s to
1.98 s. The class profiles match the truth to a maximum absolute error
of 0.083 over the 4 by 10 table under the best relabeling, and the four
class shares are 0.146, 0.203, 0.238 and 0.413. `lca_probs()` is
0.011 s and `lca_profiles()` 0.005 s. Nothing here constrains any
phase.

### ode: 66 minutes, and the closed form is a prerequisite

This is the first question the plan deferred to Phase 0. The answer is
that the segmented sensitivity solve is NOT usable at population scale.

100 subjects x 8 samples, the depot and central system, fourteen doses
over seven days written as one steady-state row plus `ii`/`addl`, the
default `n_ss = 20`:

| quantity | value |
|---|---|
| whole `frm(se = TRUE)` call | 3948 s, 65.8 minutes |
| tape build | 40.1 s |
| one gradient at the starting values | 68.7 s |
| one gradient at the optimum | 14.2 s |
| peak process working set | 1.09 GB |
| solves that failed | 0 of 100 groups |

Six and a half times the ten-minute rule, on 800 rows. The row is not
slow because the data are large; it is slow because every gradient
replays one adjoint solve per group per segment, and this design has
about 35 segments per group once the twenty steady-state cycles are
counted. That is 39.5 s of wall clock per subject, against the 0.4 s
per subject the package's own vignette quotes for the same system with
ONE dose.

It is not poor conditioning. The fit lands on the truth: `ka` 0.984
against 1.0, `ke` 0.1483 with an interval of 0.1398 to 0.1573 covering
0.15, `V` 19.63 against 20, and the two subject standard deviations
0.269 and 0.261 against 0.3 and 0.25. `nlminb` returns code 1, its
false-convergence flag, with a maximum absolute gradient of 0.0061, a
positive definite Hessian and no bad standard errors, which is the
usual shape of that code rather than a failure. The model is right and
the solver is the cost.

So item 5.4 is a PREREQUISITE and not an option. Its own check already
says so: "the Phase 0 design's wall clock, before and after". The
before is 3948 s.

### spline: the feature search is not the cost, and memory is

The plan expected the cost to be in `frm_curve_feature()`, "which calls
`predict()` per Newton step". It is not, but the first version of this
document put the wrong number on that and the correction is worth
following, because it is a limitation of the instrument as much as a
fact about the package.

`frm_curve()` MEMOIZES the joint-precision solve on the fit object.
The solve is paid ONCE per fit, by whichever curve call comes first,
and every later curve call on that fit is cheap. So there is no single
cost for any of these functions; there is a cost for the solve and a
cost for each call after it.

| call | seconds |
|---|---|
| the solve, paid once per fit by the first curve call | 0.48 |
| `frm_curve(simultaneous = TRUE)` after that | 0.071 |
| `frm_curve_feature()` after that | 0.0098 |
| `frm_curve(simultaneous = FALSE)` after that | 0.0031 |

The fit itself is 3.83 s at 8,000 rows and 258 coefficients, so the
whole curve surface is under a fifth of the fit however it is called.

**The number this document first printed for the feature search, 0.0096
s, was a warm number wearing a cold label.** `scale_interleave()`
reports a "first" per arm, and the feature is arm 3, so by the time it
ran the two band arms had already warmed the cache. Measured properly,
three fresh fits with `frm_curve_feature()` called BEFORE anything else
touches the curve machinery, its cold cost is 0.398 s (0.398, 0.418,
0.432), against 0.0097 s on the second call: a factor of 41. The
reviewer measured the same quantity independently at 0.616 s. The tier
now measures the solve explicitly, on a cold fit, through that call.

The conclusion is unchanged and the shape of it is clearer. Nothing in
the curve surface is the cost, and the cold half of it is not the
feature search or the band but the solve underneath both: a fit whose
feature search warms the cache then serves a simultaneous band in
0.074 s, and a fit whose band warms it then serves a feature search in
0.0097 s. Whichever is called first pays.

The memoized solve makes exactly ONE `predict()` call, and the
covariance it assembles agrees with core's to 2.2e-16, one ulp. Note
that `n_predict` is a constant in the source rather than a measurement
of the feature search; what is measured is the covariance agreement.

The number worth carrying forward is the peak process working set:
1.32 GB, against an R heap peak of 346 MB.

The peak location comes back at 0.44925 with a standard error of
0.00355 against a truth of 0.45, and its interval covers. The residual
standard deviation is 0.0607 against 0.06.

### coupling: the whole ladder is seconds, and the contrast recovers

FIVE models on one 4,800-row design, 40 subjects x 2 conditions x 60
frequencies, each the minimum over its passes:

| model | outer par | `frm()` with `sdreport()` | passes |
|---|---|---|---|
| `coh ~ 1` | 4 | 0.110 s | 3 |
| `coh ~ cond` | 5 | 0.150 s | 3 |
| `coh ~ cond + s(freq, by = cond)` | 9 | 1.39 s | 3 |
| the same `+ (1 \| id)`, the survey's 52.5 s model | 10 | 8.20 s | 1 |
| the same `+ (1 \| id:cond)`, the model item 2.6 names | 11 | 11.5 s | 3 |

**The survey's four benchmarks ARE recoverable, and the first draft of
this document said they were not.** The scripts are gone, but the
constructions are recorded in the survey's own session transcript,
which is where the plan's four numbers came from, and two of them are
named there by formula on 4,800 rows:

| survey | model | here |
|---|---|---|
| 52.5 s | `coh ~ cond + s(freq, by = cond) + (1 \| id)` | the `coh-id` rung |
| 28.9 s | the same `+ (1 \| id:cond)` | the `coh-full` rung |

So the promotion the plan asked for is done: the top rung is the
survey's fourth fit and the rung below it is the survey's third. The
two remaining numbers, 2.4 s at 600 rows and 14.6 s at 2,400 rows, put
`s(freq)` on three dpars and `(1 | id)` on four rather than on
coherence alone, so they are a different specification and are not
matched here.

No speedup is claimed from any of these pairs. The survey ran frmtmb
0.53.0 under `load_all()` and this runs 0.55.0 installed, which are not
comparable builds. The defensible statement is that the model that took
28.9 s in the survey takes 11.5 s here.

The top rung answers item 2.6's open question. The condition contrast
comes back at 0.488 with a 95 percent interval of 0.415 to 0.561
against a truth of 0.5, and the components read
`coh: 1 | id = 0.313` and `coh: 1 | id:cond = 0.146` against 0.35 and
0.20. On one seed, at the realistic design, with the standard error the
survey did not record, the contrast is where it should be. That is not
evidence about the survey's 0.77, which was a different construction;
it establishes that the model as item 2.6 spells it recovers its own
truth at this scale.

What the lower rungs get wrong is the INTERVAL, and not the one it
would be natural to guess. The first draft of this document said their
point estimate was "the attenuated marginal one"; it is barely
attenuated at all. Every rung that carries a condition term recovers
the contrast, and the fifth rung says which term the width depends on:

| rung | contrast, truth 0.5 | interval | width | covers |
|---|---|---|---|---|
| `coh ~ cond` | 0.471860 | (0.437641, 0.506080) | 0.0684 | yes |
| `+ s(freq, by = cond)` | 0.480248 | (0.445750, 0.514745) | 0.0690 | yes |
| `+ (1 \| id)` | 0.485773 | (0.451426, 0.520121) | 0.0687 | yes |
| `+ (1 \| id:cond)` | 0.487742 | (0.414908, 0.560576) | 0.1457 | yes |

The attenuation is real, monotone in the right direction, and small:
0.016 across the whole ladder, about a fifth of a standard error. The
INTERVAL is where the damage is, and adding `(1 | id)` does not touch
it: 0.0687 against 0.0684 for a model with no random effect at all.
That is not a defect, it is the design. The condition contrast is a
WITHIN-subject comparison, so a subject main effect is orthogonal to it
and contributes nothing to its standard error. `(1 | id:cond)` is the
whole of the difference, and it doubles the width.

So the sharper statement is not "the cheap rungs miss the subject
variance". It is that a coherence model can carry a subject random
effect, look properly hierarchical, and still report an interval on its
condition effect that is 53 percent too narrow, because the term that
contrast needs is the subject-by-condition one.

The only variance components the two rungs without a subject term
carry are the mgcv SMOOTHING penalties, `coh: s(freq):conda = 1.69`
and `coh: s(freq):condb = 2.03`, not subject variance. Reading those by
position out of `VarCorr()` would report a smoothing penalty as a
subject standard deviation, which is why every row in this tier records
its components by name.

### sample: the per-draw loops are not the problem

This is the second question the plan deferred to Phase 0. The answer is
no: item 4.6's caching is not a prerequisite for anything else in that
package.

A 2000-row binomial GLMM with two crossed factors, 4 chains x 2000
iterations run SERIALLY (`cores` defaults to `getOption("mc.cores", 1)`
and the tier's process sets neither), 4000 post-warmup draws. Sampling
is 119 s over four passes. `posterior_epred()` over all 4000 draws and
2000 rows is 1.01 s, which is 0.25 ms per draw. The five post-fit calls
together are 11.8 s, and 7.3 s of that is `loo()`, which item 4.6 does
not touch. The calls item 4.6 DOES rewrite are `posterior_epred()`,
`posterior_predict()` and `ranef()`: 1.01 s, 1.33 s and 0.28 s, that is
2.6 s in total, or 2.2 percent of the time it took to produce the
draws they run on.

**The superlinearity this document claimed in its first version is not
there, and the claim was the exact error the standing rules were
rewritten about.** "Ten times the draws cost 12.75 times the time" was
`1.59 / 0.1249`: a ratio of two single unreplicated timings. Measured
properly, one call of each arm per round for three rounds, the ratio
for a draw ratio of exactly 10 is:

| measurement | ratio | spread on the full arm | spread on the tenth arm |
|---|---|---|---|
| this tier, run 1 | 10.50 | 2.59 | 1.99 |
| this tier, run 2 | 11.72 | 1.71 | 1.83 |
| the review's independent three rounds | 9.985 | | 1.02 |

The loop is LINEAR in draws. The per-arm spreads of 1.7 to 2.6 are why
a single pair could produce anything between 10 and 13, and they are
also why no claim finer than "linear" is supported by this row. The
conclusion it was attached to is unaffected and in fact strengthened.

The slowest post-fit call by a wide margin is `loo()` at 7.27 s, of
which 1.87 s is `log_lik()` and the rest is loo's own smoothing on a
4000 by 2000 matrix, which no change in this package touches.

Recovery on one seed: the slope is 0.710 with a 95 percent interval of
0.600 to 0.821 covering the truth of 0.8; the two variance components
come back at exp(-0.4377) = 0.645 and exp(-1.113) = 0.329 against 0.7
and 0.4; `looic` is 2455.9 with no Pareto k above 0.7.

## The decision rule, applied

The plan's rule: a fit past ten minutes on its design moves that
package's cost item to the front of Phase 1.

Every row is the MINIMUM over its passes, because contention can only
make a wall clock longer.

| row | whole call | passes | past ten minutes |
|---|---|---|---|
| `eam` | 269 s | 4 | no |
| `eam-sv` | 117 s | 4 | no |
| `eam-unbounded` | 128 s | 3 | no |
| `learn`, delta | 5.7 s | 4 | no |
| `learn-rlddm` | 738 s | 3 | YES, by 1.23x |
| `hmm` | 118 s | 5 | no |
| `lca` | 0.54 s | 5 | no |
| `ode` | 3948 s | 1 | YES, by 6.6x |
| `spline` | 3.8 s | 6 | no |
| `coupling`, top rung | 11.5 s | 4 | no |
| `sample` | 119 s of sampling | 4 | no |

Two packages breach: **frmtmb.ode**, and **frmtmb.learn** through
`rlddm()`.

**The rlddm margin is thin and the first version of this document
overstated it.** That row was reported at 1175 s from a single pass,
which the notes' own procedure forbids: "a row measured OVER ten
minutes is checked against its own reproduction before the rule is
applied to it". Reproduced it is 738 s. Counting the review's independent run of the
same row at 979 s, four measurements of it read 738, 979, 1073 and
1175 s, a spread of 1.59x, and the smallest is 1.23x the rule rather
than 2.0x. Every estimate was identical across
all three and the peak working set agreed to the megabyte, so this is
the machine and not the model.

The breach holds on all three, but 23 percent of headroom is not much
to hang a phase move on, and it does not have to carry one: item 1.0b
moves with item 1.0a under Rule 3 regardless, because it is the same
parameterization and the same defect. The clock is a second reason and
a thin one, and this document says so rather than leaning on it.

### What moves

**Item 5.4, ode closed forms, moves from Phase 5 to Phase 1.** The plan
already said "Whether 5.4 is Phase 5 or Phase 1 is what Phase 0
decides", and 3948 s decides it. Little else in frmtmb.ode pays off
while a population fit takes an hour: item 3.1's NONMEM-shaped reader
hands a user a dataset that then takes an hour to fit. Item 1.2 is a
separate matter and stays where it is. The maximum-likelihood fit IS
reachable, as this row demonstrates by completing it and recovering the
truth; what item 1.2 refuses is `frm_sample()` on top of it.

**The per-subject non-decision-time bound moves to Phase 1**, in both
packages that need it. Only HALF of that is the ten-minute rule: the
rlddm half is, because frmtmb.learn breaches at 738 s, and the eam
half is not, because 269 s is under ten minutes. The eam half moves
under Rule 1, "Any item whose cost the survey could only guess at is
re-ranked after Phase 0's numbers", which is the promotion rule; Rule
3, "a silent wrong answer outranks a refusal", is an ORDERING rule and
is what puts it first once it is here. The silent wrong answer is
measured: a converged fit, no warning, nothing from `diagnose()`, and a
population non-decision time pinned 0.42 of a standard error below the
bound with a standard error of 7.2e-06 on it. The two halves are one
parameterization shared by two packages, so they would have traveled
together anyway.

It is not a new item either: the plan already carries it twice.

- Item 2.1 says "If recovery fails here, item 3.4 gains a per-subject
  bound." Recovery failed. So the per-subject bound is now a Phase 1
  item in frmtmb.eam, separated from the rest of item 3.4, which is
  censoring and is unrelated. The trigger is pulled EARLY and on one
  seed: "here" in item 2.1 meant Phase 2's 60 replicates. That is
  defensible only because the failure is structural rather than
  statistical, and it is the structure that was measured: 20 of 30
  subjects outside the representable set, 0 of 30 inconsistent with
  their own data, and the same log-likelihood to eleven digits in two
  independent processes. A second seed would move which subjects, not
  whether.
- Item 3.3 already contains "rlddm's non-decision bound per subject,
  with `ndt` estimated as a fraction of that subject's floor inside
  `lpdf`". That clause moves to Phase 1 in frmtmb.learn; the rest of
  3.3, which is `session =`, stays in Phase 3.

The plan's Rule 3 is what puts this ahead of the ode item in the
ordering within Phase 1: "a silent wrong answer outranks a refusal, and
a refusal outranks a missing feature". The `eam-sv` arm is a silent
wrong answer with a 7.2e-06 standard error on it. The ode row is slow,
which is a cost, not a wrong answer.

### What does NOT move, and why the rule alone would have got it wrong

**Item 4.6, sample throughput, stays in Phase 4.** frmtmb.sample does
not breach: sampling is 213 s. But the rule would not have caught the
real point either way, which is that item 4.6 addresses the per-draw
loops and those are 1.59 s. Even a package that DID breach on sampling
time would not have been helped by 4.6, because the breach would be in
Stan and not in R. Applying the rule mechanically to a package whose
cost item does not address its cost is how a measurement gets
misapplied, so the rule is applied here with the number that says which
part is slow.

**Item 4.2, hmm on new data, stays in Phase 4.** The plan's Phase 0
entry for hmm asks about the post-fit passes because they are R loops.
They are 0.254 s and 0.088 s.

**Item 1.3 and item 4.1, the spline curve items, stay where they are.**
The feature search the plan worried about is 0.0098 s once the memoized
solve is paid, and 0.40 s when it is the call that pays it. Either way
it is a fraction of a 6.0 s fit.

**The coupling items stay where they are.** The most expensive model in
that package's realistic design is 11.5 s.

## What this tier does not settle

**One seed.** Every estimate here is a single draw. A row whose
estimate lands near its truth says only that nothing is obviously
broken. Phase 2's 60 replicates are what recovery means.

**One machine, and a busy one.** Every number is a minimum over the
passes each row got, one or two, on a machine that had sibling lanes on
it for most of the queue. The rows that got two passes agreed on every
estimate and disagreed on wall clock by up to 2.7x. Nothing here is a
cross-platform claim.

**The `ss` cost is not separated from the dosing cost.** The ode row
carries both fourteen dose events and twenty steady-state cycles, so
its 3948 s is not attributed between them. An `ss = FALSE` arm would
have said which, and it was not run. It does not change the conclusion,
because a closed form removes both.

**Nothing was profiled.** The tier reports wall clock and one
allocation counter. Where the eam row's 4.78 GB of process working set
goes inside RTMB is not something this measures.

## Two things to file, not to wait on

Both are worth writing before Phase 1 lands, because a user fitting a
hierarchical DDM today gets a confident wrong answer with no warning.

**A NEWS bullet for frmtmb.eam**, drafted verbatim in
`extensions/frmtmb.eam/dev/scale-row.md` under "Disclosure to ship".
This lane did not edit `NEWS.md`, because that file is organized by
released version and picking the version is not a measurement lane's
call.

**A core backlog entry** for the diagnostic, lift as is:

> `diagnose()`'s `unbounded_dpar` check cannot see a distributional
> parameter pinned at the edge of a BOUNDED link. Its evidence pair,
> `abs(estimate) > 10` and a standard error larger than the estimate,
> is calibrated for `student()`'s `nu` on `logm1`, where running off
> makes the standard error explode. On a scaled logit the derivative
> vanishes at the edge, so the standard error COLLAPSES and the linear
> predictor saturates below 10. Measured on the Phase 0 eam rows:
> `ndt: (Intercept)` at 11.20 with a link-scale standard error of 2.37
> (the check needs it above 11.20), and 7.17 in the non-convergent arm,
> below the magnitude threshold entirely. Both fits sit at 0.9992 and
> 0.99999 of the bound. `extreme_theta` does look at random-effect log
> standard deviations, but its threshold of 8 is link-agnostic and the
> `ndt` component stands at 2.05. Four bounded `ndt` links ship today
> (`wiener()`, `lba()`, `rdm()`, `gddm()`), so the check needs to know
> whether a link is bounded and, if it is, to test proximity to the
> bound rather than magnitude and a large standard error.

It is also in the plan's Core seams table, needed by items 1.0a and
1.0b.


