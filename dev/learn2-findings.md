# frmtmb.learn, round two: factorization, RLDDM, and the two models left out

Lane `wt-learn2`, worktree `C:/Users/adf44/source/r/frmtmb-wt-learn2`,
branch `wt-learn2` at main's head (54d4d92, frmtmb 0.53.0). Private
library `.../scratchpad/l3-lib`, holding core, `frmtmb.eam` and
`frmtmb.learn` built from this worktree. Every number here was produced
by a script in this lane.

## The three items, and how they came out

1. The families now declare a factorization, `frm(importance =)` is
   admitted, and what the correction does to a learning fit is measured.
   It is usable at 100 trials on a dataset whose variance component the
   Laplace fit kept, moving it +0.11 to +0.18 log units toward the
   truth with a Monte Carlo error a third of that; it is not usable at
   20 trials, and the shift it reports there is its own step cap.
2. `rlddm()`, a delta rule driving the drift rate of a Wiener
   first-passage density, importing `frmtmb.eam` through one new
   export.
3. `igt_orl()` and the Kalman filter's exploration bonus, both shipped
   with recovery tables, and both tables contradict the reason the
   first release gave for leaving them out.

THREE CLAIMS THIS LANE WROTE DOWN BEFORE MEASURING AND THEN HAD TO
CORRECT, all in the same direction: an expected trade-off that was not
there.

| written first | measured |
|---|---|
| `rlddm()`'s `drift` and `bs` trade off | 0.17. `bs` with `ndt` is -0.57 and `drift` with `bias` is +0.56 |
| `igt_orl()` is weakly identified | every parameter recovers; what is weak is `betaF` against `betaP`, at -0.77 |
| the Kalman bonus does not separate from `tau` | 0.12. It does not separate from `sigmaD`, at -0.81 |

Each help page now says the measured thing.

## Reading, before any edit

`extensions/frmtmb.learn` (2178 lines of R across 12 files),
`dev/learn-findings.md` (812 lines), `dev/structured-family-protocol.md`
(567), `dev/protocol-findings.md` (565), `inst/rl/rw-delta.R` (the
reference consumer of the two new slots), `R/importance.R:100-320` and
`:660-780`, `R/structure.R:740-840`, `R/predict.R:2395-2460`,
`extensions/frmtmb.eam/R/wiener-density.R`.

## What landed, with file and line

Paths relative to the worktree root.

| what | where |
|---|---|
| the recursion, now yielding per-trial factors | `extensions/frmtmb.learn/R/engine.R:173` |
| the stacking the correction calls the slots with | `.../R/engine.R:141` |
| `loglik`, the sum of the pieces | `.../R/engine.R:331` |
| `loglik_group`, one value per subject | `.../R/engine.R:345` |
| `loglik_row`, one value per trial, with the saturated attribute | `.../R/engine.R:361` |
| the trace's summary column name, `p` or `dens` | `.../R/engine.R:117` |
| the three new spec slots (`logp`, `draw`, `choice_col`/`choice_map`) | `.../R/engine.R:429` |
| the block's reserved `group` | `.../R/family.R:84` |
| the two factorization slots declared | `.../R/family.R:262` and `:266` |
| the nominal-versus-density refusal split | `.../R/family.R:219` |
| the `importance` compat row, `refused` to `works` | `.../R/zzz.R:97` |
| `rlddm()`'s density choice rule | `.../R/rlddm.R:178` |
| `rlddm()`'s bounded non-decision-time link | `.../R/rlddm.R:277` |
| the `sim_needs` declaration and its guard | `.../R/rlddm.R:227`, `.../R/task.R` |
| `igt_orl()`'s three-store update with the rate swap | `.../R/igt-orl.R:123` |
| the Kalman exploration bonus | `.../R/bandit4arm2-kalman.R:134` |
| `wiener_lpdf()`, the one new export | `extensions/frmtmb.eam/R/extension-api.R:79` |

## ITEM 1. What the families can honestly declare, and what they now do

### The two statements, established before either slot was written

The protocol asks how finely a family's likelihood factorizes, and this
one factorizes twice. Both are properties of the recursion rather than
of any one learning rule, so both hold for all eight families.

**Over subjects.** Two subjects share the parameters and nothing else:
no value store crosses a subject boundary, which `ln_pack()` enforces by
laying the rows out one subject per matrix row. So the response's
likelihood is a product over subjects and `loglik_group` is exact, with
the subject as the group.

**Over trials.** Given the parameters, a subject's likelihood is

    L(subject) = prod_t P(choice_t | history_t, theta)

because the value store at trial `t` is a function of that subject's own
rows before `t` and of the parameters. Each factor is therefore the
conditional log-density of one ROW given earlier rows, which is exactly
what `loglik_row` is defined to be, and each depends only on its own
subject's random effects. So the per-row slot is honest here, and it is
honest for a reason `frmtmb.latent::hmm()` does not have: a row's
emission density there is not its contribution to the likelihood,
because the state that emitted it was reached through every earlier row
and is summed over. A learning rule conditions on the history rather
than integrating it out, so the factor exists.

`unit` is untouched and keeps its own answer, "one subject's trial
sequence". Dropping a trial changes every later trial's value store, so
the leave-one-out unit is the whole sequence even though the
factorization is per trial. That is the distinction the protocol makes
and this family is the counterexample it was made for.

### How they are computed: one recursion, three slots

`extensions/frmtmb.learn/R/engine.R:167` `ln_recurse()` no longer
accumulates a running total. Its `"terms"` mode returns one masked
vector per trial, over subjects, and the three slots are one line of
arithmetic each:

| slot | file:line | expression |
|---|---|---|
| `loglik` | `R/engine.R:281` | `sum(Reduce(+, terms))` |
| `loglik_group` | `R/engine.R:295` | `Reduce(+, terms)` |
| `loglik_row` | `R/engine.R:313` | the terms scattered back to their rows |

Nothing extra is computed for them. The walk always formed one factor
per trial per subject and the total it used to accumulate was the sum of
exactly these, so the three cannot drift: there is one definition.

The row slot scatters with ONE sub-assignment rather than one per trial,
because `[<-` on a taped vector copies the whole vector each time it is
called, and this package measured that copy costing 3.6x to 14.2x at
tape build last round. Padded cells are dropped rather than written: a
pad repeats its subject's FIRST row number, so writing it would
overwrite that trial's own factor with a masked zero.

### Did folding the total into the pieces cost the hot path anything?

Measured on this machine, both shapes in ONE process, interleaved, five
rounds, `bandit2arm_delta()` at 40 subjects x 100 trials, gradients over
batches of 20. The machine was running three sibling lanes, so read the
columns against each other and not against the timing table in
`dev/learn-findings.md`.

| `loglik` shape | tape build | one gradient | objective |
|---|---|---|---|
| new: sum of the per-trial pieces | 0.250 s | 98.5 ms | 2327.9311591473 |
| old: scalar accumulated per trial | 0.260 s | 103.0 ms | 2327.9311591473 |

Identical objective to full printed precision, and the new shape is not
slower on either column. It has to be: the recursion produced one vector
per trial either way, and the old shape summed each vector as it went
where the new one sums them at the end. Same additions, different order.
This reproduces what the core lane measured on `inst/rl/rw-delta.R` when
it made the same change there.

### The stacking contract

`R/engine.R:129` `ln_stack()` is the whole of what the protocol's
stacking costs a family whose loop is already vectorized across
subjects: row `i` of replicate `k` is at `i + (k - 1) * n`, so the
block's row NUMBERS shift by that and nothing else does. The walk runs
over subject-crossed-with-replicate, which is the same loop over a
longer vector.

`tests/testthat/test-factorization.R` pins it three ways: the tiled
design gives tiled values in replicate-major order, a replicate whose
parameters differ gives different values (without which the first
assertion would pass on a family that ignored the stacking entirely),
and the per-subject values are the sums of that subject's own rows in
the level order `structure_group_codes()` fixes.

### The block's reserved `group`

`R/family.R:76`. `ln_pack()` now returns `group = gv` beside the
`subject = gv` this package already had. They are one column with two
readers rather than two columns that could drift, because both are the
same object: the core reads `group` to align its own grouping against
the family's before it admits the correction and to order what
`loglik_group` returns, and `frm_value_trace()` reads `subject`.

### The by-name refusal is gone

`R/zzz.R`, the `importance` compat row, is now `works` with the
measurements below in its note. `test-surface.R` asserted the old
refusal and now asserts the new status; the `frm(importance = 50)`
error assertion is deleted rather than inverted, because a fit that
takes 30 seconds does not belong in the ungated surface file. The run
that proves the seam is in `test-factorization.R` on an 8-subject
fixture.

TWO REFUSAL SENTENCES CHANGED because declaring the slots made them
false. The deviance refusal said "the loglik slot returns one total, so
the core never sees a saturated per-row comparison"; the core sees one
now, and what is missing is the SIGN, since the response is a nominal
option code with no mean to depart from. `R/family.R:150`. The mixture
refusal said "same seam as importance", which no longer names anything;
it now says what is actually left, a `mixture()` that reads a declared
factorization instead of a rowwise density, and that this is a change
under core's `R/`.

### What the correction actually does, measured

`extensions/frmtmb.learn/dev/learn-importance.R`.
`bandit2arm_delta()` on a reversal task, 40 subjects, learning rate 0.30
before the reversal and a +1.0 logit effect after it, softmax
sensitivity 3, subject-level sd 0.5 on the learning rate's logit, one
scalar random intercept. Six datasets per cell, seeds 4301 to 4306.

PER DATASET, NOT AVERAGED, and that is a decision the first version of
this table got wrong. A variance component estimated from short binary
sessions collapses on some datasets and not others, and one collapsed
replicate moves a mean of six by more than the whole effect being
measured: the first run reported a mean Laplace `log sd(alpha)` of
-7.824 at 20 trials, which is not a number about any dataset.

#### 40 subjects x 100 trials, importance = 100

Four of six completed. The two that did not are the correction's own
convergence guard, not the seam: seeds 4303 and 4306 refuse with "the
corrected negative log-likelihood ROSE ... so the iteration moved away
from the answer".

| seed | sd Laplace | sd corrected | log shift | mcse | min ESS/draw | seconds |
|---|---|---|---|---|---|---|
| 4301 | 0.0004 | 0.0025 | +1.822 | 0.000 | 1.000 | 40.5 |
| 4302 | 0.4170 | 0.4841 | +0.149 | 0.048 | 0.925 | 31.2 |
| 4304 | 0.4543 | 0.5053 | +0.106 | 0.064 | 0.874 | 77.6 |
| 4305 | 0.3798 | 0.4541 | +0.179 | 0.049 | 0.924 | 76.9 |

Fixed effects over the four, truth in brackets: `alpha_(Intercept)`
-0.004 [-0.847], `alpha_after` +0.010 [1.000], `tau_(Intercept)` +0.001
[1.099]; largest single shift 0.017. The Laplace fits take 1.8 s and the
corrected ones 31 to 78 s.

Read the three well-behaved rows together. Truth is sd = 0.5; the
Laplace fits give 0.38 to 0.45 and the correction moves them to 0.45 to
0.51, that is UP and TOWARD the truth, by 0.11 to 0.18 log units, with a
Monte Carlo standard error of 0.05 to 0.06 and a smallest effective
sample size per draw of 0.87 to 0.93. Those diagnostics are good: near 1
means the weights inside every group are nearly uniform.

That is a smaller and better-diagnosed shift than the core lane measured
on the reinforcement-learning vignette's model (about +0.36 to +0.48 at
100 draws, min ESS/draw down to 0.03). The difference is the MODEL, not
the family: their figure is for a correlated pair of subject effects on
two parameters, `(1 | p | id)` twice, which is a weakly identified 2x2
covariance; this is one scalar variance component. Their own
scalar-model table gives +0.12 to +0.18 at 200 draws, which is what this
reproduces on a different family.

#### 40 subjects x 20 trials, importance = 100 and 400

All six completed at both draw counts, and five of the six say nothing,
because the quantity the correction exists to correct is already gone.

| seed | sd Laplace | sd at 100 draws | sd at 400 draws |
|---|---|---|---|
| 4301 | 0.3707 | 0.7292 | 0.6226 |
| 4302 | 0.0001 | 0.0006 | 0.0005 |
| 4303 | 0.0001 | 0.0008 | 0.0006 |
| 4304 | 0.0001 | 0.0006 | 0.0005 |
| 4305 | 0.0001 | 0.0007 | 0.0005 |
| 4306 | 0.0001 | 0.0005 | 0.0004 |

**The +1.822 those five rows report as a log shift is not an estimate;
it is the correction's own step cap times its round count, and the run
says so.** Checked directly on seeds 4302 and 4303: `rounds = 5`,
`capped = TRUE`, `moves = 0.3645 0.3645 0.3645 0.3645 0.3645` and
`grads = 3.02` five times over, with the warning "used all 5 of its
rounds and the estimates were still moving by 0.364 at the last one".
Five capped steps of 0.3645 is 1.8225. At 400 draws the cap resolves to
0.3004 and the reported shift is 1.502, five times that. The iteration
never converges; it walks at its step limit until the rounds run out.

So the honest statement about 20-trial data is not "the correction
overstates the variance component". It is that on five of six datasets
the Laplace fit has already collapsed sd(alpha) to 1e-4, there is
nothing at that point for the correction to reweight, and it reports its
own step cap while warning that it did not converge. Seed 4301, the one
dataset whose Laplace fit kept a real variance component, is the only
informative row: 0.371 goes to 0.729 at 100 draws (mcse 0.103, min
ESS/draw 0.743) and to 0.623 at 400 (mcse 0.036, min ESS/draw 0.928),
against a truth of 0.5. It OVERSHOOTS at 100 draws, and more draws bring
it back down, which is what a Monte Carlo error of 0.10 predicts.

The fixed effects at 20 trials move by 0.003 to 0.007 on average with a
largest single shift of 0.058, so the conclusion the recovery study
already reached, that a condition effect on a learning rate survives a
short session, is unchanged by the correction.

#### Where the correction is usable, stated plainly

* At a realistic trial count with a well-identified scalar variance
  component: YES. The shift is +0.11 to +0.18 log units, its Monte Carlo
  error is a third of that, the effective sample sizes are near 1, and
  it moves the estimate toward the truth. It costs 30 to 80 seconds
  against 2 for the Laplace fit.
* At a realistic trial count on a dataset whose Laplace fit collapsed
  the component: NO, and the run returns without complaint. Check
  `sqrt(VarCorr(fit))` before believing a shift.
* At 20 trials: NO, on this design. Either the component has already
  collapsed, in which case the reported shift is the step cap, or it has
  not and the correction overshoots by more than the bias it is
  correcting. The warning about capped rounds is the signal, and
  `fit$importance$capped` with `moves` all equal is what it looks like
  in the object.
* A model grouped on anything but the family's own subject is refused by
  name, by `check_importance_scope()`, before any of this. That refusal
  is pinned in `test-factorization.R`.

## ITEM 2. RLDDM, and the one export it needed from frmtmb.eam

### LOUD, for whoever else is in frmtmb.eam this round

**This lane added one exported function to `frmtmb.eam` and bumped that
package to 0.4.1.** The file is new and nothing existing was edited
except `NAMESPACE`, `NEWS.md` and `DESCRIPTION`'s `Version:` line, which
is deliberate: a sibling lane working in `wiener-density.R` or
`wiener-family.R` will not collide with it.

* `extensions/frmtmb.eam/R/extension-api.R`, new, one export.
* `wiener_lpdf(dt, drift, bs, bias, upper)`, the Wiener first-passage
  log density with `wiener()`'s parameterization and `wiener()`'s tape
  safety. It is `ddm_lpdf_both()` with a guard on `upper` and a help
  page; the guard runs only when `upper` is not an advector, because
  `upper` is data and the other four are parameters.
* `extensions/frmtmb.eam/NAMESPACE`: one line, `export(wiener_lpdf)`.
* `extensions/frmtmb.eam/man/wiener_lpdf.Rd`, generated.

WHY AN EXPORT RATHER THAN THREE COLONS. Every Wiener function in that
package is `@noRd`, so the only route to the density was
`frmtmb.eam:::ddm_lpdf_both()`, which is a promise nobody made. The
brief allowed one small seam and this is it. Nothing else was opened:
the series truncations, the blend between them, the across-trial
variability integrals and the CDF all stay internal.

WHAT WAS NOT EXPORTED, and the duplication that decision cost. The
non-decision time needs a logit link scaled onto `(0, min(rt))`, because
the density is zero at and below `ndt` and a log link lets an optimizer
step over that edge. `frmtmb.eam` builds exactly that link in
`ddm_ndt_finalize()`, and asking for it would have been a second export
whose subject is link construction rather than density. So
`frmtmb.learn` writes its own, `R/rlddm.R` `ln_ndt_link()`, ten lines of
the same standard construction. That is a duplication and it is recorded
here rather than hidden; the comment above the function says so too.

A SECOND SEAM IS RECORDED AND NOT TAKEN. `rlddm()` declares no
`fitted_mean`, and unlike the seven softmax families the reason is not
that the response is nominal: a response time is ordered and its
conditional mean exists. What is missing is a route to it. The mean of a
Wiener first-passage time conditional on the boundary reached is a
closed form that belongs to `frmtmb.eam` (`ddm_mean_rt()` is there,
internal), and deriving it again here would be worse than the ndt link,
because it is arithmetic a reader cannot check by eye. So `fitted()`
refuses, and the compat row says which seam would close it.

### Verified before the family was written

The convention question was settled by measurement rather than by
reading, because getting a boundary the wrong way round produces a model
that fits and is wrong.

| check | result |
|---|---|
| `wiener_lpdf(t, v, a, w, 1)` integrated over t, at v = 1.2, a = 1.5, w = 0.5 | 0.85814894 |
| analytic P(upper) from the diffusion's exit probability | 0.85814894 |
| the same at w = 0.3 | 0.67895608 against 0.67895608 |
| against `RWiener::dwiener(..., resp = "upper", give_log = TRUE)` | max abs diff 8.9e-16 |
| the lower boundary, against `resp = "lower"` | 1.8e-15 |

So `upper = 1` is the upper boundary, `bias` is the relative start
toward it, and positive drift favors it. That is also Stan's
`wiener_lpdf` convention, which is what makes the identity program a
direct transcription rather than a reparameterization.

### What rlddm() is

`extensions/frmtmb.learn/R/rlddm.R`. Pedersen, Frank and Biele (2017): a
delta rule whose value difference drives the drift rate,

    v_t = drift * (Q[upper] - Q[lower])

with the response time the first passage of a diffusion with boundary
`bs`, start `bias` and non-decision time `ndt`. Five parameters, all
ordinary distributional parameters: `alpha` (logit), `drift` (identity),
`bs` (log), `ndt` (a logit scaled onto the interval up to the fastest
response) and `bias` (logit).

**hBayesDM does not carry this model, and the help says so rather than
naming one it is not.** hBayesDM ships the two halves separately:
`choiceRT_ddm` is this choice rule with no learning, and
`bandit2arm_delta` is this learning rule with a softmax. The parameter
map is a rename where it exists at all: `choiceRT_ddm`'s `alpha` is the
boundary separation, which is `bs` here, its `beta` is `bias`, its
`delta` is `drift` times the value difference, and its `tau` is `ndt`.
`alpha` here is the LEARNING RATE, as in every other family in this
package, and that collision is why the diffusion parameters carry
`frmtmb.eam`'s names rather than hBayesDM's.

### The engine seam it needed

Three new `ln_spec()` slots, all optional and all NULL for the seven
softmax families, so nothing about them changed:

* `logp(state, d, ch)`, a log density per subject, replacing `choice`
  and the engine's softmax. There are no utilities to normalize when a
  trial's contribution is a first-passage density.
* `choice_col` and `choice_map`, because the option taken arrives in
  `dec()` rather than in `y`, coded 0 and 1 rather than 1 and 2.
* `draw(state, d)`, the simulator's counterpart, which writes its own
  answer back into `d` in the coding the addition term uses, because
  only the family knows that coding.

`frm_value_trace()`'s summary column is named `dens` rather than `p` for
such a family, and that is not cosmetic: `exp(lp)` is a probability in
the softmax case and a DENSITY in this one, which may exceed 1. Calling
it `p` would invite a reader to check it against a probability.

`simulate()` is refused, for ts_par7's reason: one trial's draw is two
numbers, the boundary and the time, and a response vector holds one.
Drawing the time alone against the observed choice is a draw from a
different model. `frm_task_simulate()` returns whole data frames with
both columns and is the route.

### Recovery, 30 subjects by 100 trials, 60 replicates

`extensions/frmtmb.learn/dev/learn-recovery-round2.R`. One random
intercept on the learning rate, sd 0.4. Everything on the NATURAL scale,
with the Wald interval endpoints pushed through the same monotone link;
`ndt`'s bound is that replicate's own fastest response, so its link
scale is not comparable across replicates and this is not a choice.

| parameter | truth | bias | mc se | sd of est | coverage | fits |
|---|---|---|---|---|---|---|
| `alpha` | 0.35 | -0.002 | 0.003 | 0.026 | 0.87 | 60 |
| `drift` | 3.00 | +0.001 | 0.009 | 0.069 | 0.98 | 60 |
| `bs` | 1.60 | 0.000 | 0.002 | 0.019 | 0.95 | 60 |
| `ndt` | 0.20 | +0.001 | 0.000 | 0.003 | 0.98 | 60 |
| `bias` | 0.50 | -0.002 | 0.001 | 0.007 | 0.95 | 60 |

Every bias is inside its Monte Carlo error and 60 of 60 replicates
produced a usable interval. Coverage is nominal for four of the five and
0.87 for `alpha`, so a learning rate's Wald interval is slightly
optimistic on this design.

**The pair that trades off is not the one to expect, and this lane wrote
the wrong one down before measuring it.** A first-passage density is
driven largely by the RATIO of drift to boundary, so `drift` and `bs`
look like the pair at risk. Measured, their estimates correlate at 0.17
across replicates, which is nothing.

| | alpha | drift | bs | ndt | bias |
|---|---|---|---|---|---|
| alpha | 1.00 | -0.28 | 0.18 | -0.01 | -0.29 |
| drift | -0.28 | 1.00 | 0.17 | -0.31 | 0.56 |
| bs | 0.18 | 0.17 | 1.00 | -0.57 | 0.04 |
| ndt | -0.01 | -0.31 | -0.57 | 1.00 | -0.49 |
| bias | -0.29 | 0.56 | 0.04 | -0.49 | 1.00 |

What co-varies is `bs` with `ndt` at -0.57 and `drift` with `bias` at
+0.56, and both make sense after the fact: the boundary and the
non-decision time are the two ways to make responses slower, and the
drift and the start point are the two ways to favor a boundary.
`?rlddm` now says that instead of the guess it opened with.

### Identity against Stan, the three new rows

Same method as last round: each program declares exactly the parameters
frmtmb estimates on the scales frmtmb estimates them on, so the map is
the identity and the comparison carries no Jacobian; the random-effect
standard deviation rides along as data. Run with
`FRMTMB_BRMS_FIT_TESTS=true`, `NOT_CRAN=true`, a private copy of the
Stan cache and `R_MAKEVARS_USER` at a copy of
`dev/stan-cache/makevars-cxx17.mk`.

| family | subjects, rows | Stan log_prob | frmtmb joint | difference | max grad on subject effects | population grad vs frmtmb's |
|---|---|---|---|---|---|---|
| `bandit4arm2_kalman_filter (bonus)` | 10, 500 | -441.0315488439 | -441.0315488439 | -2.842e-13 | 5.284e-10 | 8.793e-14 |
| `igt_orl` | 12, 960 | -935.0557464773 | -935.0557464773 | -6.821e-13 | 3.055e-11 | 6.484e-14 |
| `rlddm` | 12, 720 | -224.7348379010 | -224.7348379042 | **3.176e-09** | 1.669e-11 | 1.666e-07 |

The eight rows from last round were re-run unchanged on this worktree
and reproduce: -1.7e-13 to 5.7e-13, the same numbers to the digit where
the fixture is unchanged.

**`rlddm` is the one row that is not exact, and the reason makes it the
strongest row rather than the weakest.** Every other program in the file
evaluates the SAME arithmetic as the engine in a different order, so the
two agree to the level of double accumulation. Stan implements the
Wiener first-passage density ITSELF, with its own truncation, where
`frmtmb.eam` evaluates Navarro and Fuss's two series and blends them
smoothly in `log(u)`. So this row compares two independent
implementations of the density as well as two of the recursion.

WITHDRAWN, AND REPLACED. The first version of this section argued that
the residual is a per-row implementation difference rather than
accumulated rounding because "the residual GROWS with the row count
while the per-row figure stays flat". Its own table says neither thing.

| rows | Stan | frmtmb | difference | per row | relative |
|---|---|---|---|---|---|
| 120 | -45.2481141123 | -45.2481141135 | 1.188e-09 | 9.90e-12 | 2.6e-11 |
| 240 | -6.2154828182 | -6.2154828185 | 2.500e-10 | 1.04e-12 | 4.0e-11 |
| 720 | -224.7348379010 | -224.7348379042 | 3.176e-09 | 4.41e-12 | 1.4e-11 |
| 1440 | -290.8491641722 | -290.8491641780 | 5.837e-09 | 4.05e-12 | 2.0e-11 |

240 rows has a residual five times SMALLER than 120 rows, so the total
does not grow with the row count; and the per-row column spans a factor
of 9.5, so it is not flat. The relative column is the one that is
stable, at 1.4e-11 to 4.0e-11 against about 1e-15 on the exact rows, and
even that cannot separate a per-row bias from cancellation: four points
of a SIGNED difference between two totals never could. The argument was
wrong and the reviewer was right to say so.

THE ARGUMENT THAT DOES WORK asks a third implementation which side the
difference is on, and it needs one measurement rather than four.
Measured here over 4000 rows spanning the region the fixture visits,
`dt` in (0.05, 3), drift in (-4, 4), boundary in (0.6, 2.5) and bias in
(0.2, 0.8), `frmtmb.eam`'s density against `RWiener`'s:

| statistic | per-row absolute difference in the log density |
|---|---|
| median | 0 |
| mean | 3.38e-16 |
| 90th percentile | 8.88e-16 |
| 99th percentile | 3.55e-15 |
| max | 1.95e-14 |

The reviewer ran the same measurement on their own draw and got a max of
1.78e-14 with the same median and quantiles, so the two independent
draws agree on the conclusion as well as the order of magnitude.

So `frmtmb.eam`'s density is right to about 1e-15 per row against an
implementation neither it nor Stan shares. A per-row difference of about
4e-12 against Stan therefore sits on **Stan's** side of the comparison,
which is a stronger and more specific claim than the symmetric "the two
series agree to eleven significant figures" the first version made.

The test's own tolerance is 1e-6 relative, so the row passes with four
orders of magnitude to spare, and what it buys that no other row buys is
a check on the DENSITY. The population-gradient identity is 1.67e-07 for
the same reason, against 6.5e-14 and 8.8e-14 on the two exact rows.

The four-point check that opened this section was also run before the
family was written, at a single set of parameters: 8.9e-16 at the upper
boundary and 1.8e-15 at the lower. That was right and far too narrow to
carry the conclusion it was asked to carry.

Gated tier: **55 assertions, 0 failures**, 11 identity rows. Two
warnings, both pre-existing and both on the `bandit2arm_dual (pe)`
fixture, whose `sign(pe)` selector puts a kink in the joint density and
stops TMB's inner Newton solve short of the mode; `dev/learn-findings.md`
records that and the test relaxes only the mode assertion, never the
value.

## ITEM 3. The two models left out for time

### igt_orl, and the claim it was left out on

`extensions/frmtmb.learn/R/igt-orl.R`. Haines, Vassileva and Ahn (2018).
Three value stores per deck rather than one: the expected value `EV`,
the expected frequency of a gain `EF`, and perseverance `PS`. The two
learning rates SWAP between the played deck and the three that were not,
on the sign of what the played deck returned, which is the asymmetry the
model is named for. The choice utility is
`EV + betaF * EF + betaP * PS` with the softmax sensitivity FIXED at 1,
which is hBayesDM's parameterization and is kept.

The sign branch is on DATA (the payoffs are addition-term columns and
the indicator that picks one is the observed choice), so it reaches the
tape as a multiplication by zeros and ones rather than as a comparison,
the same way `bandit2arm_dual(split = "outcome")` does.

ONE TRANSFORM against hBayesDM, and it is documented rather than
hidden: hBayesDM puts `K` on `(0, 5)` through a probit and decays
perseverance by `3^K`; this family estimates `k` on `(0, Inf)` with a
log link and decays by the same `3^k`. Same number, same scale, wider
support.

**The first release left this family out because it is "known to be
weakly identified", and it ships with a recovery table that says the
claim is half right.** 30 subjects by 100 trials, one random intercept
on `Arew` with sd 0.4, 60 replicates, everything on the natural scale.

| parameter | truth | bias | mc se | sd of est | coverage | fits |
|---|---|---|---|---|---|---|
| `Arew` | 0.3 | -0.008 | 0.004 | 0.032 | 0.97 | 60 |
| `Apun` | 0.1 | 0.000 | 0.001 | 0.007 | 0.98 | 60 |
| `k` | 0.5 | -0.004 | 0.014 | 0.107 | 0.92 | 60 |
| `betaF` | 1.0 | -0.015 | 0.013 | 0.097 | 0.95 | 60 |
| `betaP` | 1.0 | +0.012 | 0.010 | 0.075 | 0.97 | 60 |

Every parameter recovers. Biases are at or inside their Monte Carlo
error, coverage runs 0.92 to 0.98, and 60 of 60 replicates produced a
usable interval. So a point estimate from this family is not the hazard
that a variance component from a twenty-trial binary session is, and
"weakly identified" as a blanket statement about the family is wrong.

What IS weak is the separation of the three choice-rule parameters, and
it shows up in the correlation of the estimates rather than in their
bias, which is why a bias table alone would have missed it:

| | Arew | Apun | k | betaF | betaP |
|---|---|---|---|---|---|
| Arew | 1.00 | 0.38 | -0.21 | -0.21 | -0.04 |
| Apun | 0.38 | 1.00 | -0.15 | -0.20 | -0.01 |
| k | -0.21 | -0.15 | 1.00 | 0.45 | -0.53 |
| betaF | -0.21 | -0.20 | 0.45 | 1.00 | -0.77 |
| betaP | -0.04 | -0.01 | -0.53 | -0.77 | 1.00 |

`betaF` against `betaP` at -0.77 is the pair. `k` correlates with both,
-0.53 and +0.45, which is what a decay rate on perseverance should do.
The two learning rates are comparatively clean and neither correlates
with any of the three above 0.21.

The practical statement in `?igt_orl` follows the measurement: a study
comparing groups on `Arew` or `Apun` is on solid ground, and a study
comparing them on `betaF` alone risks attributing to outcome frequency
what belongs to perseverance. Report the two together or test them
jointly. Three of the 60 fits warned about a maximum gradient of about
0.001, which is small and is recorded rather than swept up.

I had written "`k` and `betaP` do not separate well from each other"
into the help before the study ran, on the reasoning that both act on
the tendency to repeat a deck. Measured, that pair is -0.53 and the
strongest pair is `betaF` with `betaP`. The help now says the measured
thing.

### The Kalman filter's exploration bonus

`bandit4arm2_kalman_filter(bonus = TRUE)`. Daw and others' softmax with
an exploration bonus: an arm's utility becomes
`tau * (mu + phi * sqrt(s))`, so an arm the subject is uncertain about
is worth more than its posterior mean by an amount the fit estimates.
`sqrt()` of a posterior variance is safe on the tape because the
variance cannot reach zero: the diffusion adds `sigmaD^2` at every
trial.

IT IS OFF BY DEFAULT, and that is what keeps every 0.1.0 fit
bit-identical. `test-reference.R` pins it the strict way rather than by
inspection: `bonus = TRUE` with `phi` held at zero through
`bf(phi = 0)` gives the same log-likelihood as `bonus = FALSE`, so the
two objectives are the same function and the family gained a parameter
rather than a different model.

RECOVERY, 30 subjects by 100 trials, 60 replicates, with `center`,
`mu0` and `sigma0` HELD at the task's own values through
`bf(name = value)`. The family's help already records that those three
do not separate from the rest on a session this size; leaving them free
would have measured their collapse a second time and put its noise into
`phi`'s column.

| parameter | truth | bias | mc se | sd of est | coverage | fits |
|---|---|---|---|---|---|---|
| `tau` | 0.15 | +0.005 | 0.001 | 0.011 | 0.85 | 60 |
| `lambda` | 0.98 | 0.000 | 0.000 | 0.003 | 1.00 | 60 |
| `sigmaD` | 3.00 | -0.045 | 0.026 | 0.203 | 0.95 | 60 |
| `phi` | 1.50 | +0.059 | 0.018 | 0.141 | 0.93 | 60 |

| | tau | lambda | sigmaD | phi |
|---|---|---|---|---|
| tau | 1.00 | -0.25 | 0.07 | 0.12 |
| lambda | -0.25 | 1.00 | -0.10 | 0.00 |
| sigmaD | 0.07 | -0.10 | 1.00 | -0.81 |
| phi | 0.12 | 0.00 | -0.81 | 1.00 |

**The first release's stated reason for leaving the bonus out is wrong,
and the measurement says so.** It said "a bonus and `tau` are hard to
separate at the data sizes these studies collect". Measured, `phi` and
`tau` correlate at 0.12, which is nothing, and `phi` recovers with a
spread of 0.14 against a truth of 1.5, a tenth of the effect rather than
the same order. What `phi` trades off against is `sigmaD`, at -0.81, and
after the fact that is the obvious pair: `phi` multiplies the posterior
standard deviation and `sigmaD` sets how fast that standard deviation
grows, so the two scale the same term.

`phi`'s bias of +0.059 is about three times its Monte Carlo error, so it
is a small real upward bias rather than noise. `tau`'s interval
undercovers, 0.85 against the nominal 0.95; every other coverage in the
table is at or above nominal. Five of the 60 fits warned about a maximum
gradient between 0.0013 and 0.0038, which is small and is recorded here
rather than swept up.

`?bandit4arm2_kalman_filter` now says the measured thing rather than the
guess it opened with.

## What I refused, and why

* **A `fitted_mean` for `rlddm()`.** The mean of a Wiener first-passage
  time conditional on the boundary reached is a closed form, it exists
  in `frmtmb.eam` as an internal, and the response here is ordered so
  nothing about the model forbids `fitted()`. Deriving it again in this
  package would be a second copy of arithmetic a reader cannot check by
  eye, and asking for a second export was outside the one seam the
  brief allowed. Recorded in the compat row and in `?frmtmb.learn`.

* **`accepts_aterms` on these families.** frmtmb's allow-list landed
  this round and these families do not declare one, so a term none of
  them reads is still carried without effect. Every term that would
  reshape a per-row contribution (`weights`, `cens`, `trunc`, `se`) is
  already refused by name in `check_spec` with a better message than
  the allow-list gives, so what a declaration would add is the narrow
  case of a term the family neither reads nor refuses. It would also
  change the compat table's `mi()` row from `untested` to `refused`
  without anyone having tested `mi()`. Left, and named in
  `?frmtmb.learn`.

* **Deviance residuals.** The magnitude is available now, and the sign
  is not. Seven of the eight families have a nominal response with no
  mean to depart from; the eighth has no exported route to the mean of
  a first-passage time. Declaring `deviance = TRUE` without a
  `fitted_mean` would fall through core's `fitted_mean` branch and be
  refused with a WORSE message than the family's own, which the core
  lane records as a known misleading fall-through.

* **A saturated value for `rlddm()`'s rows.** A row there contributes a
  DENSITY, whose supremum over the parameters at a fixed response time
  is unbounded, so there is no constant to subtract. The seven softmax
  families attach `attr(x, "saturated") <- 0`, which is true of a
  nominal choice: a saturated fit puts probability one on the option
  taken. The distinction is a family argument on `ln_structure()`
  rather than a blanket, precisely so that adding a continuous-response
  family could not silently inherit a wrong constant.

* **`mixture()` over learning strategies.** Wanted, and now one step
  closer: core's `mixture()` combines per-ROW densities, the families
  supply per-sequence ones, and joining the two is a change under
  core's `R/` that this lane has no mandate for.

* **Changing anything in core.** `git status` shows no modification
  under the repository root's `R/`, `tests/`, `inst/`, `man/` or
  `vignettes/`. The root `README.md` gains one edited table row, as the
  first release did. Core's `DESCRIPTION` names no extension in
  `Imports` or `Suggests` and this lane added none.

* **A cross-check against hBayesDM.** Still blocked for the reason
  `dev/learn-findings.md` records at length: hBayesDM 2.x needs
  cmdstanr and CmdStan 2.39's vendored TBB does not compile under this
  machine's GCC 14, and 1.2.1's Stan programs use pre-2.33 array syntax
  that this rstan rejects. Nothing this round changes that, and
  `igt_orl()` and `rlddm()` are checked the same way the other six are:
  a Stan program of my own, a longhand reference in the opposite
  spelling, and parameter recovery.

## The defect R CMD check found that the suite did not

Worth its own section, because it is the one real bug this round shipped
and then fixed, and because the reason the ungated suite missed it is
structural.

`rlddm()`'s simulator drew through `frmtmb.eam::ddm_simulate()`, which
uses `RWiener` when it is installed and falls back to its own rejection
sampler when it is not. Run from sources with the user library on the
path, `RWiener` is there and everything passed. `R CMD check --as-cran`
sets `_R_CHECK_SUGGESTS_ONLY_`, which makes only the package's OWN
declared dependencies visible, and `RWiener` is a Suggest of
`frmtmb.eam` and was not one of `frmtmb.learn`. So the check ran the
fallback, and the fallback errored:

    wiener: the fallback simulator could not produce a draw at the
    requested boundary for 1 of 6 rows in 50 passes.

That is not an artifact of the check environment; it is what a user
without `RWiener` would get. The fallback draws the boundary from its
marginal and then rejects forward paths until one reaches THAT boundary,
which fails outright on a row where one boundary is strongly favored,
and a learning model's whole point is to make one boundary strongly
favored. Three tests failed on it.

Fixed three ways rather than papered over:

* `RWiener` joins `frmtmb.learn`'s Suggests, which is the honest
  declaration: the draw needs it.
* `frm_task_simulate()` checks it ONCE per call, through a new
  `sim_needs` constant on the family, and refuses by name if it is
  missing. An up-front refusal naming the package beats a draw that
  works on most rows and errors on the rest.
* `?rlddm` gains a section saying that drawing needs `RWiener` and
  FITTING does not: the likelihood is `wiener_lpdf()` and is exact
  either way.

The tests that draw from `rlddm()` now `skip_if_not_installed()`, which
is what makes the suite tell the truth about a machine without it.

### And a second, in the same run

`test-families.R`'s "the value store reaches the drift rate" assertion
compared the POOLED choice share at two drift scalings. It passed from
sources and failed under the fallback simulator, and it deserved to
fail: the statistic is wrong. A learner that locks onto its
first-rewarded arm locks onto a DIFFERENT arm in different subjects, so
the pooled share stays near 0.5 however decisive each subject is.
Measured over five seeds at 5 subjects by 30 trials:

| drift | pooled statistic | per-subject statistic |
|---|---|---|
| 0 | 0.013 to 0.040 | 0.047 to 0.087 |
| 3 | 0.020 to 0.313 | 0.233 to 0.313 |
| 12 | 0.027 to 0.340 | 0.320 to 0.480 |

The pooled column does not separate drift 3 from drift 0 at all; the
per-subject column separates all three cleanly. The test now asserts two
anchored claims on the per-subject statistic (below 0.15 at drift 0,
above 0.25 at drift 12) instead of comparing two random draws to each
other.

## Verification

SUPERSEDED BY THE PUNCH ROUND. The counts below are the pre-punch run;
the current ones are in "Punch-round verification" at the end of this
file. They are kept because the punch round changed the numbers and not
the conclusions, and because the `Status: OK` row is the one the review
corrected.

| check | result |
|---|---|
| `frmtmb.learn` suite, one process, ungated | **253 passing, 0 failing, 9 skipped**, 60 blocks, 7 files |
| the 9 skips | the whole of `test-stan-identity.R`, gated on `FRMTMB_BRMS_FIT_TESTS` and `NOT_CRAN` |
| gated Stan identity tier | **55 assertions, 0 failing**, 11 identity rows |
| `frmtmb.eam` suite, one process | **1359 passing, 0 failing, 0 skipped**, 18 files; 14 of those are the new `test-extension-api.R` |
| roxygen idempotence, both packages | `man/` and `NAMESPACE` byte identical after a second pass, verified by md5 over every file |
| `R CMD check --as-cran`, `frmtmb.learn`, with the manual | **Status: OK**. 0 errors, 0 warnings, 0 notes. Tests 115 s OK, examples OK, vignette rebuild OK, `checking PDF version of manual ... OK` and `checking HTML version of manual ... OK` |
| `R CMD check --as-cran`, `frmtmb.eam`, with the manual | **Status: 1 WARNING, 1 NOTE**, both pre-existing and both named below. Tests 18 min OK (including the new `test-extension-api.R`), examples OK, vignette rebuild OK, `checking PDF version of manual ... OK` |
| core `R/`, `tests/`, `inst/`, `man/`, `vignettes/` changed | none, verified by `git status` |
| core `DESCRIPTION` naming an extension | none in `Imports` or `Suggests` |

Suite counts audited by file name against `ls tests/testthat`:

| file | passing | failing | skipped |
|---|---|---|---|
| `test-engine.R` | 19 | 0 | 0 |
| `test-factorization.R` | 27 | 0 | 0 |
| `test-families.R` | 45 | 0 | 0 |
| `test-message-uniqueness.R` | 4 | 0 | 0 |
| `test-reference.R` | 12 | 0 | 0 |
| `test-stan-identity.R` | 0 | 0 | 9 |
| `test-surface.R` | 146 | 0 | 0 |
| **total, ungated** | **253** | **0** | **9** |

Against 0.1.0's 184 / 0 / 6: `test-factorization.R` is new (27),
`test-families.R` rose 36 to 45, `test-reference.R` 6 to 12,
`test-surface.R` 119 to 146, and the gated tier grew from 6 skips to 9.
Run with the tier on, `test-stan-identity.R` contributes 55 and the
total is 308 with 0 skips.

The compatibility table this package registers grew from 154 rows to
**205**: 66 `works`, 106 `refused`, 24 `untested`, 9 `conditional`. The
growth is two more families at 25 shared rows each plus their own, and
the one status that CHANGED for every family is `importance`, from
`refused` to `works`.


### The two frmtmb.eam check findings, named with their measured cause

Neither is this lane's, and the evidence is that the FIRST eam check ran
a tarball that did not yet contain `test-extension-api.R` and produced
the identical pair.

* **WARNING, `checking for unstated dependencies in 'tests'`.** The
  printed content of the warning is three Bioconductor index downloads
  failing and nothing else; no package is named. The names it would
  have printed are `frmtmb.sample` and `tmbstan`, used by
  `extensions/frmtmb.eam/tests/testthat/test-sampling.R:13-21` and
  declared in neither that package's `Imports` nor its `Suggests`. Both
  ARE installed on this machine, but `--as-cran` sets
  `_R_CHECK_SUGGESTS_ONLY_`, so the check tries to resolve them against
  the configured repositories to report them, and the machine's
  Bioconductor mirrors are unreachable. The step is then labelled
  WARNING for the download failure rather than for the dependency.

  Not fixed. The fix is two entries in `frmtmb.eam`'s `Suggests`, and
  one of them is `frmtmb.sample`, which makes an extension Suggest
  another extension. That is a monorepo dependency-shape decision that
  belongs to whoever owns that package, not to a lane that was allowed
  one exported function there. `frmtmb.learn`'s own check is `Status:
  OK`, so nothing here is a pattern this lane introduced.

* **NOTE, `checking HTML version of manual`:** "Skipping checking math
  rendering: package 'V8' unavailable". An absent optional package;
  `V8` is the only one of `frmtmb.eam`'s dependency closure not
  installed here. The core lane recorded the same NOTE for the same
  reason. `frmtmb.learn`'s manual carries no math, so its HTML step is
  OK.

## Scope

`git diff --numstat` against the branch point plus untracked files. 31
files changed, 1692 insertions, 283 deletions, and 11 files added.

### frmtmb.eam, and this is the part a sibling lane should read

| file | change |
|---|---|
| `R/extension-api.R` | ADDED. One exported function, `wiener_lpdf()`. |
| `man/wiener_lpdf.Rd` | ADDED, generated. |
| `tests/testthat/test-extension-api.R` | ADDED. 14 assertions. |
| `NAMESPACE` | +1 line, `export(wiener_lpdf)`. |
| `NEWS.md` | +14 lines, a 0.4.1 section. |
| `DESCRIPTION` | one line, `Version: 0.4.0` to `0.4.1`. |

No existing R file in that package was edited. That is deliberate: a
lane working in `wiener-density.R`, `wiener-family.R` or anywhere else
under `frmtmb.eam/R/` will not meet this change except in `NAMESPACE`,
`NEWS.md` and the version line.

### frmtmb.learn

Added: `R/rlddm.R`, `R/igt-orl.R`,
`tests/testthat/test-factorization.R`, `dev/learn-importance.R`,
`dev/learn-recovery-round2.R`, and the two generated help pages.

Changed, R sources: `R/engine.R` (+193 -15, the terms mode, the stack,
the three slots, the density-rule seam), `R/family.R` (+98 -22, the two
slot declarations, the block's `group`, the nominal-aware refusals,
`valid_y` and `family_finalize` pass-through), `R/zzz.R` (+126 -50, the
compat rows), `R/frmtmb.learn-package.R` (+155 -29, docs),
`R/bandit4arm2-kalman.R` (+75 -20, the bonus), `R/trace.R` (+33 -11,
the eight-family table and the `dens` column), `R/task.R` (+30 -3, the
continuous response and the `sim_needs` guard),
`R/bandit2arm-delta.R` (+13 -3, the shared Laplace section).

Changed, tests: `helper-stan-programs.R` (+99, three programs),
`test-stan-identity.R` (+100, three rows), `test-reference.R` (+114,
two longhand references and the bonus identity), `test-families.R`
(+112 -1), `test-surface.R` (+31 -19).

Changed, other: `DESCRIPTION` (version, the `frmtmb.eam` import,
`RWiener` in Suggests, the Description text), `NAMESPACE` (+2),
`NEWS.md` (+136), `_pkgdown.yml` (the two new families in the reference
index), `vignettes/learning.Rmd` (+40 -25).

### Repository root

`README.md`, one edited table row naming the two new families.
`dev/learn2-findings.md`, this file. Nothing else: `R/`, `tests/`,
`inst/`, `man/` and `vignettes/` under the root are untouched, and
core's `DESCRIPTION` names no extension in `Imports` or `Suggests`.

### Left in the worktree

No `dev/*.log`, no scratch scripts. Every script this lane wrote either
lives in `extensions/frmtmb.learn/dev/` as a deliverable
(`learn-importance.R`, `learn-recovery-round2.R`, both with headers
saying what they answer) or under the scratchpad, outside the worktree.
The private library, the Stan cache copy and every check directory are
under the scratchpad prefix as well.

## For the coordinator

* **`frmtmb.eam` gained an export and a version bump.** Say it twice
  because a sibling lane may be in that package: 0.4.0 to 0.4.1,
  `wiener_lpdf()`, in a new file. See the scope table above.
* **`frmtmb.learn` now Imports `frmtmb.eam`.** This is the first
  extension-to-extension dependency in the monorepo. Anything that
  installs `frmtmb.learn` must install `frmtmb.eam` first, which
  includes `dev/build-docs.R` and the CI workflow. The learn workflow
  file was not edited by this lane and will need that line; it is named
  here rather than changed, because the workflow was in another lane's
  hands last round.
* **`frmtmb.learn` Suggests `RWiener`**, for `rlddm()`'s simulator only.
  The CI job's dependency step skips `rstan` and `hBayesDM` by name;
  `RWiener` should NOT be skipped, or the rlddm tests skip themselves
  and the family ships with its simulator unexercised in CI.
* The counts in `dev/build-docs.R` ("installed five packages", "the
  four extension subsites") were already wrong before the first release
  and are still wrong. Not touched.

# Punch round, 2026-09-08

Review: `dev/reviews/2026-09-08-learn2.md`, verdict PUNCH, seven
findings, two blocking. All seven addressed, none disputed. The
reviewer's own verification of the factorization work went further than
this lane's and is quoted where it does.

## Blocker 1. `rlddm()` failed from a clean session

REPRODUCED before it was fixed, on this worktree against the correct
`frmtmb.eam` 0.4.1:

    eam namespace loaded after library(frmtmb.learn)? FALSE
    registered aterms: cens mi payoff reward se stage2 trials trunc
                       vint vreal weights
    clean session, user's own data: FAILS: Addition term `dec()` is not
    supported ...

The cause is exactly as the review states. `frmtmb.eam` was named in
DESCRIPTION `Imports:` and nowhere in `NAMESPACE`, and a DESCRIPTION
entry alone only requires a package to be INSTALLED. `dec()` reaches
frmtmb's registry from `frmtmb.eam`'s `.onLoad()`, and
`frmtmb.eam::wiener_lpdf()` loads that namespace when the function is
CALLED, which is after `frm()` has parsed the formula and refused the
term. The false premise was written into
`extensions/frmtmb.learn/R/zzz.R` in a comment crediting the DESCRIPTION
entry.

WHY NO TEST CAUGHT IT, and this is the part worth keeping. Every route
in the package that fits an `rlddm()` model called `frm_task_simulate()`
first, which reaches `ddm_simulate()` and loads the namespace as a side
effect. The help page's example, `test-factorization.R`'s rlddm block
and all four `rlddm` blocks in `test-families.R` do it. The one thing no
test did is the one thing a user does: fit data they already have.

FIXED, not worked around:

* `extensions/frmtmb.learn/R/frmtmb.learn-package.R` gains
  `@importFrom frmtmb.eam wiener_lpdf ddm_simulate`, so `NAMESPACE`
  carries a real import and R loads the namespace when this package is
  loaded. R loads a package's imports BEFORE running its own
  `.onLoad()`, so `dec()` is in the registry before this package's
  compat rules name it.
* `R/rlddm.R` calls both functions unqualified now, so the import is
  load-bearing rather than decorative.
* The comment in `R/zzz.R` says what actually registers the term, and
  says that an earlier version credited the DESCRIPTION entry and that
  `rlddm()` was unusable from a clean session for exactly that reason.

VERIFIED after the fix, same script:

    eam namespace loaded after library(frmtmb.learn)? TRUE
    registered aterms: cens dec mi payoff reward se stage2 trials trunc
                       vint vreal weights
    clean session, user's own data: works

AND PINNED, by `tests/testthat/test-clean-session.R`, in its own file so
that no earlier block can load the namespace for it. It cannot be
written in-process: by the time testthat runs, this package and its
imports are loaded, which is the state the bug hides in. So it starts a
fresh R process, builds a data frame BY HAND with no
`frm_task_simulate()` anywhere, and asserts three things the broken
version fails: the eam namespace is loaded after `library()`, `dec` is
in the aterm registry, and the fit returns a finite log-likelihood.

Two mechanical notes for whoever reads that test. `system2(env =)` is
documented as unsupported on Windows and there produced no output at
all, so the child is told the library paths by a `.libPaths()` line
written into the script instead. And the child inherits nothing else, so
it is a genuinely clean session rather than a fork.

## Blocker 2. The `frmtmb.eam` version floor

`extensions/frmtmb.learn/DESCRIPTION:39` said `frmtmb.eam (>= 0.4.0)`
while `wiener_lpdf()` is new in 0.4.1, which this package's own NEWS
already said. Raised to `frmtmb.eam (>= 0.4.1)`. The core bound had been
raised in the same diff and this one was missed.

## 3. The package doc denied a capability this lane delivered

The page said the slots "do NOT deliver `loo()` or `waic()`". That is
true of a `frmtmb_fit`, which has no draws, and false of the route
core's own `?loo` sends users to. VERIFIED here rather than accepted:

    admission: ADMITTED
    length = 6   rows in data = 180   subjects = 6
    attr unit = one subject's trial sequence
    sum = -63.4328142046   logLik = -63.4328142046

`frmtmb.sample/R/loo.R` admits any structure declaring `loglik` with
either factorization slot and reads `loglik_group %||% loglik_row`, so a
column is a SUBJECT, not a trial, and the matrix carries
`attr(x, "unit")`. The page now says that, says what a column is, and
says why a trial could not be one: every later trial's value store
depends on it.

## 4. The row-scaling argument, withdrawn and replaced

See the ITEM 2 identity section above, which now opens with WITHDRAWN
and gives the reason. The short version: 240 rows had a residual five
times smaller than 120 rows, so the total did not grow with the row
count, and the per-row column spanned a factor of 9.5, so it was not
flat. Four points of a signed difference between two totals cannot
separate a per-row bias from cancellation.

Replaced with the measurement the reviewer proposed, run here
independently: over 4000 rows spanning the region the fixture visits,
`frmtmb.eam`'s density against `RWiener`'s has a median absolute
difference of 0, a mean of 3.38e-16 and a 99th percentile of 3.55e-15 in
the log density. So eam's side is right to about 1e-15 per row against a
third implementation, and the 4e-12 per-row difference against Stan sits
on Stan's side. That is a stronger claim than the symmetric one it
replaces. The test comment and NEWS carry the same correction.

## 5. The core capped-correction warning, measured and referred upward

NOT this lane's code and not fixable from here: `R/fit.R:2055-2064`
under the repository root. The warning says

> The importance correction used all 5 of its rounds and the estimates
> were still moving by 0.364 at the last one. Raise
> frmtmb_control(importance_rounds =), or raise the draw count so each
> round lands in the same place

Following the first half of that advice makes the number worse in exact
proportion. Measured here, seed 4302 at 40 subjects by 20 trials, whose
Laplace `sd(alpha)` is 1.008e-04:

| rounds asked | rounds used | capped | move | log shift | rounds x move |
|---|---|---|---|---|---|
| 5 | 5 | TRUE | 0.3645 | 1.822 | 1.822 |
| 10 | 10 | TRUE | 0.3645 | 3.645 | 3.645 |

The shift is the step cap times the round count and nothing else, so
more rounds buy a bigger artifact rather than a sharper estimate. FOR
THE COORDINATOR: the warning should say that a capped run's move is a
cap rather than a residual, and should not advise raising the round
count when every entry of `moves` is equal. This package cannot fix it
and now documents around it: `?frmtmb.learn` and
`vignettes/learning.Rmd` both say plainly not to follow that advice on a
capped run, with these two numbers.

## 6. The `R CMD check` claim, restated

The lane reported `Status: OK` for `frmtmb.learn`. That is what this
machine produces and it is not the whole truth, because the run sets
`_R_CHECK_CRAN_INCOMING_=false` (there is no reliable network here; the
Bioconductor mirrors this machine is configured with are unreachable and
the eam check's one WARNING is that failure). With incoming feasibility
ON, the reviewer measures 1 WARNING and 1 NOTE:

* WARNING, `checking CRAN incoming feasibility`: new submission, strong
  dependencies not on CRAN, and a 301 on a DESCRIPTION URL missing a
  trailing slash. The reviewer built and checked `main`'s 0.1.0 the same
  way and got the identical WARNING with the same items, so it is
  pre-existing and not a regression.
* NOTE, `checking examples`, `bandit2arm_delta` at 4.64 user seconds.
  The example code is not touched by this diff, and the reviewer timed
  the recursion change directly against 0.1.0 installed beside it: tape
  build 0.36 s either way, gradient 85 ms against 100 ms, objective
  identical to fifteen digits. Contention from five lanes checking at
  once, not a regression, and it independently confirms this lane's own
  "cost nothing" measurement.

The honest statement is therefore: **`Status: OK` with
`_R_CHECK_CRAN_INCOMING_=false`, and 1 WARNING plus 1 NOTE with it on,
both environmental and the WARNING identical on main.** The table below
says which run each row is.

## 7. The collapsed dataset in the 100-trial cell

The per-dataset table in this file always carried it, seed 4301 at 100
trials with a Laplace `sd(alpha)` of 0.0004 and the same spurious
+1.822. The SUMMARIES elided it, and that mattered because it makes the
advice sound like a short-session caveat when it is not. `?frmtmb.learn`
now carries a five-row table naming all six datasets at each trial count
by what they are, and both it and `NEWS.md` say that a collapsed
component is not a short-session problem and that `sqrt(VarCorr(fit))`
is worth checking at every trial count.

## 8. The observation about `sqrt()` on a posterior variance

No action, and the reviewer asked for none. Recorded because the
reasoning in `R/bandit4arm2-kalman.R`'s comment is right from trial two
and rests on `sigma0^2` at trial one, which is log-linked and therefore
positive. An optimizer would have to reach `exp(-380)` to make `sqrt()`
differentiate badly there, and no fit in the recovery study or the
identity tier went near it.

## What the review verified that this lane had not

Kept here because it is stronger than what this lane claimed, and
because the claims are now load-bearing in the docs.

* **Tape Jacobians.** One `RTMB::MakeTape()` per family returning
  `c(loglik, sum(loglik_row), sum(loglik_group))`, compared at four
  parameter vectors for five families: values agree to 3.4e-13 or
  better and **Jacobians agree to exactly 0** in all 20 comparisons,
  with the gradient itself between 0.46 and 3.0 so the check is not
  vacuous. Exact gradient equality is what "one recursion, three slots"
  is supposed to mean.
* **Padding and level order**, which sum-to-total cannot check. With
  unequal trial counts (11, 20, 24, 20, 18, 30) every row matches an
  independent longhand reference to 4.4e-16 and no row is left at zero;
  with `levels(id)` reversed and the rows shuffled, rows still match to
  8.9e-16 and `loglik_group` follows `levels(group)` rather than order
  of appearance, the two orderings differing by 12.8 log units so the
  check is not vacuous. `importance = 24` completes on both layouts,
  which means core's `imp_verify()` passed on the permuted one.

Both were gaps in this lane's `test-factorization.R`, whose fixtures all
have equal trial counts and natural level order, and both are now
CLOSED IN THE SUITE rather than only recorded. Two blocks were added,
each comparing rows against a longhand delta-rule reference written one
subject and one trial at a time in the same file:

* a PADDED block, trial counts 11, 20, 24, 20, 18, 30, which is the case
  `ln_loglik_row()`'s `keep` mask exists for. It also asserts that no
  row is left at exactly zero, which is what a scatter that wrote the
  pads would produce, since a pad repeats its subject's first row
  number.
* a PERMUTED layout, `levels(id)` reversed and the rows shuffled, so
  order of appearance and level order disagree. It asserts that
  `loglik_group` follows `levels(group)`, and that the two orderings
  really differ, so the assertion is not vacuous.

`test-factorization.R` goes from 27 assertions to 34.

## Punch-round verification

Every suite re-run one process per package on the final sources, counts
audited by file name.

| file | passing | failing | skipped |
|---|---|---|---|
| `test-clean-session.R` | 4 | 0 | 0 |
| `test-engine.R` | 19 | 0 | 0 |
| `test-factorization.R` | 34 | 0 | 0 |
| `test-families.R` | 45 | 0 | 0 |
| `test-message-uniqueness.R` | 4 | 0 | 0 |
| `test-reference.R` | 12 | 0 | 0 |
| `test-stan-identity.R` | 0 | 0 | 9 |
| `test-surface.R` | 146 | 0 | 0 |
| **total, ungated** | **264** | **0** | **9** |

63 blocks over 8 files. Against the pre-punch 253 / 0 / 9 over 7 files:
`test-clean-session.R` is new (4) and `test-factorization.R` rose 27 to
34, which is the padded and permuted coverage the review found missing.

| check | result |
|---|---|
| gated Stan identity tier | **55 passing, 0 failing, 0 skipped**, 9 blocks, 11 identity rows, residuals unchanged to the digit |
| `frmtmb.eam` suite, one process | **1359 passing, 0 failing, 0 skipped**, 18 files; 14 of them the new `test-extension-api.R` |
| roxygen idempotence, both packages | `man/` and `NAMESPACE` byte identical after a second pass, md5 over every file |
| `R CMD check --as-cran`, `frmtmb.learn` | **Status: OK** with `_R_CHECK_CRAN_INCOMING_=false`; tests 112 s OK (the clean-session child runs inside the check and passes there), vignette rebuild OK, PDF manual OK, HTML manual OK. See finding 6 above for what incoming feasibility adds. |
| `R CMD check --as-cran`, `frmtmb.eam` | **Status: 1 NOTE**, and the NOTE is `checking HTML version of manual ... Skipping checking math rendering: package 'V8' unavailable`. See below for the WARNING that appeared on two earlier runs and not this one. |
| core `R/`, `tests/`, `man/`, `vignettes/`, `DESCRIPTION` | untouched. The only file changed outside the two extension directories is the repository-root `README.md`, one line; it is `.Rbuildignore`d at the root, so it is a monorepo document rather than part of the core build. The review asked for that distinction and this is it. |

The gated identity rows are unchanged by the punch round, which is the
point: the fixes were to loading and to prose, not to arithmetic.

    bandit4arm2_kalman_filter (bonus)  -441.0315488439  diff -2.842e-13
    igt_orl                            -935.0557464773  diff -6.821e-13
    rlddm                              -224.7348379010  diff  3.176e-09

### The version floor, verified by accident and then on purpose

While the floor was raised and `frmtmb.eam` 0.4.0 was still the
installed version, roxygen2 refused to load the package with

    The package "frmtmb.eam" (>= 0.4.1) is required.

which is the bound doing exactly what the review asked for: a
`frmtmb.learn` 0.2.0 sitting next to an eam 0.4.0 is now refused at
load rather than failing later at `frm()`. Installing eam 0.4.1 cleared
it.

### The eam check WARNING is non-deterministic, which settles what it was

Two earlier runs of this check gave `1 WARNING, 1 NOTE`; this one gives
`1 NOTE`, with `checking for unstated dependencies in 'tests' ... OK`.
The package did not change between them. The WARNING's printed content
was three Bioconductor index downloads failing and no package name, so
what it reports is a repository lookup that succeeds or fails with the
network rather than anything about `frmtmb.eam`. The names it would have
printed are `frmtmb.sample` and `tmbstan`, used by that package's own
pre-existing `test-sampling.R` and declared in neither its `Imports` nor
its `Suggests`; both are installed here, and `--as-cran` hides them
because they are undeclared, which is what sends the check to a
repository.

The coordinator reports the same pair on `main` verbatim, with `main`
carrying one warning and TWO notes against this worktree's one and one,
so nothing this lane did added a finding to that package and one of
main's notes is absent here.

### And the learn Status claim, corrected

The verification table above says `Status: OK` and names the flag that
makes it so. That is this machine's answer with
`_R_CHECK_CRAN_INCOMING_=false`, and it is not the whole truth: with
incoming feasibility on, the reviewer measures 1 WARNING and 1 NOTE for
`frmtmb.learn` as well, the WARNING identical on `main`'s 0.1.0 and the
NOTE an example timing that tracks machine load. Finding 6 above carries
both. The honest one-line summary of the check is therefore "clean on
everything this environment can test, with two environmental findings
this environment cannot", not "Status: OK".


## Three things now confirmed rather than claimed

The coordinator asked for these to be stated as settled, and they are
settled by measurements this lane did not make.

* **The three slots are the SAME TAPE, not merely equal numbers.**
  Jacobians of `loglik`, `sum(loglik_row)` and `sum(loglik_group)`
  differ by exactly zero across five families at four parameter
  vectors, with the gradient itself between 0.46 and 3.0 so the check
  is not vacuous. Rows match a longhand reference under padding and
  under a reversed level order with shuffled rows, and group values
  follow `levels(group)` rather than row order, which a deliberate test
  separated by 12.8 log units. The padded and permuted cases are now in
  this package's own `test-factorization.R` as well.
* **The `frmtmb.eam` footprint is exactly as this file describes it.**
  Three new files plus one line each in `NAMESPACE`, `NEWS.md` and the
  version, and no existing R file in that package touched.
* **Neither `frmtmb.eam` check finding is this lane's.** Both appear on
  `main` verbatim, and `main` carries one warning and two notes against
  this worktree's one and one, so the worktree is if anything cleaner.
  The example-timing class of note tracks machine load rather than
  code, which the reviewer confirmed by timing 0.2.0 against 0.1.0 side
  by side: tape build 0.36 s either way, gradient 85 ms against 100 ms,
  objective identical to fifteen digits. That is an independent
  confirmation of this lane's own "the shape change cost nothing"
  measurement, made with a different method on a different machine
  state.
