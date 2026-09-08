# remlopt: what `reml-ar1-se-two-optima` actually was

Lane `remlopt`, round 0.55, worktree `frmtmb-wt-remlopt` on `a2deb65`.

## Answer first

It was not two optima, and it was not a REML, `se()` or `ar1` problem.
It was one unidentified parameter sitting on a likelihood plateau that
the AD student-t density could no longer represent. `RTMB::dt()` has
two branches. The double branch is `stats::dt()` and is accurate. The
AD branch, which is the only branch a fit ever uses, forms
`lgamma((nu + 1) / 2) - lgamma(nu / 2)` as written, and those two
values agree in every leading digit once `nu` is large. The objective
past about `nu = 1e7` was noise, its gradient changed sign where the
true likelihood is monotone, and the optimizer stopped at whichever
sign change the row order sent it to.

The density is now formed without the cancellation. The two row orders
now return a BITWISE identical log likelihood, and parameter vectors
agreeing to 1.8e-11 with the whole residue in `nu`, the parameter that
has no maximum. The tier is green with the pending entry removed.

## 1. Is either point stationary? No. Both are stopped optimizers.

The gradient at both points was already small in the units the fit
judges itself in (8.7e-06 and 2.1e-05 against `grad_tol = 1e-3`), so
the gradient alone says nothing. The parameter vector does.

| | fit 1 (as given) | fit 2 (rows permuted) |
| --- | --- | --- |
| objective | 143.1657474169 | 143.1608851298 |
| `betad` = `log(nu - 1)` | 18.216655 | 24.154751 |
| `theta_1` | 0.083332 | 0.083332 |
| `theta_2` | -0.019303 | -0.019302 |
| max abs gradient | 8.7e-06 | 2.1e-05 |
| numeric Hessian eigenvals | 131.8, 57.2, **-0.10** | 131.8, 57.2, **-2.76** |

The outer vector has three entries and the whole 5.938 of distance is
in ONE of them. `sum|p1 - p2| = 5.938096`, `max|p1 - p2| = 5.938095`,
and the two `theta` differ in the seventh decimal. The Hessian is not
positive definite at either point: the third eigenvalue, the `betad`
direction, is negative at both, which is a numerical Hessian of a
direction the objective does not resolve.

Profiling the objective along `betad` with the two `theta` held at fit
1's values settles it. The values below are `obj$fn` on fit 1's own
tape:

| `betad` | `nu` | objective | `d/d betad` |
| --- | --- | --- | --- |
| 10 | 2.2e+04 | 143.1680847363 | -2.33e-03 |
| 14 | 1.2e+06 | 143.1658018133 | -4.25e-05 |
| 17 | 2.4e+07 | 143.1657589999 | -1.00e-06 |
| 18 | 6.6e+07 | 143.1657582442 | -3.68e-07 |
| 23 | 9.7e+09 | 143.1649964612 | -2.48e-09 |
| 24 | 2.6e+10 | 143.1718464035 | -6.02e-03 |
| 28 | 1.4e+12 | 143.4301802495 | -1.67e-11 |
| 30 | 1.1e+13 | 142.0235516200 | -2.27e-12 |
| 31 | 2.9e+13 | 144.8386286059 | -8.51e-13 |
| 32 | 7.9e+13 | 127.9048333421 | -1.80e+01 |

A likelihood falling smoothly to a plateau and then jumping by 17 log
units at one step of `betad` is not a likelihood. Neither reported
point is a stationary point of the real one.

## 2. Could a user get a different scientific answer? No.

`betad` is `nu_(Intercept)` on the `logm1` link, so the two "optima"
are `nu = 8.15e+07` and `nu = 3.09e+10`. Both mean "gaussian". Nothing
a user reports moved:

| | fit 1 | fit 2 |
| --- | --- | --- |
| `mu` (Intercept) | 1.07196 (se 0.12101) | 1.07196 (se 0.12101) |
| `mu` x | 0.37714 (se 0.14000) | 0.37714 (se 0.14000) |
| z, p for x | 2.6940, 0.007061 | 2.6940, 0.007061 |
| `nu` (Intercept) | 18.217 (se 1837.6) | 24.155 (se 35897.4) |
| random-effect sd | 1.0869 | 1.0869 |

`max|d mu|` between the two fits was **3.2e-09** before the fix and is
**5.6e-17** after it. The defect was confined to the one parameter
nobody reports. What it cost was reproducibility, not correctness: the
same model on the same data in a different row order gave a different
`logLik`, a different AIC and a different `nu`.

**One correction to my own first draft, which said the two fits land
"bitwise" outright.** They do on the log likelihood: `identical()` on
the two values is `TRUE` and the gap is exactly 0. They do not on the
parameter vector, where `sum|p1 - p2|` is **1.796620e-11**. My probe
printed that with `%.6f`, which renders 1.8e-11 as `0.000000`, and I
read the format as the result. Measured at full precision the residue
is entirely in `betad`: 21.85715401282329 against 21.857154012805324,
a difference of 1.797e-11, while the two `theta` differ by 1.249e-16
and 0 and `mu` by 5.551e-17. Same class of error as the timing below,
and worth stating plainly: a printed zero is not a measured zero.

## 3. Which ingredient? `student`, and only `student`.

16 cells, 5 seeds each, 2 fits per seed (160 fits). `over` counts the
seeds whose `logLik` gap under a row permutation exceeded
`fuzz_permutation_tol()`. `max_nu` is the larger of the two fits'
degrees of freedom.

| family | aterm | re | mode | over/5 was | max gap was | max `nu` | now |
| --- | --- | --- | --- | --- | --- | --- | --- |
| student | se | ar1 | reml | **2** | 3.98e+00 | 5.8e+13 | 0 |
| student | none | ar1 | reml | 0 | 8.58e-02 | 6.9e+11 | 0 |
| student | se | ar1 | ml | **1** | 7.36e-05 | 6.1e+08 | 0 |
| student | none | ar1 | ml | **2** | 1.25e+02 | 4.7e+36 | 0 |
| student | se | ri | reml | 0 | 7.14e-01 | 1.3 | 0 |
| student | none | ri | reml | 0 | 8.53e-14 | 19.1 | 0 |
| student | se | ri | ml | 0 | 5.71e-02 | 1.5 | 0 |
| student | none | ri | ml | 0 | 1.42e-13 | 16.6 | 0 |
| gaussian | (all 8 cells) | | | 0 | <= 7.3e-08 | n/a | 0 |

Read the table by `max_nu`, not by the design axes. Every gaussian cell
is at machine precision. Every student cell whose `nu` stayed small
(1.3 to 19.1) is too. The four cells over tolerance are exactly the
four whose `nu` ran past 1e8. `se()`, `ar1` and REML are not
ingredients: they only make the surface flat enough that `nu` escapes.
The worst cell in the whole table is not the recorded spec: `student`,
no aterm, `ar1`, ML reached `nu = 4.66e+36` and a `logLik` gap of
124.9, which the recorded entry's four-axis matcher would have called a
new finding.

## 4. The measurement that names the cause

`stats::dt()` is not the problem. Against a 300-bit `Rmpfr` reference
it holds 1e-13 at every `df` up to 1e50. The AD tape is. Taped
`-sum(dt(z, df = 1 + exp(e), log = TRUE))` on 90 standard normal
points, against the same reference:

| `e` | AD error | AD `d/de` | true `d/de` |
| --- | --- | --- | --- |
| 10 | 2.4e-10 | -5.95e-04 | -5.95e-04 |
| 14 | 5.7e-10 | -1.08e-05 | -1.09e-05 |
| 16 | 1.4e-06 | -1.44e-06 | -1.47e-06 |
| 18 | -1.7e-06 | **+1.70e-06** | -1.99e-07 |
| 20 | 9.3e-06 | **+2.61e-05** | -2.70e-08 |
| 24 | 6.1e-03 | **+2.53e-03** | -4.94e-10 |
| 28 | 2.6e-01 | **+1.33e-01** | -9.07e-12 |
| 32 | -1.5e+01 | -1.96e-01 | -1.59e-13 |

Fit 1 stopped at `e = 18.22` and fit 2 at `e = 24.15`. Those are two of
the sign changes in that column. The AD error at `e = 24` is 6.1e-03
and fit 2's "gain" over fit 1 was 4.86e-03: the same number.

The mechanism is one subtraction. `lgamma(a + 1/2) - lgamma(a)` with
`a = nu/2`, against a 300-bit reference:

| `a` | reference | naive error |
| --- | --- | --- |
| 1e+02 | 2.301335098 | -4.9e-14 |
| 1e+06 | 6.907755154 | 1.1e-09 |
| 1e+10 | 11.512925465 | 1.4e-05 |
| 1e+14 | 16.118095651 | 3.8e-01 |
| 1e+50 | 57.564627325 | -5.8e+01 |

At `a = 1e10` each `lgamma` is about 2.3e11, where a double is spaced
3e-5 apart, and the answer they must produce is 11.5.

## 5. What was changed

### The density (`R/families.R`)

`lgamma_shift_diff(a, s)` computes `log Gamma(a + s) - log Gamma(a)`
without the cancellation, and `dt_stable(z, nu)` uses it.
`fam_student()`'s `lpdf` calls `dt_stable()` in place of `RTMB::dt()`.

Two steps, no branch. Push `a` up 25 times with
`Gamma(x + 1) = x Gamma(x)`, one exact `log1p` per step, until Binet's
remainder series is exact; then take the Stirling difference as
`(b - 1/2) log1p(s/b) + s log(b + s) - s + mu(b + s) - mu(b)`, whose
pieces all stay the size of the answer.

Branchless, but NOT because a branch is impossible. My first version of
this paragraph said RTMB has no `CondExp` so a comparison "would freeze
at whatever value the tape was built with". The first half is true and
the second is false. Re-measured on RTMB 1.9:

| probe | result |
| --- | --- |
| `grep("^CondExp", getNamespaceExports("RTMB"))` | 0 matches |
| comparison, default `TapeConfig` | ERRORS, "generally unsafe for AD types" |
| comparison under `TapeConfig(comparison = "tape")` | TAPES correctly |
| `ifelse` on an advector test | errors under both settings |

Under `comparison = "tape"` I taped `b <- x > 25; b * (x * 2) +
(1 - b) * (x + 100)` and re-evaluated it at 1, 10, 24, 25, 26 and 100:
the value switches between 25 and 26 and the jacobian is 1 on the low
branch and 2 on the high one. Nothing freezes.

Branchless is still right, for two reasons that are true. `TapeConfig`
is PROCESS-GLOBAL and defaults to forbidding comparisons, and a package
must not flip a global setting under its user. And the one conditional
that needs no setting, a multiplicative blend, evaluates both arms.

The blend's threshold was also wrong, in my favour. I wrote that
`lgamma(a)` overflows above `a = 1.8e308`; `1.8e308` is not a
representable double (`is.finite(1.8e308)` is `FALSE`). Bisecting,
`lgamma` is finite up to `a = 2.5327e+305` and `Inf` from there. That
is `nu = 5.0655e+305`, reached at `log(nu - 1) = 703.91`, INSIDE the
709.78 the `logm1` link can produce. So the blend returns `NaN` at a
reachable `nu`, three orders of magnitude sooner than I claimed, where
`dt_stable(0, 5.0655e305)` is -0.918938533204653 against
`dnorm(0, log = TRUE)` of -0.918938533204673. The rejection stands for
a stronger reason than I gave it.

Accuracy against the 300-bit reference, worst absolute error over `a`
in `[0.5, 1e50]`:

| form | `s = 1/2` | `s = 3` |
| --- | --- | --- |
| naive `lgamma` difference | 5.8e+01 | 1.5e-03 |
| sigmoid blend (rejected) | 3.6e-13 | 6.8e-14 |
| shift + Stirling (shipped) | **1.8e-15** | **7.1e-15** |

### The cost. My first measurement was wrong, and here is how

**The claim I made was "it costs nothing, 0.75 to 0.84 times the old
cost". Do not believe it. It was a clock artifact.**

`proc.time()[["elapsed"]]` on this machine ticks at exactly 0.010000 s,
which I measured by polling it until it changed. My benchmark timed
blocks of 200 calls. At the ~200 us per call I reported, a block is
0.040 s, which is FOUR TICKS. The pair I reported, 200.0 us against
150.0 us, is four ticks against three, and 3/4 is exactly the 0.75 I
published. I was reading the quantizer, not the code. The project's
standing instruction to profile rather than guess exists for this, and
a timer below its own resolution is the failure it is meant to prevent.
`Sys.time()` resolves to 1.9e-06 s here and would have been the better
clock; sizing the block to 1.2 s makes either one adequate.

Redone: every call at a DISTINCT parameter vector (`obj$fn` caches its
last evaluation point), both tapes built in one process, visiting order
rotated every round, minimum of 9 rounds, each block grown until it
takes at least 1.2 s. A CONTROL is included: a second tape built from
the SAME density, so `old2/old` and `new2/new` are what this method
reports when the true ratio is 1.

| case | `n` | op | old (us) | new (us) | new/old | old2/old | new2/new |
| --- | --- | --- | --- | --- | --- | --- | --- |
| scalar `nu` | 5000 | `fn` | 226.7 | 188.8 | 0.833 | 0.980 | 0.984 |
| scalar `nu` | 5000 | `gr` | 434.8 | 336.3 | 0.773 | 0.977 | 1.122 |
| scalar `nu` | 50000 | `fn` | 2522.4 | 2114.3 | 0.838 | 1.003 | 0.963 |
| scalar `nu` | 50000 | `gr` | 4315.1 | 3233.1 | 0.749 | 0.964 | 1.119 |
| `nu ~ x` | 5000 | `fn` | 1291.7 | 4801.3 | **3.717** | 1.073 | 1.009 |
| `nu ~ x` | 5000 | `gr` | 3480.0 | 8520.2 | **2.448** | 0.899 | 0.852 |
| `nu ~ x` | 50000 | `fn` | 17595.7 | 62747.1 | **3.566** | 1.007 | 1.063 |
| `nu ~ x` | 50000 | `gr` | 36554.5 | 90436.5 | **2.474** | 1.051 | 0.957 |

The controls put this method's noise floor at 0.90 to 1.12, so both
effects clear it.

**The `nu ~ x` result is the one that matters and it is a real cost.**
With a vector `nu`, `a = nu/2` is a vector of length `n`, so the 25
fixed recurrence steps become 25 VECTOR `log1p` nodes where the old
code had two `lgamma` nodes. `bf(y ~ x, nu ~ z) + student()` fits
today. The reviewer measured the same direction at 2.6 to 5.4 times;
mine is 2.4 to 3.7. **Take the span, 2.4 to 5.4 times.**

**The scalar result does not agree with the reviewer's and I could not
make it.** I measure 0.75 to 0.84 (new faster); the reviewer measured
1.06 to 1.17 (new slower) on the same machine with the same
methodology. The disagreement is larger than either control band. What
I can show is the mechanism, which says the effect is small and has two
parts that nearly cancel. Timing three forms with a scalar `nu`, same
process, interleaved:

| `n` | form | `fn` (us) | vs `RTMB::dt` | ns per observation |
| --- | --- | --- | --- | --- |
| 5000 | `RTMB::dt` | 229.3 | 1.000 | 45.86 |
| 5000 | naive, lgamma pair kept scalar | 187.2 | 0.816 | 37.44 |
| 5000 | `dt_stable` (shipped) | 188.9 | 0.824 | 37.77 |
| 50000 | `RTMB::dt` | 1968.0 | 1.000 | 39.36 |
| 50000 | naive, lgamma pair kept scalar | 1569.8 | 0.798 | 31.40 |
| 50000 | `dt_stable` (shipped) | 1749.4 | 0.889 | 34.99 |

Writing the SAME formula out by hand is already 0.80 to 0.82 times
`RTMB::dt` at a scalar `df`, which points at `RTMB::dt` broadcasting
`df` per element rather than once. The recurrence then adds 1 to 11
percent back. So on the scalar path the algorithm itself is close to
free and the measured total is dominated by not calling `RTMB::dt`.
The honest published range is 0.75 to 1.17, spanning both runs, with
the effect under a quarter either way.

### The same subtraction in two more places

`student_lpdf_core()` in `R/covstruct.R` (the `gr(dist = "student")`
blocks) and the multivariate-t branch of `R/autocor.R` (`student()`
with `ar(cov = TRUE)`) both form `lgamma((nu + d)/2) - lgamma(nu/2)`.
Both now use `lgamma_shift_diff()`, and the `log(1 + qv/nu)` in
`autocor.R` is now `log1p()`, which `covstruct.R`'s own comment eight
lines away already argued for.

`covstruct.R` said "A large `nu` is how the gaussian limit is reached,
and how a user checks that a t block reduces to one." Measured on that
check, `autocor_loglik()` with `K = 4`, `G = 6`, against the gaussian
of the same scale matrix:

| `nu` | gap to the gaussian, BEFORE | AFTER |
| --- | --- | --- |
| 1e+08 | 1.2e-08 | 8.4e-08 |
| 1e+10 | 1.3e-05 | 8.4e-10 |
| 1e+12 | 7.3e-03 | 8.4e-12 |
| 1e+14 | 5.2e-01 | 1.1e-13 |
| 1e+16 | **4.9e+01** | 3.6e-14 |
| 1e+18 | **4.8e+02** | 2.8e-14 |

The AFTER column falls as `1/nu` until it hits the double floor, which
is what the limit is supposed to do.

### `diagnose()` (`R/confint.R`)

The fix makes the objective honest, and an honest objective has NO
maximum in `nu` on data with no heavy tails: the likelihood rises all
the way to the gaussian limit. So the fix makes the reported `nu` more
extreme, not less. In the ingredient table above `student`/none/`ar1`
under ML went from `nu = 4.7e+36` to `nu = 1.4e+238`. That is the
correct answer to a question with no answer, and 0.54.0's `diagnose()`
replied "No convergence problems detected" to it.

`diagnose()` gains `unbounded_dpar`. It names a DISTRIBUTIONAL
parameter whose estimate is far out on its link AND whose standard
error is larger than the estimate itself, prints the natural-scale
value, and says what to do. Primary dpars are left to
`diagnose_separation()`, whose message (separation, collinearity) does
not describe this.

False-alarm rate, measured before shipping it: 84 `student()` fits,
seven error laws (`t(3)`, `t(5)`, `t(10)`, `t(30)`, normal, uniform,
logistic), `n` of 60 and 200, six replicates each.

| group | fits | `log(nu - 1)` | `se(log(nu - 1))` |
| --- | --- | --- | --- |
| `nu` ran off (> 1e4) | 35 | [18.2, 21.6] | [4878, 1.2e+04] |
| `nu` identified | 49 | [0.081, 5.06] | [0.35, 21.8] |

The two groups do not overlap on either axis, with a factor of 3.6
between them on the estimate and 224 on the standard error. The check
requires `|est| > 10`, which excludes all 49 identified fits on its
own. **False alarms: 0 of 49.**

`?frmtmb-families` gains a "Degrees of freedom that run off" section,
which is where a user of `student()` will look, and
`vignette("diagnostics")` gains "A distributional parameter with no
maximum" next to its "NaN standard errors" section, which is where
`diagnose()`'s own message already sends people. The vignette chunk was
run by hand because pandoc is not installed on this machine, so
`R CMD check` cannot build vignettes here; it prints
`nu: (Intercept)  21.53311  7702.791  2247558151`.

## 6. What `frm_allfit()` said. It is not defective.

`frm_allfit()` saw it, loudly, and needs no change.

| | BEFORE | AFTER |
| --- | --- | --- |
| logLik spread (fit 1) | **90.1** | 4.2e-06 |
| worst optimizer (fit 1) | `nloptr_lbfgs` at **-53.05** | -143.1658 |
| logLik spread (fit 2) | 0.00491 | (same fit) |
| max fixed-effect spread | 24.8 | 5.54 |

The -53.05 was itself the cancellation: that optimizer walked further
out into the noise than the others. All four now agree to 4.2e-06. The
5.54 of remaining fixed-effect spread is the unidentified `nu`, and it
is the honest report of a parameter with no maximum, which is why the
`diagnose()` check above is the thing that had to change and
`frm_allfit()` is not.

## 7. What `diagnose()` said before. Silence, and that was the gap.

On BOTH points, 0.54.0's `diagnose()` printed four lines and the
all-clear. Fit 1's, verbatim (fit 2's is the same with a maximum
gradient of 2.056e-05):

```
Optimizer convergence code: 0 (relative convergence (4))
Max |gradient|: 8.733e-06 at theta_1
Hessian positive definite: TRUE
No convergence problems detected
```

`flat` was empty in both, because that check is gated on the covariance
having already failed and this covariance succeeds. `separation` was
`NULL` in both, because it is gated on binomial-type families and on
the primary dpar. So a fit reporting `nu = 8.2e+07` with a standard
error of 1838 got a clean bill of health, twice, at two different
answers. That is the defect item 5 addresses.

## 8. Tests, each seen failing first

| file | test | before | after |
| --- | --- | --- | --- |
| `test-families.R` | density keeps its digits as nu runs off | 6 | 0 |
| `test-families.R` | no false optimum in the large-nu tail | 2 | 0 |
| `test-diagnostics-ux.R` | a dpar with no maximum is named | 2 + 1 err | 0 |
| `test-autocor.R` | dmvnorm/dmvt, plus the gaussian limit | 1 | 0 |
| `test-tre.R` | us_t and diag_t, plus the gaussian limit | 3 | 0 |

Both tolerances in the two new `test-families.R` blocks are ratios to
something the run measures, `64 * .Machine$double.eps * max(abs(val))`,
not numbers chosen in advance. The gradient assertion needed the same
treatment: past `e = 33` the objective is flat to the last bit of a
double, the sign of a 1e-14 gradient is meaningless, and the first
draft of the test demanded more than double precision can give. It now
asserts that no ladder step reports an UPHILL gradient bigger than the
objective's own noise, which the pre-fix gradients (1.7e-06 up to
0.13) exceed by five to ten orders of magnitude.

## 8b. The fuzz tier, full plan, before and after

300 rows, seed `20260901`, `FRMTMB_FUZZ=true`, 312.7 s wall. Both runs
use the CURRENT `helper-fuzz.R`, with `reml-ar1-se-two-optima` already
deleted, so the before run is the tier asking the question with no
entry to hide behind.

| class | BEFORE | AFTER |
| --- | --- | --- |
| **real_new** | **1** | **0** |
| generator | 0 | 0 |
| known_pending | 58 | 58 |
| known_refusal | 3 | 3 |
| known_divergence | 6 | 6 |
| unconverged | 13 | 15 |

The one REAL-NEW before is the recorded spec, to the digit:

```
[REAL-NEW] row_permutation
    spec: family=student aterm=se re=ar1 special=none dpar=none
          mode=reml op=confint (seed 20331245, cover)
    logLik changed under a row permutation
    values: logLik=-143.16575; logLik_permuted=-143.16089;
            diff=0.0048622871; tol=0.00017391926
```

After the fix that spec produces no finding at all. The remaining
`known_pending` ids are `re-factor-call` (36), `double-bar-factor` (15)
and `trunc-postfit` (7). None is a `row_permutation`.

Six findings still carry the `row_permutation` invariant in both runs,
all classed `unconverged`: fits with maximum gradients of 9.0e+12,
1.2e+13, 3.2e+09, 1.7e+11, 25 and 1464, "false convergence (8)", and
NaN objectives. A fit that did not converge cannot be asked to be
permutation invariant, which is what `FUZZ_CONVERGENCE_SENSITIVE`
already says. None of them is the phenomenon fixed here: this one had a
maximum gradient of 8.7e-06.

The set of six is not identical between the runs. The recorded spec
left it, and `family=student aterm=se re=nested special=mo_int
dpar=none mode=sparse_x op=confint` entered it. That one was measured
on both libraries and is broken either way: convergence code 1, "false
convergence (8)", maximum gradient 4.6e+09 BEFORE and 1.7e+11 AFTER,
`nu = 1.52`, so it is not the runaway case. The unconverged count moves
by two for the same reason: an optimum that never arrived lands
somewhere slightly different when the objective changes in its
fifteenth digit.

## 8c. Verification

`roxygen2::roxygenise()` is clean. It emitted two link warnings on the
first pass, for `[4878, 1.2e4]` and `[0.35, 22]` in a `@noRd` block
that roxygen read as Rd link targets; both ranges are reworded and the
warnings are gone. `man/diagnose.Rd` and `man/frmtmb-families.Rd` are
regenerated.

Core test files, one per R process, `NOT_CRAN=true`, against the
private library:

| file | pass | fail | error | skip |
| --- | --- | --- | --- | --- |
| `test-families.R` | 227 | 0 | 0 | 0 |
| `test-diagnostics-ux.R` | 127 | 0 | 0 | 0 |
| `test-autocor.R` | 191 | 0 | 0 | 0 |
| `test-tre.R` | 101 | 0 | 0 | 0 |
| `test-simulate-density.R` | 435 | 0 | 0 | 0 |
| `test-numerical-robustness.R` | 679 | 0 | 0 | 0 |
| `test-custom-family.R` | 99 | 0 | 0 | 0 |
| `test-review-v28.R` | 73 | 0 | 0 | 0 |
| `test-method-residue.R` | 64 | 0 | 0 | 0 |
| `test-input-validation.R` | 43 | 0 | 0 | 0 |
| `test-confint-anova.R` | 35 | 0 | 0 | 0 |
| `test-distributional.R` | 28 | 0 | 0 | 0 |
| `test-bcm-mpt.R` | 15 | 0 | 0 | 2 |
| `test-map.R` | 7 | 0 | 0 | 0 |
| `test-fuzz.R` (`FRMTMB_FUZZ=true`) | 2 | 0 | 0 | 0 |

`R CMD check --as-cran` on the final tree: **the full suite is OK**
(`Running 'testthat.R' [305s] OK`), examples OK, examples with
`--run-donttest` OK. Status is 1 ERROR, 3 WARNINGs, 4 NOTEs, and every
one of them is a missing tool on this machine rather than anything in
the diff:

| finding | cause |
| --- | --- |
| ERROR, PDF manual without index | `pdflatex is not available` |
| WARNING, PDF version of manual | the same `pdflatex` |
| WARNING x2, vignettes, `inst/doc` | `--no-build-vignettes`; no pandoc |
| NOTE, HTML manual | `package 'V8' unavailable`, the expected one |
| NOTE, top-level files | `README.md`/`NEWS.md` need pandoc |
| NOTE, CRAN incoming | new submission, no vignette index |
| NOTE, non-standard things | `frmtmb-manual.tex`, from the failed `pdflatex` |

## 8d. Does the recurrence have to be 25 steps? Asked and answered

**As posed, no shorter recurrence is available: `m` cannot be chosen at
tape time from anything that is not a parameter.** The shift needed
does fall with `a`, but `a` is `nu / 2` and `nu` is a fitted parameter
that the optimizer moves across the whole range AFTER the tape is
built. The only facts fixed at tape time are the link's:

| link on `nu` | what it fixes | lower bound on `a` |
| --- | --- | --- |
| `logm1` | `linkinv = 1 + exp(eta) > 1` | `a > 0.5` |
| `identity` | nothing; `nu` is the predictor | unbounded below |

`logm1` buys exactly ONE step out of 25, and `identity` buys none. A
conditional on the tape-time VALUE of `nu` would be a conditional on a
parameter and is not on the table.

**The real lever is the series length, and it depends on no parameter
at all.** `m` is set by the smallest argument at which the Binet series
is accurate, and that falls fast with the number of terms. Measured
against a 400-bit reference, smallest `x` with absolute error at or
below 1e-16:

| Binet terms | smallest usable `x` | error at `x = 25` |
| --- | --- | --- |
| 4 (shipped) | 28 | 2.20e-16 |
| 5 | 17 | 4.34e-19 |
| 6 | 12 | 4.34e-19 |
| 7 | 10 | 4.34e-19 |
| 8 | 8 | 4.34e-19 |
| 10 | 7 | 4.34e-19 |

Two things follow. First, the shipped `m = 25` with four terms is a
shade inside where four terms are exact by this criterion, which is why
the roxygen on `lgamma_binet()` now states the measured 2.2e-16 rather
than claiming exactness. Second, a longer series would let `m` drop and
buy back most of the vector-path cost. Priced on the same interleaved
benchmark, `nu ~ x` at `n = 5000`, and swept for accuracy over the
review's 37 `a`-points crossed with `s` in 0.5, 1, 1.5, 3, 50:

| variant | vector cost vs 0.54.0 | worst relative error |
| --- | --- | --- |
| `RTMB::dt` (0.54.0) | 1.00 | 2.9e-02 at `a = 1e12` |
| `m = 25, k = 4` (shipped) | 3.52 | 2.938e-14 |
| `m = 12, k = 6` | **2.13** | **1.276e-14** |
| `m = 8, k = 8` | 1.90 | 2.876e-14 |

`m = 12, k = 6` is both cheaper and more accurate than what ships.

**I did not make the change.** It rewrites the constants of a density
that has just been verified at those constants over 210 `(nu, z)`
points, 148 `(a, s)` points and a full derivative sweep, and my
re-running of those sweeps would not be the independent check the
current numbers rest on. It is a clean, parameter-free, well-posed
follow-up: pick `(m, k)` on the table above, extend `BINET` to `k`
terms (the coefficients are `B_{2n} / (2n(2n-1))`, and terms 5 to 8 are
`1/1188`, `-691/360360`, `1/156`, `-3617/122400`), and redo the
accuracy sweep and the derivative sweep before shipping. The scratch
implementation and both sweeps are in `remlopt-shift.R` and
`remlopt-shift2.R` in this round's scratchpad.

## 8e. After review

The review at `dev/reviews/2026-09-08-remlopt.md` found the density
sound (6 ulp over 210 `(nu, z)` points from `nu = 0.1` to `1e300`, 3.3
ulp at the moderate `nu` where the old code was already fine, correct
gradient sign everywhere the old one flipped) and named three things to
settle. What changed here afterwards, none of it functional:

* **F1**, the message read backwards at the lower end. Fixed by the
  reviewer. `?frmtmb-families`, `vignette("diagnostics")` and NEWS
  described only the upward case; all three now name both directions.
* **F2**, the performance claim. Reproduced with a control and
  rewritten; section 5 above.
* **F3**, "a comparison would freeze". False. I verified the taped
  comparison myself and rewrote the reason.
* **F4**, `a = 1.8e308`. Not a representable double. Bisected: the
  real threshold is 2.533e305, and it is reachable.
* **F6**, the "bitwise" overclaim. Verified: the log likelihood is
  bitwise, the parameters are 1.8e-11 apart.
* **F7**, no brms coverage of `student()`. Recorded in section 10 as
  an unowned gap.
* **The `m = 25` question.** Answered in section 8d: not shortenable
  from anything that is not a parameter, and the series-length
  alternative is priced and left as a follow-up.

Two of my six errors were the same mistake: reading a printed or
quantized zero as a measured zero. `%.6f` turned 1.8e-11 into
`0.000000`, and a 10 ms clock turned a 3-to-4-tick block ratio into a
speedup. Both survived into a document full of other measurements
because nothing in the write-up asked what the instrument could
resolve. The check that would have caught both is cheap: print a
suspected zero in `%e`, and measure the timer before trusting it.

Tests rerun after these changes, one file per R process: `test-families.R`
227, `test-diagnostics-ux.R` 127, `test-autocor.R` 191, `test-tre.R`
101, all with 0 failures and 0 errors. The new Cauchy chunk in
`vignette("diagnostics")` was run by hand (pandoc is absent here) and
prints `nu: (Intercept)  -19.06957  4260.962  1`.

## 9. What I did NOT do, and why

* **No bound on `nu`.** The only source of parameter bounds in
  `frmtmb` is a prior, deliberately. A hard ceiling would be an
  arbitrary constant that changes results, and it would hide the thing
  `diagnose()` should say out loud.
* **No change to the `RTMB::dt()` call in `R/priors.R`.** A student-t
  PRIOR takes its `df` as a fixed constant from the user, typically 3
  to 7, so there is nothing to cancel. Reachable only if a user writes
  `set_prior(student_t(df = 1e10))`, at which point they have asked for
  a normal prior. Recorded, not fixed.
* **No change to the LKJ term in `R/priors.R:2366`.** Same reason: its
  `eta` is a user-supplied prior constant.
* **The fuzz invariant was not weakened.** `fuzz_permutation_tol()` is
  untouched. The `reml-ar1-se-two-optima` entry is DELETED from
  `FUZZ_KNOWN_PENDING`, replaced by a comment saying what it was, and
  the tier is green without it.

## 10. Defects found and not fixed

* **`nu = Inf` gives `NaN`.** `dt_stable()` is exact up to
  `.Machine$double.xmax` (it returns the gaussian log density to the
  digit there) and returns `NaN` at `nu = Inf` exactly, because
  `(nu + 1)/2 * log1p(z^2/nu)` is `Inf * 0`. `RTMB::dt()` on the tape
  did the same, so this is not new. It needs `log(nu - 1) > 709.78`,
  and the furthest any fit measured here reached was 548. Guarding it
  needs a branch, which is the thing RTMB cannot give.
* **`diagnose()`'s `flat` check cannot reach a case like this.** It is
  gated on `!pdHess || any(!is.finite(se))`, and a plateau this flat
  still yields a positive-definite Hessian and finite (enormous)
  standard errors. `unbounded_dpar` covers the student case measured
  here; the general question of a flat direction on a fit whose
  covariance succeeds is open.
* **`student` with `nu` near the LOWER end of `logm1` is also flat.**
  In the ingredient table, `student`/`se`/`ri` reached `nu = 1.31` and
  `nu = 1.48`, with `logLik` gaps of 1.5e-01 and 1.8e-01 that stay
  under the derived tolerance, so the tier is quiet. `nu -> 1` is the
  other end of the same identifiability problem.

  **I wrote here that `unbounded_dpar` "does not see it (the estimate
  is small, not large)". That is wrong,** and the review disproved it.
  It is true only of the two cells I happened to measure, which stopped
  short of the boundary. When `nu` actually REACHES one the estimate is
  large and NEGATIVE: on Cauchy errors `log(nu - 1)` lands between
  -20.5 and -17.7 with a standard error of 4003 to 8533, and the check
  fired on 8 of 8 such fits. The trigger was always right; only the
  printed message was, which the reviewer fixed. My docs described only
  the upward case and now name both.
* **`unbounded_dpar` misses an unidentified dpar whose estimate is
  SMALL.** `skew_normal()` on symmetric data leaves `alpha` at 0.002
  with a standard error of 14.6, unidentified and not flagged, because
  the check requires `|est| > 10`. `alpha` also takes an IDENTITY link,
  so that threshold is not scale free there. Recorded by the review;
  not addressed here, because a check that fires on a small estimate
  with a large standard error would fire on every genuinely-zero
  coefficient in the package.
* **`unbounded_dpar`'s message names `student()` but the check is
  general.** It correctly fires on `negbinomial()` whose `shape` runs
  to the Poisson limit (log shape 16.0 to 19.0 on Poisson-generated
  data, 4 of 4 seeds), where the advice to "refit with `gaussian()`" is
  the wrong family to name. Naming `poisson()` there too is a small
  addition someone should make.
* **F7: the gated brms tier does not test this density at all.**
  `student()` appears in no brms comparison anywhere in the suite. The
  only match for `student` under `test-brms-*.R`, `test-bcm-*.R` and
  `test-rl-example.R` is the `student_t` PRIOR in
  `test-brms-priors.R`, which lives in `R/priors.R` and which this lane
  deliberately did not touch; `test-brms-likelihood.R`'s family roster
  is poisson, Gamma, negbinomial and bernoulli. Every gated count
  matches the 0.54.0 release figure, but that is evidence of no
  collateral damage and nothing else. The new tests use `stats::dt`,
  which is the right reference and which the review verified
  independently, so the density is not unvalidated. But if a future
  change moves it, brms will not notice. **One row for `student()` in
  `test-brms-likelihood.R` would close this.** Coverage gap, unowned.
* **A pre-existing `sparse_x` defect this lane made visible.** `family
  =student aterm=se re=nested special=mo_int mode=sparse_x`, seed
  20354693, does not converge on either library: code 1, "false
  convergence (8)", no positive-definite Hessian, 3 non-finite standard
  errors, `nu = 1.52312` on BOTH so nowhere near the cancellation
  regime. It became a `row_permutation` finding only because one tape
  now returns `NaN` at the midpoint where it returned a number. Not
  caused here and not worsened here; file it on its own.
