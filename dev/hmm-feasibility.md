# Hidden Markov models in frmtmb: feasibility probe

Date: 2026-09-02. Branch `wt-hmmprobe`, frmtmb 0.32.0, RTMB 1.9,
depmixS4 1.5.4, hmmTMB 1.1.2, R 4.6.1, Windows 11.

**No package code was changed.** Everything below runs through
`custom_family()` + `vint()`, which is already arbitrary R evaluated on
the AD tape, or through a bare `RTMB::MakeADFun`. Probe scripts are in
`dev/hmm/`; they run from the worktree root and the frmtmb ones want
`FRMTMB_LIB` pointing at a library holding the worktree build.

## Verdict

The thesis holds, without qualification.

1. The forward algorithm tapes in RTMB, gives the exact marginal
   likelihood, and its gradient is correct to 1e-8 against finite
   differences. Cost is linear in the number of time steps, and it is
   the number of *rows*, not the number of sequences, that costs.
2. It composes with Laplace over Gaussian random effects exactly as the
   v0.19 group-mixture insight predicts: sum the discrete states inside,
   integrate `b` outside. Measured against adaptive Gauss-Hermite
   quadrature over the per-group scalar `b`, the Laplace error is 0.126
   in the log-likelihood over 40 groups of 25 observations
   (8.9e-5 relative) and 4.4e-4 in the parameters.
3. **Rung 1 works today.** A gaussian HMM with covariate-dependent
   transitions and per-state random effects fits through
   `frm(bf(y | vint(g, t) ~ ...) + custom_family(...))` and reaches the
   same optimum as a hand-rolled `MakeADFun`, as depmixS4, and as
   hmmTMB. `check_custom_family()` passes. Categorical emissions - the
   covid19retrospective shape - work the same way.

Two things spoil the ergonomics rather than the capability, and both are
what rung 2 exists to fix:

- The post-fit surface **silently lies**. `fitted()`, `residuals()` and
  `predict(type = "response")` all return state 1's mean at every row,
  with no warning, because `response_mean()` falls back to `dpars$mu`
  when a family supplies no `mean_fn`. There is no per-row `mean_fn`
  that could be right: `E[y_t]` needs the state-occupancy probability,
  which only forward-backward produces.
- The cold start reaches a **converged local optimum** 8.1 log-likelihood
  units below the global one on the random-effect model (probe D4),
  with `convergence == 0`, `max|grad| == 3.5e-4` and a positive-definite
  Hessian. `diagnose()` reports no problem. Multimodality is inherent to
  HMM ML, as it is to mixtures, but nothing in the rung-1 recipe warns.

Recommendation: build `hmm(K, family)` as a first-class family, on the
mixture machinery. Sizing at the end: 7-9 days.

## 1. Probe A: the bare forward algorithm

`dev/hmm/probeA1-forward-rtmb.R`. 2-state gaussian HMM, one sequence of
T = 200, free initial distribution, transition matrix a row-wise
multinomial logit.

Two formulations were taped and compared: log-space via
`RTMB::logspace_add`, and Zucchini scaling (rescale alpha, accumulate
`log(sum)`).

| check | value |
| --- | --- |
| numeric scaled vs numeric log-space forward, at the truth | 3.7e-13 |
| taped log-space vs taped scaled, at their optima | 2.8e-13 |
| taped optimum vs the numeric forward at the same parameters | 3.1e-13 |
| AD vs central-difference gradient, log-space | 8.1e-8 |
| AD vs central-difference gradient, scaled | 1.3e-8 |
| depmixS4 (best of 5 random EM starts) | -235.204363002 |
| RTMB | -235.204362981 |
| difference | 2.1e-8 |

Every depmixS4 parameter matches to five decimals (transition matrix
0.94377/0.05623/0.30748/0.69252 both; means -1.00209 vs -1.00208; SDs
0.57349/0.78621 vs 0.573494/0.786188).

### Which formulation

Both are correct and both tape. **Log-space is what to build on**, but
the scaled form is measurably cheaper on the gradient (0.061 ms vs
0.111 ms per `gr()` at T = 200 - it avoids `logspace_add`'s branch-free
`log1p(exp(-|d|))` in every step) and the tape build is a wash
(21.7 ms log, 24.6 ms scaled). The scaled form's exponentials underflow
on long sequences of low-density observations, which is exactly what a
K-state model with a poorly-fitting state does in early optimizer
iterations, so the robustness is worth the 1.8x on `gr`.

### Cost vs T

Log-space, ms per call, one sequence, K = 2. (`obj$fn(obj$par)` in a
loop times TMB's *cache*, not the tape - every number here moves the
parameter on each call.)

| T | tape | fn | gr | full nlminb |
| --- | --- | --- | --- | --- |
| 200 | 17.5 | 0.043 | 0.102 | 5.9 |
| 1 000 | 77.9 | 0.181 | 0.489 | 38.8 |
| 5 000 | 438.2 | 0.993 | 3.201 | 212.7 |
| 20 000 | 1 862.7 | 4.720 | 8.380 | 754.3 |
| 50 000 | 5 275.1 | 10.677 | 17.695 | 2 454.2 |

Evaluation is linear as predicted (250x the data costs 248x on `fn`).
The **tape build** is the part that grows slightly faster than linear
(301x for 250x) and it is what becomes annoying: 5.3 s to tape a single
50 000-step chain, against 10 ms to evaluate it. For an interactive
model-comparison workflow the wall hits somewhere around T = 20 000
(1.9 s per tape build, so ~2 s of dead time on every `frm()` call);
below T = 5 000 nothing is noticeable.

### Number of sequences is free

`dev/hmm/probeA2-group-scaling.R` holds the row count at 5 000 and
varies how it is cut into sequences.

| N | T | tape (ms) | fn (ms) | gr (ms) | nlminb (ms) |
| --- | --- | --- | --- | --- | --- |
| 1 | 5 000 | 543.4 | 0.697 | 1.893 | 122.9 |
| 10 | 500 | 629.1 | 0.857 | 2.392 | 140.4 |
| 50 | 100 | 607.8 | 0.953 | 2.721 | 172.5 |
| 250 | 20 | 797.1 | 0.918 | 2.442 | 143.0 |
| 1 000 | 5 | 652.8 | 0.877 | 3.353 | 358.2 |

Rows are the unit of cost. A thousand short sequences cost about the
same as one long one; the per-sequence loop only adds that sequence's
own initial and terminal reductions. This is the opposite of the ODE
finding (`dev/ode-feasibility.md` section 7), where the per-group solve
was the whole cost.

### The RTMB gotchas, measured

- `[[i, j]]` was never needed. The trick that avoids advector-matrix
  indexing entirely is to hold the emission log-densities as a **list of
  K length-n advector columns** and the transition matrix as a list of
  K rows. Then one forward step is `K - 1` calls to `logspace_add`,
  each vectorized over the K target states, plus a `c()` of K scalars
  for the emission row.
- `ADoverload("c")` is required and is **lexical**: every helper that
  builds a vector with `c()` needs its own `"c" <- RTMB::ADoverload("c")`
  line, including the ones called from inside an `lpdf`.
- The `[<-` gotcha is much milder here than the >1000x folklore. Writing
  the recursion with per-cell `new[j] <- ...` instead of the vectorized
  fold costs **1.54x on the tape build and nothing at all on `fn`/`gr`**
  (T = 200, K = 2), and gives a bit-identical value. The >1000x figure
  is about assigning into n-length vectors; here the assigned vector is
  length K, so the cost scales with `K * T` extra tape nodes, not
  `n * T`. It would still bite at K = 10.
- `RTMB::solve`, not `base::solve`, for the stationary distribution
  (section 6). `RTMB::matrix`, not `base::matrix`, to assemble it.
- No branching on parameter values was ever needed. Missing
  observations, which are the obvious place to want one, are handled
  with a **data multiplier** (section 6, F4).

## 2. Probe B: many sequences, covariate-dependent transitions

`dev/hmm/probeB1-multiseq-covtrans.R`. N = 30 sequences of T = 20,
K = 2, each transition-matrix row a multinomial logit in a covariate
`x`, plus a free initial distribution. 9 parameters.

Tape 0.12 s, optimization 0.01 s, convergence 0.

| check | value |
| --- | --- |
| taped logLik | -891.019360499 |
| numeric forward at the taped optimum | -891.019360499 |
| difference | 0 |
| AD vs FD gradient | 9.7e-9 |
| depmixS4 (`transition = ~ x`, `ntimes = rep(20, 30)`) | -891.019361782 |
| difference | 1.3e-6 (EM tolerance) |

Every depmixS4 coefficient matches to four decimals, which pins the
**convention mapping** exactly, and it is worth writing down because
depmixS4 and hmmTMB disagree about it:

- Each *row* of the transition matrix is its own multinomial logit.
- depmixS4's reference cell is **state 1** for every row. hmmTMB's is
  the **diagonal** (`ref = 1:n_states`). Both span the same model; only
  the reported coefficients differ. The probes use depmixS4's, because
  it is the one a `tr12 ~ x` formula spelling reads naturally.
- The covariate value at time `t` drives the transition **from t to
  t + 1**. The covariate on a sequence's last row is never used.

## 3. Probe C: the rung-1 test

`dev/hmm/probeC1-frm-custom-family.R`, family in
`dev/hmm/hmm-family.R`. The probe B1 model, verbatim, through `frm()`.

```r
hmm2_lpdf <- function(y, dpars, aterms, extra = list(hmm_ldel = 0)) {
  "c" <- RTMB::ADoverload("c")
  rows_by_g <- hmm_seq_index(aterms$vint1, aterms$vint2)
  lp <- list(RTMB::dnorm(y, dpars$mu,  dpars$sigma1, log = TRUE),
             RTMB::dnorm(y, dpars$mu2, dpars$sigma2, log = TRUE))
  lg <- tpm_logs_ad(list(list(0, dpars$tr12),
                         list(0, dpars$tr22)), 2L)
  ld <- log(softmax0_ad(extra$hmm_ldel, 2L))
  llv <- NULL
  for (gi in seq_along(rows_by_g)) {
    v <- fwd_ad_log_tv(lp, lg, rows_by_g[[gi]], ld, 2L)
    llv <- if (is.null(llv)) v else c(llv, v)
  }
  first <- vapply(rows_by_g, function(r) r[1], integer(1))
  S <- Matrix::sparseMatrix(i = first, j = seq_along(first), x = 1,
                            dims = c(length(y), length(first)))
  as.vector(S %*% llv)
}

fit <- frm(bf(y | vint(g, t) ~ 1, mu2 ~ 1, sigma1 ~ 1, sigma2 ~ 1,
              tr12 ~ x, tr22 ~ x) + fam, data = dat)
```

### Does the contract permit it

Yes, and nothing had to bend. Three details make it work:

1. `R/objective.R:250` forms `nll - sum(w * ll)` and nothing else, so an
   `lpdf` may return **anything that sums to the right total**. Putting
   each sequence's log-likelihood on its first row and zero elsewhere is
   legal by construction. Weights would break it (`weights()` must be
   refused), and so would `cens()`/`trunc()`, which index into `ll` per
   row.
2. `vint()` reaches the `lpdf` as `aterms$vint1`, `vint2`, ... verbatim
   (`R/parse.R:58-67`), integer-checked at frame time
   (`R/frame.R:838`). Group id and time order fit exactly. A third
   payload carries the missing-data flag (section 6, F4).
3. `extra_pars` gives the initial-distribution logit a home outside the
   dpar system, reaching the `lpdf` as its fourth argument
   (`R/objective.R:188-192`). It is counted in `df`, appears in
   `confint()` as `hmm_ldel_1`, and is read back with
   `fit$estimates[["hmm_ldel"]]`.

The **scatter should be a sparse matrix multiply**, not `[<-`: `S %*%
llv` with a constant `n x N` indicator puts no assignment nodes on the
tape at all, and it is the same construction `mix_g$Gt` already uses.

### Does it agree

| check | value |
| --- | --- |
| `frm()` logLik | -891.019360495 |
| hand-rolled `MakeADFun` (probe B1) | -891.019360499 |
| difference | 4.2e-9 |
| max coefficient difference (8 betas) | 8.4e-6 |
| `hmm_ldel` | matches to 8 digits |
| `df` | 9 (8 betas + 1 extra) |
| fit time | 0.22 s |

### check_custom_family()

**Passes, unmodified.** No skip, no adaptation, no loosened tolerance.
Its finite-difference check is over the `dpars` only, so a
non-factorizing `lpdf` is no obstacle - the objective is smooth in every
dpar. Two notes:

- `check_custom_family()` calls `family$lpdf(y, p, aterms)` with **three
  arguments**. A family with `extra_pars` therefore needs a default for
  the fourth (`extra = list(hmm_ldel = 0)`), or the check errors on a
  missing argument. That is a one-word fix in the family and a candidate
  one-line fix in `R/interop.R` (pass `family$extra_pars(y, aterms)`
  through when it exists).
- The test `dpars` must be full-length vectors, not scalars, because the
  recursion indexes them by row.

### Post-fit methods

This is the honest weak spot.

| method | behavior |
| --- | --- |
| `summary()` | ok |
| `confint(method = "wald")` | ok, all 9 rows including `hmm_ldel_1` |
| `predict(dpar = "tr12")` | **meaningful**: the per-row transition logit |
| `predict(type = "link")` | state 1's mean at every row |
| `fitted()` | **silently wrong**: state 1's mean at every row |
| `residuals()` | **silently wrong**: `y - mu1` |
| `predict(type = "response")` | **silently wrong**: same |
| `simulate()` | refuses cleanly: "family 'hmm2_gaussian' has no simulator yet" |
| `ranef()` / `VarCorr()` | ok (they never touch the family) |

`fitted()` runs even with `post = list()`, because `response_mean()`
falls back to `dpars$mu`. Supplying `mean_fn = function(dpars, aterms)
dpars$mu` changes nothing except making the intent explicit: on probe
B1's data `fitted()` is constant at 0.0382 while `y` ranges over
[-2.11, 5.81], and `cor(fitted, y)` is undefined.

There is no fix at the family level. `E[y_t]` under an HMM is
`sum_k P(S_t = k | y, theta) mu_k`, and the occupancy probability comes
from forward-backward over the whole sequence, which a per-row
`mean_fn(dpars, aterms)` cannot compute. Computing it post hoc from the
fit is easy and cheap (probe C section 7): local decoding recovers
98.7% of the simulated states, and the occupancy-weighted mean
correlates 0.914 with `y`. **That is rung 2's most valuable single
deliverable.**

## 4. Probe D: random effects, and how good the Laplace is

`dev/hmm/probeD1-random-effects.R`. N = 40 sequences of T = 25, a random
intercept on state 1's mean, true `sd(b) = 0.8`.

`frm()` fits it in 1.1 s and agrees with a hand-rolled
`MakeADFun(random = "b")` over the same forward algorithm to **6.2e-9**.
That only says the plumbing is right, since it is the same Laplace. The
real question needs a reference that is not a Laplace at all.

Because `b` is a per-group scalar and the group likelihood given `b` is
*exact* (the forward algorithm), the marginal is a one-dimensional
integral. **Adaptive** Gauss-Hermite - recentred on the conditional mode
and rescaled by its curvature - converges immediately:

| nodes | logLik |
| --- | --- |
| 3 | -1419.1586 |
| 5 | -1419.1437 |
| 11 | -1419.14156 |
| 21 | -1419.14150 |
| 41 | -1419.141497 |

Naive prior-weighted Gauss-Hermite does **not**: -1433.09, -1421.68,
-1419.70, -1419.04 at 11, 21, 41, 81 nodes. A sequence of 25
observations makes the conditional posterior of `b` far sharper than its
prior, which is the textbook failure mode; anyone repeating this
measurement must adapt.

| quantity | Laplace | quadrature (nq = 41 / 15) |
| --- | --- | --- |
| logLik | -1419.26702 | -1419.14150 |
| mu1 | -0.010693 | -0.011133 |
| mu2 | 2.969338 | 2.969014 |
| sigma1 | 0.615161 | 0.615018 |
| sigma2 | 0.591292 | 0.591528 |
| tr12 | -1.636458 | -1.636049 |
| tr22 | 1.411608 | 1.411772 |
| sd(b) | 0.870877 | 0.870490 |

Laplace bias: **-0.126 in the log-likelihood** (8.9e-5 relative), **4.4e-4
absolute in the parameters**, 0.04% relative on `sd(b)`. The requested
1e-4 tolerance on logLik is *not* met and cannot be: the integrand is a
mixture over state sequences, not a Gaussian, so the Laplace is
genuinely approximate here - the same situation as the v0.19 group
mixtures, where the recorded bias was ~0.1 in logLik and 0.01 in the
parameters. This is a **better** result than v0.19 on the parameters by
about an order of magnitude, presumably because a long sequence makes
each group's conditional tighter and more Gaussian.

### Third-party agreement: hmmTMB

`dev/hmm/probeD2-hmmtmb.R`. A model both packages express exactly:
K = 2, gaussian, constant transitions, **stationary** initial
distribution, with and without per-state random intercepts on the means.

Fixed effects only, both started at frm's own optimum:

| quantity | frm | hmmTMB |
| --- | --- | --- |
| logLik | -1216.40337453 | -1216.40337453 (8.7e-10) |
| mu1, mu2 | -0.18119852, 3.07922909 | -0.18119873, 3.07922925 |
| sd1, sd2 | 0.80618908, 0.79517108 | 0.80618881, 0.79517087 |
| tpm | 0.83056236 0.16943764 0.15862241 0.84137759 | 0.83056251 0.16943749 0.15862219 0.84137781 |

Max parameter difference 2.8e-7. That is a genuine third-party
validation of the whole rung-1 stack.

### hmmTMB sharp edge, found the hard way

`dev/hmm/probeD3-hmmtmb-knownstate.R`. **hmmTMB silently reads a data
column named `state` as KNOWN STATES** and maximizes the complete-data
likelihood instead. The obvious thing to do in a simulation study -
keep the generating state sequence in the data frame - therefore
compares against a different model, with no message:

| data columns | hmmTMB `llk()` | plain forward at the same parameters |
| --- | --- | --- |
| ID, t, state, y | -1247.248426 | -1217.775828 |
| ID, t, y | -1216.403375 | -1216.403375 |
| ID, y | -1216.403375 | -1216.403375 |

Anyone validating against hmmTMB must drop that column first. Worth an
upstream documentation issue.

### The multimodality finding

`dev/hmm/probeD4-re-multimodality.R`. On the random-effect model frm's
cold start and hmmTMB land 8.1 log-likelihood units apart. The arbiter
is a hand-rolled `MakeADFun(random = c("b1", "b2"))` over the same
model.

- The two objectives are **the same function**: evaluated at frm's own
  cold-start vector, frm gives -1449.0619967355 and the hand-rolled
  gives -1449.0619967353 (1.6e-10).
- The hand-rolled optimum, reached from *both* starting points, is
  **-1087.99646521**, which is hmmTMB's value to the last digit, with
  `sigma` 0.610455/0.602054 and `sd(b)` 0.645297/0.427590 matching
  hmmTMB's 0.610454/0.602054 and 0.64530/0.42759.
- frm from its default cold start stops at **-1096.09575602** with
  `convergence == 0`, `max|grad| == 3.5e-4`, a positive-definite
  Hessian, and `diagnose()` reporting "No convergence problems
  detected".
- frm **restarted** at the global point reaches -1087.9964652.

So the likelihood is right and the optimizer is honest; the surface has
a second mode and the quantile-spread cold start finds it. Rung 2 must
say so and must ship the remedy (`frm_allfit`, multiple starts, or a
short deterministic-annealing pass).

## 5. Probe E: categorical emissions

`dev/hmm/probeE1-categorical.R`. K = 2 states, C = 4 observed
categories, per-state category probability vectors as multinomial
logits, covariate-dependent transitions, N = 40 sequences of T = 15.
This is the covid19retrospective shape. Ten parameters, eight of them
dpars with their own formulas.

The emission term is written branch-free with data indicators: pick the
numerator out of the state's logit vector with `sum_c 1(y == c) * eta_c`
and subtract the shared log-normalizer. No `[[ ]]`, no lookup on the
tape.

| check | value |
| --- | --- |
| `check_custom_family()` | PASS |
| family `lpdf` vs numeric forward at the truth | 0 |
| `frm()` logLik | -772.154040795 (df 10, 0.29 s) |
| numeric forward at frm's estimates | -772.154040795 (1.1e-13) |
| independent BFGS on the numeric forward | -772.154040795 (1.1e-13) |
| max parameter difference vs that BFGS | 1.5e-8 |

**The covid model class is reachable.** Two honest caveats: one emission
probability went to the 0 boundary (state 2, category 2: fitted 0.0000
against a true 0.0724), which an unpenalized multinomial logit will do
whenever a category is rare within a state - a rung-2 `hmm()` should
document it and point at `set_prior()` on the emission logits; and local
decoding recovers only 80.8% of the states, against 98.7% for the
well-separated gaussian model, because overlapping category profiles are
much less identifiable.

## 6. Probe F: sharp edges

`dev/hmm/probeF1-sharp-edges.R`.

### F1. Label switching and initialization

Probe B1 data, reference logLik -891.0193605. `mixture()`'s own
convention is `quantile(y, k / (K + 1))`.

| mu, mu2 init quantiles | logLik | mu1 | mu2 |
| --- | --- | --- | --- |
| 1/3, 2/3 (mixture's own) | -891.01936051 | 0.0383 | 3.0881 |
| 0.25, 0.75 | -891.01936050 | 0.0383 | 3.0881 |
| 0.10, 0.90 | -891.01936050 | 0.0383 | 3.0881 |
| 2/3, 1/3 | -891.01936051 | 3.0881 | 0.0383 |
| 0.45, 0.55 | -891.01936049 | 0.0383 | 3.0881 |
| **0.50, 0.50** | **-1153.15064660** | 1.2569 | 1.2569 |

Label switching is benign: the reversed init reaches the identical
optimum with the labels swapped, which is the expected `K!` symmetry.
What breaks it is a start on the **symmetry axis**. With both means at
the median the objective's gradient is symmetric under a label swap, the
optimizer never leaves the diagonal, and the fit converges to the
one-state solution 262 units below. A fully symmetric start (both means
at the median *and* both transition logits 0) does the same. The
mixture()-style quantile spread is exactly the right defence, and a
rung-2 `hmm()` should additionally refuse (or nudge) a user-supplied
start whose component means coincide.

### F2. Sequences of length 1 degenerate to a mixture

With every sequence of length 1, the forward algorithm collapses to
`sum_k delta_k f_k(y)`, which is a finite mixture. n = 400.

| | logLik | mu1 | mu2 | sigma1 | sigma2 | weight |
| --- | --- | --- | --- | --- | --- | --- |
| HMM (rung 1) | -696.486410877 | 0.115657 | 2.952683 | 0.907375 | 0.585206 | 0.630406 |
| `mixture(gaussian(), gaussian())` | -696.486410877 | 0.115657 | 2.952683 | 0.907375 | 0.585206 | 0.630406 |

Difference 3.4e-13, every parameter identical to six digits. A clean
internal-consistency check and a good regression test for rung 2. The
`df` differ, 7 against 5: the HMM carries two transition logits that a
length-1 sequence never uses, so they are unidentified and the reported
`df`/AIC are wrong. A rung-2 `hmm()` should detect an all-singleton
grouping and refuse, or drop the transition block.

### F3. Stationary initial distribution

The on-tape linear solve works. `delta (I - Gamma + 1 1') = 1'` with
`RTMB::solve` on a `RTMB::matrix` of the taped transition probabilities:
`check_custom_family()` passes on the stationary family, and the taped
value matches the numeric stationary forward to **4.6e-13**. It
differentiates (that check *is* a gradient check).

On probe D1's data:

| initial distribution | logLik | df |
| --- | --- | --- |
| free (one extra logit) | -1608.76549264 | 7 |
| stationary (on-tape solve) | -1609.41007570 | 6 |

`2 * (free - stationary) = 1.29` on 1 df. Both are worth offering;
stationary is the right default for long sequences and the only sane
choice when sequences are numerous and short (hmmTMB warns about exactly
this). The constraint is that the transition matrix must be **constant**
- with covariate-dependent transitions there is no single stationary
vector - so `init = "stationary"` has to be refused at frame time
whenever any transition dpar has a non-intercept term.

### F4. Missing observations mid-sequence

The correct likelihood drops the emission term at a missing time point
and keeps the transition. Two routes, N = 30, T = 20, 90 of 600 rows
missing at t = 7, 8, 14.

- **Masked route (correct).** Keep the row, carry the missing flag as a
  third `vint()` payload, and multiply the emission log-density by
  `1 - vint3`. That is *data*, so it is a constant weight on the tape:
  no branching, no `NA` ever touching an advector. `frm()` gives
  -691.400711096 and the numeric reference at the same parameters gives
  -691.400711096. Difference **0**.
- **`na.omit` route (silently wrong).** Leaving `NA` in the response
  lets `na.action` drop the rows, which *shortens the chain*: time 6 and
  time 9 become adjacent, and one transition stands in for three.

| tpm | 1->1 | 1->2 | 2->1 | 2->2 |
| --- | --- | --- | --- | --- |
| masked | 0.88464 | 0.11536 | 0.22422 | 0.77578 |
| na.omit | 0.87158 | 0.12842 | 0.24900 | 0.75100 |
| true | 0.90 | 0.10 | 0.25 | 0.75 |

The bias is small here because only three of twenty points are dropped
and the chain is sticky; it grows with the gap length. A rung-2 `hmm()`
owns the response's `NA` handling: keep the row, mask the emission,
never let `na.omit` reach an HMM response.

### F5. REML

REML **runs without complaint and should not**. `frm(..., REML = TRUE)`
integrates only the `primary_dpars` fixed effects
(`R/families.R:32-36`), which for this family is `mu` alone. The other
five dpars' coefficients - including both transition blocks and both
state SDs - stay in the outer problem. The result is a *partial*
restricted likelihood that corresponds to no standard definition:

- fixed-effect HMM: ML -1608.76549264, "REML" -1610.95250460
- with `(1 | g)`: "REML" -1420.30322741, df 7

A rung-2 `hmm()` should **refuse REML outright**. There is no useful
restricted likelihood for a model whose location parameters are per
state, and the current silent partial answer is a trap.

## 7. Rung 1: the recipe as it stands

It works cleanly, and it is vignette-ready with three warnings attached.
The recipe:

1. `vint(g, t)` - and a third payload if the response has gaps.
2. Dpars: one per state-dependent parameter (`mu`, `mu2`, `sigma1`,
   `sigma2`, ...) plus `K * (K - 1)` transition logits, each with its
   own formula, so covariate-dependent transitions and mixture-of-experts
   gating are free.
3. `extra_pars` for the initial distribution, or the on-tape stationary
   solve if the transitions are constant.
4. `lpdf` returns each sequence's log-likelihood on its first row, zero
   elsewhere, scattered with a constant sparse indicator.
5. `init_dpars` on spread response quantiles, mixture()'s convention.
6. Run `check_custom_family()` first; give `lpdf` a default fourth
   argument so it can.

The warnings: `fitted()`/`residuals()`/`predict(type = "response")` are
wrong and silent; `weights()`, `cens()`, `trunc()` and `mi()` must not
be used; REML must not be used; and the cold start can converge to a
local optimum without any diagnostic firing.

## 8. Rung 2: the `hmm()` design

### What transfers verbatim from `mixture()`

- **Suffixed per-state dpars with the full formula grammar.**
  `mixture()` already generates `mu1..muK`, `sigma1..sigmaK` from the
  component families and gives each a linear predictor, random effects
  included. `hmm(K, gaussian())` is the same construction with the same
  component-family walk.
- **The multinomial-logit weight block.** `mixture()`'s `theta` dpars
  are a `K`-way multinomial logit with full linear predictors, which is
  what makes mixture-of-experts free. An HMM's transition block is `K`
  copies of that same object, one per source state.
- **The logsumexp primitive.** `RTMB::logspace_add` folded over states.
- **Quantile-spread mean initialization** (`R/families.R:1979`,
  `quantile(y, k / (K + 1))`) - probe F1 confirms it is the right
  default and shows the one start it must refuse.
- **The group-structure frame object.** `frame$mix_g` already carries
  `G`, `Gt`, `first`, `gindex`, `levels` per response
  (`R/frame.R:936-956`). An HMM needs exactly this plus a within-group
  time order.
- **The sum-inside-integral objective shape.** `R/objective.R:154-171`
  is the template: a family-level branch that replaces the row-wise
  `sum(w * ll)` with a per-group reduction, evaluated *inside* the
  Laplace. The HMM branch is that block with the logsumexp-over-classes
  replaced by the forward recursion.
- **`mixture_probs()`** (`R/families.R:2106`) is the shape of
  `hmm_probs()`.

### The genuinely new pieces

**Spec.** `hmm(K, family, time = , group = , init = c("stationary",
"estimated", "uniform"), trans = ~ 1)`. `time` and `group` are captured
expressions, like `mixture(..., groups = g)` already does
(`R/frame.R:537, 660, 936`). `trans` is the default formula for every
transition cell; per-cell overrides come through the ordinary dpar
formula slots in `bf()`, named `tr12`, `tr13`, ..., `trK,K-1` on the
depmixS4 (state-1 reference) convention.

**Frame.** One new structure per response, alongside `mix_g`:

```
hmm_g = list(rows   = <list of N integer vectors, time-ordered>,
             first  = <N indices>,
             S      = <n x N sparse scatter>,
             const_trans = <logical: every transition dpar
                            intercept-only>,
             miss   = <n logical: response NA>)
```

Checks at frame time, not on the tape: `time` unique within `group`; no
all-singleton grouping (F2); `init = "stationary"` only when
`const_trans` (F3); `weights()`/`cens()`/`trunc()`/`mi()`/`rescor`
refused; `REML` refused (F5). NA responses are **kept** and masked
(F4), which is a deliberate departure from the package-wide
`na.action` and needs its own message.

**Objective.** A third branch beside `rescor` and `mix_g`, about 40
lines: build the K emission log-density columns from the per-state
dpars, the `K x K` log-transition columns from the transition dpars via
one shared normalizer per source state, the initial log-distribution
(uniform, an `extra_par` logit, or `RTMB::solve` on the constant
transition matrix), then loop groups running the log-space recursion and
`sum` the results. Random effects need no special handling at all -
they are already in the dpar linear predictors, and the Laplace wraps
the whole thing.

**Post-processing.** The largest single piece and the most valuable.
Numeric, off-tape, from the fitted dpars:

- `hmm_probs(fit)` - `n x K` smoothed state probabilities
  (forward-backward), the `mixture_probs()` analog.
- `hmm_states(fit)` - Viterbi global decoding.
- `fitted()`/`residuals()`/`predict(type = "response")` route through
  `hmm_probs()` so `E[y_t] = sum_k P(S_t = k | y) mu_k(x_t)` is what
  they report. This is what makes rung 2 worth building.
- `simulate()` - forward-simulate the state path per group, then the
  emissions. Needs a family-level simulator hook that can see the
  transition dpars, which `sim(dpars, aterms, n)` cannot; the same gap
  `mixture_mvn` logged.

**Cost of the post-processing.** Forward-backward is `O(n K^2)`,
identical in shape to the forward pass, and it runs once. Probe C's
naive R implementation takes 68.5 ms for n = 600, which is 100x the
taped forward evaluation and would be visible at n = 1e5; it should be
written with the same list-of-columns vectorization as the tape version,
or in a `Reduce` over pre-multiplied matrices, not a per-row R loop.

### What to refuse initially

`weights()`, `cens()`, `trunc()`, `mi()` on the response, `rescor` /
`mvbf` components, quadrature (`frm(quadrature = )` - the integrand is
not the one it assumes), OSA residuals, `REML`, `predict(se.fit =
TRUE)` on the response scale (the delta method would need the gradient
of the occupancy probabilities). Continuous-time / irregular-interval
transition matrices are a separate model and stay out. Higher-order
Markov chains stay out. `hmm()` components that are themselves mixtures
stay out.

### Sizing

Mixture-wave scale, as the gaps entry guessed.

| piece | days |
| --- | --- |
| `hmm()` family constructor, dpar and transition naming, links, inits | 1 - 1.5 |
| frame: time/group contract, ordering, the six refusals, NA masking | 1 |
| objective branch, plus RE composition (no work expected) | 0.5 - 1 |
| post-processing: forward-backward, Viterbi, hmm-aware fitted/predict/residuals | 1.5 - 2 |
| `simulate()` and the family-level simulator hook | 0.5 - 1 |
| multi-start / label-switching guidance and diagnostics | 0.5 |
| tests: vs depmixS4, hmmTMB, hand-rolled, the `mixture()` degeneracy, the masked-NA reference | 1 - 1.5 |
| vignette | 1 |
| **total** | **7 - 9.5** |

The two risks are both measured rather than guessed: the post-processing
leg is the one that is genuinely new work (nothing in the package
computes a backward pass), and the multimodality of probe D4 means the
test suite needs deterministic starts or it will be flaky.

## Probe scripts

Run from the worktree root. The frmtmb ones want
`FRMTMB_LIB=<library holding the 0.32.0 build>`.

| script | what it establishes |
| --- | --- |
| `dev/hmm/hmm-common.R` | simulation, numeric forward (scaled and log-space), adaptive GH, the taped recursions |
| `dev/hmm/hmm-family.R` | the rung-1 `custom_family()`, free and stationary initial distribution |
| `dev/hmm/probeA1-forward-rtmb.R` | the bare recursion, both formulations, depmixS4, cost vs T, the `[<-` cost |
| `dev/hmm/probeA2-group-scaling.R` | cost is rows, not sequences |
| `dev/hmm/probeB1-multiseq-covtrans.R` | many sequences, covariate transitions, depmixS4, the convention mapping |
| `dev/hmm/probeC1-frm-custom-family.R` | the rung-1 test through `frm()`, the post-fit survey, post-hoc forward-backward |
| `dev/hmm/probeD1-random-effects.R` | Laplace vs adaptive quadrature over the per-group `b` |
| `dev/hmm/probeD2-hmmtmb.R` | the on-tape stationary solve; exact agreement with hmmTMB |
| `dev/hmm/probeD3-hmmtmb-knownstate.R` | hmmTMB's silent `state` column |
| `dev/hmm/probeD4-re-multimodality.R` | frm's cold start converges to a local optimum 8.1 below the global one |
| `dev/hmm/probeE1-categorical.R` | categorical emissions: the covid19retrospective shape |
| `dev/hmm/probeF1-sharp-edges.R` | label switching, T = 1 vs `mixture()`, stationary vs free, missing data, REML |

`probeC1` and `probeF1` read the `.rds` files `probeB1` and `probeD1`
write, so run them in order.
