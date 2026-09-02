# Cluster-robust (sandwich) covariance: design record

Delivered in `R/sandwich.R` (`vcov_cluster()`, `cluster_scores()`),
with a gated weights hook in `R/objective.R`, a `cluster =` surface on
`vcov()`, and `vcov =` acceptance in `confint()`, `hypothesis()` and
`summary()`. Probes live in `dev/sandwich/`.

## Why per-cluster scores exist and per-observation ones do not

The optimized objective is a MARGINAL likelihood. Once a random effect
is integrated out, an observation has no separable contribution to it,
so `sandwich::estfun()` has nothing to return, and the package still
ships no `estfun()` method (the older note in `feature-gaps.md` stands).

A CLUSTER does have a contribution, whenever the marginal likelihood
factors: `nll = sum_g nll_g`. That holds exactly when every random
effect is nested in, or equal to, the clustering factor, and when every
likelihood term is a product over rows. The Laplace approximation does
not break the factorization - the joint inner Hessian is block diagonal
over clusters, so `log det H = sum_g log det H_g` and the approximated
objective factors as exactly as the true one does.

## The score route: two candidates, both probed

`dev/sandwich/probe-scores.R` implements and times both. They rest on
the same identity: set cluster `g`'s per-row data weight to 1 and every
other cluster's to 0, and the resulting objective IS `nll_g`. The
zero-weight clusters keep their random-effect prior at weight 1, so
their Laplace contribution is the Gaussian prior integrated over the
whole space - exactly 1, so exactly 0 on the log scale, with exactly
zero gradient in `theta`. Nothing is discarded and nothing is
approximated.

* **route A, re-tape.** Bake the 0/1 mask into
  `frame$aterm_values[[r]]$weights` and call `RTMB::MakeADFun()` once
  per cluster. Needs no change to `objective.R` at all.
* **route B, weight parameter.** Give the objective one extra
  parameter per cluster, `clw`, multiplying that cluster's rows. ONE
  tape; read `obj$gr(c(theta_hat, e_g))` on the `theta` block. Needs
  the (gated) hook in `build_objective()`.

Probe numbers, R 4.6.1:

| model | n | G | route A | route B | max abs A - B |
|---|---|---|---|---|---|
| poisson glm, one cluster per row | 120 | 120 | 0.33 s | 0.010 s | 0 |
| gaussian LMM, cluster = g | 240 | 30 | 0.20 s | 0.020 s | 0 |
| bernoulli GLMM | 600 | 100 | 0.92 s | 0.160 s | 4.3e-14 |
| nested REs, cluster = school | 240 | 20 | 0.19 s | 0.030 s | 7.1e-14 |

The two agree to machine precision everywhere, and `sum_g s_g` matches
`obj$gr()` at the optimum to between 3e-15 and 1.6e-8. Route B wins on
cost by 5 to 30 times, entirely because it tapes once rather than G
times; both are O(G) full-model gradient evaluations otherwise. Route B
ships.

The hook it needs is one line, gated on a frame field that is only ever
set on a private copy of the frame inside `cluster_scores_at()`:

```r
clw_idx <- frame[["cluster_w"]]           # NULL for every ordinary fit
...
if (!is.null(clw_idx)) w <- w * pars[["clw"]][clw_idx[[r]]]
```

It sits at the one place a per-row weight already enters, `sum(w * ll)`.
Every likelihood whose contribution is NOT a product over rows -
`autocor()`, `hmm()`, `rescor`, `mi()`/`me()` - is refused by the
guards rather than reached by the hook.

## What is shipped, and what is refused

`V = a_G B (sum_g s_g s_g') B` with `B` the inverse observed
information of the marginal likelihood (`vcov(full = TRUE)`), over the
WHOLE outer parameter vector. `a_G` follows clubSandwich exactly: CR0
= 1, CR1 = `G/(G-1)`, CR1p = `G/(G-p)`, CR1S = `G(N-1)/((G-1)(N-p))`,
with `p` the number of estimated outer parameters.

CR2/CR3 are refused. Both are built from the hat matrix of a linear or
GLS model; clubSandwich supports `lmerMod` by representing it as
weighted least squares with a working covariance, which a nonlinear
Laplace-marginal likelihood has no analogue of. No derivation, no
adjustment.

The nesting guard is structural rather than name based: it takes the
row / coefficient-column incidence of every `Z`, partitions the
coefficient columns into prior-independent units (one per level, except
for `gr_cov`, `gr_prec`, `car` and `spde`, whose fixed across-level
matrices make the whole block one unit), and requires each unit's rows
to lie in one cluster. That one test covers crossed effects, `mm()`
pooled levels, global smooths, `gp()`, `hsgp()`, `car()` and the SPDE
without a special case for any of them.

## Validation

1. Poisson GLM, one cluster per row: `cluster_scores()` equals
   `sandwich::estfun()` to 6.8e-6 absolute, and CR0/CR1p equal
   `sandwich::vcovHC()` HC0/HC1 to 2.2e-6 relative. What is left is the
   gap between frmtmb's and `glm()`'s optima, not the estimator.
2. Gaussian LMM against `clubSandwich::vcovCR(lmer, type = "CR0"/"CR1")`
   on the matched `REML = FALSE` fit (fixed effects agree to 9.3e-7).
   clubSandwich CONDITIONS on the variance parameters, so the exact
   comparison is against the conditional form rebuilt from
   `cluster_scores()` and the fixed-effect block of the information:
   2.2e-4 relative, again the two optimizers' disagreement about
   `theta`. The shipped joint form differs by the cost of estimating
   `theta` - 2.4% on the SEs here, against the 20% by which the
   model-based SEs fall short of the robust ones on the same fit - and
   that difference is documented on the manual page rather than hidden.
3. Coverage under misspecification (`dev/sandwich/probe-coverage.R`,
   60 replicates, G = 40, m = 12): model-based 0.700, cluster-robust
   0.967 for a nominal 0.95. Tuning note - plain cluster
   heteroskedasticity barely dents the model-based interval (0.967 vs
   0.950), because `sigma` absorbs its average; what breaks it is an
   omitted random slope on a predictor that varies mostly BETWEEN
   clusters.
4. `sum_g s_g` equals the objective gradient to 1e-7 and that gradient
   is below 1e-3 at the optimum (the optimizer's own tolerance, not the
   score machinery's).
5. Group-level mixtures (`dev/sandwich/probe-mixture.R`): a masked-out
   mixture group contributes `log(sum_k pi_k) = 0` rather than dropping
   off the tape, a different route to the same exact zero, so it was
   checked separately. `sum_g s_g` matches the gradient to 1.4e-13 at
   the mixture grouping and to 7.1e-14 at a coarser one; a finer
   clustering is refused.

## Suggests added

`sandwich` and `clubSandwich`, both test-only. They are the reference
implementations gate 1 and gate 2 are measured against, and calling
them live is what keeps those comparisons from going stale (the
`helper-reference.R` G5.4b rule).
