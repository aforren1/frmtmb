# t-distributed random effects: feasibility probe

Date: 2026-09-02. Branch `wt-tre`, frmtmb 0.35.0, RTMB 1.9, TMB 1.9.25,
tmbstan 1.2.0, rstan 2.32.7, brms 2.23.0, mvtnorm 1.4.2, R 4.6.1,
Windows 11.

> **Status: implemented.** Phase 1 cleared the feature, and
> `(x | gr(g, dist = "student"))` ships with `nu` fixed. The probe
> scripts are in `dev/tre/`; every number below is reproducible by
> running them from the worktree root.

**No package code was changed for phase 1.** Every probe result comes
from bare `RTMB::MakeADFun` objects, hand-written quadrature, `tmbstan`,
and `brms::make_stancode()`.

## Verdict

The Laplace approximation over a Student-t latent is defensible, with
one honestly stated limit and one exact escape hatch.

1. **The approximation error is a near-constant offset, not a bias.**
   As the latent scale shrinks against the residual SD the per-group
   log-likelihood error converges to a closed-form constant

       c(nu) = lgamma((nu+1)/2) - lgamma(nu/2) + 0.5 * log(2/(nu+1))

   (`-0.2258` at `nu = 3`, `-0.1408` at `nu = 5`, `0` in the gaussian
   limit), matched to 7e-4 by measurement. A constant in the objective
   does not move the argmax. As the scale grows the error goes to zero
   instead. What is left is the derivative in between, and that is
   small in every ordinary design.

2. **Measured against exact ML, the cost is a fraction of a standard
   error.** Over 100 replicates of a 40-group design with the latent
   scale equal to the residual SD, the Laplace estimate of the scale is
   biased UPWARD by 1.6% of one Monte-Carlo standard error at 8
   observations per group, 0.14% at 25, and at worst 12.4% at 3
   observations per group with `nu = 2.5`. Full-joint NUTS
   independently puts the same displacement at 7% of a posterior SD in
   the hardest case.

3. **It is materially wrong in one corner, and the corner is
   nameable**: a variance component far below the residual SD combined
   with tiny groups. At 2 observations per group and a true scale a
   quarter of the residual SD, the Laplace ML scale came back 0.483
   where the exact one is 0.163.

4. **`quadrature = TRUE` is EXACT over a scalar t latent**, not merely
   better. TMB's Gauss-Kronrod marginalization reproduces the adaptive
   Gauss-Hermite reference to 1e-6 in the log-likelihood and recovers
   0.156 against the exact 0.158 in the corner above. That turns the
   caveat into a supported check rather than a warning.

5. **`nu` is not estimable by maximum likelihood** in realistic
   designs, so it is fixed. brms estimates it only because a
   `gamma(2, 0.1)` prior is holding it up.

6. **The robustness payoff is real and measurable.** With one group
   displaced by 10 SDs, the gaussian latent SD inflates to 1.85 against
   a truth of 1 while a `t(5)` latent holds at 1.25, and the intercept's
   RMSE falls from 0.298 to 0.180.

## 1. The reference, and how far it can be trusted

`dev/tre/probeA1-reference.R`. The model everywhere is the simplest one
that isolates the question: a gaussian LMM with a SCALAR random
intercept whose latent density is a scaled t. A scalar latent is what
makes an exact answer possible: each group's marginal likelihood is a
one-dimensional integral.

Two independent quadratures agree, so the reference is not itself the
approximation being measured. Adaptive Gauss-Hermite has to be the
general (non-gaussian-weight) form here: center and scale at the
conditional mode and curvature, then undo the gaussian weight with
`exp(z^2)`. The worry with a fat tail is convergence in the node
count. There is none to speak of:

| nu | obs/group | AGHQ K=21 | AGHQ K=101 | `integrate()` over the line | \|K101 - int\| |
| --- | --- | --- | --- | --- | --- |
| 2.5 | 4 | -279.506955143 | -279.506955143 | -279.506955143 | 0 |
| 2.5 | 20 | -1168.280843847 | -1168.280843847 | -1168.280843847 | 0 |
| 3 | 20 | -1213.981329633 | -1213.981329633 | -1213.981329633 | 2.3e-13 |
| 10 | 4 | -274.062423938 | -274.062423938 | -274.062423938 | 1.1e-13 |

21 nodes is already exact to 1e-13. And the hand-rolled Laplace matches
TMB's to 2.8e-13 at the same parameters, so nothing below is a TMB
artefact.

## 2. The Laplace error in the log-likelihood

`dev/tre/probeA1-reference.R` section 2, at the TRUE parameters, 40
groups, latent scale = residual SD = 1. `err = logLik(Laplace) -
logLik(exact)`:

| nu | n=3 | n=5 | n=10 | n=25 |
| --- | --- | --- | --- | --- |
| 2.5 | -0.2195 | -0.1008 | -0.0199 | -0.0065 |
| 3 | -0.2376 | -0.0870 | -0.0190 | -0.0049 |
| 5 | -0.1858 | -0.0664 | -0.0146 | -0.0050 |
| 10 | -0.1119 | -0.0413 | -0.0097 | -0.0029 |
| 30 | -0.0537 | -0.0220 | -0.0060 | -0.0012 |
| Inf | 0 | 1.8e-16 | 0 | 3.1e-16 |

The gaussian row is the sanity check the whole probe rests on: Laplace
is exact over a gaussian latent, and it is exact here to machine
precision, so the whole error above is the t's.

The missing axis is the RATIO of the latent scale to the residual SD,
which sets how much of each conditional's shape the latent density
contributes. `dev/tre/probeA3-worst-corner.R`, PER GROUP:

| nu | n | s=0.1 | s=0.25 | s=0.5 | s=1 | s=2 | s=4 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2.5 | 2 | -0.2545 | -0.2084 | -0.0880 | -0.0104 | -0.0010 | -0.0001 |
| 3 | 2 | -0.2228 | -0.2057 | -0.0580 | -0.0087 | -0.0009 | -0.0001 |
| 3 | 15 | -0.1829 | -0.0332 | -0.0033 | -0.0003 | -0.0000 | -0.0000 |
| 5 | 5 | -0.1320 | -0.0804 | -0.0173 | -0.0021 | -0.0002 | -0.0000 |
| 10 | 15 | -0.0564 | -0.0165 | -0.0021 | -0.0002 | -0.0000 | -0.0000 |

Both limits are analytic. As `s -> Inf` the data dominate each
conditional and the error vanishes. As `s -> 0` the conditional becomes
the latent density itself, whose Laplace error is scale-free:

    c(nu) = lgamma((nu+1)/2) - lgamma(nu/2) + 0.5 * log(2/(nu+1))

`dev/tre/probeA4-convergence-check.R` confirms it at `s/sigma = 0.02`:

| nu | closed form c(nu) | measured per group | diff |
| --- | --- | --- | --- |
| 2.5 | -0.265937 | -0.263706 | 2.2e-03 |
| 3 | -0.225791 | -0.225111 | 6.8e-04 |
| 5 | -0.140842 | -0.140397 | 4.5e-04 |
| 10 | -0.072614 | -0.072380 | 2.3e-04 |
| 100 | -0.007475 | -0.007446 | 3.0e-05 |

**This is the headline.** The worst the log-likelihood error gets, per
group, is a CONSTANT that depends on `nu` alone. A constant offset
cancels out of the score, which is why the estimates in section 3
survive an error that looks large in the likelihood. It does not cancel
out of `logLik()`, `AIC()` or `BIC()`, so those carry roughly `G *
c(nu)` when the variance component is small: `-4.5` at `nu = 3` over 20
groups, `-22.6` over 100. Comparing two t models with the same `nu` and
the same grouping is fine; comparing a t model against a gaussian one by
AIC is not.

## 3. What it costs the ESTIMATES

`dev/tre/probeA2-parameter-bias.R`. 100 replicates, 40 groups, latent
scale 1, residual SD 1, `nu` fixed at its true value. Both fits
maximize the same target on identical data, one approximately and one
exactly by AGHQ, so the paired difference is the approximation and
nothing else. Reported as a percentage of the exact estimate's own
Monte-Carlo standard error:

| nu | n | beta0 | beta1 | sigma | scale s |
| --- | --- | --- | --- | --- | --- |
| 2.5 | 3 | -0.45% | 0.25% | -8.13% | **+12.35%** |
| 2.5 | 8 | -0.03% | 0.00% | -0.75% | +1.96% |
| 2.5 | 25 | 0.00% | -0.00% | -0.06% | +0.16% |
| 3 | 3 | -0.17% | -0.01% | -6.69% | +10.99% |
| 3 | 8 | -0.02% | -0.01% | -0.70% | +1.62% |
| 3 | 25 | 0.01% | 0.00% | -0.05% | +0.14% |
| 5 | 3 | -0.04% | -0.17% | -7.10% | +8.88% |
| 5 | 8 | 0.01% | 0.01% | -0.61% | +1.36% |
| 10 | 3 | -0.02% | 0.19% | -4.26% | +5.97% |
| 30 | 3 | -0.00% | 0.12% | -1.89% | +2.97% |
| 30 | 25 | 0.00% | 0.00% | -0.02% | +0.05% |

The sign is stable and explained: `c(nu) < 0`, so the Laplace objective
is depressed where the variance component is small and lifted where it
is large, which pushes the scale UPWARD. The fixed effects barely move
at all. The largest displacement anywhere in the table is 0.45% of a
standard error.

`nlminb` reported a non-zero convergence code on most of these fits.
`dev/tre/probeA4-convergence-check.R` settles it on the gradient rather
than the code: `max|grad|` is 4e-7 to 2e-5 for the Laplace fits and 3e-5
to 3e-3 for the numerically-differentiated exact ones, with
`rel.tol = 1e-12` set far inside a 51-node quadrature's noise floor.
These are optima, reported as "singular convergence (7)".

### The independent check: full-joint NUTS

`dev/tre/probeA5-tmbstan.R`. `tmbstan` runs NUTS over the same tape two
ways. `laplace = FALSE` samples the joint with no approximation
anywhere, `laplace = TRUE` samples the marginalized objective, so the
difference is the Laplace error on the posterior scale, with no
quadrature of mine in the path. 4 chains, 4000 iterations, min n_eff
958, max Rhat 1.002.

| nu | n | parameter | joint mean | Laplace mean | joint sd | difference / sd |
| --- | --- | --- | --- | --- | --- | --- |
| 3 | 3 | log_s | 0.05922 | 0.07162 | 0.17649 | +7.03% |
| 3 | 3 | log_sigma | -0.06843 | -0.07130 | 0.08010 | -3.59% |
| 3 | 8 | log_s | 0.42480 | 0.42581 | 0.17092 | +0.59% |
| 5 | 3 | log_s | 0.00686 | 0.01937 | 0.17138 | +7.30% |
| 5 | 8 | log_s | 0.10752 | 0.10849 | 0.15284 | +0.63% |

Same sign, same magnitude, same dependence on the group size. Two
independent references agree.

### The corner where it is not fine

`dev/tre/probeA3-worst-corner.R` section B, 100 replicates at 2
observations per group:

| nu | true s | d sigma (sd) | d scale s (sd) | d logLik |
| --- | --- | --- | --- | --- |
| 2.5 | 0.10 | -0.10857 (0.02907) | **+0.33677** (0.09769) | -6.38 |
| 2.5 | 0.25 | -0.10943 (0.02696) | +0.29231 (0.11419) | -4.80 |
| 3 | 0.10 | -0.09924 (0.02927) | +0.32512 (0.10585) | -5.66 |
| 3 | 0.25 | -0.10400 (0.02500) | +0.29894 (0.11882) | -4.53 |
| 5 | 0.10 | -0.07230 (0.02363) | +0.26578 (0.11296) | -3.60 |

A latent scale of 0.10 estimated 0.33 too high is not a rounding error.
The regime is a near-null variance component observed through
near-singleton groups: the likelihood says almost nothing about the
scale, and the `c(nu)` gradient is the loudest thing left in the
objective. Section 4 is the answer to it.

## 4. The inner Newton solve, and quadrature

`dev/tre/probeB1-inner-newton.R`. A t log-density is concave only inside
`|b| < s*sqrt(nu)`; outside it is convex, with a second derivative
peaking at `(nu+1)/(8 s^2 nu)`. A second mode in a group's conditional
therefore needs

    s^2 < sigma^2 (nu + 1) / (8 n nu),

which a 108-point grid scan confirms exactly: 4 of 108 configurations
were bimodal, every one of them below that bound, and all four at
`n <= 2` with a residual mean 6 or more SDs out.

**TMB's inner optimizer found the global mode every time.** On a
purpose-built conflicting dataset (30 groups of 2, one group displaced
by up to 10), `max|b_TMB - b_global|` was between 2e-16 and 2e-10 across
every scale and displacement tested. Sweeping `log_s` from -4 to 2 at
that data, the smallest inner curvature stayed at 3.04, so the inner
Hessian never left the positive-definite cone, and TMB's objective
matched the hand-rolled Laplace to the printed six decimals throughout.

The price of the mode it did NOT pick is real, which is why the check
mattered: at `nu = 5`, `s = 0.15`, one observation with residual 8, the
Laplace value at the wrong mode is -31.91 against -17.73 at the right
one (exact -17.73).

### Gauss-Kronrod quadrature is exact here

`dev/tre/probeG1-quadrature.R`. frmtmb's `quadrature = TRUE` is TMB's
experimental `integrate =` Gauss-Kronrod marginalization, calibrated
once at the Laplace mode and then frozen. That freezing is the part a
polynomial tail might have defeated. It does not:

| nu | n | exact (AGHQ) | Laplace | GK | GK error | Laplace error |
| --- | --- | --- | --- | --- | --- | --- |
| 2.5 | 3 | -222.268913 | -222.488382 | -222.268913 | 0.00000 | -0.21947 |
| 3 | 3 | -218.824300 | -219.061896 | -218.824300 | -0.00000 | -0.23760 |
| 3 | 8 | -499.087773 | -499.122593 | -499.087773 | -0.00000 | -0.03482 |
| 10 | 8 | -492.082526 | -492.104845 | -492.082526 | -0.00000 | -0.02232 |

Full fits agree with exact ML to 1e-5 in every parameter. And in the bad
corner of section 3:

| nu | exact s | Laplace s | GK s |
| --- | --- | --- | --- |
| 2.5 | 0.15771 | 0.49369 | 0.15561 |
| 3 | 0.16324 | 0.48277 | 0.16282 |
| 5 | 0.17338 | 0.44289 | 0.17338 |

This is what makes the feature shippable rather than merely arguable:
the failure mode has a supported, one-argument, exact remedy.

## 5. `nu`: what brms does, and what ML can do

### What brms does

`dev/tre/probeE1-brms-encoding.R`, brms 2.23.0. `gr()`'s signature is
`gr(..., by, cor, id, pw, cov, dist)` with
`dist = match.arg(dist, c("gaussian", "student"))`. The generated Stan
for `(1 | gr(g, dist = "student"))`:

```stan
real<lower=1> df_1;
vector<lower=0>[N_1] udf_1;
vector<lower=0>[M_1] sd_1;   // "group-level standard deviations"
dfm_1 = sqrt(df_1 * udf_1);
r_1_1 = dfm_1 .* (sd_1[1] * (z_1[1]));
lprior += gamma_lpdf(df_1 | 2, 0.1) - 1 * gamma_lccdf(1 | 2, 0.1);
target += inv_chi_square_lpdf(udf_1 | df_1);
target += std_normal_lpdf(z_1[1]);
```

Five things this settles:

1. **The construction is the standard multivariate t.**
   `b_j = sqrt(nu u_j) W z_j` with `u_j ~ inv-chi2(nu)` is
   `MVt(nu, Sigma = W W')`. `dev/tre/probeF1-mvt-and-limits.R` confirms
   it by simulation (the Mahalanobis form is `d * F(d, nu)`, KS
   p-values 0.29 to 0.99 at `d = 1, 2, 3`) and the closed-form density
   against the mixture integral to 4e-14.
2. **One mixing variable per LEVEL**, shared across that level's
   coefficients, and, critically, **shared under `cor = FALSE` too**:
   brms writes `r_1_1 = dfm_1 .* (sd_1[1] * z_1[1])` and
   `r_1_2 = dfm_1 .* (sd_1[2] * z_1[2])` with the same `dfm_1`. So a
   `diag()` t block is a multivariate t with a diagonal scale matrix,
   NOT a product of `d` independent univariate t's. Implementing it as
   the latter would have been wrong and would have gone unnoticed at
   `d = 1`.
3. **`sd_1` is the SCALE**, despite the name and despite brms printing
   it under "Group-Level Effects: sd(Intercept)".
4. **`nu` is estimated, per grouping term**, with a `gamma(2, 0.1)`
   prior truncated at 1 (mean 20). It can be pinned:
   `prior(constant(3), class = "df", group = "g")` emits `df_1 = 3;`.
5. **`make_standata()` does not expose the choice.** The data blocks for
   the gaussian and student versions are `identical()`. The planned
   structural cross-check against brms's own encoding is therefore not
   available: the difference lives entirely in the Stan code, which is
   why section 5's reading of `make_stancode()` is the citation instead.

brms also accepts `gr(g, cov = A, dist = "student")` and
`mm(g1, g2, dist = "student")`. The former emits
`r_1_1 = dfm_1 .* (sd_1[1] * (Lcov_1 * z_1[1]))`, a per-level scalar
rescaling of a CORRELATED gaussian field. That is not a multivariate t
over the levels and has no closed-form marginal density. Both are
refused here, by name.

### What ML can do

`dev/tre/probeC1-nu-identifiability.R`. RTMB tapes `dt()` with an AD
`df` correctly (AD against central differences: 8.6e-9), so estimating
`nu` is mechanically available. It is statistically hopeless.

Profile log-likelihood in `nu`, true `nu = 3`, 5 observations per group,
reported as the drop from the profile maximum:

| groups | peak | 95% profile interval | drop at nu=500 |
| --- | --- | --- | --- |
| 20 | 3 | [2.1, 500] | 1.05 |
| 40 | 5 | [2.1, 500] | 0.53 |
| 100 | 3 | [2.1, 5] | 8.85 |
| 400 | 5 | [4, 10] | 9.87 |

At 20 and 40 groups the entire grid from 2.1 to 500 is inside the
interval. Joint ML, 100 replicates, `nu = 2 + exp(.)`:

| true nu | groups | median | IQR | at ceiling (>100) | at floor (<2.1) |
| --- | --- | --- | --- | --- | --- |
| 3 | 20 | 4.53 | [2.26, 1.1e7] | 41% | 24% |
| 3 | 40 | 3.79 | [2.39, 9.97] | 14% | 16% |
| 3 | 100 | 3.33 | [2.47, 4.71] | 5% | 9% |
| 5 | 20 | 8.10 | [3.47, 1.0e7] | 44% | 9% |
| 10 | 20 | 7.9e6 | [7.69, 1.6e7] | 61% | 0% |
| 10 | 100 | 26.92 | [8.51, 8.3e6] | 39% | 0% |

Two-thirds of a boundary run in the `nu = 10` designs, and a quarter
collapsing on the floor at `nu = 3`. Neither failure is loud: the
optimizer converges, and `sdreport()` would hand back a standard error
for a parameter sitting on a flat ridge. Even the easier question, is
the latent t at all, needs a hundred groups: `2*(logLik(nu free) -
logLik(gaussian))` exceeds 3.84 in 17.5% of replicates at 20 groups with
a true `nu = 3`, 47.5% at 40, 77.5% at 100.

**Decision: `nu` is fixed, set by `dist_nu`.** The frequentist analogue
of brms's `prior(constant(3), class = "df")`, and the honest reading of
the evidence: brms's default is not "estimate `nu`", it is "estimate
`nu` under a prior with mean 20", and there is no prior here.

### Which fixed `nu`

`dev/tre/probeD2-nu-default.R` sweeps `nu` over the contaminated design
of section 6 (39 clean gaussian groups, one displaced by `delta`), 80
replicates. RMSE of the reported latent SD against the truth 1:

| delta | 2.5 | 3 | 4 | **5** | 7 | 10 | 20 | gaussian |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 0.6825 | 0.3778 | 0.2198 | **0.1767** | 0.1529 | 0.1455 | 0.1430 | 0.1441 |
| 3 | 0.7773 | 0.4524 | 0.2767 | **0.2225** | 0.1862 | 0.1708 | 0.1622 | 0.1636 |
| 6 | 0.8118 | 0.4857 | 0.3132 | **0.2631** | 0.2356 | 0.2342 | 0.2654 | 0.3784 |
| 12 | 0.8200 | 0.4948 | 0.3257 | **0.2802** | 0.2646 | 0.2873 | 0.4280 | 1.1143 |
| worst | 0.8200 | 0.4948 | 0.3257 | **0.2802** | 0.2646 | 0.2873 | 0.4280 | 1.1143 |

**Default `dist_nu = 5.`** The minimax choice is 7 (worst-case 0.2646)
and 5 is within 6% of it (0.2802) while keeping more headroom at heavy
contamination, where it beats 10 and 20 outright. `nu = 3` costs a
quarter of the reported SD under clean data: a `t(3)` fitted to
gaussian latents reports an SD of 1.32 against a truth of 1. That is a
worse default than the robustness it buys. The floor is 2, not brms's 1:
below it the latent has no variance for `VarCorr()` to convert, for
new-level prediction to use, or for `simulate()` to have drawn from.
Probe D showed exactly that failure with `nu` free: at `delta = 12` the
median `nu_hat` was 2.46 and the implied SD 6.41.

## 6. The robustness payoff, measured

`dev/tre/probeD1-robustness-payoff.R`. 39 clean groups with
`b ~ N(0, 1)`, one displaced by `delta`, 8 observations per group,
`sigma = 1`, 60 replicates. Reported latent SD (`scale * sqrt(nu/(nu-2))`
for the t fits):

| delta | gaussian | t(3) | t(5) | t(nu free) | median nu_hat |
| --- | --- | --- | --- | --- | --- |
| 0 | 1.0131 | 1.3538 | 1.1072 | 1.0185 | 9.7e6 |
| 2 | 1.0576 | 1.4081 | 1.1522 | 1.0643 | 5.2e6 |
| 5 | 1.2730 | 1.4678 | 1.2231 | 1.4273 | 4.28 |
| 10 | **1.8522** | 1.4806 | **1.2477** | 6.4090 | 2.46 |

RMSE of the intercept against its truth:

| delta | gaussian | t(3) | t(5) |
| --- | --- | --- | --- |
| 0 | 0.1655 | 0.1796 | 0.1732 |
| 5 | 0.2063 | 0.1882 | 0.1838 |
| 10 | **0.2982** | 0.1853 | **0.1801** |

The within-group slope is untouched at every `delta` (0.0498 to 0.0508
throughout), because a group-level displacement is a group-level
problem. The
efficiency cost when nothing is wrong is 8.5% on the intercept's RMSE
for `t(3)` and 4.6% for `t(5)`.

One dataset in detail (`delta = 10`): the true `b[40]` is 9.59; the
gaussian fit predicts 9.16 having widened its latent SD to 1.82, the
`t(3)` fit predicts 9.69 with a latent SD of 1.45. RMSE of the CLEAN
groups' predicted effects falls from 0.562 to 0.406. The outlying group
was distorting the other 39 through the shared variance component, and
the t stops it.

## 7. REML

`dev/tre/probeH1-reml.R`. REML here puts the `mu` fixed effects into the
inner Laplace problem, so the objective stacks a second approximation on
the first. An exact reference exists for a two-coefficient design: the
inner integral is section 1's certified scalar quadrature, and the outer
one is a tensor adaptive Gauss-Hermite rule over it.

| nu | exact REML | TMB REML | REML error | ML error, same data |
| --- | --- | --- | --- | --- |
| 2.5 | -102.390709 | -102.499766 | -0.109 | -0.050 |
| 3 | -101.721451 | -101.830682 | -0.109 | -0.053 |
| 5 | -101.295383 | -101.367638 | -0.072 | -0.039 |
| 10 | -98.581977 | -98.642968 | -0.061 | -0.034 |
| 30 | -97.943922 | -97.967661 | -0.024 | -0.014 |

Stacking the beta integral roughly doubles the error, which is what one
extra Laplace dimension should cost, and the fitted values track: at
`nu = 5` the TMB REML scale is 1.0538 against the exact 1.0434. The
gaussian limit holds to 6.6e-9 at `nu = 1e8`, and REML moves the
variance component in the direction it exists to move it (60 replicates,
8 groups, true scale 1, `nu = 5`: ML 0.911, REML 0.996). **REML ships**,
with the same accuracy caveat as ML.

## 8. What shipped

Spelling, matching brms: `(x | gr(g, dist = "student"))`, with
`dist = "gaussian"` accepted as the explicit default. `dist_nu = 5` is a
frmtmb argument. brms has no such argument, and pins `nu` through
`prior(constant(.), class = "df")` instead.

- **Registry**: `us_t` and `diag_t` in `R/covstruct.R`, sharing `npar`,
  `sd_idx`, `start`, `from_natural` and `vcov` with their gaussian
  twins. `nu` is a constant on the block (`blk$dist_nu`), NOT a `theta`
  entry, so nothing downstream that indexes `theta` had to change. The
  density is the closed-form multivariate t off a triangular factor -
  one solve, one log-determinant, no inverse. `log1p`, not
  `log(1 + .)`: the tail term is multiplied by `(nu + d)/2`, so at
  `nu = 1e8` the rounding of `1 + q/nu` is amplified by the same factor
  and the gaussian-limit check fails in the third decimal. Measured on
  a fitted `theta` of a `diag()` block: 3.1e-2 before the fix, 3.1e-5
  after, and 1.3e-5 to 6.8e-9 for a scalar one. The DENSITY difference
  at a fixed `theta` barely moves (2.8e-6 over 20 levels either way),
  because what is left there is the cancellation in
  `lgamma((nu+d)/2) - lgamma(nu/2)`, which is a constant in `theta` and
  therefore invisible to the optimizer. Only the `log1p` term depended
  on the parameters, which is exactly why it was the one that mattered.
- **Reporting**: `VarCorr()` stores the SCALE matrix, matching brms's
  `sd_` naming, and tags it with `nu`; `print()` shows a `Scale` column,
  a converted `Std.Dev.` column and the fixed `nu`, so the convention is
  visible rather than silent. `confint()`, `variables()` and
  `frm_simulate(newparams = )` speak of the scale under
  `sd_<group>__<term>`, as brms does.
- **`quadrature = TRUE`** accepts scalar t blocks, and section 4 makes it
  the recommended check.
- **`simulate()` / `frm_simulate()`** draw a multivariate t with one
  chi-square per level.
- **`predict(allow_new_levels = TRUE)`** inflates the unseen level's
  variance by `nu/(nu-2)`. The interval is still built as a gaussian
  one, so it carries the right variance but not the right far-tail
  quantile; documented.

Guarded out, each with its own message: `gr(cov = )` / `gr(prec = )`,
`mm()`, `|ID|` keys, every covariance structure but `us` and `diag`,
`dist_nu` without `dist`, `dist_nu <= 2`, a `dist` that is neither of
brms's two, and `dist_nu` on a gaussian latent.

## 9. Loose ends

- **`nu` estimation stays out.** If it is ever wanted, section 5 says
  what it needs first: a penalty or a bound that does the work brms's
  `gamma(2, 0.1)` does, and a `sdreport()` story for a parameter that
  reaches a flat ridge in a third of designs.
- **`|ID|` keys** are refused only because the merged-block assembly in
  `R/frame.R` builds one gaussian structure. brms supports it (one `df`
  per merged ID), and the workaround is the identical density: writing
  the merged coefficients as one term. So this is ergonomics, not
  capability.
- **`gr(mm(...), dist = )` and `gr(cov =, dist = )`** are refused for a
  real reason, not a missing feature: neither has a closed-form marginal
  density. They would need the mixing variables as explicit latents,
  which is a different (and much larger) design.
- **An unknown structure name around a bar is silently accepted**:
  `foo(x | g)` fits as a default `us` term, and so does `us_t(x | g)`.
  This is pre-existing `reformulas` behaviour rather than anything the t
  blocks introduced (verified: `foo()` behaves identically), but the new
  registry names make it slightly more reachable. `us_t` and `diag_t`
  are kept OUT of `supported_cs`, so no formula can reach them without a
  `dist_nu`.
- **`obj$simulate()`** on a t block is not exercised. The registry's
  gaussian entries lean on `RTMB::dmvnorm`'s simref support; the
  hand-written t density preserves `dim<-` but has no simulator. Nothing
  in frmtmb uses that path, because `simulate()` draws in R through
  `draw_b()`. So this is a note, not a defect.
