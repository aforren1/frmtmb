# Extension gaps: the plan

Date: 2026-09-08, re-read against 0.55.0 the same day. frmtmb 0.55.0;
frmtmb.eam 0.5.1, frmtmb.sample 0.3.2, frmtmb.latent 0.2.1,
frmtmb.learn 0.2.1, frmtmb.ode 0.1.1, frmtmb.spline 0.4.0,
frmtmb.coupling 0.2.0.

Source: a survey of `extensions/` against what a practitioner in each
field needs for a paper. The survey found that every extension is a
correct likelihood with a thin workflow around it, and that the same
four things block real use in every package: no evidence at realistic
scale, hierarchy exercised only at the one-intercept level, no path for
real data in, and a post-fit surface that refuses more than it answers.
This document orders the work to close that, phase by phase, with a
check and a size for each item.

Scope: the extensions. An item that cannot be built without a change
under core's `R/` is not built in an extension; it is listed once in
"Core seams" with the extension item that waits on it. No extension
reaches core with `:::`.

## Rules every item is held to

1. **Measure before building.** Phase 0 is a gate. Any item whose
   cost the survey could only guess at is re-ranked after Phase 0's
   numbers, and this document is edited to say so.
2. **The shipping bar for a family or a method** is the one the
   extensions already hold themselves to: an identity against an
   independent reference at the same estimates, a recovery table at
   the realistic scale in the table below, refusals by name for what
   it does not do, compat rows for what it composes with, and a timing
   taken in a fresh process with both arms warmed and the seed count
   stated (`extensions/frmtmb.eam/dev/dev-findings.md`, "How this
   package times things"). Since 0.55.0 it also includes the hazard-container guard: a new package or a new `R/` file is covered by `frm_hazard_reads()` in that package's own `test-bracket-access.R`, because a `$` read on a hazard container used to reach the release tally before anything caught it. And a timing is not a timing until the instrument has been checked: two claims in the 0.55.0 round measured below `proc.time()`'s 10 ms tick or outside the hot path, so interleave the arms in one process and carry a control that must report 1.0.
3. **A silent wrong answer outranks a refusal, and a refusal outranks
   a missing feature.** Where a phase has both, the silent one goes
   first.
4. **A data-shaping helper is not model machinery**, and that is a
   reason to write it, not to leave it out: it is where a user with
   real data stops.

### Realistic scale, per field

These are the designs Phase 0 times and Phase 2 recovers at. They are
what a typical paper in the field collects, not a stress test.

| package | design | rows |
|---|---|---|
| eam | 30 subjects x 400 trials, 2 conditions, random effects on drift, boundary and non-decision time | 12,000 |
| learn | 100 subjects x 200 trials, correlated random effects on every parameter | 20,000 |
| latent, hmm | 50 sequences x 500 steps, K = 3, random effect on a transition | 25,000 |
| latent, lca | n = 2000, 10 items, K = 4, two covariates on membership | 2,000 |
| ode | 100 subjects x 8 samples, a depot and a central state, twice-daily dosing for 7 days | 800; 14 dose events per subject and, with one `ss` row at the default `n_ss = 20`, about 35 solves per group per gradient |
| spline, curves | 40 subjects, per-subject `fs` curves, 200 points each | 8,000 |
| spline, RP | 2000 subjects, 40 percent censored, one frailty | 2,000 |
| coupling | 40 subjects x 2 conditions x 60 frequencies, smooth in frequency plus random effects | 4,800 |
| sample | a 2000-row GLMM with two crossed factors, 4 chains x 2000 | 2,000 |

## Phase 0. The measurement tier

DONE, 2026-09-08. One fit per package at the scale above, in a gated
test tier (`FRMTMB_SCALE_TESTS=true`) that lives in each package's
`tests/testthat/test-scale.R`, recording tape build, one gradient, the
whole `frm()` call including `sdreport()`, peak memory and the estimate
against the truth. The numbers are in `dev/scale-findings.md` as one
table, how the run was made is in `dev/scale-lane-notes.md`, and each
package's `dev/scale-row.md` points at its own row.

| package | fit | decided |
|---|---|---|
| eam | wiener with the three random effects; a second run with `variability = "sv"`; a third on the same data with the ndt bound lifted | MINUTES, not hours: 309 s over three passes. The bounded ndt link does NOT tolerate a random effect. 20 of the 30 subjects have a true `ndt` above `min(rt)`, which is the bound; the plain arm does not converge (gradient 1.25e11, seven NaN standard errors) and the `sv` arm converges to the link's edge, coming back PINNED 0.42 of a standard error below the bound with a standard error of 7.2e-06 on a population `ndt` that is wrong by ten percent, with nothing from `diagnose()`. Lifting the bound on the same data recovers everything and finds a log-likelihood 121.4 units higher at the same parameter count |
| learn | bandit2arm_delta with `(1 \| p \| id)` on both parameters, then rlddm the same way | the recursion scales: the delta learner is 7.0 s at 20,000 rows and recovers. rlddm is 12.3 minutes over two passes, 738 s against a 600 s rule, and fails the eam failure because it takes its diffusion from frmtmb.eam |
| hmm | K = 3 gaussian with `tr12 ~ (1 \| id)` | the post-fit passes are NOT a cost: over three interleaved rounds `hmm_probs()` is 0.10 s and `hmm_viterbi()` 0.06 s against a 118 s fit. The tape build is 0.84 s at 25,000 rows, near where `dev/hmm-feasibility.md` put it |
| lca | K = 4 | cheap, as expected: 0.54 s |
| ode | the design above with `ii`/`addl` and one `ss` row | NOT usable at population scale: 3948 s, 65.8 minutes, on 800 rows, with the truth recovered. The closed-form path is a prerequisite, not an option |
| spline | `s(t) + s(t, subject, bs = "fs")` at 40 subjects, then `frm_curve()` and `frm_curve_feature()` on it | the curve surface is NOT the cost, and it decomposes into a memoized joint-precision solve paid ONCE per fit by whichever curve call comes first, 0.48 s, and cheap calls after it: the simultaneous band 0.071 s, the feature search 0.0098 s, the pointwise band 0.0031 s, against a 6.0 s fit. The feature search's own cold cost is 0.40 s when it is the call that pays the solve. The number to carry forward is memory: 1.32 GB of process working set |
| coupling | five models on one 4,800-row design | seconds, not minutes, across the whole ladder. The survey's benchmarks ARE recoverable, from its own session transcript: its 52.5 s and 28.9 s fits are the ladder's top two rungs, so the promotion this row asked for is done, though the survey ran 0.53.0 under `load_all()` and this runs 0.55.0 installed so no speedup is claimed. Both rungs recover the contrast, 0.486 (0.451, 0.520) and 0.488 (0.415, 0.561) against 0.5, so the survey's 0.77 was noise. What the cheap rungs get wrong is the interval WIDTH, and adding `(1 \| id)` does not fix it: the contrast is within-subject, so only `(1 \| id:cond)` widens it, and without that term the interval is 53 percent too narrow |
| sample | 4 x 2000 on the GLMM, then `posterior_epred()` and `loo()` on the draws | the per-draw loops are NOT the cost: `posterior_epred()` over 4000 draws and 2000 rows is 1.01 s against 120 s of sampling, and the three calls item 4.6 rewrites are 2.2 percent of the time it took to produce the draws. Phase 4's caching is not needed before anything else |

Decision rule: a fit past ten minutes on its design moves that
package's cost item to the front of Phase 1. TWO packages breached:
frmtmb.ode at 3948 s, comfortably, and frmtmb.learn at 738 s through `rlddm()`, by only 1.23x once that row was reproduced as the procedure requires.
What moved is in Phase 1 below, and why is in
`dev/scale-findings.md`, "The decision rule, applied". Actual size: one
day of work and about four hours of runs.

## Phase 1. Small, high leverage

Each of these is under two days and removes a wall a user hits on day
one, EXCEPT the three Phase 0 moved here, which are larger and are at
the top because nothing behind them is worth doing first.

| item | what | where | check | days |
|---|---|---|---|---|
| 1.0a | The per-subject non-decision-time bound in eam, moved here from item 3.4 by Phase 0. `ndt`'s scaled logit is bounded by `min(y)`, the GLOBAL fastest response, so a subject deviation on it is a deviation on a fraction of one subject's floor: at the plan's own eam design 20 of 30 subjects have a true `ndt` above that bound and none is above its OWN fastest response. `ndt` becomes a fraction of the group's floor inside `lpdf`, and the help page says that is what `ndt`'s linear predictor means | `frmtmb.eam/R/ddm-shared.R`, `wiener-family.R` | the eam scale row converges with a positive definite Hessian and recovers `sd(ndt)`; the log-likelihood is at least the -7027.4 the bound-lifted arm reaches, against -7148.8 today. That target is reachable under a per-subject bound: at that optimum every subject's `ndt` is below its own fastest response with 22 ms to spare | 1.5 |
| 1.0b | The same bound in rlddm, moved here from item 3.3 by Phase 0. `rlddm()` takes its diffusion parameterization from frmtmb.eam and inherits the failure: at 100 x 200 the fit runs 12.3 minutes to a gradient of 6.36e9 with `sd(ndt)` at 2.24 on its link. It moves with 1.0a under Rule 3 rather than on its clock, whose margin over the ten-minute rule is only 1.23x. The rest of item 3.3, which is `session =`, stays in Phase 3 | `frmtmb.learn/R/rlddm.R` | the learn-rlddm scale row converges with a positive definite Hessian and no NaN standard errors | 0.5 |
| 1.0c | `frm_lincmt()`, moved here from item 5.4 by Phase 0. The plan already said Phase 0 would decide whether 5.4 is Phase 5 or Phase 1; the design's `frm()` call is 3948 s, so it is Phase 1. Details and check are unchanged from 5.4 | `frmtmb.ode/R/lincmt.R` | identity with `frm_ode()` on the same schedule to 1e-8; the Phase 0 design's wall clock, before and after. The before is 3948 s | 3 |
| 1.1 | `frm_sample(stan_control = list(adapt_delta =, max_treedepth =))`, passed to tmbstan's `control`. `control` stays `frmtmb_control()`; the new name is explicit so a brms user's `control = list(adapt_delta = )` keeps refusing by name and the refusal names `stan_control` | `frmtmb.sample/R/sample.R` | divergence count on a centered funnel falls between 0.8 and 0.99; `nuts_params()` reads the tightened step size back | 0.5 |
| 1.2 | A pre-flight refusal for a fit the sampler cannot run. `frm_sample()` reads the compat registry for rows whose feature is `frm_sample` and whose status is `refused`, and stops before taping. `frmtmb.ode` registers that row until a patched RTMBode is released, and the three patches in `frmtmb.ode/dev/upstream/` get filed | `frmtmb.sample/R/sample.R`, `frmtmb.ode/R/zzz.R` | an ODE fit refuses in under a second with no Stan call; the eam and latent rows stay `conditional` | 1 |
| 1.3 | A difference curve. `frm_curve(object, newdata, contrast = newdata2)` returns `A1 - A2` with `(A1 - A2) V (A1 - A2)'`, pointwise and simultaneous, and `frm_curve_feature()` on it locates where the difference crosses zero. Everything it needs is already assembled in `curve-cov.R` | `frmtmb.spline/R/curve.R`, `curve-cov.R`, `curve-feature.R` | identity with `gratia::difference_smooths()` on the same mgcv fit; simultaneous-band coverage as a rate over seeds, the way `test-gratia.R` is posed | 1.5 |
| 1.4 | A units guard in eam. Every default assumes seconds; `valid_y` warns when the smallest response is above 20, names milliseconds as the likely cause, and says what to divide by | `frmtmb.eam/R/ddm-shared.R` | a millisecond fixture warns; a seconds fixture does not | 0.25 |
| 1.5 | `reward(pay)` with one column, for data that records the received outcome alone. The arity-1 spelling duplicates the column and `simulate()` refuses on it by name, since it cannot know what the other arm paid | `frmtmb.learn/R/zzz.R`, `bandit2arm-delta.R` | the two spellings give one log-likelihood; `simulate()` names the missing column | 0.5 |
| 1.6 | `frm_cross_spectrum()` splits at NA runs instead of refusing them, so an artifact-rejected record is one call; a `window = "hann"` option beside the sine tapers, with the degrees of freedom it costs stated in `n` | `frmtmb.coupling/R/cross-spectrum.R` | a record with two rejected spans gives the sum of the three clean pieces; Hann against `stats::spec.pgram()` on one clean segment. Note that core's `whittle()` stopped refusing Hann-tapered responses in 0.55.0, so the option no longer needs a caveat about the sibling refusal | 1 |

Phase 1 total: about 10 days. Items 1.1 through 1.6 are independent and
can run as parallel lanes. Items 1.0a, 1.0b and 1.0c are the Phase 0
moves; 1.0b waits on 1.0a, because the fix is one parameterization
shared by two packages, and 1.0c is independent of both.

The order WITHIN the three is the plan's own Rule 3, "a silent wrong
answer outranks a refusal, and a refusal outranks a missing feature".
Item 1.0a is a silent wrong answer: the `sv` arm of the eam row
converges, reports no warning, and gives a population non-decision time
that is wrong by ten percent with a standard error of 7.2e-06 on it,
because the scaled logit's derivative vanishes at the edge it has been
pushed to. Item 1.0c is a cost. So the eam and rlddm bound goes first.

## Phase 2. Hierarchy: validate what is claimed

No new features. Each row is a recovery table and a third-party
identity at the realistic scale, and the result is written into the
family's help page whatever it says.

| item | what | reference | expected trouble | days |
|---|---|---|---|---|
| 2.1 | eam: `mu ~ cond + (1 \| s)`, `bs ~ (1 \| s)`, `ndt ~ (1 \| s)`, 30 x 400, 60 replicates; `check_laplace()` on one dataset, and `diagnose()` on every replicate. WAITS ON ITEM 1.0a: at this design today the fit does not converge, so 60 replicates of it would measure the bound and not the model | hierarchical Wiener in brms on one dataset (gated tier), estimates within Monte Carlo error | Phase 0 settled the expected trouble and it is worse than expected. The ndt bound is the GLOBAL minimum response time and 20 of 30 subjects have a true `ndt` above it; the plain fit does not converge and the `sv` fit converges to the edge with a 7.2e-06 standard error on a wrong answer. 0.55.0's `unbounded_dpar` check does not fire, and NOT because it fails to look: it reads the fixed `ndt` intercept, which stands at 11.20 on the link with a link-scale standard error of 2.37, PASSES its `abs(est) > 10` half and fails the `se > abs(est)` half. Its evidence pair is calibrated for an unbounded link where running away explodes the standard error; on a bounded link the derivative vanishes at the edge and the standard error COLLAPSES instead, so both halves point the wrong way. `extreme_theta` does read variance components, but its threshold of 8 is link-agnostic and the `ndt` block's log sd is 2.05. So a replicate that recovers badly and diagnoses clean is exactly what this design produces today; the seam is filed | 1.5 |
| 2.2 | learn: `(1 \| p \| id)` on every parameter of bandit2arm_delta and rlddm, 100 x 200, 60 replicates; the importance correction's success rate at that scale, counted with 0.55.0's `imp_stalled()` rather than by eye. The rlddm half WAITS ON ITEM 1.0b; the delta half does not, and Phase 0 already showed it fits in 7 s and recovers on one seed | the family's own Stan program with the correlated block added | variance-component collapse is documented at 20 to 60 trials; this measures 200, which is where the field's designs sit. 60 replicates of rlddm at 12.3 minutes each is 12 hours, which is a second reason 1.0b comes first. A stalled correction is now distinguishable from a slow one, and the threshold that separates them has margins of only 1.54 and 2.10, so a rate near either edge is a finding about the threshold as much as about the model | 1.5 |
| 2.3 | hmm: K = 3 gaussian against depmixS4; `tr12 ~ (1 \| id)` against hmmTMB with the same random effect | depmixS4, hmmTMB | multimodality. The probe found a cold start 8.1 log-likelihood units under the optimum with every diagnostic clean. Ships with `hmm_starts(fit, n, jitter)`: n refits from jittered starting values, the best returned, the spread reported. Generalize to core later (see Core seams) | 2 |
| 2.4 | lca: K = 4 at n = 2000 against poLCA | poLCA | none expected | 0.5 |
| 2.5 | spline RP: `(1 \| centre)` frailty against rstpm2's log-normal frailty on the same data; a random effect on `gamma1` | rstpm2 | the current frailty test is a smoke test for finiteness | 1 |
| 2.6 | coupling: `coh ~ cond + s(freq, by = cond) + (1 \| id) + (1 \| id:cond)`, promoted from the scratchpad benchmark to a recovery assertion on the condition contrast, with its standard error captured this time | the simulator's truth | the survey's constructions ARE recorded, in its own session transcript, and its 0.77 came from `coh ~ cond + s(freq, by = cond) + (1 \| id)`, a subject effect WITHOUT the subject-by-condition effect the simulator has. Phase 0 runs both on one seed at the realistic design: this row's own model gives 0.488 (0.415, 0.561) with components 0.313 and 0.146 against 0.35 and 0.20, and the survey's model gives 0.486 (0.451, 0.520). Both cover 0.5, so the 0.77 was noise. What is left for this row is the replicate count, and the interval WIDTH, which is where the misspecified rungs actually go wrong: they are 53 percent narrower than the correct model's | 0.5 |

Phase 2 total: about 7 days, parallel by package.

## Phase 3. Real data in

| item | what | where | check | days |
|---|---|---|---|---|
| 3.1 | `frm_ode_records(d, id, time, evid, amt, cmt, rate, ii, addl, ss)`: one NONMEM-shaped table in, `list(data, events)` out. Pure reshaping; every unknown column refused by name | `frmtmb.ode/R/records.R` | round trip against `rxode2::et()` on the same schedule; the vignette's Theoph example rewritten through it | 1.5 |
| 3.2 | coupling ingestion: `frm_cross_spectrum()` takes a list of epochs of unequal length with a `group` vector, the way core's `frm_periodogram()` does; `frm_cross_pairs(X, pairs)` stacks channel pairs with a `pair` factor, so `(1 \| pair)` and `s(freq, by = pair)` are the multiple-comparison story, and the vignette says so | `frmtmb.coupling/R/cross-spectrum.R`, `pairs.R` | an epoch list equals the concatenated-and-split call when lengths agree; a 4-channel pairwise frame fits with `(1 \| pair)` and the per-pair coherences match four two-channel fits | 2.5 |
| 3.3 | learn: `session =` on every family, so `init` runs again at each session boundary and the importance block stays the subject. The per-subject non-decision bound that used to be the second half of this item is now item 1.0b, moved by Phase 0 | `frmtmb.learn/R/engine.R`, `family.R` | a two-session dataset equals two single-session fits with shared parameters | 1 |
| 3.4 | eam censoring, and censoring ONLY: the per-subject non-decision bound that item 2.1 said this item would gain is now item 1.0a, moved by Phase 0. `wiener()` gains `lcdf` as the log-sum of the lower-boundary distribution function in `wiener-cdf.R` and its reflection for the upper, blended the way the density is; `gddm()` reads the absorbed mass off its own grid; `lba()` uses the race identity, one minus the product of survivors. Then `cens()` and `trunc()` work on all three, and a deadline design is `rt \| dec(r) + cens(c)` | `frmtmb.eam/R/wiener-family.R`, `wiener-cdf.R`, `gddm.R`, `lba.R` | `RWiener::pwiener()` to 1e-10 over the same grid the density is pinned on; right-censored log-likelihood against the likelihood written by hand, as `test-rdm-gng.R` already does for rdm | 3 |
| 3.5 | eam contaminant. `wiener(contaminant = TRUE)` mixes the density with a uniform over the observed response range at a mixing dpar `lambda` on a logit link, inside the family, so it needs no core mixture and no uniform family. Same option on `lba()` and `rdm()` | `frmtmb.eam/R/wiener-family.R` | recovery of `lambda` at 5 percent contamination, 30 x 400; the log-likelihood at `lambda = 0` equals the plain family's | 1.5 |

Phase 3 total: about 9.5 days, one half-day having moved to Phase 1.
Item 3.4 keeps its three days: item 2.1 said 3.4 would GAIN a
per-subject bound rather than already carry one, and that gain became
item 1.0a's day and a half.

## Phase 4. The post-fit surface

| item | what | where | check | days |
|---|---|---|---|---|
| 4.1 | RP predictions. `rp_predict(fit, newdata, times, type = c("survival", "hazard", "cumhaz", "rmst", "median"), level, simultaneous)`, on the basis rebuilt from an exported `rp_knots(fit)`, with delta-method bands through the same covariance path `frm_curve()` uses; delayed entry and interval censoring get tests against flexsurvspline's counting-process and `interval2` forms | `frmtmb.spline/R/rp-predict.R` | point estimates identical to `flexsurv::summary()` at the same parameters; bands compared loosely, since flexsurv bootstraps; delayed-entry log-likelihood identical to flexsurv's | 3 |
| 4.2 | hmm on new data. `hmm_probs(fit, newdata)` and `hmm_viterbi(fit, newdata)` under the same group and time contract; `residuals(fit, type = "pseudo")`, the forward pseudo-residuals of Zucchini and others, chapter 6 | `frmtmb.latent/R/hmm.R` | decoding on the training frame passed as newdata equals the fit's own; pseudo-residuals against hmmTMB's on the same fit. Phase 0 measured the existing passes at 0.25 s and 0.09 s on 25,000 rows, so this item is about the newdata contract and not about speed | 2.5 |
| 4.3 | lca tooling. `lca_profiles(se = TRUE)` by the delta method; `lca_enumerate(formula, data, K = 2:6)` returning one row per K with log-likelihood, AIC, BIC, entropy and the smallest class share; `lca_blrt(small, large, nsim)` by simulate-and-refit, gated | `frmtmb.latent/R/lca-tools.R` | profile SEs against poLCA's `probs.se`; the enumeration table's BIC minimum on carcinoma agrees with the literature | 2.5 |
| 4.4 | learn's two missing pieces. frmtmb.eam exports the conditional mean of a Wiener first-passage time, the one export `dev/learn2-findings.md` asked for, and rlddm gains `fitted()`; `learn_family(init, update, choose, ...)` becomes a public constructor with the rule contract on its help page and `check_learn_family()` beside it, mirroring core's `check_custom_family()` | `frmtmb.eam/R/wiener-post.R`, `frmtmb.learn/R/family.R` | a user-written delta learner through `learn_family()` equals `bandit2arm_delta()` to 1e-12; `check_learn_family()` passes it | 3 |
| 4.5 | eam model checking. `ddm_qp(fit, newdata, quantiles)`: the quantile-probability table per condition and boundary from `simulate()`, with a plot method, which is the figure every DDM paper shows | `frmtmb.eam/R/wiener-post.R` | observed and simulated tables agree on data the fit simulated | 1 |
| 4.6 | sample throughput and vocabulary. Build the design once and multiply per draw in `posterior_epred()`, `posterior_linpred()`, `posterior_predict()`, `ranef()` and `coef()`; bulk and tail ESS from `posterior::ess_bulk()` in `summary()`; `constant()` as a prior spelling that maps the parameter; a `class = "Intercept"` route to ordinal thresholds; `laplace = TRUE` draws refused by `posterior_epred()`, `ranef()` and `hypothesis()` the way `loo()` already refuses them | `frmtmb.sample/R/methods-draws.R`, `sample.R` | Phase 0's `posterior_epred()` time falls by the factor the tier measures; identical values before and after. Phase 0 measured that time at 1.59 s for 4000 draws over 2000 rows against 213 s of sampling, so this item is worth doing for the vocabulary half and NOT for the throughput half, and it stays in Phase 4 | 3 |

Phase 4 total: about 15 days.

## Phase 5. The model menu

These are new families and measures. Each ships with the full bar in
"Rules".

| item | what | where | check | days |
|---|---|---|---|---|
| 5.1 | learn: `bandit_delta(n_option = K, lapse = FALSE, perseveration = FALSE, decay = FALSE, q0 = FALSE)`, one family covering hBayesDM's `bandit4arm_lapse`, `bandit4arm_2par_lapse`, `bandit4arm_4par`, `bandit4arm_lapse_decay` and `banditNarm_*` by option, with those names as aliases in `frm_learn_families()`; the engine already vectorizes over K for the Kalman filter | `frmtmb.learn/R/bandit-delta.R` | a Stan program per option set, identity to 1e-12; recovery at 100 x 200 per option | 4 |
| 5.2 | hmm multi-stream emissions. A matrix response, following `lca()`'s precedent rather than `mvbf()`, with `family = list(gamma(), von_mises())` one per column, per-state dpars per stream, and the product of stream densities in the forward pass; circular quantile starts for angles | `frmtmb.latent/R/hmm.R`, `hmm-streams.R` | step and angle on the elk tracks moveHMM ships, against moveHMM and hmmTMB to 1e-6 in log-likelihood | 5 |
| 5.3 | coupling inference. `frm_imag_coherence(fit)`, the volume-conduction-robust quantity `sqrt(C) sin(phase)`, with a delta-method interval; `frm_coherence_test(fit)`, a likelihood-ratio test of zero coherence against two independent Whittle fits, with the boundary correction, and a max-over-frequency version from the simultaneous band on the coherence curve | `frmtmb.coupling/R/coupling.R`, `test.R` | size of the test at the null over 500 simulated pairs; imaginary coherence against FieldTrip's on one exported dataset | 2.5 |

Phase 5 total: about 12 days. Item 5.4, the ode closed forms, is no
longer here: Phase 0 measured the design at 3948 s and moved it to
Phase 1 as item 1.0c, which is what this line used to say Phase 0 would
decide.

## Core seams

Filed for core, in the order the phases above need them. Each names the
extension item that is waiting.

| seam | needed by | what |
|---|---|---|
| `frm_sample()` consults compat refusals | 1.2 | the sampler is in an extension, so this is really a frmtmb.sample change; listed here because the registry contract in core's `?frmtmb_register_compat` has to promise that a `refused` row is consulted, not only printed |
| `frm_multistart(fit, n, jitter)` | 2.3 | the general form of `hmm_starts()`: refit from jittered starts, return the best, report the spread. Every mixture-shaped model in core wants it too |
| log-scale left and interval censoring, and the truncation normalizer | 4.1 | `dev/spline-core-findings.md`, "What it does not close". Delayed entry in a survival model is a left truncation, and it round-trips through a probability today |
| `required_aterms` as a disjunction | closes eam dev-findings #2 | "one of `dec()`, `vint()`" cannot be declared; the check is hand-rolled in `valid_y` |
| a `constraints` slot on `frmtmb_family()` | closes eam dev-findings #9 | functions of the whole dpar vector, checked at the start and reported at the optimum; the joint constraints on `sz`/`bias` and `st`/`ndt` have no other home |
| a prior-only objective | 4.6, deferred | `sample_prior = "only"` needs the likelihood term dropped at tape build, which is a core control, not a sampler option |
| `diagnose()`'s `unbounded_dpar` inverts on a BOUNDED link | 1.0a, 1.0b | filed by Phase 0. The check's evidence pair is `abs(est) > 10` and a standard error LARGER than the estimate, calibrated on `student()`'s `nu` where running off an unbounded link explodes the standard error. On a scaled logit the derivative vanishes at the edge, so the standard error collapses and the linear predictor saturates below 10. Measured on the Phase 0 eam rows: `ndt: (Intercept)` at 11.20 with a link-scale se of 2.37, which is 0.99999 of the bound and does not fire; and 7.17 in the non-convergent arm, below the magnitude threshold entirely at 0.9992 of the bound. `extreme_theta` does look at random-effect log standard deviations but its threshold of 8 is link-agnostic and the `ndt` component stands at 2.05. Four bounded `ndt` links ship today (`wiener()`, `lba()`, `rdm()`, `gddm()`), so the check needs to know whether a link is bounded and, if it is, to test proximity to the bound rather than magnitude plus a large standard error |
| `mixture()` over families whose primary dpar is not `mu` | eam races | `mixture(rdm(), ...)` is refused because the dpars are `v1..vn`; a structure-level opt-in would unblock it |

## Not doing, and why

- **More than two choices in `gddm()`.** `dev/gddm-nchoice-feasibility.md`:
  the only route that keeps state-dependent drift costs 12x the tape
  and 72x the gradient for an answer an order of magnitude worse.
  `lba()` and `rdm()` cover the choice count.
- **An hBayesDM cross-check.** Blocked on its toolchain, as
  `dev/rl-findings.md` records, and it would be an
  estimate-against-estimate comparison where a purpose-written Stan
  program gives an identity.
- **Estimated lag times, event times and dose intervals in ODE.** They
  decide where the solve is split, which is settled before the tape
  is built.
- **Time-frequency coherence and directed measures.** Different
  models, said so in the coupling vignette; the per-frequency Wishart
  cannot give them.
- **Horseshoe and R2D2 priors.** They carry auxiliary parameters the
  tape would have to hold; a real design, not a vocabulary entry.
- **A functional-response object like `refund::pffr`.** A scope
  difference, recorded in `dev/feature-gaps.md`.
- **Semi-Markov and continuous-time hidden Markov models.** Separate
  models; the discrete-time chain is the one the protocol fits.

## Sequence and size

| phase | days | runs in parallel with | gate |
|---|---|---|---|
| 0, measurement | 1, done | nothing; it was first | none |
| 1, small and high leverage | 10 | 2 | none; Phase 0's table has landed |
| 2, hierarchy validation | 7 | 1 | 2.1 waits on 1.0a, 2.2's rlddm half on 1.0b |
| 3, real data in | 9.5 | 4 | none left: 2.1 decided 1.0a and it moved |
| 4, post-fit surface | 15 | 3 | 1.3 lands before 4.1 |
| 5, model menu | 12 | nothing | none left: Phase 0 moved 5.4 to 1.0c |
| total | about 54.5 | | |

Days are single-lane estimates at the sizing convention of
`dev/hmm-feasibility.md`, which came in on the low side of its own
estimate. Phase 0 came in at one day of work against an estimate of
three, with about four hours of runs on top. Phases 1 and 2 are
per-package lanes and can run as seven lanes at once; the sequence
within a package is what the gate column constrains.

## What this document is not

It is not a promise of order beyond Phase 0. The survey's ranking rested
on doc figures and one afternoon of benchmarks, and Phase 0 existed to
replace it with measurements. It has: three items moved, two phase
totals changed, and two of the survey's cost worries turned out to be
nothing while a defect it had not looked for turned out to be the
worst thing in the table. The items moved, and the
table above is edited rather than appended to, so that the plan reads as
one document and not as a changelog.
