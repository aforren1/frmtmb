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
| 1.0a | DONE, and it beat its target. `ndt_group()` names the grouping and the bound becomes that group's fastest response; without the term the scalar bound stays in the LINK exactly where 0.6.0 put it, so an ungrouped model is bitwise unchanged, 0 ulp over 50 fixed-parameter probes across all five families. At the eam scale row, seed 20260908: log-likelihood -7003.01 against a target of -7027.4 and a baseline of -7148.8, convergence 0 at a maximum gradient of 9.87e-04, positive definite Hessian, 0 of 7 NaN standard errors, `sd(ndt)` 0.02543 against a truth of 0.02629, and 30 of 30 subjects below their own fastest response with 27.3 ms to spare. The bound-lifted arm reproduces Phase 0's -7027.37 and measures the row's own "22 ms" claim at 21.9 ms. `ndt_time()` reads the value in seconds. `gddm()` keeps the scalar bound and refuses `ndt_group()` by name, on SCOPE: it reads every dpar at the first row of its condition | `frmtmb.eam/R/ddm-shared.R`, `wiener-family.R` | met | 1.5 |
| 1.0a-note | The acceptance criterion "recovers `sd(ndt)`" is met IN APPEARANCE and Phase 2 should not lean on it. An estimator with the random effect on `ndt` switched off returns 0.02748 against a truth of 0.02629 where the shipped model returns 0.02543, and at 100 trials per subject the two agree in every digit: the observed floors carry the spread. What DOES separate the model from its floors is the per-subject error, and it is the argument for this item: the per-group arm halves with more data, 20.91 to 14.08 to 7.67 ms at 100, 200 and 400 trials, while the global bound's does not, 27.23 to 20.51 to 30.37. The global bound does not converge on the truth as data accumulates. Phase 2 should assert the per-subject error and the log-likelihood, not `sd(ndt)` | | | |
| 1.0b | DONE. `rlddm()` takes `ndt_group()` through a six-function seam exported from frmtmb.eam, whose contract stops at the bound: `ndt_bound()`, `ndt_bound_pending()`, `ndt_bound_attach()`, `ndt_bound_of()`, `ndt_bound_key()` and `ndt_apply()`. `ddm_ndt_install()` stays internal deliberately, because it rewrites a family's slots by name. The learn-rlddm scale row goes from convergence code 1 at a maximum gradient of 6.36e+09 with a non positive definite Hessian and four NaN standard errors, to code 0 at 3.12e-03 with a positive definite Hessian, no NaN, and 100 of 100 learners below their own fastest response with 20.9 ms to spare. The floors do not hand it that: per-learner error is 13.996 ms at a correlation of 0.940, against an analytic bound of 21.075 ms and 0.848 for the best constant fraction of each learner's own floor, and 22.107 ms for the same model with the component off. BREAKING: a constant `ndt` outside the link's range is refused at parse rather than reaching the optimizer as a NaN; the capability is respelled, not lost. Two measurements matter more than the feature. `sd(ndt)` is WORSE than doing nothing, 0.0334 against 0.0380 with the component off and a truth of 0.0393, so `?rlddm` says to read the per-learner error instead. And the tier's `cor_true` was 0 when the truth in the fitted parameterization is 0.776847 at that seed, computable per replicate with no fit; see the warning added to item 2.2 | `frmtmb.learn/R/rlddm.R`, `frmtmb.eam/R/ndt-seam.R` | met | 0.5 |
| 1.0c | DONE. `frm_lincmt()` is the analytic one- to three-compartment solution with the same `events` grammar, refusing by name what superposition cannot carry. On the Phase 0 design, whose `frm_ode()` arm reproduced at 3716.6 s against this row's recorded 3948 s, the closed form fits in 55 s at its default and 119 s with the run-in written out: at least 67x and at least 31x, taking each arm's slowest measurement over two passes against the other's fastest, which is what the design supports because its arms run one after another rather than interleaved. At 5 and 25 subjects, where `nlminb` takes identical iteration counts in every arm, the factors are 95x and 47x, and 89x and 39x. The load-independent count is 33 `RTMBode::ode()` calls per subject per objective evaluation against 0. Every estimate matches `dev/scale-findings.md` to seven significant digits, and `nlminb` returns code 0 at 4.7e-04 where the solver arm returns code 1 at 6.1e-03. The identity holds at 2.6e-12 of the trajectory's scale over 118 schedules against `frm_ode()` at `atol = rtol = 1e-12`, with the review's own 104 random schedules and 18 corners agreeing at worst 1.8e-11; against a 240-bit reference the error relative to the trajectory's maximum is 3.6e-10 on a seven-decade rate box, every coalescence swept to exact equality. It found no defect in `frm_ode()`'s arithmetic and one in its DEFAULT, which is item 1.0d. It also unblocks sampling: `frm_sample()` refuses `frm_ode()` and runs on `frm_lincmt()` | `frmtmb.ode/R/lincmt.R` | met | 3 |
| 1.0d | DONE, and not by raising the default. `frm_ode(ss_extrapolate = TRUE)` sums the geometric tail from the run-in's own last differences, which is arithmetic on tape variables and so runs DURING a fit at whatever parameters the fit has reached, at no extra solve. Worst error over a dosing interval against the exact limit, at the shipped `n_ss = 20`: 3.9e-03 to 9.9e-09 at a 23 h terminal half-life, 1.2e-02 to 2.6e-09 at 107 h, 2.1e-01 to 1.4e-08 at 670 h. The reason this was a wrong ANSWER and not a slow one is the gradient: at the 107 h design `d/dlog(k21)` from the truncated tape is +10.37 where the exact limit is -1.30, the wrong SIGN, and truncation manufactures precision on `lk21` by shrinking its standard error up to 41.6x while inflating `lke`'s up to 8.3x. The reviewer's `"auto"` proposal was measured and REJECTED: at the point `frm()` actually starts a nonlinear model it picks 7 cycles where the worst subject needs 169, and floored at 20 it is still dominated, reaching the integrator's own floor in 135 solves where this reaches it in 21 and truncation needs 201. The correction stands down smoothly because a fit walks there; both junctions read the curvature `2 / (1 - r)^3` to five figures, the same law the shape reads where there is no junction at all. It costs 258 tape nodes; `ss_extrapolate = FALSE` is 0.3.0 node for node. The warning now reports a distance rather than the cycle-to-cycle movement, and `?frm_ode` says plainly that it DETECTS rather than measures. Found on the way: past about `n_ss` = 1000, chained-solve error at `atol = 1e-8` contributes more than the extra cycles remove, so the docs say to tighten tolerances alongside it | `frmtmb.ode/R/ode.R` | met | 1 |
| 1.0e | DONE. `gddm()` refuses a distributional parameter that varies inside its condition. `gd_densities()` read every dpar at the FIRST ROW of its condition, so adding 100 to a covariate on 118 of 120 rows left the objective BITWISE identical for `mu`, `ndt`, `bs` and `bias`. It reached the random effects too: with a subject factor crossing the index, 2 of 4 subject deviations came back bitwise 0 and the fit reported `sd(s) = 7.6e-11` with no error and no warning, which a user reads as a finding about their data. Refusing rather than reading, because making the density read it is one solve per row instead of one per condition, and the count needs no clock since solves per evaluation IS the condition count: 60x. 0 false alarms over 21 designs the field writes, 26 of 26 on the paired drop arm, each naming the parameter, the variable and the user's own condition label. The tolerance is per pair plus a floor, because an exact comparison false-alarms on `poly()`, which returns rows 391 ulp apart for bitwise identical inputs; scaling the whole tolerance by the column maximum was tried and rejected in review because it made the answer depend on whether the user had centered. Residual band, stated as measured: a difference `s` on a column of maximum `M` is missed once `M / s` passes about 1e11, which neither search reached from a real covariate. Applying the floor only to computed columns closes that band and was built, measured and DECLINED: it false-alarms on 3 of 21 correct designs, because a precomputed `poly()` stored under a bare name carries exactly the noise the floor exists for. `gddm_simulate()` had the same defect and now refuses too, 5 of 5 caught where the previous build caught 0 | `frmtmb.eam/R/gddm.R` | met | 1 |
| 1.1 | DONE, and not as written. `frm_sample(control = list(adapt_delta =, max_treedepth =))` takes brms's spelling AND brms's meaning, and the fit-time options move to `fit_control =`. The `stan_control` of the original row was rejected: brms compatibility is a locked goal (SPEC.md section 5, "brms code ports mechanically"), `control` IS the Stan control in brms, and a second name would have left `vignette("brms-posterior")`'s "no route" row standing with a footnote instead of removing it. Putting the sampler options in `frmtmb_control()` was also considered and rejected: core has no Stan dependency by design (SPEC.md section 3) and `frmtmb_control()` validates eagerly, so core would validate options for a sampler it does not know, and `frm()` would accept `adapt_delta` and ignore it. The rename is breaking; the old shape is refused by name, keyed on the `frmtmb_control()` field names, which share none of rstan's fifteen | `frmtmb.sample/R/sample.R` | divergence count on a centered funnel falls between 0.8 and 0.99; `nuts_params()` reads the tightened step size back | 0.5 |
| 1.2 | DONE as written, except that the patches are prepared for filing rather than filed. `frm_sample()` reads the registry rows that name it with status `refused` and stops before taping; `frmtmb.ode` registers `frm_ode()` as a feature and `frm_ode() x frm_sample` as `refused`. The breakage was re-confirmed on 0.55.0 against RTMBode at commit 5242257, which is still the only build available. What the pre-flight can match against a model is a formula call and a family name; a refused row it cannot match warns rather than passing. `?frmtmb_register_compat` now promises that a `refused` row is consulted, which is a documentation change in core. The three patches are current and unapplied; the filing text is in `dev/stanctl-findings.md` | `frmtmb.sample/R/sample.R`, `frmtmb.ode/R/zzz.R` | an ODE fit refuses in under a second with no Stan call; the eam and latent rows stay `conditional` | 1 |
| 1.3 | DONE. `frm_curve(object, newdata, contrast = newdata2)` returns the difference with pointwise and simultaneous bands, and `frm_curve_feature()` finds its zero crossing. The gratia identity holds and the estimate gap is provably the fit rather than the arithmetic (5.55e-16); simultaneous coverage 0.970 over 200 seeds against a binomial mcse of 0.0154. An exact `gp()` is answered where the kriging residual is the same random variable, decided on the DESIGN rather than the numbers, and refused otherwise. It uncovered a pre-existing one-grid defect: `sp_sim_crit()` standardizes by a variance it did not draw from, so a simultaneous band over an exact `gp()` off the observed positions is at least 17 percent too narrow. That waits on the same core seam as the general cross-covariance | `frmtmb.spline/R/curve.R`, `curve-cov.R`, `curve-feature.R` | done | 1.5 |
| 1.4 | DONE. `valid_y` warns above 20 s across all five eam families and names milliseconds. 0 false alarms over 192 designs a two-choice task produces. Where it CAN fire on a correct model the rate is reported with the fastest response beside it, because the rate is not a function of the median: at about 50 s it spans 0.185 to 0.935 depending on drift | `frmtmb.eam/R/ddm-shared.R` | done | 0.25 |
| 1.5 | HALF DONE, and the half that is missing is a core seam. `reward(pay1)` is refused by CORE: `frmtmb_register_aterm(arity =)` holds one integer and the parser compares for equality, so registering at 1 would refuse the two-column spelling every example writes. Filed as a min-and-max arity, whose shape is settled by a measurement: under a plain range the derivation makes one group of three and the check skips for the TWO-column spelling as well, undoing this item's own fix. What WAS buildable closed a silent wrong answer: a duplicated schedule fits identically and simulates a different experiment, and is now refused across all eight families from a derivation rather than a declaration | `frmtmb.learn` | done | 0.5 |
| 1.6 | DONE. `frm_cross_spectrum()` splits at NA runs, so an artifact-rejected record is one call, and takes `window = "hann"`, measured to cost no degrees of freedom (8.021 against 8.008 at nominal 8) and identical to `stats::spec.pgram()` to 4.8e-16. Hann with `smooth` is refused, measured, because it delivers 4.67 to 5.66 where it claims 8. Recorded on the way: the cross-row correlation the docs understated is 0.397 at the default and 0.753 at `tapers = 4`, one independent frequency in four, though shipped `n` itself is right | `frmtmb.coupling/R/cross-spectrum.R` | done | 1 |

**Phase 1 is complete.** About 13 days of work against an estimate of
10, the overrun being the three items Phase 1 did not start with: 1.0d
and 1.0e were found by 1.0c and 1.0a rather than planned, and 1.0b's
export seam turned out to be a design decision rather than a rename.

Rule 3 ordered this phase from beginning to end, and it kept being
right. Every one of the three unplanned items was a SILENT WRONG
ANSWER found while building something else, and each outranked the
feature work that was already queued:

- 1.0a's `sv` arm converged, reported no warning, and gave a population
  non-decision time wrong by ten percent with a standard error of
  7.2e-06 on it. That arm now converges to 0.2479 against a truth of
  0.25 with 30 of 30 subjects inside their own floor.
- 1.0d's steady-state truncation was never reported to a fit at all,
  and biased every group in the same direction. Its gradient had the
  wrong SIGN, and it manufactured precision on the parameter a
  pharmacokineticist cares most about.
- 1.0e's first-row read reached every dpar of `gddm()`, including the
  random effects, where it could collapse a variance component to
  7.6e-11 with nothing printed.

None of the three was in the survey that produced this plan. All three
were found by measuring something adjacent.

The 0.55.1 disclosure in `frmtmb.eam/NEWS.md` telling users not to put
a random effect on `ndt` is discharged, and so is the `gddm()` note the
Phase 1 prose used to carry as an unnumbered item: it is 1.0e.

## Phase 2. Hierarchy: validate what is claimed

No new features. Each row is a recovery table and a third-party
identity at the realistic scale, and the result is written into the
family's help page whatever it says.

| item | what | reference | expected trouble | days |
|---|---|---|---|---|
| 2.1 | eam: `mu ~ cond + (1 \| s)`, `bs ~ (1 \| s)`, `ndt ~ (1 \| s)`, 30 x 400, 60 replicates; `check_laplace()` on one dataset, and `diagnose()` on every replicate. UNBLOCKED: item 1.0a landed, the design now converges with a positive definite Hessian, and the `diagnose()` seam Phase 0 filed is still open but no longer gates this row. Three things to carry in. First, do NOT assert `sd(ndt)`: the note under 1.0a explains why, the observed floors carry the spread and at 100 trials an estimator with no random effect on `ndt` agrees in every digit. Assert the per-subject error and the log-likelihood instead. Second, a question this round raised and could not settle at one seed: on a design where `ndt` is constant and only the boundary varies, the parameterization reports a 6.4 ms spread that is not there and a population `ndt` 23.6 ms low, but the global bound is 23.9 ms low on the same data, so the bias is the design rather than the change. Sixty replicates would say whether it matters at the designs the field produces. Third, the `eam-sv` row's condition effect no longer covers at the tier's seed, 0.822 (0.762, 0.883) against 0.9, under EITHER bound, and neither arm recovers `sv`; the tier records that rather than asserting it | hierarchical Wiener in brms on one dataset (gated tier), estimates within Monte Carlo error | Phase 0's expected trouble is discharged: the `sv` arm now converges to 0.2479 against a truth of 0.25 with a positive definite Hessian and no NaN standard errors, where it used to sit 0.42 of a standard error below the bound with a standard error of 7.2e-06. `diagnose()`'s `unbounded_dpar` still inverts on a bounded link and stays filed under Core seams, but the design it was misreading now converges | 1.5 |
| 2.2 | learn: `(1 \| p \| id)` on every parameter of bandit2arm_delta and rlddm, 100 x 200, 60 replicates; the importance correction's success rate at that scale, counted with 0.55.0's `imp_stalled()` rather than by eye. UNBLOCKED: item 1.0b landed and the rlddm row now converges. THREE warnings from 1.0b, all of which would otherwise cost this row its replicates. First, do NOT assert `sd(ndt)`: it is WORSE than switching the component off, 0.0334 against 0.0380 with a truth of 0.0393, so assert the per-learner error and the log-likelihood instead, exactly as item 2.1 was told. Second, the diagonal truth is the WRONG truth for a grouped arm: the design draws independent deviations but the block this fits is not the block it drew, so the correlation between the `bs` and `ndt` deviations in the fitted parameterization is about -0.78 before any fit, computable per replicate from the truths and the floors with no fit at all. Sixty replicates against a diagonal truth would fail on nearly every one for a reason unrelated to the estimator. Third, a per-learner floor is an observed minimum, so it moves with the trial count; say which count each figure came from | the family's own Stan program with the correlated block added | variance-component collapse is documented at 20 to 60 trials; this measures 200, which is where the field's designs sit. 60 replicates of rlddm at 12.3 minutes each is 12 hours, which is a second reason 1.0b comes first. A stalled correction is now distinguishable from a slow one, and the threshold that separates them has margins of only 1.54 and 2.10, so a rate near either edge is a finding about the threshold as much as about the model | 1.5 |
| 2.3 | DONE. hmm recovers on both arms. K = 3 gaussian over 40 replicates: coverage 95.6 percent at 459 of 480, state means within 0.021, transition matrix within 0.014, identity against depmixS4 over 3 replicates at 3.1e-14 to 4.9e-12 relative. `tr12 ~ (1 \| id)` over 15 replicates: coverage 94.4 percent at 184 of 195, `sd(tr12\|id)` 0.615 plus or minus 0.067 against 0.6, identity against hmmTMB over 6 replicates at 6.1e-12 to 6.9e-11. The expected trouble was real: the cold start sits 8.099 units under the optimum with convergence 0, a positive definite Hessian, and `diagnose()` printing "No convergence problems detected" verbatim, and `frm_allfit()` reaches that same wrong optimum on all four optimizers with a spread of 7.8e-07, so agreement between optimizers is not evidence here. `hmm_starts()` ships and recovers from it on 5 of 5 seeds at jitter 1 and above, 1 of 5 at 0.5. It costs 175.7 s per refit at this scale, so `?hmm` recommends no default `n` and says why. A comparison trap worth carrying: hmmTMB's `initial_state = "estimated"` fits one initial distribution PER SEQUENCE where frmtmb fits one shared, 38 extra parameters and 19.678 units on 20 sequences, with nothing saying so | depmixS4, hmmTMB | discharged, and the remedy ships | 2 |
| 2.4 | DONE, and the row's expectation was WRONG. "None expected" was false: on 8 of 200 replicates at K = 4, n = 2000, the 0.2.2 starting rule reached a local optimum 243 to 284 log-likelihood units below `poLCA(nrep = 10)`, and SEVEN of the eight passed the pair this plan's own scale tier asserts, a positive definite Hessian and a relative gradient under 1e-3. The cause was the START and not the surface on at least 2 of the 8: there all 10 of 10 poLCA single random starts found the global optimum and none ever visited frmtmb's mode. The old rule cut subjects on the mean of their item codes, and on this design every class endorses three of ten items at 0.85 and seven at 0.2, so that score explains 0.000814 of its own variance between the true classes: blind by construction. `lca()` now clusters on the response pattern and shrinks each class profile toward the pooled one, reaching poLCA's optimum on 200 of 200 of the block its one constant was tuned on AND on 200 of 200 of a second block of fresh seeds. The shrink weight sits at 0.9 rather than the first value that worked, 0.5, which loses a seed out of sample by 256 units; the usable interval runs about 0.5 to 0.99 and 0.9 is its interior. Identity against poLCA at 4.3e-14 relative on the log-likelihood, 5.6e-07 on profiles, 2.2e-06 on gating coefficients. Two coefficients UNDER-COVER and there is no remedy yet: pooled over 400 replicates `class3:x2` covers 91.75 percent and `class4:x2` 90.75 against a nominal 95, both binary gating slopes, while the third and the other six contain 95. A profile interval was measured and does not help, 103 of 118 both ways, because the standardized error has no heavy tail and no curvature: the shortfall is in the MAGNITUDE of the standard error rather than the shape of the likelihood, and `?lca` says so rather than recommending something unmeasured | poLCA | measured, and it found a defect | 0.5 |
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
| `frm_sample()` consults compat refusals | 1.2, DONE | the sampler is in an extension, so this was really a frmtmb.sample change; the core half was the contract, and `?frmtmb_register_compat` now promises that a `refused` row is enforced wherever the enforcing code can see the feature, with that code saying what it can see |
| `frm_features_used(spec)`, or the same answer from a `bf()` plus data | 1.2, for three kinds of the four | The pre-flight has to decide whether the model in front of it uses a refused feature, and a name cannot decide it: what identifies a feature in a model is WHERE its call sits in the formula. `se` inside the bar on the response's left is the `se()` addition term; `se` in a nonlinear body is a user's own helper, and frmtmb fits that model. ONE position is available without core: an addition term is written inside that bar and nowhere else, and the `bf()` object carries it, so frmtmb.sample matches the `aterm` kind by position in twelve lines and covers 9 of the 20 call-route rows. The seam is what the other three kinds need, `special`, `autocor` and `grammar`, where the parser's knowledge does not leave core: those are matched by name, which is justifiable but not sound, so a refusal registered on `mo()`, `ma()`, `cosy()` or `mm()` would fire on a user's function of that name. A core function returning the vocabulary rows a spec actually uses would let every consulting caller refuse exactly and warn never. Needed before a second package registers a refused row on one of those kinds. `dev/stanctl-findings.md` section 2.4 has the audit |
| A resolved-answer cache for `frm_compat()`, beside the vocabulary cache already in `compat_cache_store()` | 1.2, filed not built | `frm_compat(feature, status =)` re-resolves the whole pair table on every call, and `frmtmb.sample::frm_sample()` now calls it once per sampling call. Measured, interleaved, five rounds of twenty calls with a control arm: 0.117 s per call with two extensions loaded (111 vocabulary rows, 447 rules), 0.324 s with seven (138, 991), 0.581 s at a synthesized fifteen-extension scale (178, 1191), and twenty identical calls each pay it. The cost tracks the RULE count and the pair table it resolves, not the package count. Core already has the discipline this needs and uses it for the VOCABULARY: `compat_cache_store()` validates against `frmtmb_compat_contrib$features` rather than a dirty flag, so a caller that restores the list cannot leave a stale flag behind. Caching the resolved per-feature answer under that same key keeps ONE code path, keeps the call-time semantics that make a package's refusal appear the moment it loads, and takes the cost to near zero after the first call in a session. The alternative the sampler considered and declined, a cheap superset scan of `frm_compat_rules()` to skip the resolved read, is the wrong fix for a refusal: it is a second code path that can disagree with the first. `dev/stanctl-findings.md` section 2.6 has the table |
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
| 1, small and high leverage | 13, DONE | 2 | none |
| 2, hierarchy validation | 7, 2.5 done | 1 | none left: 1.0a and 1.0b both landed |
| 3, real data in | 9.5 | 4 | none left: 2.1 decided 1.0a and it moved |
| 4, post-fit surface | 15 | 3 | 1.3 lands before 4.1 |
| 5, model menu | 12 | nothing | none left: Phase 0 moved 5.4 to 1.0c |
| total | about 57.5, with Phase 1 done and 2.3 and 2.4 of Phase 2 |  |  |

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
