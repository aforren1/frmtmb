# brms vignette port: scorecard

## Re-measurement at frmtmb 0.67.0, 2026-10-05

Lane `vigport`, worktree `frmtmb-wt-vigport` at 9e902909. Build under
test: the installed 0.67.0 release library `rellib-r5`; 0.65.0
(`rellib-r3`) and 0.66.0 (`rellib-r4`) were measured too, to date
changes. brms 2.23.0 and the same nine vignettes as below. The full
record, with every defect and repro, is `dev/vigport-findings.md`.

The definitions are the 0.34.0 ones, unchanged: CLEAN is decided by
the raw pass alone, and `classify()` in `dev/brms-port/port-lib.R`
scores every row. The transform is unchanged too: 250 of 250
expressions have the same transformed source as at 0.34.0. Two things
in the instrument changed. The runner loads the installed build
(`PORT_LIB`) where it used `pkgload::load_all()`, and every expression
is seeded from its id so that a brms run builds the same simulated
data; the 0.67.0 counts were identical before and after seeding.

### Headline

<!-- BEGIN generated: mechanical headline -->
| calls | outcome | 0.34.0 | 0.67.0 |
|---|---|---|---|
| model (42) | CLEAN | 16 (38%) | 36 (86%) |
| model (42) | SPELLING | 21 (50%) | 3 (7%) |
| model (42) | FAIL | 4 (10%) | 2 (5%) |
| model (42) | CASCADE | 1 (2%) | 1 (2%) |
| post (102) | CLEAN | 35 (34%) | 63 (62%) |
| post (102) | SPELLING | 28 (27%) | 8 (8%) |
| post (102) | FAIL | 31 (30%) | 25 (25%) |
| post (102) | CASCADE | 8 (8%) | 6 (6%) |
<!-- END generated -->

36 of 42 model calls now run on the mechanical transform alone. Of the
remaining six: two nonlinear fits need `start =` because the transform
deletes the vignette's priors (one more is their cascade), group-level
`mo()` is refused, the JSS paper's `brm(threshold = )` fails in brms
2.23.0 as well, and the Stan-side `custom_family()` takes down its
model call.

### Per vignette

Model calls, then post-processing calls, as
`total | clean spelling fail cascade`:

<!-- BEGIN generated: mechanical per vignette -->
| vignette | model 0.34.0 | model 0.67.0 | post 0.34.0 | post 0.67.0 |
|---|---|---|---|---|
| brms_overview | 4 \| 1 2 1 0 | 4 \| 3 0 1 0 | 4 \| 1 1 1 1 | 4 \| 1 0 2 1 |
| brms_multilevel | 8 \| 2 4 2 0 | 8 \| 6 2 0 0 | 13 \| 3 2 6 2 | 13 \| 7 4 2 0 |
| brms_distreg | 4 \| 4 0 0 0 | 4 \| 4 0 0 0 | 12 \| 10 2 0 0 | 12 \| 11 0 1 0 |
| brms_nonlinear | 6 \| 3 3 0 0 | 6 \| 5 1 0 0 | 20 \| 5 7 7 1 | 20 \| 13 4 3 0 |
| brms_phylogenetics | 6 \| 5 1 0 0 | 6 \| 6 0 0 0 | 18 \| 15 2 1 0 | 18 \| 16 0 2 0 |
| brms_monotonic | 6 \| 0 5 1 0 | 6 \| 5 0 1 0 | 11 \| 0 8 2 1 | 11 \| 7 0 3 1 |
| brms_multivariate | 3 \| 0 3 0 0 | 3 \| 3 0 0 0 | 11 \| 0 4 7 0 | 11 \| 4 0 7 0 |
| brms_missings | 3 \| 0 3 0 0 | 3 \| 3 0 0 0 | 7 \| 0 2 5 0 | 7 \| 3 0 4 0 |
| brms_customfamilies | 2 \| 1 0 0 1 | 2 \| 1 0 0 1 | 6 \| 1 0 2 3 | 6 \| 1 0 1 4 |
<!-- END generated -->

<!-- BEGIN generated: mechanical moves -->

model calls: rows are the 0.34.0 class, columns the 0.67.0 class

| 0.34.0 \\ 0.67.0 | CLEAN | SPELLING | FAIL | CASCADE |
|---|---|---|---|---|
| CLEAN | 16 | 0 | 0 | 0 |
| SPELLING | 18 | 3 | 0 | 0 |
| FAIL | 2 | 0 | 2 | 0 |
| CASCADE | 0 | 0 | 0 | 1 |

post calls: rows are the 0.34.0 class, columns the 0.67.0 class

| 0.34.0 \\ 0.67.0 | CLEAN | SPELLING | FAIL | CASCADE |
|---|---|---|---|---|
| CLEAN | 32 | 0 | 3 | 0 |
| SPELLING | 23 | 2 | 3 | 0 |
| FAIL | 7 | 4 | 19 | 1 |
| CASCADE | 1 | 2 | 0 | 5 |
<!-- END generated -->

### The patch table at 0.67.0

The 4 retry rules (default family, bare constructor, `formula.`,
`newdata`) fire on no expression. Of the 9 per-expression patches,
only `brms_nonlinear.9.1` and `brms_multilevel.15.3` (`start =` for
the two loss-triangle models) are still needed, measured by a pass
that applies those two and nothing else. Two more patches,
`brms_overview.6.1` and `brms_distreg.5.1`, were wrong at 0.67.0: they
strip the relation from a `hypothesis()` string, which 0.61.0 refuses
as brms 2.23.0 does. Both are deleted from `patches.R` (punch round 1
of the review); the spell pass rerun without them changes no class and
runs 3 more post rows. The table, row by row, is in
`dev/vigport-findings.md`.

### Keep-prior transform

A second transform, beside the first and not in its place: keep
`prior =` and run the standalone prior code. Since 0.43.0 `frm()`
translates a brms prior and since 0.44.0 a located nonlinear prior
places the start.

<!-- BEGIN generated: keep-prior transform -->
model calls: raw transform 36 of 42 run; keep-prior transform 37 of 42 run
post calls: raw transform 63 of 102 run; keep-prior transform 70 of 102 run

rows whose outcome differs between the two transforms:

- brms_monotonic.10.1 (post): raw OK, keep-prior ERROR: object 'fit4' not found
- brms_monotonic.9.2 (model): raw OK, keep-prior ERROR: object 'prior4' not found
- brms_multilevel.15.3 (model): raw ERROR, keep-prior OK
- brms_multilevel.16.1 (post): raw ERROR, keep-prior OK
- brms_multilevel.17.1 (post): raw ERROR, keep-prior OK
- brms_multilevel.18.3 (post): raw ERROR, keep-prior OK
- brms_multilevel.18.4 (post): raw ERROR, keep-prior OK
- brms_multilevel.19.2 (model): raw ERROR, keep-prior OK
- brms_nonlinear.10.1 (post): raw ERROR, keep-prior OK
- brms_nonlinear.10.3 (post): raw ERROR, keep-prior OK
- brms_nonlinear.11.3 (post): raw ERROR, keep-prior OK
- brms_nonlinear.11.4 (post): raw ERROR, keep-prior OK
- brms_nonlinear.9.1 (model): raw ERROR, keep-prior OK
- brms_overview.7.1 (model): raw OK, keep-prior ERROR: No random-effect correlations match class=cor. These blocks have no correlation parameter: 1 | patient [us]
<!-- END generated -->

### Estimate plausibility, against brms fits

The table at 0.34.0 compared with the posterior means the vignettes
print, which for the simulated data sets belong to rows nobody can
rebuild. This one fits every model with brms 2.23.0 on the same rows
(`dev/brms-port/brms-fit.R`): 2 chains of 1000 iterations, 1000 draws,
the vignette's priors and `control`, and the vignette's seed where it
gives one. `brms_overview` was refitted at the vignette's own budget
(4000 draws); its `fit1` does not mix at either, so that row is no
reference. `z` is `(frmtmb - brms posterior mean) / brms posterior
SD`. Effective sample sizes were not recorded; at 1000 draws the Monte
Carlo error of a posterior mean is about 0.03 to 0.045 posterior SD,
so a median `|z|` near that size is agreement within sampling noise.
`frmtmb conv` is frmtmb's `opt$convergence`, and the summary is
given over the fits where both sides converged.

<!-- BEGIN generated: plausibility -->
**provenance**

frmtmb: frmtmb 0.67.0 from C:/Users/adf44/source/r/rellib-r5, mode spell, cap 120 s 
brms:   brms 2.23.0, rstan 2.32.7; 2 chains x 1000 iter (500 warmup) 
refitted at the vignette's own budget: brms_overview.3.1, brms_overview.7.1, brms_overview.10.1 

**per model (generated)**

| model | frmtmb | frmtmb conv | brms | coefs | median abs z | max abs z (coef) | sds | max abs z, sd | unmatched names | brms max Rhat | divergent | draws |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `overview.3.1` | OK | 0 | OK | 7 | 0.540 | 0.583 (`age`) | 3 | 1.74 | 0 | 1.59 | 1 | 4000 |
| `overview.7.1` | OK | 0 | OK | 7 | 0.0435 | 0.0568 (`diseasePKD`) | 2 | 1.68 | 0 | 1.00 | 0 | 4000 |
| `overview.10.1` | OK | 0 | OK | 6 | 0.169 | 0.520 (`Intercept[3]`) | 1 | 0.855 | 0 | 1.01 | 0 | 4000 |
| `overview.11.1` | ERROR |  | ERROR | 0 |  |  | 0 |  | 0 |  |  |  |
| `multilevel.2.1` | OK | 0 | OK | 4 | 0.0384 | 0.0558 (`child`) | 0 |  | 0 | 1.01 | 0 | 1000 |
| `multilevel.5.1` | OK | 0 | OK | 6 | 0.110 | 0.175 (`zi_Intercept`) | 0 |  | 0 | 1.01 | 0 | 1000 |
| `multilevel.9.1` | OK | 0 | OK | 4 | 0.0775 | 0.179 (`t2areayearc_3`) | 2 | 0.160 | 0 | 1.02 | 0 | 1000 |
| `multilevel.12.2` | OK | 0 | OK | 8 | 0.136 | 0.914 (`sigma_Intercept`) | 2 | 0.122 | 0 | 1.03 | 0 | 1000 |
| `multilevel.15.3` | OK | 0 | OK | 3 | 0.0732 | 0.195 (`theta_Intercept`) | 2 | 0.611 | 0 | 1.01 | 0 | 1000 |
| `multilevel.19.2` | OK | 1 | OK | 3 | 0.216 | 0.302 (`omega_Intercept`) | 4 | 1.07 | 0 | 1.01 | 0 | 1000 |
| `multilevel.23.1` | OK | 0 | OK | 1 | 0.0722 | 0.0722 (`Intercept`) | 2 | 0.607 | 0 | 1.03 | 0 | 1000 |
| `multilevel.26.1` | OK | 0 | OK | 1 | 0.191 | 0.191 (`Intercept`) | 1 | 0.0690 | 0 | 1.02 | 0 | 1000 |
| `distreg.2.1` | OK | 0 | OK | 4 | 0.0280 | 0.291 (`sigma_Intercept`) | 0 |  | 0 | 1.01 | 0 | 1000 |
| `distreg.7.1` | OK | 0 | OK | 4 | 0.0266 | 0.0531 (`camper`) | 0 |  | 0 | 1.00 | 0 | 1000 |
| `distreg.9.1` | OK | 0 | OK | 6 | 0.0689 | 0.151 (`zi_child`) | 0 |  | 0 | 1.01 | 0 | 1000 |
| `distreg.12.1` | OK | 1 | OK | 5 | 0.899 | 2.30 (`sx1_1`) | 2 | 0.896 | 0 | 1.02 | 0 | 1000 |
| `nonlinear.2.2` | OK | 0 | OK | 2 | 0.00794 | 0.0140 (`b2_Intercept`) | 1 | 0.349 | 0 | 1.01 | 0 | 1000 |
| `nonlinear.4.1` | OK | 0 | OK | 2 | 0.00967 | 0.0187 (`x`) | 1 | 0.252 | 0 | 1.00 | 0 | 1000 |
| `nonlinear.9.1` | OK | 0 | OK | 3 | 0.0447 | 0.151 (`theta_Intercept`) | 2 | 0.583 | 0 | 1.01 | 0 | 1000 |
| `nonlinear.13.1` | OK | 0 | OK | 2 | 0.0292 | 0.0364 (`ability`) | 0 |  | 0 | 1.00 | 0 | 1000 |
| `nonlinear.15.1` | OK | 0 | OK | 2 | 0.167 | 0.196 (`eta_ability`) | 0 |  | 0 | 1.00 | 0 | 1000 |
| `nonlinear.18.1` | OK | 0 | OK | 2 | 0.0560 | 0.0810 (`eta_ability`) | 0 |  | 0 | 1.00 | 0 | 1000 |
| `phylogenetics.3.1` | OK | 0 | OK | 2 | 0.122 | 0.199 (`Intercept`) | 2 | 0.164 | 0 | 1.01 | 0 | 1000 |
| `phylogenetics.7.1` | OK | 0 | OK | 2 | 0.146 | 0.292 (`Intercept`) | 3 | 0.197 | 0 | 1.01 | 0 | 1000 |
| `phylogenetics.11.1` | OK | 0 | OK | 3 | 0.0118 | 0.292 (`Intercept`) | 3 | 0.162 | 0 | 1.01 | 0 | 1000 |
| `phylogenetics.15.1` | OK | 0 | OK | 1 | 0.0384 | 0.0384 (`Intercept`) | 2 | 1.41 | 0 | 1.02 | 1 | 1000 |
| `phylogenetics.18.1` | OK | 0 | OK | 2 | 0.0581 | 0.0655 (`Intercept`) | 2 | 0.731 | 0 | 1.02 | 0 | 1000 |
| `phylogenetics.20.1` | OK | 0 | OK | 2 | 0.0574 | 0.112 (`Intercept`) | 2 | 1.22 | 0 | 1.02 | 1 | 1000 |
| `monotonic.2.1` | OK | 0 | OK | 2 | 0.204 | 0.207 (`moincome`) | 1 | 0.440 | 0 | 1.01 | 0 | 1000 |
| `monotonic.4.2` | OK | 0 | OK | 2 | 0.0525 | 0.0689 (`Intercept`) | 1 | 0.237 | 0 | 1.01 | 0 | 1000 |
| `monotonic.6.2` | OK | 0 | OK | 4 | 0.00366 | 0.0639 (`income4`) | 1 | 0.404 | 0 | 1.01 | 0 | 1000 |
| `monotonic.9.2` | OK | 0 | OK | 2 | 0.218 | 0.228 (`moincome`) | 1 | 0.408 | 0 | 1.01 | 0 | 1000 |
| `monotonic.13.1` | OK | 0 | OK | 4 | 0.0222 | 0.0326 (`age`) | 1 | 0.576 | 0 | 1.02 | 0 | 1000 |
| `monotonic.16.1` | ERROR |  | OK | 0 |  |  | 0 |  | 4 | 1.01 | 0 | 1000 |
| `multivariate.2.2` | OK | 0 | OK | 8 | 0.0183 | 0.0457 (`tarsus_Intercept`) | 6 | 0.375 | 0 | 1.02 | 0 | 1000 |
| `multivariate.6.3` | OK | 0 | OK | 5 | 0.00921 | 0.0236 (`tarsus_sexMale`) | 6 | 0.261 | 0 | 1.03 | 0 | 1000 |
| `multivariate.9.3` | OK | 0 | OK | 8 | 0.0392 | 0.181 (`sigma_tarsus_sexFem`) | 4 | 0.228 | 0 | 1.02 | 0 | 1000 |
| `missings.3.1` | OK |  | OK | 4 | 0.0487 | 0.0524 (`age:chl`) | 0 |  | 0 | 1.87 | 0 | 5000 |
| `missings.8.2` | OK | 0 | OK | 6 | 0.154 | 0.783 (`bmi_age`) | 2 | 1.05 | 0 | 1.01 | 0 | 1000 |
| `missings.11.2` | OK | 0 | OK | 6 | 0.167 | 0.824 (`bmi_age`) | 2 | 1.05 | 0 | 1.01 | 2 | 1000 |
| `customfamilies.2.1` | OK | 0 | OK | 4 | 0.0411 | 0.125 (`period4`) | 1 | 0.451 | 0 | 1.03 | 0 | 1000 |
| `customfamilies.7.1` | ERROR |  | OK | 0 |  |  | 0 |  | 4 | 1.02 | 0 | 1000 |

models fitted by both with matched coefficients: 39 of 42 model calls
matched coefficients: 149; abs z: median 0.0558, 90th pct 0.453, max 2.30
coefficients with abs z > 1: 2; > 2: 1
brms fits with max Rhat > 1.05: 2 (brms_overview.3.1, brms_missings.3.1); with divergences: 4
converged brms fits only (max Rhat <= 1.05): 37 models, 138 coefficients; abs z median 0.0521, 90th pct 0.291, max 2.30; > 1: 2
frmtmb fits that did not converge (opt$convergence != 0): 2
  brms_multilevel.19.2: code 1, Optimizer did not report convergence: false convergence (8). This is a nonlinear model; un
  brms_distreg.12.1: code 1, Optimizer did not report convergence: singular convergence (7)
both converged: 35 models, 130 coefficients; abs z median 0.0488, 90th pct 0.229, max 0.914; > 1: 0
variance components, brms converged: 62; abs z median 0.256, max 1.68; > 1: 6
variance components, both converged: 56; abs z median 0.217, max 1.68; > 1: 5
variance components past abs z 1, both converged:
  brms_overview.7.1      sd_patient__Intercept          frmtmb 0.000251  brms    0.408 (sd 0.243)  z -1.68  within one posterior SD of zero
  brms_phylogenetics.15.1 sd_phylo__Intercept            frmtmb   0.0107  brms   0.0677 (sd 0.0405)  z -1.41  within one posterior SD of zero
  brms_phylogenetics.20.1 sd_phylo__Intercept            frmtmb   0.0573  brms    0.714 (sd 0.539)  z -1.22  within one posterior SD of zero
  brms_missings.8.2      sd_residual____bmi             frmtmb     2.53  brms     3.34 (sd 0.767)  z -1.05  residual SD
  brms_missings.11.2     sd_residual____bmi             frmtmb     2.53  brms     3.30 (sd 0.737)  z -1.05  residual SD
<!-- END generated -->

On the 35 models where both sides converged, no coefficient is past
`|z| = 1`. Over all models the two past 1 are `sx1_1` and `sx2_1` in
`fit_smooth1`, a fit frmtmb flags with singular convergence at an
optimum four optimizers confirm. They are a defect, not a
disagreement: frmtmb's `s()` null-space column is brms's column times
a data-dependent factor (+5.372 and +5.431 on this model's rows;
magnitude 5.03 to 7.68, and negative in 3 of 12 columns, over
`gamSim` seeds 1 to 6), under brms's name; `sigma_sx0_1` has it too
(+6.084). Carried onto brms's column the two read z +0.03 and -0.04.
The cause is brms's `smoothCon(diagonal.penalty = TRUE)` against
frmtmb's default at `R/frame.R:3110`; a prior or a one-sided
hypothesis on `sx1_1` therefore targets a different parameter than in
brms. Details and logs: `dev/vigport-findings.md`, defect 1.

Variance components: on the both-converged fits, 5 of 56 are past
`|z| = 1`, all with frmtmb below brms. Three are within one posterior
SD of zero, the boundary case recorded at 0.34.0 (`sd(patient)` in the
kidney model is 2.5e-4 against brms's 0.41). Two are not: the `bmi`
residual SD of the two `mi()` models, 2.53 against 3.34 and 3.30,
estimated from 16 observed rows with 4 coefficients; the ML
degrees-of-freedom factor, an approximation here because `chl` is
missing in 3 of those rows and enters the `bmi` formula, explains
about half of that gap. The two
frmtmb fits that did not converge are `fit_smooth1` (above) and
`fit_loss2`, whose optimum eight refits leave weakly determined.

### Status of FN-1 to FN-15

`dev/brms-port/repros.R` on 0.67.0 (the FN items at their original
spelling):

| item | 0.67.0 |
|---|---|
| FN-1 default family | fixed 0.35.0 |
| FN-2 `mm()` | fixed 0.35.0 |
| FN-3 effects on nonlinear and `mo()` fits | fixed 0.35.0 |
| FN-4 `surface = TRUE` | fixed 0.66.0 |
| FN-5 one-sided `hypothesis()` | fixed 0.35.0 |
| FN-6 bare family constructor | fixed 0.35.0 |
| FN-7 `threshold =` | in the family since 0.67.0; as an `frm()` argument still `unused argument`, and brms 2.23.0 fails it too |
| FN-8 `a + b ~ 1` | fixed 0.35.0 |
| FN-9 `update()` spellings | fixed 0.35.0 |
| FN-10 `lf()` | fixed 0.35.0 |
| FN-11 `frm_multiple` post-processing | refusals that name the workaround (0.43.0) |
| FN-12 `pp_check()` on a multivariate fit | still refused on the ML fit; the draws method takes `resp =` |
| FN-13 `conditional_effects()` on an `mi()` fit | fixed |
| FN-14 `conditional_smooths()` | fixed 0.66.0 |
| FN-15 `add_criterion()` | still missing; `bayes_R2()` on a fit is a refusal naming `frm_sample()` |

### Families

<!-- BEGIN generated: families -->
```text
# frmtmb 0.67.0 from C:/Users/adf44/source/r/rellib-r5 
                    family frmtmb
                      acat   TRUE
              asym_laplace   TRUE
                 bernoulli   TRUE
                      Beta   TRUE
                  binomial   TRUE
               categorical   TRUE
                       cox   TRUE
                    cratio   TRUE
                cumulative   TRUE
                exgaussian   TRUE
               exponential   TRUE
                     Gamma   TRUE
                  gaussian   TRUE
                 geometric   TRUE
              hurdle_gamma   TRUE
          hurdle_lognormal   TRUE
        hurdle_negbinomial   TRUE
            hurdle_poisson   TRUE
          inverse.gaussian   TRUE
                 lognormal   TRUE
               multinomial   TRUE
               negbinomial   TRUE
                   poisson   TRUE
         shifted_lognormal   TRUE
               skew_normal   TRUE
                    sratio   TRUE
                   student   TRUE
                 von_mises   TRUE
                   weibull   TRUE
        zero_inflated_beta   TRUE
    zero_inflated_binomial   TRUE
 zero_inflated_negbinomial   TRUE
     zero_inflated_poisson   TRUE
    zero_one_inflated_beta   TRUE
                 dirichlet  FALSE
          discrete_weibull  FALSE
                   frechet  FALSE
         gen_extreme_value  FALSE
           logistic_normal  FALSE
                    wiener  FALSE

covered: 34 of 40 
missing: dirichlet, discrete_weibull, frechet, gen_extreme_value, logistic_normal, wiener 
```
<!-- END generated -->

34 of 40, against 29 of 40 below.

### What blocks the remaining calls

<!-- BEGIN generated: gap ranking -->
**mechanical port (spell pass, 0.67.0): calls blocked per cause**

| rank | cause | model | post | other | of which cascade |
|---|---|---|---|---|---|
| 1 | loo(), LOO(), waic() on a maximum-likelihood fit | 0 | 9 | 0 | 0 |
| 2 | brms's Stan-side custom_family()/stanvar() spelling | 1 | 4 | 2 | 5 |
| 3 | plot() of a fit lacks plot.brmsfit()'s N, variable, regex | 0 | 5 | 0 | 0 |
| 4 | frm_multiple() result has no pooled post-processing | 0 | 4 | 2 | 2 |
| 5 | add_criterion() does not exist | 0 | 3 | 0 | 0 |
| 6 | group-level mo() is refused | 1 | 1 | 0 | 1 |
| 7 | pp_check() refuses a multivariate fit | 0 | 2 | 0 | 0 |
| 8 | threshold = as an frm() argument (brms 2.23 drops it silently) | 1 | 1 | 0 | 1 |
| 9 | bayes_R2() on a maximum-likelihood fit | 0 | 1 | 0 | 0 |
| 10 | summary() refuses brms's swallowed waic = | 0 | 1 | 0 | 0 |

rows: 38 failing (3 model, 31 post, 4 other), 0 unclassified

**hand translation (0.67.0): open rows per cause**

| cause | ML | SAMPLE | of which cascade |
|---|---|---|---|
| loo(), LOO(), waic() on a maximum-likelihood fit | 11 | 6 | 0 |
| brms-only names with no frmtmb analog (threading, shinystan, ...) | 7 | 0 | 0 |
| frm_multiple() result has no pooled post-processing | 5 | 1 | 0 |
| brms's Stan-side custom_family()/stanvar() spelling | 3 | 2 | 0 |
| conditional_effects(method = 'predict') on draws | 0 | 5 | 0 |
| dirichlet / simo prior on a mo() simplex | 4 | 1 | 0 |
| families brms has and frmtmb lacks | 5 | 0 | 0 |
| add_criterion() does not exist | 3 | 1 | 0 |
| class = sd without nlpar on a nonlinear model (brms refuses it too) | 0 | 4 | 3 |
| summary() refuses brms's swallowed waic = | 2 | 2 | 0 |
| hypothesis() with no relation (brms refuses it too) | 3 | 0 | 0 |
| plot() of a fit lacks plot.brmsfit()'s N, variable, regex | 3 | 0 | 0 |
| pp_check() refuses a multivariate fit | 2 | 1 | 0 |
| pp_mixture(), stancode(), standata() on a fit: no refusal | 3 | 0 | 0 |
| group-level mo() is refused | 1 | 1 | 0 |
| plot() of draws | 0 | 2 | 0 |
| threshold = as an frm() argument (brms 2.23 drops it silently) | 2 | 0 | 0 |
| bayes_R2() on a maximum-likelihood fit | 1 | 0 | 0 |
| brms attached after frmtmb (bf() masked) | 1 | 0 | 0 |
| class Intercept, dpar sigma without a sigma formula (brms refuses it too) | 1 | 0 | 0 |
| custom family without a simulator | 1 | 0 | 0 |
| fixef() returns brms's matrix since 0.61.0 (translation stale) | 1 | 0 | 0 |
| horseshoe prior | 1 | 0 | 0 |
| log_lik() on an mi() model | 0 | 1 | 0 |
| mixture: brms's mu1 ~ formula beside the main formula | 1 | 0 | 0 |
| multinomial() needs K, which brms reads from the data | 1 | 0 | 0 |
| rescor_matrix() on draws (refused, names the fit) | 0 | 1 | 0 |
| update() on draws | 0 | 1 | 0 |

rows: 91 open, 0 unclassified
<!-- END generated -->

### Reproduce

```sh
export PORT_LIB=C:/Users/adf44/source/r/rellib-r5
export PORT_OUT=$PWD/dev/vigport-port-out/r5
Rscript dev/brms-port/run-all.R 120 raw
Rscript dev/brms-port/run-all.R 120 spell
Rscript dev/brms-port/summarize.R
```

The runs above write to `PORT_OUT`, which is gitignored. The files
this directory tracks from 0.34.0 (`results/`, `results-spell/`,
`results-v035/`, `results-merged.*`, `summary.txt`, `tables.txt`) are
the 0.34.0 record and were not overwritten; `summarize.R` reads
`results-merged.rds` from here as its baseline.

---

## The 0.34.0 measurement, 2026-09-02, kept as the baseline

Everything below describes 0.34.0 and is not relabeled. Its scripts now
read `PORT_LIB` and `PORT_OUT` (`dev/brms-port/env.R`) in place of the
removed `wt-audit` worktree.

Audit date 2026-09-02, frmtmb v0.34.0 (worktree `wt-audit`), brms 2.23.0,
R 4.6.1 on Windows.

The claim under test: *brms vignette code runs mostly unmodified after
`brm` becomes `frm` and the MCMC-specific arguments are dropped.*

The measurement is scripted end to end. `dev/brms-port/` holds the
extractor, the transform, the runner and the patch table; every number
below comes out of `dev/brms-port/summarize.R` and is reproducible with
two commands. Nothing in `R/` was changed: the gaps are recorded, not
fixed.

## Method

1. **Extract.** Every R chunk of the nine audited vignettes, in source
   order, from the `.Rmd` files brms ships. `brms_overview` and
   `brms_multilevel` ship as pre-rendered JSS LaTeX, so their code comes
   out of the `Sinput` environments; those environments also hold printed
   output, separated by the `R>` prompt in one file and by "does it parse
   as R" in the other. Chunks with `eval = FALSE` are extracted anyway:
   they are exactly the code a reader copies.
2. **Transform**, mechanically, over the parse tree rather than by regex:
   `brm(` becomes `frm(`, `brm_multiple(` becomes `frm_multiple(`, and
   the arguments `prior priors chains iter warmup cores backend threads
   refresh seed control init file silent save_pars sample_prior algorithm
   future normalize stanvars stan_funs` are deleted from `frm()`,
   `frm_multiple()` and `update()` calls. Everything else is verbatim:
   `bf()`, families, formulas, data code, post-processing. A standalone
   `set_prior()` / `prior()` / `get_prior()` expression goes to a
   brms-only bucket and is not run; a `prior =` argument *inside* `brm()`
   is simply dropped, so the model still runs.
3. **Run**, one OS process per vignette, chunks in order in one
   environment, each expression in `tryCatch` with a 120 s elapsed cap and
   a checkpoint written after every expression. `brms` is deliberately
   **not** attached: a post-processing function must resolve to frmtmb or
   fail, which is the thing being measured. Only data-side conveniences
   are shimmed (brms datasets, brms's unexported multi-membership
   simulator, a mirror for the dead UCLA URL in `brms_multilevel`).
4. **Second pass** with `dev/brms-port/patches.R`: the same code plus the
   documented spelling changes, so the "how much more works with a
   deliberate edit" number is measured, not estimated.
5. **Third pass** with only the two patches standing in for the gaps
   slated for v0.35, so the projected tally below is measured too.

CLEAN is decided by pass 1 alone. An expression that errored under the
mechanical transform is not clean even when the edit that rescued it was
made to a *different* expression: that still costs the porter an edit, so
those rows are `SPELLING: upstream`. One model call is in that position,
`brms_multivariate` `fit3`, which fails pass 1 because the `lf()` line
above it does.

Reproduce:

```sh
Rscript dev/brms-port/run-all.R 120 raw     # pass 1
Rscript dev/brms-port/run-all.R 120 spell   # pass 2
Rscript dev/brms-port/run-all.R 120 v035    # projection
Rscript dev/brms-port/summarize.R
```

## Headline

**Model calls: 42 across the eight runnable vignettes.**

| outcome | n | share |
|---|---|---|
| CLEAN: runs on `brm` -> `frm` plus argument removal alone | 16 | 38% |
| SPELLING: runs after a one-line documented change | 21 | 50% |
| FAILS | 4 | 10% |
| blocked by an upstream non-model failure | 1 | 2% |

The 21 spelling changes are not equal in weight:

| change | n | what it is |
|---|---|---|
| add `family = gaussian()` | 13 | brms defaults to gaussian; `frm()` has no default. See FN-1. |
| write `family = cumulative()` not `family = cumulative` | 1 | See FN-6. |
| a real, deliberate edit | 6 | nonlinear `start` values (3), `update()` formula spelling (2), one nlpar formula split |
| unblocked by the `lf()` rewrite on the line above it | 1 | See FN-10. |

So the honest reading of the claim is: **16 of 42 model calls run on the
mechanical transform alone; 14 more run after a purely clerical
family-argument fix; 7 more after an edit a reader of the migration guide
would make; 4 fail outright.**

`lf()` (FN-10) and the missing family default (FN-1) are both slated for
v0.35, and they land together. Measured on the projection pass, with only
those two fixed and nothing else: **model calls go from 16 to 30 of 42
(38% to 71%) and post-processing from 35 to 52 of 102**, leaving the
remaining spelling changes at 7 and the failures unchanged at 4.

**Post-processing calls: 102.** 63 run (35 clean; 28 after a spelling
change, of which 27 were unblocked by an edit to the model call above
them and 1 is its own rewrite), 8 are blocked by an upstream failure, and
31 fail. Of those 31, 9 are the known `loo` / `LOO` gap and 8 are one bug
(FN-3).

**Families (`brms_families`, a parameterization reference with no
runnable chunks, audited by name coverage): 29 of 40 named families are
accepted by `frm()`.** Missing: `categorical` (spelled `multinomial()`
with a matrix response in frmtmb), `cox`, `dirichlet`,
`discrete_weibull`, `frechet`, `gen_extreme_value`,
`hurdle_negbinomial`, `logistic_normal`, `von_mises`, `wiener`,
`zero_one_inflated_beta`. `hurdle_negbinomial` and
`zero_one_inflated_beta` stand out as cheap: both are compositions of
families frmtmb already has.

### Per vignette

Model calls, then post-processing calls, as
`total | clean | spelling | fail | cascade`:

| vignette | model | post |
|---|---|---|
| brms_overview | 4 \| 1 2 1 0 | 4 \| 1 1 1 1 |
| brms_multilevel | 8 \| 2 4 2 0 | 13 \| 3 2 6 2 |
| brms_distreg | 4 \| 4 0 0 0 | 12 \| 10 2 0 0 |
| brms_nonlinear | 6 \| 3 3 0 0 | 20 \| 5 7 7 1 |
| brms_phylogenetics | 6 \| 5 1 0 0 | 18 \| 15 2 1 0 |
| brms_monotonic | 6 \| 0 5 1 0 | 11 \| 0 8 2 1 |
| brms_multivariate | 3 \| 0 3 0 0 | 11 \| 0 4 7 0 |
| brms_missings | 3 \| 0 3 0 0 | 7 \| 0 2 5 0 |
| brms_customfamilies | 2 \| 1 0 0 1 | 6 \| 1 0 2 3 |

Notes per vignette:

- **brms_distreg** is the clean sweep: all four models and all twelve
  post-processing calls run, the last two only after rewriting one
  directional hypothesis. Distributional regression, zero-inflation with
  a `zi` predictor, and smooths in both `mu` and `sigma` all port
  verbatim.
- **brms_phylogenetics** is next best: five of six models clean,
  including `gr(phylo, cov = A)` with `data2`, `se()` meta-analysis, and
  the Poisson observation-level model. Only `update()` needed rewriting,
  and only `loo()` fails at the end.
- **brms_monotonic** is entirely blocked in pass 1 by the missing family
  default and then almost entirely clean: five of six models, including
  `mo(income) * age`. Group-level `mo()` is the one refusal, and it is a
  documented gap.
- **brms_multivariate** ports at the model level once edited (`mvbind`,
  `set_rescor`, per-response families, a smooth in one response) but no
  model call survives pass 1: two need the family default and the third
  is blocked by `lf()`. It also loses most of its post-processing:
  `add_criterion`, `bayes_R2`, `loo`, and `pp_check` on a multivariate
  fit.
- **brms_missings** ports all three models, including the two `mi()`
  models, but `frm_multiple` objects support almost no post-processing.
- **brms_customfamilies** stops at `custom_family()`: the brms call is
  Stan-shaped (`lb`, `ub`, `type`, `vars`, a Stan `stanvars` block) and
  frmtmb's takes an R `lpdf`. This is documented and expected; the
  cascade costs the vignette's remaining five calls.
- **brms_overview** and **brms_multilevel** are the JSS papers and the
  oldest code; they carry `LOO()`, `threshold =`, `formula. =` and
  `mm()`, which is why they score lowest.

Full per-expression tables are in `dev/brms-port/tables.txt` (generated,
not hand-maintained) and the merged record in
`dev/brms-port/results-merged.csv`.

## Estimate plausibility

Where a fit is CLEAN the point estimates sit on top of the posterior
means the vignettes print. Spot checks:

| fit | brms (posterior mean) | frmtmb (ML) |
|---|---|---|
| `fit_zinb1`, zero-inflated Poisson | -1.01, 0.87, -1.36, 0.80; zi 0.41 | -0.998, 0.872, -1.361, 0.796; zi 0.409 |
| `fit_zinb2`, zi ~ child | -1.07, 0.89, -1.17, 0.78; -0.95, 1.21 | -1.057, 0.889, -1.168, 0.771; -0.915, 1.186 |
| `fit_rent1`, t2 smooth + RE | 7.80, -1.00, 0.75, -0.07; sigma 1.95 | 7.800, -1.005, 0.736, -0.073; sigma 1.953 |
| `fit_loss1`, nonlinear | ult 5273.7, omega 1.34, theta 46.07 | 5292.5, 1.337, 45.90 |
| `model_simple`, phylogenetic | 38.38, 5.17; sigma 9.24 | 39.80, 5.176; sigma 9.171 |
| `fit1`, kidney lognormal + cens | 2.73, 0.01, 2.42, -0.40, -0.52, 0.60, -0.02; sigma 1.15 | 2.681, 0.016, 2.473, -0.429, -0.529, 0.615, -0.023; sigma 1.121 |

One expected divergence worth stating in any future vignette section: in
the kidney model brms reports `sd(patient) = 0.40` and frmtmb returns
3e-4. The variance component is at the boundary; brms's half-Cauchy prior
holds it away from zero and maximum likelihood does not. Nothing is
wrong, but a reader comparing the two summaries will notice it first.

## FAILS-NEW

Each item below is reproduced by `dev/brms-port/repros.R`, which prints
one line per item and can be rerun after a fix. Setup for the snippets:

```r
d <- data.frame(x = rnorm(60), g = gl(6, 10))
d$y <- 1 + 2 * d$x + rnorm(60)
d$o <- factor(sample(1:4, 60, TRUE), ordered = TRUE)
d$k <- as.integer(cut(d$y, 3))
```

### FN-1. `frm()` has no default family (13 model calls)

The single highest-leverage finding. brms defaults to `gaussian()`;
`frm()` refuses, so the most ordinary line in every vignette fails.

```r
frm(y ~ x, data = d)
#> Error: No family specified. Attach one with `bf(...) + gaussian()` or
#>   the `family` argument of frm()
```

The error message is also now stale: it names the `+` spelling first,
where the house style is `frm(bf(y ~ x), family = gaussian(), data = d)`.

### FN-2. `mm()` multi-membership grouping is unsupported (2 model calls)

`brms_multilevel`'s fourth worked example, both the plain and the
weighted form. Not mentioned anywhere in `?frmtmb-compat`,
`frm_compat_rules()` or the migration vignette, so a user meets it as a
missing-function error rather than a refusal.

```r
frm(y ~ 1 + (1 | mm(g, g)), data = d, family = gaussian())
#> Error: could not find function "mm"
```

### FN-3. `conditional_effects()` sees no predictors on nonlinear or `mo()`-only fits (8 post calls)

The largest single post-processing failure. `conditional_effects()`
enumerates plottable terms from the `mu` linear predictor only, so a
nonlinear model (whose covariates live in the nlpar formulas) and a model
whose only fixed term is `mo()` both come up empty. Passing the effect
explicitly works, which is what makes this a lookup bug rather than a
missing feature: `conditional_effects(fit5, "income:age")` succeeds on
the same class of fit.

```r
fnl <- frm(bf(y ~ a * x + b, a ~ 1, b ~ 1, nl = TRUE), data = d,
           family = gaussian(), start = list(beta = c(1, 1)))
conditional_effects(fnl)
#> Error: No plottable predictors found for dpar 'mu'

fmo <- frm(y ~ mo(o), data = d, family = gaussian())
conditional_effects(fmo)
#> Error: No plottable predictors found for dpar 'mu'
```

### FN-4. `conditional_effects(surface = TRUE)` errors on a 2-D smooth

```r
d2 <- d; d2$z <- rnorm(60)
fsurf <- frm(y ~ t2(x, z), data = d2, family = gaussian())
conditional_effects(fsurf, surface = TRUE)
#> Error: requires numeric/complex matrix/vector arguments
```

### FN-5. `hypothesis()` rejects brms's one-sided form

brms writes directional hypotheses with `>` or `<` and reports an
evidence ratio. frmtmb parses only `=`, so the brms idiom errors with a
message that does not say what to write instead.

```r
flin <- frm(y ~ x, data = d, family = gaussian())
hypothesis(flin, "x > 0")
#> Error: Hypothesis 'x > 0' must evaluate to a single number at the
#>   fitted estimates
```

The frequentist reading (`hypothesis(flin, "x")` and read the sign and
interval) is available; nothing points the user at it.

### FN-6. A family passed as a bare constructor is rejected

`family = cumulative` without parentheses is legal brms and legal
`stats::glm`. The error names an internal class and misleads.

```r
frm(k ~ x, data = d, family = cumulative)
#> Error: Cannot interpret `family` of class frmtmb_family
```

### FN-7. `threshold = "equidistant"` is unsupported

Ordinal threshold parameterizations (`"equidistant"`, `"sum_to_zero"`)
have no frmtmb equivalent and no named alternative.

```r
frm(k ~ x, data = d, family = cumulative(), threshold = "equidistant")
#> Error: unused argument (threshold = "equidistant")
```

### FN-8. A multi-parameter nlpar formula is rejected

brms lets one formula name several nonlinear parameters; the error
frmtmb gives is about identifier syntax and does not suggest the split.

```r
frm(bf(y ~ a * x + b, a + b ~ 1, nl = TRUE), data = d,
    family = gaussian(), start = list(beta = c(1, 1)))
#> Error: Invalid parameter name 'a + b': names must be alphanumeric
#>   without dots or underscores (they collide with coefficient naming)
```

### FN-9. `update()` does not speak brms's argument names

Three separate gaps in one method. frmtmb's `update()` splices any named
argument straight into the stored `frm()` call, so brms's `formula.` and
`newdata` both land as unused arguments, and a one-sided update formula
(`~ . + z`) is never expanded against the original.

```r
update(flin, formula. = ~ . + I(x^2))
#> Error: unused argument (formula. = ~. + I(x^2))
update(flin, formula = y ~ x, newdata = d)
#> Error: unused argument (newdata = d)
update(flin, formula = ~ . + I(x^2))
#> Error: The model formula needs a response (left-hand side)
```

`formula.` and `newdata` are one-line aliases; the one-sided update is
`stats::update.formula` against the stored formula.

### FN-10. `lf()` is missing

`bf(...) + lf(sigma ~ 0 + sex)` is the brms spelling for a dpar formula
added to an existing `bf()`, used in `brms_multivariate` where each
response carries its own dpar.

```r
lf(sigma ~ x)
#> Error: could not find function "lf"
```

### FN-11. `frm_multiple` fits support almost no post-processing (4 post calls)

`summary()` works. `plot()`, `conditional_effects()`, and the posterior
package entry points do not.

```r
fm <- frm_multiple(y ~ x, data = imps, family = gaussian())
plot(fm, variable = "^b", regex = TRUE)
#> Error: 'x' is a list, but does not have components 'x' and 'y'
conditional_effects(fm, "x")
#> Error: no applicable method for 'conditional_effects' applied to an
#>   object of class "frmtmb_multiple"
as_draws_array(fm)       #> Error: arg 'fits' is non-atomic
nchains(fm)              #> Error: no applicable method for 'nchains'
```

`plot()` failing with a base-graphics message rather than a refusal is
the worst of these: it looks like a bug in the user's call.

### FN-12. `pp_check()` refuses multivariate fits (2 post calls)

A clean refusal, not a crash, but `pp_check(fit, resp = "tarsus")` is the
standard multivariate check and there is no alternative named.

```r
pp_check(fmv, resp = "y")
#> Error: pp_check() is not supported yet for multivariate fits
```

### FN-13. `conditional_effects()` on an `mi()` fit cannot build its grid

The imputed predictor has no complete values to form a grid from, and the
prediction path demands them.

```r
conditional_effects(fit_imp2, "age:chl", resp = "bmi")
#> Error: mi(chl): newdata must supply complete values
```

(The reduced repro in `repros.R` fails one step earlier, with `'from'
must be a finite number`; same cause.)

### FN-14. `conditional_smooths()` is missing

The companion of `conditional_effects()` for smooth terms, used in
`brms_multilevel`.

```r
conditional_smooths(fsurf)
#> Error: could not find function "conditional_smooths"
```

### FN-15. `add_criterion()` and `bayes_R2()` are missing (4 post calls)

`add_criterion()` is worth a stub even in a frequentist package: it is
how brms users attach *any* criterion to a fit, and AIC/BIC are the
natural payload.

```r
add_criterion(flin, "loo")   #> Error: could not find function "add_criterion"
bayes_R2(flin)               #> Error: could not find function "bayes_R2"
```

## FAILS-KNOWN

Confirmed present, no action implied:

- `loo()` (6), `LOO()` (3) and `expose_functions()` (1): 10
  post-processing calls, plus one `stanvar()` assignment that is not a
  post-processing call. Sampling-based or Stan-code-specific by
  construction.
- `custom_family()` takes an R `lpdf` and not brms's `lb` / `ub` /
  `type` / `vars` plus a Stan `stanvars` block. Documented in the
  migration vignette; costs `brms_customfamilies` five cascaded calls.
- Group-level `mo()`: `(mo(income) | city)` is refused with a clear
  message naming the supported forms. Recorded in `dev/feature-gaps.md`.
- `set_prior()` / `prior()` / `get_prior()` pipelines: 4 expressions,
  diverted to the brms-only bucket by design. frmtmb's `set_prior()` is a
  penalty, not a prior, so auto-conversion would fabricate agreement.

## Porting guide

Distilled from the transform that produced the numbers above. This can
seed a vignette section.

**The mechanical part**, which a user can do with find-and-replace:

1. `brm(` -> `frm(`, `brm_multiple(` -> `frm_multiple(`.
2. Delete these arguments: `prior`, `priors`, `chains`, `iter`, `warmup`,
   `cores`, `backend`, `threads`, `refresh`, `seed`, `control`, `init`,
   `file`, `silent`, `save_pars`, `sample_prior`, `algorithm`, `future`,
   `normalize`, `stanvars`, `stan_funs`. Note that `control` and `priors`
   are also real `frm()` arguments with different meanings, so they must
   be deleted rather than kept.
3. Delete standalone `set_prior()` / `prior()` / `get_prior()`
   expressions. frmtmb has `set_prior()` but it is a penalized-likelihood
   term, not a prior; converting mechanically would be wrong.
4. Everything else stays: `bf()`, `mvbind()`, `set_rescor()`, addition
   terms, `gr(..., cov = )`, `data2 = `, families and links, smooths,
   `mo()`, `mi()`, `gp()`, `vint()`.

**The one clerical fix**, needed on 31% of the vignettes' model calls:

5. Give every model a family. brms defaults to gaussian; `frm()` does
   not. Add `family = gaussian()` as an argument:

   ```r
   # brms
   brm(y ~ x, data = d)
   # frmtmb
   frm(bf(y ~ x), family = gaussian(), data = d)
   ```

   Use the `family =` argument rather than `bf(y ~ x) + gaussian()`. The
   argument spelling maps to brms position for position and matches
   lme4/glmmTMB. The `+` form is still valid, and is the only way to give
   a *per-response* family in a multivariate model. Call the constructor:
   `family = cumulative()`, not `family = cumulative`.

**The deliberate edits**, which need a person:

6. Nonlinear models need `start`. In brms the priors put the sampler in
   the right region; frmtmb starts from zero and the gradient is `NaN`.
   Read the vignette's own prior means straight across:

   ```r
   frm(bf(cum ~ ult * (1 - exp(-(dev/theta)^omega)),
          ult ~ 1 + (1 | AY), omega ~ 1, theta ~ 1, nl = TRUE),
       data = loss, family = gaussian(),
       start = list(beta = c(5000, 1, 45)))
   ```

7. One formula per nonlinear parameter: `b1 + b2 ~ 1` becomes
   `b1 ~ 1, b2 ~ 1`.
8. `update()` takes `formula` (not `formula.`), `data` (not `newdata`),
   and a complete two-sided formula rather than a `~ . + z` delta.
9. Directional hypotheses lose the operator: `"a > b"` becomes `"a - b"`,
   and the sign of the estimate with its confidence interval replaces the
   evidence ratio.
10. `bf(...) + lf(sigma ~ z)` becomes `bf(..., sigma ~ z)`.
11. `ar()` / `ma()` / `arma()` need `cov = TRUE`, and `sigma` is then the
    marginal residual SD, not brms's innovation SD. (Documented in
    `?frmtmb-autocor`; no audited vignette exercises it.)

**What has no port**: `loo()` / `waic()` / `add_criterion()` become
`AIC()` / `BIC()` / `anova()`; `stanvar()` and `expose_functions()` have
no analog; a `custom_family()` must be rewritten with an R `lpdf`.

## Follow-up worth filing

Ordered by vignette calls unblocked per unit of work:

1. A gaussian default for `frm()` (FN-1) and `lf()` (FN-10), the two
   already slated for v0.35. Together they take model calls from 16 to 30
   of 42 and post-processing from 35 to 52 of 102, measured on the
   projection pass.
2. `conditional_effects()` term enumeration for nlpar and `mo()` fits
   (FN-3, 8 calls).
3. `update()` argument aliases and one-sided formulas (FN-9).
4. A `hypothesis()` path for `>` / `<`, or an error that names the
   replacement (FN-5).
5. `mm()`, or a refusal that says so (FN-2).
6. `add_criterion()` returning AIC/BIC, `conditional_smooths()`,
   `bayes_R2()` -> an R2 that is not Bayesian (FN-14, FN-15).
7. `frm_multiple` post-processing methods (FN-11).
8. `hurdle_negbinomial()` and `zero_one_inflated_beta()`: both compose
   from families already present.

Once item 1 lands the audit should be rerun: the numbers above are the
pass-1 baseline for v0.34.0, and the projection pass is a stand-in, not a
substitute for measuring the shipped fix.

## Files

- `dev/brms-vignette-port.md` (this file)
- `dev/brms-port/port-lib.R`: chunk extraction and the AST transform
- `dev/brms-port/patches.R`: the documented spelling changes
- `dev/brms-port/run-vignette.R`, `run-all.R`: the runner (one process per
  vignette, checkpointed per expression; modes `raw`, `spell`, `v035`)
- `dev/brms-port/summarize.R`, `tables.R`, `sanity.R`: reduction. The
  CLEAN / SPELLING / FAIL / CASCADE rule lives in `classify()` in
  `port-lib.R` so both reducers score rows the same way.
- `dev/brms-port/families-coverage.R`: `brms_families` name coverage
- `dev/brms-port/repros.R`: minimal repro for every FAILS-NEW item
- `dev/brms-port/results/`, `results-spell/`, `results-v035/`:
  per-vignette logs and RDS
- `dev/brms-port/results-merged.csv`, `summary.txt`, `tables.txt`:
  generated tables
- `tests/testthat/test-brms-port.R`: six tests pinning seven ported fits,
  under 3 s

Not audited: `brms_threading` (pure Stan performance).
