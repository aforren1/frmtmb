# Lane remeasure: the brms vignette port and bin 2 of brms's suite at 0.70.0

Date 2026-10-13. Worktree `frmtmb-wt-remeasure`, base d2ef6312
(frmtmb 0.70.0). Build judged: `C:/Users/adf44/source/r/rellib-r8`
(frmtmb 0.70.0, frmtmb.sample 0.18.0, frmtmb.eam 0.11.3). brms 2.23.0,
rstan 2.32.7, loo 2.10.1, R 4.6.1 on Windows.

This is a measurement lane. Nothing under `R/`, `extensions/*/R`,
`tests/`, `man/` or `vignettes/` changed. Every count in a generated
block is pasted from a reducer's output by `dev/remeasure-fill.R`;
the prose quotes those blocks.

STATUS: final after fix 2b (`dev/reviews/2026-10-13-remeasure.md`).
Every figure in this paragraph is the generated bin-2 summary's
(section "Task B", block "The totals"). Bins 1 and 2 of brms's suite
pass 534 of 823 (64.9%), against "389 of 823" with bin 2 unmeasured:
bin 2T 23 of 59, bin 2P 122 of 270 (178 run through the adapter, 122
pass). Over every assertion except the 114 that frmtmb has no
counterpart for by decision, 534 of 709 pass (75.3%); the 65 gaps
(functions, families or CDFs frmtmb lacks), 31 rows not ported (cost)
and 5 port artifacts stay in that denominator. The denominator excludes
exactly the "cannot transfer" rows of the classes stan (26),
brms-internal (63), mcmc (4), no-draws (3), no-density (12),
no-quantile (3) and brms-backends (3), each row with its reason. Of the
145 bin-2 passes, 29 compare a value and 116 check only a shape,
length, names or type.

### Fix 2b: what changed

1. **G2's `inv_link_<ordinal>()` array maps are ported** through
   `brms_adapt_ilord()` over the existing `brms_adapt_dord()` (brms's
   `x` is `thres - eta` for cumulative and sratio, `eta - thres` for
   cratio and acat). `distributions:300` (one loop of 8 assertion
   lines) has a verdict per line (`dev/brmsport-bin2-split.tsv`): 7
   pass in 32 of 32 iterations, and the cumulative identity, 16 of 32,
   is a divergence (domain), the ordered-threshold limit of `:238`.
   `distributions:479` (2 assertions) holds as a whole. Evidence: the
   review's `dev/remeasure-rev3-g2.R` (seed 300), rerun by the lane,
   `dev/remeasure-out/probe/g2-inv-link.txt`.
2. **`dcategorical()` (`:527`) and `inv_link_categorical()` (`:549`)**
   are not ported (cost), the class of the categorical `log_lik` and
   `posterior_predict` rows. The other 48 brms-internal 2P rows stay
   decisions: the `link_*()` maps and the "inverts" blocks (22), G5
   (18), G6 (5) and G7 (3).
3. **Wording.** "168 ported" now reads "run through the adapter" with
   its pass count; categorical reads "not ported (cost); the existing
   adapter refuses brms's mu1/mu2".
4. **New defect R8**, open, low: the `$lcdf` slot returns the CDF.

### Fix round 2: what changed

Superseded where fix 2b says otherwise (the headline, the rate, the
G2 rows).

1. **G1 (brms's exported d/p/q/r functions) is not a decision.** The 24
   d and r rows of families frmtmb has are ported through shims with
   brms's signatures over `row_lpdf()` and `sim_draw()`
   (`brms_g1_shims` in `dev/brmsport-bin2-adapters.R`, the review's
   construction `dev/remeasure-rev2-g1.R`): all 24 hold. So does
   `pinv_gaussian()` (`distributions:111`), through `row_lpdf()` of a
   left-censored row, since inverse.gaussian is the one family here
   with a CDF. The other 15 p rows are gaps (no CDF, class `no-cdf`),
   the 9 frechet and gen_extreme_value rows absent (gap; famuni is
   adding them, so the consolidation rerun may flip them), the 3 q
   rows a decision (`no-quantile`) and the 12 multivariate normal and
   t rows a decision (`no-density`), both citing `SPEC.md:153-154`
   (numeric d/p/q/r functions come from gamlss.dist or extraDistr; no
   family object carries a quantile, and no multivariate density
   exists outside a fit). The 3 wiener rows compare brms's two wiener
   backends with each other; frmtmb.eam has one implementation, which
   agrees with both to 5.7e-14 on the log scale at the block's inputs
   (`dev/remeasure-bin2-wiener-g1.R`): a decision, `brms-backends`.
2. **`distributions:238` has a verdict per assertion line**
   (`dev/brmsport-bin2-split.tsv`, from the review's
   `dev/remeasure-rev2-d238.R`, rerun: `dev/remeasure-out/probe/d238.txt`):
   7 pass (the sratio, cratio and acat identities and the four
   `dim()` lines, 32 of 32 iterations each) and 1 divergence, domain
   (the cumulative identity, 16 of 32, which needs ordered
   thresholds). The reason text now says that frmtmb's NaN falls on
   either side of a decreasing pair: at draw 5 of `bin2-probe3.txt`
   categories 1 and 2 are NaN where brms's category 1 is 0.805.
3. **`log_lik:302` (com_poisson)** is a divergence, class
   parameterization (R6), in the denominator.
4. **The unported G3/G4 rows** are "not ported (cost)", in the
   denominator: frmtmb has each of these models. In the categorical
   blocks, categorical: not ported (cost); the existing adapter refuses
   brms's mu1/mu2. Multinomial is absent (brms's bare `multinomial()`,
   gap rank 26), and dirichlet,
   dirichlet_multinomial, dirichlet2 and logistic_normal absent
   (famsimplex is adding dirichlet and logistic_normal). G5 to G7 and
   G2's link maps stay "brms-internal", a decision.
5. **The headline and the rate** above are restated from the
   regenerated ledger, with the line naming the excluded classes.
6. **R6:** 7.30 is `exp()` of the posterior mean of the Intercept.

### Fix round 1: what changed

Superseded where fix round 2 says otherwise (the G1 classes, the rate,
`distributions:238`, `log_lik:302`, the unported G3/G4 rows).

1. **The rate.** "Cannot transfer" is split into no counterpart
   (decision), absent (gap) and port artifact, in bins 1 and 2 alike;
   the rate excludes the decisions only. The first record's "87.1% of
   551" left the 28 gaps out. Bin 2P's 67 brms d/p/q/r rows are now
   class `no-export`, a decision.
2. **Ordinal rows ported.** The 12 ordinal rows of `tests.log_lik.R`
   (`:461` to `:478`) and `tests.posterior_predict.R` (`:301` to
   `:327`) go through each ordinal family object's exported
   `ord_lpdf_make(identity)` and `ord_sim_make(identity)`: all 12 pass.
   `distributions:446` (2 assertions) passes; `distributions:238` (8)
   is a divergence, class domain (section "Task B").
3. **What the passes assert** is a column of the bin-2 summary.
   `log_lik:32` tests the adapter's own `ll * w[i]`, not frmtmb, and is
   now cannot transfer, class adapter (a port artifact). The "can say
   no" demonstration is redone on frmtmb's density (section "Task B").
4. **R6 restated** as a brms upstream bug: brms 2.23.0 fits `mu` as the
   COM-Poisson rate and post-processes with `mu = rate^(1/shape)`.
5. **R1:** the octopus-merge wording is corrected, and the review's
   mechanism (`predictor_par_units()`, the `optimHess()` step) is
   recorded. Owner lane optstart.
6. **Median `|z|`** is compared on the common coefficient set: 0.047
   over 124 coefficients on every build. The own-set figures are not
   comparable, and the claim that `brms_missings.8.2` is the whole
   change is withdrawn.
7. **Families** are fitted, not only named: 34 of 40 fit in core, 35
   with frmtmb.eam. `multinomial` fits through `multinomial(K = 3)`.
8. **Gap list:** rank 7 (`threshold =`) is "deliberate: as brms", like
   ranks 13, 15, 19 and 20; rank 18 and the six 2T rows on several
   models are open work for postfix, as `dev/brms-api-diff.md` (c)
   lists them; R1 to R4 are ranked with owners in the gap block.
9. The reviewer's scratch trees `dev/remeasure-rev-wt/`, `-b2/`, `-b3/`
   and `-tables/` are deleted; the reviewer's scripts are kept.

## Ranked gaps at 0.70.0

What still blocks a vignette call on 0.70.0, one root cause per failing
row, a cascade charged to the row that should have made its object.
"deliberate" is a refusal frmtmb makes on purpose for a stated reason
(a maximum-likelihood fit has no posterior; there is no Stan program;
brms refuses the same call). "open" is a brms call frmtmb could serve
and does not yet. "instrument" is a row the measurement scores wrong.
The lane column names the lane of this round that works on the gap.

Read with these notes:

- **Lane postfix** owns the open ranks 6, 8, 11, 18 and 23: 3
  mechanical calls and 15 hand rows, and ranks 2 to 5 of the defect
  table at the end of the block. Rank 3 (`plot()`, `as_draws_array()`,
  `nchains()` on a `frm_multiple()` result) is scored deliberate because
  each needs posterior draws that m pooled ML fits do not have; the
  pooling of DRAWS is rank 23, `frm_sample()` of a `frm_multiple()`
  result.
- **Rank 7** (`threshold =` as an `frm()` argument) is "deliberate: as
  brms", like ranks 13, 15, 19 and 20: brms 2.23.0 fails the same JSS
  call ("passing unknown arguments: threshold") and drops it from
  `stancode()` silently (vigport repro R7). The family argument
  `sratio(threshold = "equidistant")` works since 0.67.0. Postfix's
  brief names it; the most it calls for is a refusal that names the
  family argument.
- **Rank 18** (`loo()` of several models on draws) is open:
  `dev/brms-api-diff.md` (c) lists `compare =` and `model_names =` on
  `loo()` and `waic()` as Small open work, not a decision.
- **Rank 24 is the instrument**, not a gap: the hand row
  `brms_families` "ML: wiener" probes core alone, and frmtmb.eam
  exports `wiener()`, which fits brms's spelling (section "Families").
- Group-level `mo()` (rank 6 at 0.67.0) is gone: 0.70.0 fits
  `(mo(income) | city)`.
- Not in the vignette ranking, because the call still runs: two `mi()`
  fits of `brms_missings` lost their standard errors at 0.70.0 (defect
  R1 below, commit ea711734, owner optstart). It heads the defect table
  at the end of the block, with R2 to R4 (owner postfix).

<!-- BEGIN generated: remeasure-out/r8/gaps.txt -->
### frmtmb 0.70.0, generated by dev/remeasure-gaps.R r8

Mechanical port, spell pass: 35 rows fail (2 model, 29 post, 4 other; 8 of them cascades). Hand translation: 81 open rows (57 ML, 24 SAMPLE).

Ranked by mechanical calls blocked (model + post + other), then by all calls.

| rank | cause | status | lane | mech model/post/other (cascade) | hand ML/SAMPLE (cascade) |
|---|---|---|---|---|---|
| 1 | loo(), LOO(), waic() on a maximum-likelihood fit | deliberate: no posterior | no lane | 0/9/0 (0) | 11/0 (0) |
| 2 | brms's Stan-side custom_family()/stanvar()/expose_functions() | deliberate: Stan-side | no lane | 1/4/2 (5) | 3/1 (0) |
| 3 | frm_multiple() result has no draws (plot, as_draws_array, nchains) | deliberate: no posterior | no lane | 0/3/2 (2) | 4/0 (0) |
| 4 | plot() of a fit refuses plot.brmsfit()'s N, variable, regex (posterior panels) | deliberate: no posterior | no lane | 0/5/0 (0) | 3/0 (0) |
| 5 | add_criterion(fit, 'loo') on a maximum-likelihood fit | deliberate: no posterior | no lane | 0/3/0 (0) | 3/0 (0) |
| 6 | pp_check() on a multivariate model (resp =) | open | postfix | 0/2/0 (0) | 2/1 (0) |
| 7 | threshold = as an frm() argument (brms 2.23.0 fails the same call) | deliberate: as brms | postfix (at most a refusal naming the family argument) | 1/1/0 (1) | 2/0 (0) |
| 8 | summary(waic = TRUE), which brms's summary() swallows | open | postfix | 0/1/0 (0) | 2/2 (0) |
| 9 | bayes_R2() on a maximum-likelihood fit | deliberate: no posterior | no lane | 0/1/0 (0) | 1/0 (0) |
| 10 | brms-only names with no analog (threading, backend, shiny, marginal_effects) | deliberate: Stan-side or retired in brms | no lane | 0/0/0 (0) | 7/0 (0) |
| 11 | conditional_effects(method = 'predict') on draws | open | postfix | 0/0/0 (0) | 0/5 (0) |
| 12 | LOO(), WAIC() on draws (spellings brms deprecated) | deliberate: deprecated in brms | no lane | 0/0/0 (0) | 0/4 (0) |
| 13 | class = sd without nlpar on a nonlinear model (brms refuses it too) | deliberate: as brms | no lane | 0/0/0 (0) | 0/4 (3) |
| 14 | stancode(), standata() on frmtmb objects | deliberate: Stan-side | no lane | 0/0/0 (0) | 2/1 (0) |
| 15 | hypothesis() with no relation (brms refuses it too) | deliberate: as brms | no lane | 0/0/0 (0) | 3/0 (0) |
| 16 | families brms has and core lacks: discrete_weibull, frechet, gen_extreme_value | open | famuni | 0/0/0 (0) | 3/0 (0) |
| 17 | families brms has and core lacks: dirichlet, logistic_normal | open | famsimplex | 0/0/0 (0) | 2/0 (0) |
| 18 | loo() of several models on draws (brms compares them in one call; dev/brms-api-diff.md (c) lists it as Small open work) | open | postfix | 0/0/0 (0) | 0/2 (0) |
| 19 | class Intercept, dpar sigma without a sigma formula (brms refuses it too) | deliberate: as brms | no lane | 0/0/0 (0) | 1/0 (0) |
| 20 | a simo prior other than dirichlet (brms refuses it too) | deliberate: as brms | no lane | 0/0/0 (0) | 1/0 (0) |
| 21 | brms attached after frmtmb masks bf() | deliberate: coexistence, the refusal says how | no lane | 0/0/0 (0) | 1/0 (0) |
| 22 | rescor_matrix() on draws (a frmtmb function; the refusal names the fit) | deliberate: frmtmb API | no lane | 0/0/0 (0) | 0/1 (0) |
| 23 | frm_sample() of a frm_multiple() result (pooled draws) | open | postfix | 0/0/0 (0) | 0/1 (0) |
| 24 | wiener scored missing: the hand row probes core, frmtmb.eam exports wiener() | instrument: covered by frmtmb.eam | famuni (parity) | 0/0/0 (0) | 1/0 (0) |
| 25 | mixture: brms's mu1 ~ formula beside the main formula | open | no lane | 0/0/0 (0) | 1/0 (0) |
| 26 | multinomial() needs K, which brms reads from the data | open | no lane | 0/0/0 (0) | 1/0 (0) |
| 27 | update() on draws | open | no lane | 0/0/0 (0) | 0/1 (0) |
| 28 | log_lik() on an mi() model | open | no lane | 0/0/0 (0) | 0/1 (0) |
| 29 | horseshoe prior | open | no lane | 0/0/0 (0) | 1/0 (0) |
| 30 | custom family without a simulator (pp_check) | open | no lane | 0/0/0 (0) | 1/0 (0) |
| 31 | fixef() returns brms's matrix since 0.61.0 (hand translation stale) | instrument: stale translation | no lane | 0/0/0 (0) | 1/0 (0) |

Row ids per cause (mechanical: vignette.chunk.expr without `brms_`; hand: vignette:n, the n-th non-data bv() row of the script):

1. loo(), LOO(), waic() on a maximum-likelihood fit
   mechanical: overview.8.1, multilevel.7.1, multilevel.20.1, nonlinear.7.1, nonlinear.17.1, phylogenetics.23.1, monotonic.8.1, multivariate.8.1, customfamilies.11.1 
   hand: coexistence:21 (ML: loo(gaussian) with brms attached); coexistence:22 (ML: loo(binomial) with brms attached); customfamilies:12 (ML: loo(fit1, fit2)); monotonic:10 (ML: loo(fit1, fit2, fit3)); multilevel:6 (ML: LOO(fit_zinb1, fit_zinb2)); multilevel:26 (ML: LOO(fit_loss1, fit_loss2)); multivariate:10 (ML: loo(fit1, fit2)); nonlinear:12 (ML: loo(fit1, fit2)); nonlinear:30 (ML: loo(fit_ir1, fit_ir2)); overview:17 (ML: LOO(fit1, fit2)); phylogenetics:25 (ML: loo(model_pois, model_normal)) 
2. brms's Stan-side custom_family()/stanvar()/expose_functions()
   mechanical: customfamilies.4.1, customfamilies.6.1, customfamilies.7.1, customfamilies.8.1, customfamilies.9.1, customfamilies.13.1, customfamilies.15.1 
   hand: customfamilies:5 (ML: stan_funs <- '<Stan source>'); customfamilies:6 (ML: stanvars <- stanvar(scode = stan_funs, block = 'function); customfamilies:10 (ML: expose_functions(fit2, vectorize = TRUE)); customfamilies:27 (SAMPLE: expose_functions(fit2 draws)) 
3. frm_multiple() result has no draws (plot, as_draws_array, nchains)
   mechanical: missings.5.1, missings.6.2, missings.6.3, missings.6.4, missings.6.5 
   hand: missings:4 (ML: plot(fit_imp1, variable = '^b', regex = TRUE)); missings:5 (ML: as_draws_array(fit_imp1)); missings:6 (ML: nchains(fit_imp1)); missings:7 (ML: posterior::nchains(fit_imp1)) 
4. plot() of a fit refuses plot.brmsfit()'s N, variable, regex (posterior panels)
   mechanical: distreg.3.2, nonlinear.10.2, phylogenetics.4.2, monotonic.3.2, monotonic.11.1 
   hand: distreg:3 (ML: plot(fit1, N = 2, ask = FALSE)); monotonic:4 (ML: plot(fit1, variable = 'simo', regex = TRUE)); monotonic:17 (ML: plot(fit4, variable = 'prior_simo', regex = TRUE, N = 3)) 
5. add_criterion(fit, 'loo') on a maximum-likelihood fit
   mechanical: multivariate.3.1, multivariate.7.1, multivariate.10.1 
   hand: multivariate:2 (ML: fit1 <- add_criterion(fit1, 'loo')); multivariate:8 (ML: fit2 <- add_criterion(fit2, 'loo')); multivariate:13 (ML: fit3 <- add_criterion(fit3, 'loo')) 
6. pp_check() on a multivariate model (resp =)
   mechanical: multivariate.4.1, multivariate.4.2 
   hand: missings:30 (SAMPLE: pp_check(s_imp2)); multivariate:4 (ML: pp_check(fit1, resp = 'tarsus')); multivariate:5 (ML: pp_check(fit1, resp = 'back')) 
7. threshold = as an frm() argument (brms 2.23.0 fails the same call)
   mechanical: overview.11.1, overview.12.1 
   hand: monotonic:24 (ML: appendix cumulative() with threshold = 'equidistant'); overview:21 (ML: fit4 (the vignette's exact call)) 
8. summary(waic = TRUE), which brms's summary() swallows
   mechanical: overview.5.1 
   hand: overview:9 (ML: summary(fit1, waic = TRUE)); overview:24 (ML: summary(fit4, waic = TRUE)); overview:28 (SAMPLE: summary(fit1, waic = TRUE)); overview:45 (SAMPLE: summary(fit4, waic = TRUE)) 
9. bayes_R2() on a maximum-likelihood fit
   mechanical: multivariate.5.1 
   hand: multivariate:6 (ML: bayes_R2(fit1)) 
10. brms-only names with no analog (threading, backend, shiny, marginal_effects)
   hand: overview:12 (ML: marginal_effects(fit1)); overview:13 (ML: launch_shiny(fit1)); threading:1 (ML: threading(threads, grainsize, static)); threading:2 (ML: brm(..., threads = threading(2))); threading:3 (ML: reduce_sum in the generated Stan code); threading:4 (ML: brm(..., backend = 'cmdstanr')); threading:5 (ML: grainsize tuning) 
11. conditional_effects(method = 'predict') on draws
   hand: distreg:46 (SAMPLE: conditional_effects(fit1, method = 'predict')); missings:29 (SAMPLE: conditional_effects(s_imp2, method = 'predict')); multivariate:42 (SAMPLE: conditional_effects(..., method = 'predict')); nonlinear:46 (SAMPLE: the vignette's exact call (method = 'predict')); overview:32 (SAMPLE: conditional_effects(fit1, method = 'predict')) 
12. LOO(), WAIC() on draws (spellings brms deprecated)
   hand: coexistence:65 (SAMPLE: LOO(draws) with brms attached); coexistence:66 (SAMPLE: WAIC(draws) with brms attached); multilevel:38 (SAMPLE: LOO(fit_zinb1, fit_zinb2)); overview:38 (SAMPLE: LOO(fit1, fit2)) 
13. class = sd without nlpar on a nonlinear model (brms refuses it too)
   hand: multilevel:48 (SAMPLE: fit_loss1); multilevel:50 (SAMPLE: conditional_effects(fit_loss1)); multilevel:51 (SAMPLE: conditional_effects(conditions, re_formula = NULL, m); multilevel:52 (SAMPLE: posterior_predict(fit_loss1) [the named substitute]) 
14. stancode(), standata() on frmtmb objects
   hand: customfamilies:29 (SAMPLE: stancode(s2)); overview:7 (ML: stancode(fit1)); overview:8 (ML: standata(fit1)) 
15. hypothesis() with no relation (brms refuses it too)
   hand: coexistence:9 (ML: hypothesis(gaussian) with brms attached); coexistence:10 (ML: hypothesis(binomial) with brms attached); missings:11 (ML: hypothesis(fit_imp1, 'age')) 
16. families brms has and core lacks: discrete_weibull, frechet, gen_extreme_value
   hand: families:9 (ML: discrete_weibull); families:17 (ML: frechet); families:18 (ML: gen_extreme_value) 
17. families brms has and core lacks: dirichlet, logistic_normal
   hand: families:24 (ML: dirichlet); families:25 (ML: logistic_normal) 
18. loo() of several models on draws (brms compares them in one call; dev/brms-api-diff.md (c) lists it as Small open work)
   hand: multivariate:37 (SAMPLE: loo(fit1, fit2)); phylogenetics:58 (SAMPLE: loo(model_pois, model_normal)) 
19. class Intercept, dpar sigma without a sigma formula (brms refuses it too)
   hand: phylogenetics:30 (ML: frm(..., prior = ) [the penalized fit]) 
20. a simo prior other than dirichlet (brms refuses it too)
   hand: monotonic:13 (ML: set_prior(class = 'simo') with a supported density) 
21. brms attached after frmtmb masks bf()
   hand: coexistence:5 (ML: frm(bf(...)) where bf() is brms's) 
22. rescor_matrix() on draws (a frmtmb function; the refusal names the fit)
   hand: multivariate:32 (SAMPLE: rescor_matrix(fit1)) 
23. frm_sample() of a frm_multiple() result (pooled draws)
   hand: missings:18 (SAMPLE: fit_imp1 (the pooled object)) 
24. wiener scored missing: the hand row probes core, frmtmb.eam exports wiener()
   hand: families:21 (ML: wiener) 
25. mixture: brms's mu1 ~ formula beside the main formula
   hand: distreg:20 (ML: appendix bf(y ~ 1, mu1 ~ x1, mu2 ~ x1)) 
26. multinomial() needs K, which brms reads from the data
   hand: families:32 (ML: multinomial [brms spelling]) 
27. update() on draws
   hand: phylogenetics:43 (SAMPLE: model_repeat2 (update)) 
28. log_lik() on an mi() model
   hand: missings:31 (SAMPLE: loo(s_imp2)) 
29. horseshoe prior
   hand: overview:5 (ML: set_prior('horseshoe(1)')) 
30. custom family without a simulator (pp_check)
   hand: customfamilies:16 (ML: pp_check() on a custom family with no sim) 
31. fixef() returns brms's matrix since 0.61.0 (hand translation stale)
   hand: missings:8 (ML: per-submodel check (the substitute for summarise_draws)) 

### totals by status

deliberate  causes 16  mechanical rows  32  hand rows  53
open        causes 13  mechanical rows   3  hand rows  26
instrument  causes  2  mechanical rows   0  hand rows   2

### open calls per lane

postfix            mechanical  3  hand 15  (pp_check() on a multivariate model (resp =); summary(waic = TRUE), which brms's summary() swallows; conditional_effects(method = 'predict') on draws; loo() of several models on draws (brms compares them in one call; dev/brms-api-diff.md (c) lists it as Small open work); frm_sample() of a frm_multiple() result (pooled draws))
famuni             mechanical  0  hand  3  (families brms has and core lacks: discrete_weibull, frechet, gen_extreme_value)
famsimplex         mechanical  0  hand  2  (families brms has and core lacks: dirichlet, logistic_normal)
no lane            mechanical  0  hand  6  (mixture: brms's mu1 ~ formula beside the main formula; multinomial() needs K, which brms reads from the data; update() on draws; log_lik() on an mi() model; horseshoe prior; custom family without a simulator (pp_check))

### Defects outside the vignette calls, ranked by severity

| rank | defect | owner | evidence |
|---|---|---|---|
| 1 | R1: two `mi()` fits of brms_missings lose their standard errors at an unchanged optimum (0.70.0, commit ea711734) | optstart | `dev/remeasure-repro-missings.R`, secmp |
| 2 | R2: loo() on draws silently ignores newdata | postfix | methods:469 (bin 2T, hollow) |
| 3 | R3: logLik() of draws has no method | postfix | methods:434 (bin 2T, method) |
| 4 | loo() and waic() of several models on draws | postfix | methods:1015, methods:1016, methods:446, methods:447, methods:448, methods:449 (bin 2T, argument) |
| 5 | R4: bayes_R2() on draws refuses newdata and re_formula | postfix | methods:121 (bin 2T, argument) |
<!-- END generated -->

## Headline: 0.67.0, 0.69.1 and 0.70.0

Every build was RUN by this lane with the same scripts, the 0.67.0 build
reinstalled from 9e902909 into the lane library (frmtmb 0.67.0,
frmtmb.sample 0.15.0), 0.68.1 and 0.69.1 from `rellib-r6` and
`rellib-r7`. The 0.67.0 run reproduces the vigport record exactly: every
model and post row's raw status equals the one the record's history
block implies, and the hand translation's 448 rows equal the record's
row for row (`dev/remeasure-seed-history.R`, block below). The 0.69.1
run reproduces the 0.69.0 release record (post CLEAN 64, post spell 72,
hand 453; `dev/round-20261007.md`). The figures the brief quotes as
0.67.0's (post CLEAN 64, spell 72, hand 453) are 0.69.x's.

<!-- BEGIN generated: remeasure-out/r8/history.txt from ^## totals by build to ^## mechanical port -->
```text
## totals by build

model calls 42, post calls 102, hand rows 539 

        build model_raw_ok model_clean model_spell_ok post_raw_ok post_clean
 0.34.0(v034)           16          16             37          35         35
 0.65.0(v065)           36          NA             NA          57         NA
 0.66.0(v066)           36          NA             NA          62         NA
 0.67.0(v067)           36          36             39          63         63
   0.68.1(r6)           36          36             39          63         63
   0.69.1(r7)           36          36             39          64         64
   0.70.0(r8)           37          37             40          65         65
 post_spell_ok hand_ok
            63      NA
            NA     437
            NA     446
            71     448
            71     448
            72     453
            73     458

```
<!-- END generated -->

Seeding the history from the surviving records, and checking it against
the fresh 0.67.0 run:

<!-- BEGIN generated: remeasure-out/seed-history.txt -->
```text
mechanical rows: 250  same transformed source as 0.34.0: 250 
rows listed in the record: 54  the record says: 54 
raw OK, model rows: 0.34.0=16  0.65.0=36  0.66.0=36  0.67.0=36 
raw OK, post rows:  0.34.0=35  0.65.0=57  0.66.0=62  0.67.0=63 
rebuilt raw OK, model 0.34.0=16 0.65.0=36 0.66.0=36 0.67.0=36
rebuilt raw OK, post  0.34.0=35 0.65.0=57 0.66.0=62 0.67.0=63
rebuilt 0.67.0 against the fresh 0.67.0 run, model and post rows: 0 rows differ
other rows (not in the record's block) where 0.34.0 and the fresh 0.67.0 run differ: 1 

hand rows: 539  ok in the fresh 0.67.0 run: 448 
hand rows listed in the record as differing between builds: 11 
rows: 539   ok: 0.65.0=437  0.66.0=446  0.67.0=448   missing keys: 0  
rebuilt ok: 0.65.0=437 0.66.0=446 0.67.0=448
rebuilt 0.67.0 against the fresh run: 0 rows differ

wrote outcomes-mech.tsv (250 rows), outcomes-hand.tsv (539 rows), outcome-builds.tsv
```
<!-- END generated -->

- **Mechanical model calls:** 37 of 42 CLEAN at 0.70.0 (36 at 0.67.0);
  40 of 42 run in the spell pass (39). The one more is group-level
  `mo()`, `brms_monotonic.16.1`.
- **Post-processing calls:** 65 of 102 CLEAN (63 at 0.67.0, 64 at
  0.69.1); spell pass 73 of 102 (71, 72).
- **Hand translation:** 458 of 539 rows run (448 at 0.67.0, 453 at
  0.69.1).
- **Keep-prior transform:** 40 of 42 model and 73 of 102 post calls run
  (37 and 70 at 0.67.0). Its two 0.67.0 losses run now: the Dirichlet
  prior on a `mo()` simplex (0.70.0) and the stale `lkj(2)` prior that
  `update()` drops (0.69.0).
- **Families named in brms_families:** 34 of 40 fit in core and 35 of
  40 with the extensions (frmtmb.eam's `wiener()`), each fitted on
  simulated data. At 0.67.0 the instrument counted names only, 34 of 40.
- **Estimates:** on the 34 models where both sides converged on every
  build, the median `|z|` is 0.047 over 124 coefficients on 0.67.0 and
  on 0.70.0 alike (the 90th percentile moves from 0.208 to 0.213 at
  0.68.1, with the `s()` null-space fix). Each build's own set gives
  0.0488 over 130 (0.67.0) and 0.0458 over 129 (0.70.0), but the sets
  differ by two models, so those two figures are not comparable. None
  past 1.
- **Regressions:** none in any call's outcome. Two `mi()` fits regress
  in their standard errors, which the outcome does not see (R1).

## The mechanical port

`brm` to `frm`, MCMC arguments deleted, scored by the 0.34.0 rules
(`classify()` in `dev/brms-port/port-lib.R`), unchanged:

<!-- BEGIN generated: remeasure-out/blocks/headline.md -->
| calls | outcome | 0.34.0 | 0.67.0 | 0.70.0 |
|---|---|---|---|---|
| model (42) | CLEAN | 16 (38%) | 36 (86%) | 37 (88%) |
| model (42) | SPELLING | 21 (50%) | 3 (7%) | 3 (7%) |
| model (42) | FAIL | 4 (10%) | 2 (5%) | 1 (2%) |
| model (42) | CASCADE | 1 (2%) | 1 (2%) | 1 (2%) |
| post (102) | CLEAN | 35 (34%) | 63 (62%) | 65 (64%) |
| post (102) | SPELLING | 28 (27%) | 8 (8%) | 8 (8%) |
| post (102) | FAIL | 31 (30%) | 25 (25%) | 24 (24%) |
| post (102) | CASCADE | 8 (8%) | 6 (6%) | 5 (5%) |
<!-- END generated -->

Raw, spell and keep-prior passes:

<!-- BEGIN generated: remeasure-out/blocks/passes.md -->
| pass | calls | 0.67.0 | 0.70.0 |
|---|---|---|---|
| raw OK | model (42) | 36 | 37 |
| spell OK | model (42) | 39 | 40 |
| raw OK | post (102) | 63 | 65 |
| spell OK | post (102) | 71 | 73 |
| raw OK | other (106) | 97 | 97 |
| spell OK | other (106) | 97 | 97 |
| keep-prior OK | model (42) | 37 | 40 |
| keep-prior OK | post (102) | 70 | 73 |
<!-- END generated -->

Per vignette, model calls then post calls, as
`total | clean spelling fail cascade`:

<!-- BEGIN generated: remeasure-out/blocks/per-vignette.md -->
| vignette | model 0.67.0 | model 0.70.0 | post 0.67.0 | post 0.70.0 |
|---|---|---|---|---|
| brms_overview | 4 \| 3 0 1 0 | 4 \| 3 0 1 0 | 4 \| 1 0 2 1 | 4 \| 1 0 2 1 |
| brms_multilevel | 8 \| 6 2 0 0 | 8 \| 6 2 0 0 | 13 \| 7 4 2 0 | 13 \| 7 4 2 0 |
| brms_distreg | 4 \| 4 0 0 0 | 4 \| 4 0 0 0 | 12 \| 11 0 1 0 | 12 \| 11 0 1 0 |
| brms_nonlinear | 6 \| 5 1 0 0 | 6 \| 5 1 0 0 | 20 \| 13 4 3 0 | 20 \| 13 4 3 0 |
| brms_phylogenetics | 6 \| 6 0 0 0 | 6 \| 6 0 0 0 | 18 \| 16 0 2 0 | 18 \| 16 0 2 0 |
| brms_monotonic | 6 \| 5 0 1 0 | 6 \| 6 0 0 0 | 11 \| 7 0 3 1 | 11 \| 8 0 3 0 |
| brms_multivariate | 3 \| 3 0 0 0 | 3 \| 3 0 0 0 | 11 \| 4 0 7 0 | 11 \| 4 0 7 0 |
| brms_missings | 3 \| 3 0 0 0 | 3 \| 3 0 0 0 | 7 \| 3 0 4 0 | 7 \| 4 0 3 0 |
| brms_customfamilies | 2 \| 1 0 0 1 | 2 \| 1 0 0 1 | 6 \| 1 0 1 4 | 6 \| 1 0 1 4 |
<!-- END generated -->

The patch table is unchanged in effect: of the explicit patches only the
two nonlinear `start =` values are still needed, and the four retry
rules fire on no expression:

<!-- BEGIN generated: remeasure-out/r8/port/summary.txt from ^## patch necessity .explicit to ^## per vignette -->
```text
## patch necessity (explicit per-expression patches)
brms_nonlinear.2.2         raw=OK      spell=OK      UNNEEDED: the raw line runs
brms_nonlinear.9.1         raw=ERROR   spell=OK      NEEDED
      raw!  The nonlinear fit failed from its default starting values (NA/NaN gradient evaluation). Nonlinear models need starting values in the right region: par_template( 
brms_multilevel.15.3       raw=ERROR   spell=OK      NEEDED
      raw!  The nonlinear fit failed from its default starting values (NA/NaN gradient evaluation). Nonlinear models need starting values in the right region: par_template( 
brms_multilevel.19.2       raw=ERROR   spell=OK      NEEDED
      raw!  object 'fit_loss1' not found 
brms_overview.7.1          raw=OK      spell=OK      UNNEEDED: the raw line runs
brms_phylogenetics.11.1    raw=OK      spell=OK      UNNEEDED: the raw line runs
brms_multivariate.9.1      raw=OK      spell=OK      UNNEEDED: the raw line runs

## patch necessity, need pass (only the patches marked * applied)
brms_nonlinear.2.2           need=OK      UNNEEDED once the patches upstream of it are applied
brms_nonlinear.9.1         * need=OK      
brms_multilevel.15.3       * need=OK      
brms_multilevel.19.2         need=OK      UNNEEDED once the patches upstream of it are applied
brms_overview.7.1            need=OK      UNNEEDED once the patches upstream of it are applied
brms_phylogenetics.11.1      need=OK      UNNEEDED once the patches upstream of it are applied
brms_multivariate.9.1        need=OK      UNNEEDED once the patches upstream of it are applied

## patch necessity (generic retry rules)
default-family       fired on 0 expressions, 0 of them then ran
update-formula-dot   fired on 0 expressions, 0 of them then ran
update-newdata       fired on 0 expressions, 0 of them then ran
family-call          fired on 0 expressions, 0 of them then ran

```
<!-- END generated -->

<!-- BEGIN generated: remeasure-out/blocks/keepprior-diff.md -->
rows of 0.70.0 whose outcome differs between the
raw and the keep-prior transform:

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
<!-- END generated -->

## The hand translation

Labels are fixed in the scripts at their 0.42.0 values, so the measure
is the drift of each row from its label (`dev/brms-vignettes/_drift.R`):

<!-- BEGIN generated: remeasure-out/blocks/hand-drift.md -->
| label at 0.42.0 | rows | ran at 0.70.0 | erred at 0.70.0 |
|---|---|---|---|
| CLEAN | 233 | 229 | 4 |
| SPELLING | 90 | 86 | 4 |
| BEHAVIOR | 117 | 109 | 8 |
| MISSING | 61 | 17 | 44 |
| REFUSAL | 38 | 17 | 21 |
| all | 539 | 458 | 81 |

| drift | ML | SAMPLE | all |
|---|---|---|---|
| AS-LABELED | 237 | 187 | 424 |
| FAILS-NOW | 11 | 5 | 16 |
| RUNS-NOW | 23 | 11 | 34 |
| STILL-OPEN | 46 | 19 | 65 |
<!-- END generated -->

## Every row that changed outcome since 0.67.0, and the regression check

By construction in `dev/remeasure-history.R`, from the tracked outcome
tables `dev/brms-port/outcomes-mech.tsv` and `outcomes-hand.tsv`, which
now carry a column per build (0.34.0, 0.65.0 and 0.66.0 rebuilt from
the vigport record; 0.67.0, 0.68.1, 0.69.1 and 0.70.0 run). A row that
ran on any build from 0.67.0 on and fails on 0.70.0 is a regression.

<!-- BEGIN generated: remeasure-out/r8/history.txt from ^## mechanical port: rows to ^## end -->
```text
## mechanical port: rows whose raw or spell status changed since 0.67.0

rows: 3 of 250 
columns: raw/spell per build, 0.67.0(v067) 0.68.1(r6) 0.69.1(r7) 0.70.0(r8) 

brms_monotonic.16.1      model ERROR/ERROR ERROR/ERROR ERROR/ERROR OK/OK  fit6 <- frm(ls ~ mo(income) * age + (mo(income) | city), data = dat)
brms_monotonic.17.1      post  ERROR/ERROR ERROR/ERROR ERROR/ERROR OK/OK  summary(fit6)
brms_missings.7.1        post  ERROR/ERROR ERROR/ERROR OK/OK OK/OK  conditional_effects(fit_imp1, "age:chl")

## hand translation: rows whose outcome changed since 0.67.0

rows: 10 of 539 
columns: ok per build, 0.67.0(v067) 0.68.1(r6) 0.69.1(r7) 0.70.0(r8) 

ERR ERR ok  ok   [MISSING] brms_distreg | 22 | ML: appendix pp_mixture(fit_mix2)
ERR ERR ok  ok   [MISSING] brms_distreg | 27 | SAMPLE: plot(fit1, N = 2, ask = FALSE)
ERR ERR ok  ok   [REFUSAL] brms_missings | 9 | ML: conditional_effects(fit_imp1, 'age:chl')
ERR ERR ERR ok   [MISSING] brms_monotonic | 12 | ML: prior4 = prior(dirichlet(c(2, 1, 1)), class = 'simo', coef = 'moincome1')
ERR ERR ERR ok   [MISSING] brms_monotonic | 14 | ML: fit4 [the vignette's exact translation]
ERR ERR ERR ok   [REFUSAL] brms_monotonic | 21 | ML: fit6
ERR ERR ok  ok   [MISSING] brms_monotonic | 28 | SAMPLE: plot(fit1, variable = 'simo', regex = TRUE)
ERR ERR ERR ok   [MISSING] brms_monotonic | 34 | SAMPLE: fit4 [the vignette's exact translation]
ERR ERR ERR ok   [REFUSAL] brms_monotonic | 40 | SAMPLE: fit6
ERR ERR ok  ok   [MISSING] brms_multivariate | 28 | SAMPLE: fit1 <- add_criterion(fit1, 'loo')

## regressions: ran on a build from 0.67.0 on, fails on 0.70.0(r8) 

mechanical (raw or spell pass): 0 
mechanical raw pass only (CLEAN lost): 0 
hand translation: 0 

## against 0.34.0: ran then (raw or behind a patch), fail now

rows: 6 
  brms_overview.5.1        CLEAN              summary(fit1, waic = TRUE)
      now: summary() cannot honor `waic`: WAIC averages the pointwise likelihood over posterior draws, and a maximum likelihood fit has none (brms 2.23.0's summa
  brms_distreg.3.2         CLEAN              plot(fit1, N = 2, ask = FALSE)
      now: plot() cannot honor `N`: brms sets how many parameters' posterior panels go on one page with it. plot() of a maximum-likelihood fit draws residual dia
  brms_nonlinear.10.2      SPELLING: upstream plot(fit_loss, N = 3, ask = FALSE)
      now: plot() cannot honor `N`: brms sets how many parameters' posterior panels go on one page with it. plot() of a maximum-likelihood fit draws residual dia
  brms_phylogenetics.4.2   CLEAN              plot(model_simple, N = 2, ask = FALSE)
      now: plot() cannot honor `N`: brms sets how many parameters' posterior panels go on one page with it. plot() of a maximum-likelihood fit draws residual dia
  brms_monotonic.3.2       SPELLING: upstream plot(fit1, variable = "simo", regex = TRUE)
      now: plot() cannot honor `variable`: brms chooses the parameters whose posterior draws it plots with it. plot() of a maximum-likelihood fit draws residual 
  brms_monotonic.11.1      SPELLING: upstream plot(fit4, variable = "prior_simo", regex = TRUE, N = 3)
      now: plot() cannot honor `variable`: brms chooses the parameters whose posterior draws it plots with it. plot() of a maximum-likelihood fit draws residual 

```
<!-- END generated -->

The 6 rows lost against 0.34.0 are the 0.58.0 dots refusal, now refused
by name since 0.69.0 (`plot()` cannot honor `N` or `variable`;
`summary()` cannot honor `waic`). They failed on 0.65.0 already; the
vigport record dates them.

What moved, read row by row:

- `brms_monotonic.16.1` and `.17.1`, and the hand rows `ML: fit6` and
  `SAMPLE: fit6`: group-level `mo()` fits (0.70.0, `mo()` grammar).
- `brms_monotonic` hand rows 12, 14 and 34: brms's `dirichlet()` prior
  on a `mo()` simplex is accepted (0.70.0).
- `brms_missings.7.1` and its hand row: the pooled
  `conditional_effects()` of a `frm_multiple()` result (0.69.0).
- hand `pp_mixture(fit_mix2)`, `SAMPLE: plot(fit1, N = 2, ask = FALSE)`,
  `SAMPLE: add_criterion(fit1, 'loo')`, `SAMPLE: plot(fit1, variable =
  'simo', regex = TRUE)`: lane surface and optima at 0.69.0.

## Estimate agreement with brms

`dev/brms-port/plausibility.R`: each model fitted by brms itself on the
same rows (`dev/brms-port/brms-fit.R`, every expression seeded by
`port_seed(id)`, 2 chains of 1000, `brms_overview` at the vignette's
own budget). The brms fits were made once by this lane
(`dev/remeasure-out/brms-ref`) and every build is compared against the
same ones, so a difference between builds is frmtmb's. `z` is
`(frmtmb - brms posterior mean) / brms posterior SD`.

<!-- BEGIN generated: remeasure-out/r8/plausibility.txt from ^models fitted by both to ^## detail -->
```text
models fitted by both with matched coefficients: 40 of 42 model calls
matched coefficients: 153; abs z: median 0.0558, 90th pct 0.292, max 0.914
coefficients with abs z > 1: 0; > 2: 0
brms fits with max Rhat > 1.05: 2 (brms_overview.3.1, brms_missings.3.1); with divergences: 4
converged brms fits only (max Rhat <= 1.05): 38 models, 142 coefficients; abs z median 0.0521, 90th pct 0.240, max 0.914; > 1: 0
frmtmb fits that did not converge (opt$convergence != 0): 3
  brms_multilevel.19.2: code 1, Optimizer did not report convergence: false convergence (8). This is a nonlinear model; un
  brms_monotonic.16.1: code 1, 
  brms_missings.8.2: code 1, Optimizer did not report convergence: false convergence (8)
both converged: 35 models, 129 coefficients; abs z median 0.0458, 90th pct 0.216, max 0.914; > 1: 0
variance components, brms converged: 65; abs z median 0.261, max 1.68; > 1: 6
variance components, both converged: 56; abs z median 0.233, max 1.68; > 1: 4
variance components past abs z 1, both converged:
  brms_overview.7.1      sd_patient__Intercept          frmtmb 0.000251  brms    0.408 (sd 0.243)  z -1.68  within one posterior SD of zero
  brms_phylogenetics.15.1 sd_phylo__Intercept            frmtmb   0.0107  brms   0.0677 (sd 0.0405)  z -1.41  within one posterior SD of zero
  brms_phylogenetics.20.1 sd_phylo__Intercept            frmtmb   0.0574  brms    0.714 (sd 0.539)  z -1.22  within one posterior SD of zero
  brms_missings.11.2     sd_residual____bmi             frmtmb     2.53  brms     3.30 (sd 0.737)  z -1.05  residual SD
the 10 largest abs z over the coefficients where both converged:
  brms_multilevel.12.2     sigma_Intercept            frmtmb      0.685  brms      0.713 (sd 0.0313)  z -0.914
  brms_missings.11.2       bmi_age                    frmtmb     -0.100  brms      -4.93 (sd 5.86)  z 0.824
  brms_distreg.12.1        sigma_sx0_1                frmtmb     -0.448  brms      -2.25 (sd 2.42)  z 0.744
  brms_multilevel.12.2     sigma_t2areayearc_2        frmtmb     0.0814  brms     0.0547 (sd 0.0418)  z 0.640
  brms_missings.11.2       bmi_michl:age              frmtmb    -0.0211  brms   -0.00548 (sd 0.0261)  z -0.600
  brms_overview.10.1       Intercept[3]               frmtmb       4.91  brms       5.15 (sd 0.450)  z -0.520
  brms_overview.10.1       Intercept[2]               frmtmb       3.68  brms       3.84 (sd 0.317)  z -0.506
  brms_phylogenetics.7.1   Intercept                  frmtmb       38.2  brms       36.1 (sd 7.32)  z 0.292
  brms_phylogenetics.11.1  Intercept                  frmtmb       38.2  brms        36. (sd 7.87)  z 0.292
  brms_distreg.2.1         sigma_Intercept            frmtmb     -0.183  brms     -0.144 (sd 0.134)  z -0.291

```
<!-- END generated -->

Build against build (`dev/remeasure-plauscmp.R`):

<!-- BEGIN generated: remeasure-out/r8/plauscmp.txt from ^models whose status to ^## end -->
```text
models whose status, optimizer code or max abs z (by more than 0.01) differ between builds: 4 

### brms_distreg.12.1 
  v067   code 1 Optimizer did not report convergence: singular convergence (7)
         Intercept=15.2(z 0.0591) sigma_Intercept=0.722(z 0.0233) sx1_1=1.83(z -2.3) sx2_1=5.47(z   -2) sigma_sx0_1=-0.0737(z 0.899) 
  r6     code 1 Optimizer did not report convergence: singular convergence (7)
         Intercept=15.2(z 0.0591) sigma_Intercept=0.722(z 0.0233) sx1_1=9.83(z 0.0301) sx2_1=29.7(z -0.0363) sigma_sx0_1=-0.448(z 0.744) 
  r7     code 1 
         Intercept=15.2(z 0.0591) sigma_Intercept=0.722(z 0.0233) sx1_1=9.83(z 0.0301) sx2_1=29.7(z -0.0363) sigma_sx0_1=-0.448(z 0.744) 
  r8     code 0 
         Intercept=15.2(z 0.0591) sigma_Intercept=0.722(z 0.0233) sx1_1=9.83(z 0.0301) sx2_1=29.7(z -0.0363) sigma_sx0_1=-0.448(z 0.744) 

### brms_monotonic.16.1 
  v067   code NA 
  r6     code NA 
  r7     code NA 
  r8     code 1 
         Intercept=47.5(z -0.056) age=-0.281(z -0.123) moincome=12.2(z -0.197) moincome:age=0.0682(z 0.196) 

### brms_multivariate.9.3 
  v067   code 0 
         tarsus_Intercept=-0.412(z 0.0273) back_Intercept=0.00129(z -0.0255) tarsus_sexMale=0.771(z 0.00511) tarsus_sexUNK=0.218(z 0.00782) sigma_tarsus_sexFem=-0.308(z -0.181) sigma_tarsus_sexMale=-0.249(z -0.0511) sigma_tarsus_sexUNK=-0.408(z -0.112) back_shatchdate_1=-0.144(z 0.0571) 
  r6     code 0 
         tarsus_Intercept=-0.412(z 0.0273) back_Intercept=0.00129(z -0.0255) tarsus_sexMale=0.771(z 0.0051) tarsus_sexUNK=0.218(z 0.00781) sigma_tarsus_sexFem=-0.308(z -0.181) sigma_tarsus_sexMale=-0.249(z -0.0511) sigma_tarsus_sexUNK=-0.408(z -0.112) back_shatchdate_1=-1.01(z -0.221) 
  r7     code 0 
         tarsus_Intercept=-0.412(z 0.0273) back_Intercept=0.00129(z -0.0255) tarsus_sexMale=0.771(z 0.0051) tarsus_sexUNK=0.218(z 0.00781) sigma_tarsus_sexFem=-0.308(z -0.181) sigma_tarsus_sexMale=-0.249(z -0.0511) sigma_tarsus_sexUNK=-0.408(z -0.112) back_shatchdate_1=-1.01(z -0.221) 
  r8     code 0 
         tarsus_Intercept=-0.412(z 0.0273) back_Intercept=0.00129(z -0.0255) tarsus_sexMale=0.771(z 0.0051) tarsus_sexUNK=0.218(z 0.00781) sigma_tarsus_sexFem=-0.308(z -0.181) sigma_tarsus_sexMale=-0.249(z -0.0511) sigma_tarsus_sexUNK=-0.408(z -0.112) back_shatchdate_1=-1.01(z -0.221) 

### brms_missings.8.2 
  v067   code 0 
         bmi_Intercept=13.4(z -0.0242) chl_Intercept= 140(z 0.0666) bmi_age=-0.099(z 0.783) chl_age=  31(z 0.00913) bmi_michl=0.109(z -0.241) bmi_michl:age=-0.0211(z -0.562) 
  r6     code 0 
         bmi_Intercept=13.4(z -0.0242) chl_Intercept= 140(z 0.0666) bmi_age=-0.099(z 0.783) chl_age=  31(z 0.00913) bmi_michl=0.109(z -0.241) bmi_michl:age=-0.0211(z -0.562) 
  r7     code 0 
         bmi_Intercept=13.4(z -0.0242) chl_Intercept= 140(z 0.0666) bmi_age=-0.099(z 0.783) chl_age=  31(z 0.00913) bmi_michl=0.109(z -0.241) bmi_michl:age=-0.0211(z -0.562) 
  r8     code 1 Optimizer did not report convergence: false convergence (8)
         bmi_Intercept=13.4(z -0.0242) chl_Intercept= 140(z 0.0666) bmi_age=-0.099(z 0.783) chl_age=  31(z 0.00913) bmi_michl=0.109(z -0.241) bmi_michl:age=-0.0211(z -0.562) 

## both converged, per build (each build's own set)

v067   35 models, 130 coefficients; abs z median 0.0488, 90th pct 0.229, max 0.914; > 1: 0
r6     35 models, 130 coefficients; abs z median 0.0488, 90th pct 0.229, max 0.914; > 1: 0
r7     35 models, 130 coefficients; abs z median 0.0488, 90th pct 0.229, max 0.914; > 1: 0
r8     35 models, 129 coefficients; abs z median 0.0458, 90th pct 0.216, max 0.914; > 1: 0

## both converged on EVERY build: the common set

v067   not in the common set: brms_missings.8.2
r6     not in the common set: brms_missings.8.2
r7     not in the common set: brms_missings.8.2
r8     not in the common set: brms_distreg.12.1
v067   34 models, 124 coefficients; abs z median 0.047, 90th pct 0.208, max 0.914; > 1: 0
r6     34 models, 124 coefficients; abs z median 0.047, 90th pct 0.213, max 0.914; > 1: 0
r7     34 models, 124 coefficients; abs z median 0.047, 90th pct 0.213, max 0.914; > 1: 0
r8     34 models, 124 coefficients; abs z median 0.047, 90th pct 0.213, max 0.914; > 1: 0

```
<!-- END generated -->

Reading it:

- **Defect 1 of the vigport record is fixed** (the `s()` null-space
  coefficient on another scale than brms's): from 0.68.1 on, `sx1_1`
  and `sx2_1` of `brms_distreg.12.1` read 9.83 and 29.7 with `z` +0.030
  and -0.036, and `back_shatchdate_1` of `brms_multivariate.9.3` moves
  from -0.144 to -1.01 (z -0.221). These are the values the vigport
  record computed on brms's column.
- `brms_distreg.12.1` converges (code 0) at 0.70.0, where 0.67.0 to
  0.69.1 reported singular convergence: it enters the both-converged
  set (5 coefficients).
- `brms_missings.8.2` leaves the set: same optimum to 9 digits, but
  0.70.0 reports false convergence (8) and loses every standard error
  (R1). With `brms_distreg.12.1` entering, the two builds' own sets
  differ by two models, which is why their medians (0.0488, 0.0458) are
  not comparable. On the common set of 124 coefficients both read
  0.047; between 0.69.1 and 0.70.0 the common-set figures are equal.
- `brms_monotonic.16.1`, new, stops with singular convergence (7) and
  frmtmb's boundary message: `cor_city__Intercept__moincome` is at the
  limit of its range (`dev/remeasure-out/probe/replay-monotonic16.txt`).
  Its coefficients agree with brms (max `|z|` 0.197).
- The median is at the Monte Carlo floor of 1000 draws (about 0.03 to
  0.045 posterior SD), as in the vigport record.

Standard errors, build against build (`dev/remeasure-secmp.R`; a fit
is listed when an estimate or a standard error moved by more than 1%,
an SE turned non-finite, or the optimizer code changed):

<!-- BEGIN generated: remeasure-out/r8/secmp.txt -->
```text
fits compared: 40 ; builds: v067 r6 r7 r8 

brms_distreg.12.1 (against v067)
  r6     code 1 -> 1; max rel change: estimate 5.08, SE 0; SE lost: Intercept sigma_Intercept sx1_1 sx2_1 sigma_sx0_1; SE gained: -
  r7     code 1 -> 1; max rel change: estimate 5.08, SE 5.08; SE lost: -; SE gained: -
  r8     code 1 -> 0; max rel change: estimate 5.08, SE 5.08; SE lost: -; SE gained: -


brms_missings.11.2 (against v067)
  r8     code 0 -> 0; max rel change: estimate 0.000114, SE 0.349; SE lost: -; SE gained: -


brms_missings.8.2 (against v067)
  r8     code 0 -> 1; max rel change: estimate 0.000262, SE 0; SE lost: bmi_Intercept chl_Intercept bmi_age chl_age bmi_michl bmi_michl:age; SE gained: -


brms_monotonic.16.1 (against v067)
  r6     does not run
  r7     does not run
  r8     runs (not in the first build)


brms_multivariate.9.3 (against v067)
  r6     code 0 -> 0; max rel change: estimate 5.98, SE 5.98; SE lost: -; SE gained: -
  r7     code 0 -> 0; max rel change: estimate 5.98, SE 5.98; SE lost: -; SE gained: -
  r8     code 0 -> 0; max rel change: estimate 5.98, SE 5.98; SE lost: -; SE gained: -


fits with a change: 5 
```
<!-- END generated -->

## Families

`dev/brms-port/families-coverage.R`, fixed twice. First, the name
probe: the 0.34.0 probe of core is kept, so its count compares across
builds, and every installed `frmtmb.*` package is probed too. A name
counts for an extension when the extension EXPORTS a function of that
name whose default call `as_frmtmb_family()` accepts; `covered_by` names
the package. Its three cases are observed: `wiener` found in frmtmb.eam,
`dirichlet` found nowhere, and an eam export that is not a family
(`ddm_simulate`) not counted (last block of the wiener probe below).
Second (fix round 1), a name is not a fit: every named family is FITTED
by `frm()`, `y ~ x` on data simulated from a model of its kind (seed
20261013, n = 300; the review's generators,
`dev/remeasure-rev-famfit.R`), the family passed by brms's name through
`frmtmb::brmsfamily()`, and counts only when a fit with a finite
log-likelihood comes back. `multinomial` is fitted through
`multinomial(K = 3)` on a matrix response with `trials()`; brms's bare
`multinomial()` fails in frmtmb (gap rank 26). `wiener` is fitted on
`frmtmb.eam::ddm_simulate()` data.

<!-- BEGIN generated: remeasure-out/r8/families.txt from ^ +family to ^frmtmb registry -->
```text
                    family  core covered_by                          fit  fits
                      acat  TRUE     frmtmb fit, code 0, logLik -396.780  TRUE
              asym_laplace  TRUE     frmtmb fit, code 1, logLik -427.363  TRUE
                 bernoulli  TRUE     frmtmb fit, code 0, logLik -188.242  TRUE
                      Beta  TRUE     frmtmb  fit, code 0, logLik 180.057  TRUE
                  binomial  TRUE     frmtmb fit, code 0, logLik -550.421  TRUE
               categorical  TRUE     frmtmb fit, code 0, logLik -328.423  TRUE
                       cox  TRUE     frmtmb fit, code 0, logLik -306.869  TRUE
                    cratio  TRUE     frmtmb fit, code 0, logLik -410.753  TRUE
                cumulative  TRUE     frmtmb fit, code 0, logLik -405.489  TRUE
                exgaussian  TRUE     frmtmb fit, code 0, logLik -482.037  TRUE
               exponential  TRUE     frmtmb fit, code 0, logLik -446.389  TRUE
                     Gamma  TRUE     frmtmb fit, code 0, logLik -334.437  TRUE
                  gaussian  TRUE     frmtmb fit, code 0, logLik -435.436  TRUE
                 geometric  TRUE     frmtmb fit, code 0, logLik -560.396  TRUE
              hurdle_gamma  TRUE     frmtmb fit, code 0, logLik -414.477  TRUE
          hurdle_lognormal  TRUE     frmtmb fit, code 0, logLik -479.065  TRUE
        hurdle_negbinomial  TRUE     frmtmb fit, code 0, logLik -491.532  TRUE
            hurdle_poisson  TRUE     frmtmb fit, code 0, logLik -493.838  TRUE
          inverse.gaussian  TRUE     frmtmb fit, code 0, logLik -370.394  TRUE
                 lognormal  TRUE     frmtmb fit, code 0, logLik -434.542  TRUE
               multinomial  TRUE     frmtmb fit, code 0, logLik -997.152  TRUE
               negbinomial  TRUE     frmtmb fit, code 0, logLik -498.707  TRUE
                   poisson  TRUE     frmtmb fit, code 0, logLik -539.529  TRUE
         shifted_lognormal  TRUE     frmtmb fit, code 0, logLik -435.667  TRUE
               skew_normal  TRUE     frmtmb fit, code 0, logLik -430.969  TRUE
                    sratio  TRUE     frmtmb fit, code 0, logLik -388.482  TRUE
                   student  TRUE     frmtmb fit, code 0, logLik -468.864  TRUE
                 von_mises  TRUE     frmtmb fit, code 0, logLik -225.127  TRUE
                   weibull  TRUE     frmtmb fit, code 0, logLik -351.072  TRUE
        zero_inflated_beta  TRUE     frmtmb  fit, code 0, logLik -18.784  TRUE
    zero_inflated_binomial  TRUE     frmtmb fit, code 0, logLik -581.061  TRUE
 zero_inflated_negbinomial  TRUE     frmtmb fit, code 0, logLik -533.357  TRUE
     zero_inflated_poisson  TRUE     frmtmb fit, code 0, logLik -505.269  TRUE
    zero_one_inflated_beta  TRUE     frmtmb  fit, code 0, logLik -70.268  TRUE
                    wiener FALSE frmtmb.eam fit, code 0, logLik -176.759  TRUE
                 dirichlet FALSE          -                            - FALSE
          discrete_weibull FALSE          -                            - FALSE
                   frechet FALSE          -                            - FALSE
         gen_extreme_value FALSE          -                            - FALSE
           logistic_normal FALSE          -                            - FALSE

extensions probed: frmtmb.coupling 0.6.0, frmtmb.eam 0.11.3, frmtmb.latent 0.6.0, frmtmb.learn 0.7.2, frmtmb.ode 0.7.2, frmtmb.sample 0.18.0, frmtmb.spline 0.10.1 
named by core (the 0.34.0 name probe): 34 of 40 
named by core plus extensions: 35 of 40 
FITTED by core: 34 of 40 
FITTED by core plus extensions: 35 of 40 
named but not fitted:  
missing from core: wiener, dirichlet, discrete_weibull, frechet, gen_extreme_value, logistic_normal 
missing from core plus extensions: dirichlet, discrete_weibull, frechet, gen_extreme_value, logistic_normal 

```
<!-- END generated -->

0.67.0 (core only, no extension installed beside it, the name probe
only): 34 of 40, the same six missing.

**Does frmtmb.eam's `wiener()` fit brms's wiener usage?** No vignette
fits a wiener model: `brms_families.Rmd` describes it (drift `mu`,
`bs`, `ndt`, `bias`; the response `rt | dec(dec)`) and
`brms_multilevel.ltx` mentions it. So the check is brms's spelling of
that description (`dev/remeasure-wiener.R`, seed 20261013, n = 400):

<!-- BEGIN generated: remeasure-out/probe/wiener-r8.txt -->
```text
frmtmb 0.70.0 frmtmb.eam 0.11.3 from C:/Users/adf44/source/r/rellib-r8 

## brms 2.23.0 accepts the spelling (stancode, no compile)
ok; 4 lines of the Stan code name wiener 

## frmtmb.eam wiener(), the same formula
                  Estimate  Est.Error       Q2.5      Q97.5
Intercept       0.64367573 0.09410635  0.4592307 0.82812079
bs_Intercept    0.33853225 0.02496676  0.2895983 0.38746619
ndt_Intercept   2.19761188 0.20291576  1.7999043 2.59531947
bias_Intercept -0.04056514 0.07121037 -0.1801349 0.09900463
x               0.47680605 0.07482577  0.3301502 0.62346187
logLik -238.528500, convergence 0

## brms's family object handed to frm(), as a brms user would
ERROR: wiener is not a supported family. Supported families are: gaussian, poisson, binomial, Gamma, lognormal, student, negbinomial, nbinom2, nbinom1, beta, Beta, tweedie, compois, zero_inflated_poisson, zero_inflated_negbinomial, hurdle_poisson, hurdle_negbinomial, multinomial, cumulative, hurdle_cumulative, beta_binomial, zero_inflated_beta_binomial, skew_normal, inverse.gaussian, exgaussian, bernoulli, geometric, exponential, weibull, shifted_lognormal, hurdle_gamma, hurdle_lognormal, zero_inflated_binomial, zero_inflated_beta, zero_one_inflated_beta, xbeta, asym_laplace, zero_inflated_asym_laplace, huber, sratio, cratio, acat, von_mises, cox, categorical 

## independent ML through rtdists::ddiffusion and optim
rtdists vs RWiener, two rows: max rel diff 2.2e-16
Intercept         x        bs       ndt      bias 
0.6436759 0.4768064 1.4028870 0.3108021 0.4898598 
logLik -238.528501, convergence 0
logLik difference frmtmb - rtdists ML: 4.86e-08
drift intercept 0.643676 vs 0.643676, slope 0.476806 vs 0.476806

## the extension probe's absent case: an eam export that is not a family
ddm_simulate exported: TRUE ; accepted as a family: FALSE 

## the scale of the ndt and bias coefficients against brms's links
brms links: bs log  ndt log  bias logit 
exp(bs_Intercept) 1.402887 vs reference bs 1.402887
exp(ndt_Intercept) 9.003486 vs reference ndt 0.310802
plogis(ndt_Intercept) * min(rt) 0.310802 (min rt 0.345322)
plogis(bias_Intercept) 0.489860 vs reference bias 0.489860
```
<!-- END generated -->

- The spelling fits: `rt | dec(dec) ~ x, bs ~ 1, ndt ~ 1, bias ~ 1`
  under `frmtmb.eam::wiener()` reaches the maximum of an independent
  likelihood built on `rtdists::ddiffusion()`, log-likelihood within
  4.9e-08, drift coefficients equal to six digits.
- **`ndt_Intercept` is not on brms's scale.** brms's `wiener()` puts
  `ndt` on a log link; frmtmb.eam puts it on a logit of `ndt / max_ndt`
  (`max_ndt` the fastest response by default): 2.198 here, which is
  `log(0.3108)` = -1.169 on brms's link. The same coefficient name
  carries a different parameter, as the `sx_1` defect did. For lane
  famuni (wiener parity).
- **brms's family object is refused without a pointer:** `frm(...,
  family = brms::wiener())` stops with "wiener is not a supported
  family" and does not name frmtmb.eam.
- RWiener's `dwiener()` reads `delta[1]` for every row (it is not
  vectorized over the drift), so a reference built on it is wrong for a
  row-varying drift; the probe uses rtdists and shows the two agree on
  scalar arguments.

## Task B: bin 2 of brms's suite

`dev/brms-suite-audit.md` bins brms 2.23.0's suite; bin 2 is 329
assertions, 59 "2T" (a posterior, restated) and 270 "2P" (a brms
internal object). Until now bin 2 was never measured, so "389 of 823"
counted all of it as failing.

### How it is ported

The bin-1 machinery, unchanged in its rules: brms's text, wrapped in
`brms_port()` with a verdict, run by the package's own harness
(`tests/testthat/helper-brms-suite.R`), so a pass is an assertion
observed to hold and never on a missing function, a stale object or a
hollow `$` read. New files only:

| file | role |
|---|---|
| `dev/brmsport-bin2-blocks.R` | finds the 2T and 2P blocks; 59 in 14 and 270 in 68, checked |
| `dev/brmsport-bin2-gen.R` | writes `dev/brmsport-bin2/<package>/test-brms-suite-bin2-*.R` |
| `dev/brmsport-bin2-run.R`, `dev/brmsport-bin2-record.sh` | stage a file beside its package's helpers and run it, recording or asserting |
| `dev/brmsport-bin2-verdicts-2t.R` | the 2T verdicts and targets, as reviewable R |
| `dev/brmsport-bin2-survey.R` | the 2P survey: groups, adapters, classes |
| `dev/brmsport-bin2-adapters.R`, `-adapters-spec.R` | the thin 2P adapter and the blocks it ports |
| `dev/brmsport-bin2-ledger.R` | the strict ledger, and the totals composed with bin 1 |
| `dev/brmsport-bin2-all.sh` | the whole bin-2 pipeline in order |

The generated tests live under `dev/`, not in `tests/`: bin 2 is a
measurement, and `dev/brmsport-bin2-run.R` copies a file with its
package's `helper-*.R` into `dev/brmsport-bin2-run/` to run it. No
bin-1 file changed. The bin-1 regeneration has no script to hook into
(it is run by hand at consolidation), so the hook is one more command
after it: `sh dev/brmsport-bin2-all.sh`.

**2T, the target rule.** A 2T assertion reads a posterior, and
frmtmb.sample's draws are frmtmb's posterior. So each of the 58 methods
rows runs twice, brms's text unchanged: on `frm_sample()` draws of the
fixture (the TARGET, which decides the row) and on the
maximum-likelihood fixture fit (the restatement for ML, recorded and
not counted). `priors:130` is core-tier and is judged on frmtmb. The
draws side binds all six fixtures (`brms_fixture_draws(1:6)`, 25 draws
each as brms's examples have; bin 1 binds four because no bin-1 block
reads fit4 or fit6 on draws).

### The totals

<!-- BEGIN generated: brmsport-bin2-summary.md -->
<!-- generated by dev/brmsport-bin2-ledger.R; do not edit by hand -->

### Bin 2 per file and bin

| file (bin) | assertions | pass (value / shape) | defect | divergence | cannot transfer: decision / gap / cost / artifact |
|---|---|---|---|---|---|
| `tests.brmsfit-helpers.R` (2P) | 8 | 0 (0 / 0) | 0 | 0 | 8 / 0 / 0 / 0 |
| `tests.brmsfit-methods.R` (2T) | 58 | 23 (1 / 22) | 9 | 7 | 0 / 17 / 0 / 2 |
| `tests.brmsterms.R` (2P) | 18 | 0 (0 / 0) | 0 | 0 | 18 / 0 / 0 / 0 |
| `tests.distributions.R` (2P) | 113 | 43 (12 / 31) | 0 | 2 | 40 / 24 / 4 / 0 |
| `tests.log_lik.R` (2P) | 65 | 38 (16 / 22) | 0 | 1 | 0 / 10 / 14 / 2 |
| `tests.posterior_predict.R` (2P) | 66 | 41 (0 / 41) | 0 | 0 | 0 / 12 / 13 / 0 |
| `tests.priors.R` (2T) | 1 | 0 (0 / 0) | 0 | 1 | 0 / 0 / 0 / 0 |
| **total** | 329 | 145 (29 / 116) | 9 | 11 | 66 / 63 / 31 / 4 |

### Bin 2 outcome by class

| bin | outcome | class | kind | assertions |
|---|---|---|---|---|
| 2T | pass | - | - | 23 |
| 2T | defect | argument | - | 7 |
| 2T | defect | hollow | - | 1 |
| 2T | defect | method | - | 1 |
| 2T | divergence | policy | - | 8 |
| 2T | cannot transfer | absent | gap | 17 |
| 2T | cannot transfer | fixture | artifact | 1 |
| 2T | cannot transfer | upstream | artifact | 1 |
| 2P | pass | - | - | 108 |
| 2P | pass | per-line | - | 14 |
| 2P | divergence | domain | - | 2 |
| 2P | divergence | parameterization | - | 1 |
| 2P | cannot transfer | absent | gap | 31 |
| 2P | cannot transfer | adapter | artifact | 1 |
| 2P | cannot transfer | brms-backends | decision | 3 |
| 2P | cannot transfer | brms-internal | decision | 48 |
| 2P | cannot transfer | no-cdf | gap | 15 |
| 2P | cannot transfer | no-density | decision | 12 |
| 2P | cannot transfer | no-quantile | decision | 3 |
| 2P | cannot transfer | not-ported | cost | 31 |
| 2P | cannot transfer | reference | artifact | 1 |

Bin 2T passes: 23 of 59 (39.0%), judged on the target object.
Bin 2P: 178 of 270 assertions run through the adapter, 122 of them pass (the rest it refused or they failed); the other 92 are classified from the survey.
Of the 145 bin-2 passes, 29 compare a value and 116 check only a shape, length, names or type.
The 58 methods rows of 2T on the maximum-likelihood fit (the restatement, not counted): divergence 52, pass 6.

### Bins 1 and 2 together

| outcome | bin 1 | bin 2T | bin 2P | bins 1 and 2 |
|---|---|---|---|---|
| pass | 389 | 23 | 122 | 534 |
| defect | 2 | 9 | 0 | 11 |
| divergence | 52 | 8 | 3 | 63 |
| cannot transfer: no counterpart (decision) | 48 | 0 | 66 | 114 |
| cannot transfer: absent (gap) | 2 | 17 | 46 | 65 |
| cannot transfer: not ported (cost) | 0 | 0 | 31 | 31 |
| cannot transfer: port artifact | 1 | 2 | 2 | 5 |
| **total** | **494** | **59** | **270** | **823** |

Bins 1 and 2 pass: 534 of 823 (64.9%); bin 1 alone 389 of 494, as dev/brmsport-ledger.tsv records it (regenerated once at consolidation).
Over every assertion but the 114 that frmtmb has no counterpart for by decision, 534 of 709 pass (75.3%). The 65 gaps, 31 rows not ported (cost) and 5 port artifacts stay in that denominator.
The rate's denominator excludes exactly the "cannot transfer" rows of these classes: brms-backends 3, brms-internal 63, mcmc 4, no-density 12, no-draws 3, no-quantile 3, stan 26.
<!-- END generated -->

### The 2T rows that do not pass, on draws

| row | verdict | evidence |
|---|---|---|
| `brmsfit-methods:121` | defect, argument | `bayes_R2(draws, newdata =, re_formula = NA)` refused by name; brms passes both to `posterior_epred()` |
| `brmsfit-methods:142` | cannot transfer, fixture | `coef(fit4 draws)$subject` is 8 x 4 x 7: brms's 8 less `disc_Intercept`, which fixture 4 drops |
| `brmsfit-methods:434` | defect, method | `logLik()` of draws: no method; brms returns the 25 x 40 log-likelihood matrix |
| `brmsfit-methods:469` | defect, hollow | holds only because `fit6$data` is NULL on draws, and `loo()` on draws swallows `newdata`: `loo(d6, newdata = 20 rows)` returns 40 pointwise rows (R2) |
| `brmsfit-methods:471` | cannot transfer, upstream | loo 2.10.1's `loo_compare()` puts the model name in column 1, so `[2, 1]` is a string; brms's own fixture 6 gives the same string |
| `brmsfit-methods:442`, `:443`, `:1005`, `:1021` | divergence, policy | `LOO()`, `WAIC()`, brms's deprecated spellings, refused by name |
| `brmsfit-methods:446` to `:449`, `:1015`, `:1016` | defect, argument | `loo(fit1, fit1)`, `waic(fit1, fit1)`: one model per call, refused by name, naming `loo_compare()`; open work in `dev/brms-api-diff.md` (c), owner postfix (fix round 1; they read "divergence, policy" before) |
| `brmsfit-methods:478` to `:480` | divergence, policy | `loo_subsample()` refused by name |
| 17 rows of `loo_R2`, `loo_epred`, `loo_linpred`, `loo_predict`, `loo_predictive_interval`, `loo_model_weights` | cannot transfer, absent | no such function in frmtmb or frmtmb.sample; a gap, not a decision |
| `priors:130` | divergence, policy | `default_prior()` describes `frm()`, flat in every slot (bin-1 `priors:85`) |

Constructions: `dev/remeasure-bin2-probe1.R`, log
`dev/remeasure-out/probe/bin2-probe1.txt`. On the ML fit (not counted),
6 rows hold as brms wrote them (the `VarCorr()` shapes) and 52 do not:
posterior quantities refused by name, and `coef()` of a fit, which is
levels x coefficients (lme4's shape, the policy of bin-1
`brmsfit-methods:811`), nested by nonlinear parameter on a nonlinear
fit, and without brms's `sigma_` columns.

### The 2P survey

Groups by the brms object the assertion reads. The names each block
calls are read from its parsed source by `dev/brmsport-bin2-survey.R`,
not typed:

<!-- BEGIN generated: remeasure-out/blocks/bin2-groups.md -->
| group | blocks | assertions | ported | reads | frmtmb counterpart | cost | class if unported |
|---|---|---|---|---|---|---|---|
| G1 | 13 | 67 | 52 | brms's exported d/p/q/r functions of its own families (dstudent_t, dskew_normal, dvon_mises, dbeta_binomial, dzero_inflated_*, dhurdle_*, dwiener, ...) and the multivariate normal and Student t | every exported family constructor's $lpdf(y, dpars, aterms) (a log density) and $sim(dpars, aterms, n); $lcdf where the family has one. No quantile function, and no multivariate density outside a fit | low for the d and r functions of the families frmtmb has: a shim per brms signature over row_lpdf() and sim_draw() (ported, fix round 2). p needs a CDF, which only inverse.gaussian has; q has no counterpart | cannot transfer: no-density |
| G2 | 13 | 46 | 20 | brms's unexported ordinal and categorical densities and link maps on arrays (brms:::dcumulative, dsratio, link_cumulative, inv_link_sratio, dcategorical, link_categorical, ...), called unqualified inside brms's own test environment | the ordinal densities: each ordinal family object's exported ord_lpdf_make(identity), a thin read (ported). The link maps on arrays and the categorical functions: none, frmtmb exports no map between linear predictors and category probabilities | low for the two density blocks; high for the link maps (the array convention, draws x observations x categories, would be new frmtmb functions) | cannot transfer: brms-internal |
| G3 | 18 | 65 | 52 | brms:::log_lik_<family>(i, prep), with prep a hand-built brmsprep (dpars as draws x observations matrices, data$Y, trials, weights, cens) | every exported family constructor's $lpdf(y, dpars, aterms), documented in ?frmtmb_family: a thin adapter evaluates it per draw at observation i; an ordinal family's ord_lpdf_make(identity) takes brms's thresholds. Autocorrelation (ARMA, SAR, FCOR), multivariate, categorical and custom-family preps have no such read | low for the univariate families (one adapter, a dpar name map and brms's reference densities bound as data); high for the rest | cannot transfer: not-ported |
| G4 | 17 | 66 | 54 | brms:::posterior_predict_<family>(i, prep), the same brmsprep | every exported family constructor's $sim(dpars, aterms, n), documented in ?frmtmb_family: a thin adapter draws per draw at observation i; an ordinal family's ord_sim_make(identity). The same exclusions as G3 | low for the univariate families; high for the rest | cannot transfer: not-ported |
| G5 | 5 | 18 | 0 | brms:::brmsterms(), check_re_formula(), update_re_terms(): brms's parsed-formula object and its formula editors | none exported: frm(dry_run = 'frame') returns frmtmb's assembled frame, whose structure is not brmsterms's (no allvars, no dpars list in brms's shape); re_form_arg() resolves a re_formula against a fit, not against a formula | high: an adapter would re-derive brms's object from frmtmb's frame | cannot transfer: brms-internal |
| G6 | 1 | 5 | 0 | brms:::get_cor_matrix_ar1(), _ma1(), _arma1(), _cosy(), _ident(): autocorrelation matrices from parameter draws | autocor_matrix(fit) is exported but reads a FITTED model's correlation, and returns a correlation matrix where brms's AR and ARMA matrices are the covariance at unit innovation variance; there is no function of the parameter | medium: an adapter would fit a model at the given parameter or call internals | cannot transfer: brms-internal |
| G7 | 1 | 3 | 0 | brms:::insert_refcat(): the reference category of a categorical linear predictor | none | high: no counterpart | cannot transfer: brms-internal |
| **total** | **68** | **270** | **178** | | | | |
<!-- END generated -->

An adapter was built only where it is a thin read of exported frmtmb
objects and unlocks a group: G3, G4 and the two ordinal density blocks
of G2 (`dev/brmsport-bin2-adapters.R`). `brms_adapt_log_lik(fun, i, prep)` and
`brms_adapt_posterior_predict(fun, i, prep)` take brms's family name to
`frmtmb::brmsfamily(fun)`, read column `i` of brms's prep, and hand it
to `frmtmb::row_lpdf()` (the row density the objective and `log_lik()`
run, `cens()` and `trunc()` folded in) or `frmtmb::sim_draw()` (every
simulator's draw, with `trunc()` rejection). The density and the draw
are frmtmb's; the adapter only moves numbers, and refuses any prep field
it does not map (autocorrelation, thresholds, multivariate) rather than
dropping it. The reference side stays brms's exported density or base
R's, bound in the generated header. Setup substitutions: brms's
internal `inv_logit` and `inv_cloglog` by copies of their one-line
bodies (bit-identical data), `brms:::rasym_laplace` by the exported
`brms::rasym_laplace`, and `skip_if_not_installed("extraDistr")` and
`("betareg")`, which only brms's own family code needs, by a no-op.
An ordinal family (cumulative, sratio, cratio, acat,
hurdle_cumulative) takes brms's thresholds through its family object's
exported `ord_lpdf_make(identity)` and `ord_sim_make(identity)`, and
brms's link through `brmsfamily(fun, link =)`; this is the review's
construction (`dev/remeasure-rev-ord.R`, seed 446, against
`brms:::log_lik_cumulative` to 6.7e-13). brms's `dcumulative()`,
`dsratio()`, `dcratio()` and `dacat()` go through `brms_adapt_dord()`
on the same objects; a category-specific (matrix) `eta` is passed as
`mu = 0` with thresholds `thres - eta`, the same numbers, since all four
families read `eta` only through `thres - eta`.
brms's exported d and r functions (G1, fix round 2) go through
`brms_g1_shims`: each keeps brms's signature and maps its arguments onto
the family's dpars (brms's `shape1, shape2` onto `mu, phi`, `shape,
scale` onto `mu, shape`, `size` onto `trials()`), then calls
`row_lpdf()` or `sim_draw()`. `dev/brmsport-bin2-adapters-spec.R` lists
the 37 blocks it reaches (11 in each of G3 and G4, 15 in G1 and G2,
the G2 `inv_link_<ordinal>()` blocks added at fix 2b); left to the
survey: the multivariate, ARMA, SAR, FCOR, wiener and custom blocks of
G3 and G4, G1's multivariate and wiener blocks, G2's link maps and
categorical maps, and G5 to G7. A family frmtmb lacks
inside a listed block (discrete_weibull, negbinomial2, frechet,
gen_extreme_value) is refused by `brmsfamily()` in frmtmb's words, so
those rows are cannot transfer, absent, and the binomial, poisson,
gamma, weibull and other comparisons beside them are kept.

The ported 2P rows that do not pass, with their verdicts, are written
by `dev/brmsport-bin2-verdicts-2p.R`; their counts per class are in the
generated "outcome by class" table above. Fix round 2 added: 15 p rows
(`no-cdf`, gap), 3 q rows (`no-quantile`, decision), 9 frechet and
gen_extreme_value rows (absent), the categorical block's rows
(categorical not ported, cost; multinomial, dirichlet,
dirichlet_multinomial, dirichlet2 and logistic_normal absent), and the
per-line split of `distributions:238`. The first rounds' rows, with
their round-2 classes:

| row | verdict | evidence |
|---|---|---|
| `log_lik:271`, `posterior_predict:159` | cannot transfer, absent | discrete_weibull |
| `log_lik:288`, `posterior_predict:169` | cannot transfer, absent | negbinomial2 (brms's `sigma = 1 / shape` spelling; frmtmb's negbinomial row `log_lik:285` holds) |
| `log_lik:324`, `posterior_predict:184` | cannot transfer, absent | frechet |
| `log_lik:349`, `posterior_predict:190` | cannot transfer, absent | gen_extreme_value |
| `log_lik:267` | cannot transfer, reference | brms's reference `dbeta_binomial()` needs extraDistr, not installed; frmtmb's side agrees with the closed form to 1.7e-14 |
| `log_lik:302` | divergence, parameterization | the reference is `brms:::dcom_poisson()` (recorded, not run), whose `mu` is `lambda^(1/shape)`; frmtmb's is the mean (R6) |
| `log_lik:392` | cannot transfer, absent | censored `von_mises`: frmtmb's has no CDF; the uncensored row holds |
| `log_lik:32` | cannot transfer, adapter | holds, but asserts brms's `log_lik_weight()`, `ll * weights[i]`, which the adapter computes itself: frmtmb weights a row outside the family (in the objective and in `log_lik()` of a fit). Not frmtmb's answer, so not a pass; a port artifact, kept in the rate's denominator |
| `distributions:238` (8 assertions, one loop) | 7 pass, 1 divergence (domain) | per assertion line (`dev/brmsport-bin2-split.tsv`): the sratio, cratio and acat identities and the four `dim()` lines hold in 32 of 32 iterations; the cumulative identity holds in 16 of 32. brms's block draws unsorted thresholds, and at a decreasing pair brms's `dcumulative()` returns a NEGATIVE probability (4 of 20 cells at seed 238), which its closed-form helper reproduces; frmtmb's cumulative returns NaN with a warning in categories on either side of the pair (at draw 5 of `bin2-probe3.txt`, categories 1 and 2, where brms's category 1 is 0.805). On sorted thresholds every cell agrees. A fit never reaches decreasing thresholds |

**What the passes assert.** The bin-2 summary's pass column is split
into "value" (the assertion compares a computed number with a reference
or another computation) and "shape" (a shape, length, names or type
only), read from brms's assertion text by `dev/brmsport-bin2-ledger.R`.
Of 136 bin-2 passes, 24 are value checks and 112 shape checks. The 23
2T passes are frmtmb.sample's, all shape checks on draws except
`brmsfit-methods:1020` (`waic2` equals `waic_pointwise`). In 2P, the
value passes are 16 `log_lik` rows that compare frmtmb's density with
brms's exported density or base R's, `distributions:446` (2,
`dsratio()` against `dcratio()` and the cloglog inequality), the three
ordinal identities of `distributions:238`, and the student integral
and `dt()` identity (`distributions:4`, `:5`); the other G1 passes are
`length()` checks, as brms wrote them; every
`posterior_predict` pass is a length or `dim()` check, and the adapter
sets the draw count itself, so those 41 say that frmtmb's simulator
answered, not what it drew. This is how brms wrote these blocks ("runs
without errors"). The review's count (73 shape, 16 2P value, before the
ordinal port) is `dev/remeasure-rev-2tpass.R`.

The passes are genuine by the harness's rules (no VACUOUS or STALE
row holds), and the reference differences behind the value passes are
0 or 4.4e-16 to 1.4e-14 absolute per family
(`dev/remeasure-bin2-probe2.R`, log
`dev/remeasure-out/probe/bin2-probe2.txt`); 0 of 450 truncated draws
fall outside their bound.

**The instrument can say no, on frmtmb's density.** The weights
construction of the first record (`dev/remeasure-bin2-absent.R`) only
showed that `log_lik:32` tests the adapter. Redone on frmtmb:
`dev/remeasure-bin2-sayno.R` stages the adapter with frmtmb's row
density, `row_lpdf()` and each ordinal `ord_lpdf_make()` row, multiplied
by `1 + 1e-6`, and records the log_lik file again. Every one of the 16
value passes stops holding and every one of the 22 shape passes still
holds, which is both the evidence that the value rows test frmtmb's
density at a relative 1e-6 and the measure of how little the shape rows
can see (`dev/remeasure-out/probe/bin2-sayno.txt`):

<!-- BEGIN generated: remeasure-out/probe/bin2-sayno.txt from ^log_lik passes to ^value rows still -->
```text
log_lik passes on the true build: 38 
       holds_with_wrong_density
kind    FALSE TRUE
  shape     0   22
  value    16    0

```
<!-- END generated -->

A pass recorded as a divergence turns the tier red, `pass=57 fail=1`,
then green again on the true verdict
(`dev/remeasure-out/probe/bin2-tier-guard.txt`).

## The whole re-measure as one command

    bash dev/remeasure-run.sh <libs> <tag> [jobs] [stages]

`<libs>` goes before the user library (`PORT_LIB`, `BV_LIB`); `<tag>`
names `dev/remeasure-out/<tag>` and the build's column in the outcome
tables. It runs the provenance check (a fit against the two reference
log-likelihoods of `dev/correct-libcheck.R`), the four mechanical passes
and the twelve hand-translation scripts (8 R processes at a time), the
brms reference fits if `dev/remeasure-out/brms-ref` is missing, then the
reducers: `summarize.R`, `_drift.R`, `families-coverage.R`,
`plausibility.R`, `remeasure-history.R`, `remeasure-gaps.R`,
`remeasure-plauscmp.R`, `remeasure-secmp.R`. Then
`Rscript dev/remeasure-blocks.R v067 <tag>` and
`Rscript dev/remeasure-fill.R` refresh this file. For the merged build
of this round:

    bash dev/remeasure-run.sh <merged lib> merged
    Rscript dev/remeasure-blocks.R r8 merged

Wall time on this machine, with six other lanes running
(`dev/remeasure-out/r8/wall.txt`):

<!-- BEGIN generated: remeasure-out/r8/wall.txt -->
```text
## 2026-10-13 15:46:16 stages=ports jobs=8 lib=C:/Users/adf44/source/r/rellib-r8
provenance 1 s
ports (4 mechanical passes and the hand translation) 86 s
total 89 s
## 2026-10-13 15:52:10 stages=brms jobs=8 lib=C:/Users/adf44/source/r/rellib-r8
provenance 2 s
brms reference fits 1418 s
total 1420 s
## 2026-10-13 16:16:42 stages=all jobs=8 lib=C:/Users/adf44/source/r/rellib-r8
provenance 1 s
ports (4 mechanical passes and the hand translation) 125 s
reducers 8 s
total 135 s
## 2026-10-13 16:23:16 stages=reduce jobs=8 lib=C:/Users/adf44/source/r/rellib-r8
provenance 3 s
reducers 8 s
total 11 s
## 2026-10-13 16:24:40 stages=reduce jobs=8 lib=C:/Users/adf44/source/r/rellib-r8
provenance 2 s
reducers 9 s
total 12 s
## 2026-10-13 17:14:32 stages=reduce jobs=8 lib=C:/Users/adf44/source/r/rellib-r8
provenance 3 s
reducers 14 s
total 18 s
```
<!-- END generated -->

The 135 s run is the command as the consolidation will run it (the brms
reference reused); the brms reference fits took 1418 s once, 3 at a
time. Bin 2 is a separate command, `sh dev/brmsport-bin2-all.sh`.

Provenance of the judged build:

<!-- BEGIN generated: remeasure-out/r8/provenance.txt -->
```text
libPaths: C:/Users/adf44/source/r/rellib-r8 ; C:/Users/adf44/AppData/Local/R/win-library/4.6 ; C:/Program Files/R/R-4.6.1/library 
C:/Users/adf44/source/r/rellib-r8: 9 package dirs, 0 hollow 
frmtmb           0.70.0   C:/Users/adf44/source/r/rellib-r8
frmtmb.sample    0.18.0   C:/Users/adf44/source/r/rellib-r8
frmtmb.eam       0.11.3   C:/Users/adf44/source/r/rellib-r8
frmtmb.latent    0.6.0    C:/Users/adf44/source/r/rellib-r8
frmtmb.spline    0.10.1   C:/Users/adf44/source/r/rellib-r8
frmtmb.coupling  0.6.0    C:/Users/adf44/source/r/rellib-r8
frmtmb.learn     0.7.2    C:/Users/adf44/source/r/rellib-r8
frmtmb.ode       0.7.2    C:/Users/adf44/source/r/rellib-r8
brms             2.23.0   C:/Users/adf44/AppData/Local/R/win-library/4.6
rstan            2.32.7   C:/Users/adf44/AppData/Local/R/win-library/4.6
tmbstan          1.2.1    C:/Users/adf44/AppData/Local/R/win-library/4.6
RTMB             2.0      C:/Users/adf44/AppData/Local/R/win-library/4.6
fit check: gaussian (1 | g) logLik -295.602189818, poisson -332.137876329
```
<!-- END generated -->

## Defects found, not fixed

**R1. 0.70.0 loses the standard errors of two `mi()` fits at an
unchanged optimum.** `brms_missings` `fit_imp2`,
`bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) +
set_rescor(FALSE)` on `mice::nhanes`: the log-likelihood is
-113.672565129 on every build, but 0.70.0 reports false convergence
(8), `diagnose()` says "Hessian positive definite: FALSE", and every
standard error is NaN. `fit_imp3` (the same with `chl | mi(se)`) keeps
code 0 but withholds both residual-SD standard errors as "flat" and
reports the others up to 35% too large (6.10 against 4.52 for
`bmi_age`). The finite-difference Hessian of the same objective at the
same point (`stats::optimHess()` on `fit$obj$fn`, frmtmb's SE code not
involved) is positive definite, eigenvalues 0.0016 to 5.45e5, and its
standard errors equal 0.67.0's to five digits. Builds 0.67.0, 0.68.1
and 0.69.1 report those standard errors.

Which commit (corrected at fix round 1): d2ef6312, the 0.70.0 release,
is an octopus merge of 2884ccbd, 58b97a2c, ea711734 and ebd12b88, and
each of those has the one parent 48c1e6c9 (0.69.0). So "58b97a2c is
right, ea711734 is wrong", as the first record put it, is not a
bisection. What holds: the common parent 48c1e6c9 is clean; ea711734
("scale-invariant fits and Newton steps from the stop") ALONE
reproduces both failures; the other three parents, 2884ccbd, 58b97a2c
and ebd12b88, are clean (my runs: ea711734 and 58b97a2c,
`dev/remeasure-out/probe/repro-missings-{ea711734,58b97a2c}.txt`; the
review's: all five, `dev/remeasure-rev-repro-<build>.txt`).
`frmtmb_control(autoscale = FALSE)` does not restore them.

The mechanism, found by the review (`dev/remeasure-rev-r1units.R`,
`dev/remeasure-rev-r1cause.R`, logs `.txt`): at 0.70.0 `fit_imp2`
carries `par_units` 8 8 8 8 128 128 1 1, the response units that
ea711734's `predictor_par_units()` gives every coefficient of a
response whatever its column's scale, under `autoscale = FALSE` as
well (0.69.1 carries none, which is why that switch does not help).
`autoscale_sdreport()` then takes `optimHess()` in `q = par / u` at the
default `ndeps = 1e-3`: a step of 0.008 in `bmi_michl:age` (SE 0.019,
column `age * chl` up to several hundred) and of 0.128 in the `chl`
coefficients. In those units the finite-difference Hessian is
indefinite from the step, not from the objective: eigenvalues -0.0115
to 3.48e7 at `ndeps` 1e-3, 1.026 to 3.49e7 and the right standard
errors at 1e-4, and in raw space at 1e-3 0.0016 to 5.45e5 with the
right standard errors. `autoscale_fit_sound()` takes the same Hessian.
Any fit with a response unit above 1 and a column of large spread can
meet this, so the reach is unknown and wants a scan.

Repro `dev/remeasure-repro-missings.R <lib>`; logs
`dev/remeasure-out/probe/repro-missings-{v067,r7,r8}.txt`. The vignette
call `summary(fit_imp2)` still "runs", so the port's outcome tables
cannot see it; `dev/remeasure-secmp.R` now does. **Owner: lane
optstart** (assigned at fix round 1, with the cause above).

<!-- BEGIN generated: remeasure-out/probe/repro-missings-r8.txt from ^## fit_imp2: to ^diagnose -->
```text
## fit_imp2: bmi | mi() ~ age * mi(chl), chl | mi() ~ age 
logLik -113.672565129  convergence 1  message false convergence (8)
warnings: Optimizer did not report convergence: false convergence (8) 
vcov() standard errors:
bmi_Intercept chl_Intercept       bmi_age       chl_age     bmi_michl 
          NaN           NaN           NaN           NaN           NaN 
bmi_michl:age 
          NaN 
max |gradient| at the optimum 4.51e-05
finite-difference Hessian: 8 parameters, eigenvalues 0.001617 to 545000
     beta      beta      beta      beta      beta      beta     betad     betad 
 6.772400  4.519000  0.034198  0.019395 22.324000 11.765000  0.183540  0.175030 
```
<!-- END generated -->

**R2. `loo()` on draws swallows `newdata`.** `loo.frmtmb_draws(x,
ndraws, resp, ...)`: `loo(d6, newdata = <20 of 40 rows>)` returns 40
pointwise rows and the same elpd as `loo(d6)`, with no message; brms
computes the elpd on `newdata`. `log_lik()` on draws refuses `newdata`
by name, so the two disagree. A wrong answer with no message, the most
severe of R2 to R4. Repro `dev/remeasure-bin2-probe1.R`. Owner: lane
postfix (assigned at fix round 1), as for R3 and R4.

**R3. `logLik()` of draws has no method** ("no applicable method"),
where brms's `logLik()` of a fit is its draws x observations
log-likelihood. Ported row `brmsfit-methods:434`. Owner: postfix.

**R4. `bayes_R2()` on draws refuses `newdata` and `re_formula`.** Row
`brmsfit-methods:121`. Owner: postfix.

**R5. frmtmb.eam's `ndt_Intercept` is on a different link than brms's
`ndt_Intercept`** (section "Families"), and `frm(family =
brms::wiener())` does not point to frmtmb.eam.

**R6. `com_poisson`: three meanings of `mu`, and brms contradicts
itself** (restated at fix round 1; the first record read brms's R-side
helper as brms's model, which it is not). frmtmb's
`brmsfamily("com_poisson")` maps onto `compois`, whose `mu` is the MEAN
(`RTMB::dcompois2`). brms 2.23.0's FITTED `mu` is the COM-Poisson rate
lambda, and its R-side helpers take `mu = lambda^(1/shape)`; the
upstream entry below gives the evidence. So a brms `com_poisson` model
and frmtmb's do not estimate the same `mu`: on the review's data
(lambda 9, nu 2, n = 400, `dev/remeasure-rev-compois.R`) frmtmb's
`exp(Intercept)` is 2.7825, the sample mean, and `exp()` of brms's
posterior mean of the Intercept, `exp(1.988)`, is 7.30, the rate at
nu = 1.80. The coefficients are
not comparable (Intercept 1.02 against 1.99). frmtmb needs a
"differences from brms" note for `com_poisson`; no lane of this round
owns it.

Upstream, in the style of `dev/upstream-bugs.md` (that file is not
edited here; the entry is for whoever files it):

> ### brms-24. `com_poisson`: the R-side pmf is not the pmf the Stan
> program fits
>
> - Seen: 2.23.0. Confirmed 2026-10-13 (source read and a fit).
> - Repro: brms's Stan function `com_poisson_log_lpmf(y | log_mu, nu)`
>   (`chunks/fun_com_poisson.stan`) is `y * log_mu - nu * lgamma(y + 1)
>   - log_Z_com_poisson(log_mu, nu)`, with series terms
>   `(k - 1) * log_mu - nu * lgamma(k)`: the pmf is proportional to
>   `mu^y / (y!)^nu`, so the fitted `mu` is Shmueli's rate lambda.
>   brms's R side, `brms:::dcom_poisson(x, mu, shape)` (used by
>   `log_lik()`, `posterior_predict()` and `loo()`), is
>   `shape * (x * log(mu) - lgamma(x + 1)) - log_Z`, with
>   `brms:::log_Z_com_poisson()` summing `shape * (k * log_mu - lfac)`:
>   proportional to `(mu^x / x!)^shape`, so there `mu = lambda^(1/shape)`.
>   At `mu = 3, shape = 2` the two pmfs have different means: 1.453549
>   on the Stan side, 2.737078 on the R side (exact sums over 0:200,
>   `dev/remeasure-compois-means.R`; the block after this entry).
> - Expected: one pmf. Actual: post-processing of a `com_poisson` fit
>   evaluates and simulates a different pmf from the one sampled, for
>   every `shape != 1`; `shape = 1` agrees (Poisson).
> - Workaround: none in frmtmb; frmtmb does not call either. Its own
>   `mu` is the mean, a third parameterization.
> - Record: `dev/remeasure-findings.md` R6; the review
>   `dev/reviews/2026-10-13-remeasure.md` section C,
>   `dev/remeasure-rev-compois.R` (brms fit: 2 chains of 1000,
>   seed 20261013).
> - Status: not reported.

<!-- BEGIN generated: remeasure-out/probe/compois-means.txt -->
```text
Stan side (rate 3, nu 2): mean 1.453549
R side (mu 3, shape 2), normalized exactly: mean 2.737078
brms:::dcom_poisson(mu 3, shape 2): mass 1.000141, mean 2.737078
```
<!-- END generated -->

The brms R-side normalizer is also approximate in its asymptotic
branch (mass 1.00014 at `mu = 3, shape = 2`), which is the helper's
property and not the fitted model's.

**R7. The ported-suite harness lets staleness lapse through a setup
line.** `brms_setup()` marks the names a line assigns fresh when it
succeeds, even when the line READ a stale name: after a failed
`a <- stop()`, `b <- a + 1` succeeds on the old `a`, and
`expect_equal(b, 2)` holds, where brms's block would have stopped
(`dev/remeasure-bin2-stale-setup.R`, log `bin2-stale-setup.txt`).
Reach, measured statically over all 177 blocks of bins 1 and 2 and
their 1372 recorded rows (`dev/remeasure-stale-reach.R`, which finds
its own constructed case and not the repaired one): 0 held assertions
read a name tainted that way. So no ledger row rests on it today; the
gap is in `tests/testthat/helper-brms-suite.R`, which this lane does not
change. The bin-1 records read are those in `dev/brmsport-log/`, from
the last regeneration.

**R8 (open, low). A frmtmb family's `$lcdf` slot returns the CDF, not
its log.** By documented convention (`R/families.R:113`: "`lcdf`
Optional vectorized AD log-safe CDF `(q, dpars, aterms)` returning
probabilities"), so the slot name is misleading:
`brmsfamily("inverse.gaussian")$lcdf(c(1.2, 3), list(mu = 1, shape =
3), list())` returns 0.7257933 and 0.990027, which are
`brms::pinv_gaussian()`'s probabilities, not their logs. Its only
caller, `fam_lcdf()` (`R/families.R:7869`), uses it correctly, and
`row_lpdf()` of a censored row gives the right log CDF and log survivor
values (`dev/remeasure-lcdf.R`, log
`dev/remeasure-out/probe/lcdf.txt`). Found while shimming
`pinv_gaussian()`; no wrong answer follows from it. No owner.

## What I did not do

- No package code changed, and no defect was fixed.
- The hand translation's labels stay at their 0.42.0 values; the drift
  table carries the change.
- The 0.65.0 and 0.66.0 columns were rebuilt from the vigport record,
  not rerun (no library of those builds exists). The rebuild matches the
  record's counts and its 0.67.0 column equals the fresh 0.67.0 run row
  for row.
- For R1 I ran two merge parents' builds (ea711734, 58b97a2c),
  installed one at a time into the lane library and replaced by 0.67.0
  again afterwards. The mechanism inside ea711734 is the review's, not
  my own trace.
- `git archive` (read-only) extracted the 0.67.0 tree and the two
  merge-parent trees into the session scratchpad.
- Bin 2P groups G5, G6 and G7, G2's link maps, and G1's multivariate
  and wiener blocks were surveyed, not ported: none has a thin read of
  an exported frmtmb object (the group table says what each would
  need). `autocor_matrix()` takes a fit, not a parameter, so G6 is not a
  thin read either. The unported G3 and G4 blocks (multivariate, ARMA,
  SAR, FCOR, wiener, custom) and the categorical rows are "not ported
  (cost)" and stay in the denominator; a fuller adapter could restate
  them.
- The posterior_predict passes check only lengths and dims, as brms
  wrote them; I did not add a distributional check of the draws.
- `dev/brmsport-bin2-run.R` stages into a fixed directory per package,
  so two runs of one package at the same time collide; the pipeline
  runs files one after another.
- The bin-1 ledger was not regenerated (the consolidation does it once);
  the composed total reads `dev/brmsport-ledger.tsv` as it is.

## Files

New: `dev/remeasure-run.sh`, `dev/remeasure-provenance.R`,
`dev/remeasure-history.R`, `dev/remeasure-seed-history.R`,
`dev/remeasure-gaps.R`, `dev/remeasure-plauscmp.R`,
`dev/remeasure-secmp.R`, `dev/remeasure-blocks.R`,
`dev/remeasure-fill.R`, `dev/remeasure-replay.R`,
`dev/remeasure-wiener.R`, `dev/remeasure-repro-missings.R`,
`dev/remeasure-bin2-probe1.R`, `dev/remeasure-bin2-probe2.R`,
`dev/remeasure-bin2-absent.R`, `dev/remeasure-bin2-stale-setup.R`,
`dev/remeasure-bin2-fams.R`, `dev/remeasure-bin2-fams2.R`,
`dev/remeasure-stale-reach.R`, and at fix round 1
`dev/remeasure-bin2-probe3.R`, `dev/remeasure-bin2-sayno.R`,
`dev/remeasure-compois-means.R`, and at fix round 2
`dev/remeasure-bin2-wiener-g1.R` and `dev/brmsport-bin2-split.tsv`, and
at fix 2b `dev/remeasure-lcdf.R`;
`dev/brmsport-bin2-*.R` (blocks, gen,
run, survey, ledger, adapters, adapters-spec, verdicts-2t,
verdicts-2p),
`dev/brmsport-bin2-record.sh`, `dev/brmsport-bin2-all.sh`, the
generated `dev/brmsport-bin2/`, the tracked tables
`dev/brms-port/outcomes-mech.tsv`, `outcomes-hand.tsv`,
`outcome-builds.tsv`, `dev/brmsport-bin2-*.tsv`,
`dev/brmsport-bin2-summary.md`, and this file. Changed:
`dev/brms-port/families-coverage.R` (the extension probe; at fix round
1, a fit per family), `dev/brms-port/plausibility.R` (`PORT_BRMS`, the
worst-rows list). Kept from the review: `dev/remeasure-rev-*.R` and
their logs.
Logs, local only: `dev/remeasure-out/`, `dev/brmsport-bin2-log/`.
