# Lane nb2dm: negbinomial2 and dirichlet_multinomial

Round opened 2026-10-20, base 4854e425 (frmtmb 0.72.0). Brief:
`next-round-briefs.md`, section "Lane nb2dm". Status: DONE 2026-10-20;
the sections below are current as of their last edit.

Every number below was pasted from the log named beside it (logs under
`dev/nb2dm-out/`, gitignored). Library: `wt-nb2dm-lib` (this tree),
then `rellib-r10`, then the user library.

## What changed

### negbinomial2 (`brmsfamily("negbinomial2")`)

- `fam_negbinomial2()` (R/families.R) builds `fam_negbinomial()` and
  hands it to a new `fam_reparam()`, which wraps EVERY slot of the
  inner family that takes `dpars` (top-level functions and every
  `post$*`) so that it first adds `shape = 1 / sigma`, and keeps every
  other slot (accepts_aterms, valid_y, type, ...) as the inner one has
  it. `reparam_slot()` keeps each slot's formals (`fam_lcdf()` reads
  their number) and stores the inner function as the attribute
  `frmtmb_inner`, which `mean_is_mu()` (R/predict.R) now unwraps.
- The translation is `shape = 1 / sigma` and nothing else. (The first
  build also set `.eta_shape = -log(sigma)`; removed in fix round 1,
  because it changed no reachable result.) `valid_y` is
  `count_y("negbinomial2")`, so a refusal names the family the user
  asked for.
- No exported constructor, as in brms ("not yet officially supported";
  `"negbinomial2" %in% getNamespaceExports("brms")` is FALSE,
  dev/nb2dm-out/brms-probe.txt). Reached by `brmsfamily()`, by
  `family = "negbinomial2"` and by a brms family object.
- brms has no zero-inflated or hurdle negbinomial2: the brms namespace
  holds `.family_negbinomial2` and no `.family_*negbinomial2*` variant
  (brms-probe.txt, "negbinomial2 variants").
- Links: brms's sets, generated (`dev/famlink-gen-links.R`, now with
  negbinomial2 and dirichlet_multinomial; `R/links-brms.R`
  regenerated).
- Least-squares start: `start_ls_rules$negbinomial2` is negbinomial's
  rule with `sigma = 1 / shape` (R/start-ls.R).
- compat: vocabulary, group `discrete`, row `rate() x negbinomial2`
  works.
- frmtmb.sample: `brms_family_prior()` gains brms's hook,
  `sigma ~ gamma(0.4, 0.3)` for a sigma without a formula.

### What negbinomial2 reuses from negbinomial (for the cdfs merge)

negbinomial2 calls `fam_negbinomial(link, "log")` at construction and
reuses, through `fam_reparam()`: `lpdf`, `sim`, `post$mean_fn`,
`post$var_fn`, `post$dev_fn`, `valid_y`, `init_dpars$mu` (and
`init_dpars$shape`, inverted, for sigma), `accepts_aterms`, `type`, and
ANY slot negbinomial gains later: `lcdf`, `lccdf`,
`post$trunc_mean_fn`, a new aterm in `accepts_aterms`. Nothing is
named, so the cdfs lane's CDF reaches negbinomial2 with no edit here.
Asserted by test-negbinomial2.R "every negbinomial slot reaches
negbinomial2, a later one too", which mocks `fam_negbinomial()` with an
added `lcdf` and `cens` term and reads both off
`brmsfamily("negbinomial2")`. At the merge, compat needs the cdfs rows
for negbinomial repeated for negbinomial2 (the same test file compares
the two families' `cens()` status and will fail until they agree).

### dirichlet_multinomial (`dirichlet_multinomial()`)

- One implementation with multinomial(): `fam_multinomial_impl()`
  takes `link_phi`; with it, the family is brms's
  `dirichlet_multinomial_logit2_lpmf()`, alpha = phi softmax(mu),
  T = rowSums(y). The response check, trials() check, categories read
  from the response, `refcat`, `K`, the mean (`post$mean_fn`, phi does
  not enter) and the deferred placeholder are multinomial's code.
  The simulator draws the probabilities from the Dirichlet and then
  the counts (brms's posterior_predict_dirichlet_multinomial()).
- The two lgamma differences of brms's formula are formed with
  `lgamma_shift_diff()`, the student-t remedy, so the density keeps
  its digits as phi grows (section "dirichlet_multinomial at a large
  phi" below). brms parity is kept: frmtmb differs from brms's Stan
  function only by brms's own error against an exact reference.
- `multinomial_families` / `fam_is_multinomial()` replace every
  `identical(fam$family, "multinomial")` in core (families.R,
  conditional-effects.R, diagnostics.R, frame.R, predict.R, priors.R)
  and frmtmb.sample (`draws_count_cats`, methods-draws.R, sample.R).
- Starting phi by moments (`dm_phi_start()`).
- compat: vocabulary, groups `matrix_response` and `trials_families`,
  rows trials() (conditional), simulate, fitted (works).
- frmtmb.latent's `hmm()` still names "multinomial" only, so a
  dirichlet_multinomial state distribution is refused with its "needs a
  'mu' parameter" message. Not changed: no brms counterpart, and the
  extension is not in this lane's scope.

## Densities against brms's Stan code

Script `dev/nb2dm-density.R`, log `dev/nb2dm-out/density.txt`. brms's
`dirichlet_multinomial_logit2_lpmf()` is cut out of the program brms
writes; negbinomial2 is brms's statement
`neg_binomial_2_log_lpmf(Y | mu, inv(sigma))` and its rate() form
`neg_binomial_2_log_lpmf(Y | mu + log_denom, inv(sigma) .* denom)`
(both read from stancode() in the log). Compiled with
`rstan::expose_stan_functions()` (rstan 2.32.7). frmtmb is evaluated on
doubles and through an RTMB tape carrying `.eta_<dpar>`, as a fit sees
it. The exact references are sums of logs of rising factorials.

```
negbinomial2: 384 grid points
  finite: stan 384 dnbinom 384 frmtmb doubles 384 tape 384
  dbl_nb vs stan: max rel 1.42e-10 (at 50/3/1e+06); vs dnbinom: max rel 1.96e-09
  tap_nb vs stan: max rel 1.42e-10 (at 50/3/1e+06); vs dnbinom: max rel 1.96e-09
  stan vs dnbinom: max rel 1.96e-09 (at 1/0/1e-08)
  vs the exact reference (max rel): frmtmb tape 5.95e-12, stan 1.38e-10, dnbinom 1.96e-09
    sigma 1e-08  max rel to exact: frmtmb 7.8e-14  stan 5e-15  dnbinom 2e-09
    sigma 0.0001 max rel to exact: frmtmb 2.7e-14  stan 1.1e-13  dnbinom 8.4e-14
    sigma 0.01   max rel to exact: frmtmb 1.7e-13  stan 1.3e-13  dnbinom 1.9e-13
    sigma 0.3    max rel to exact: frmtmb 8.1e-14  stan 4e-14  dnbinom 1.4e-13
    sigma 1      max rel to exact: frmtmb 1.1e-13  stan 1.4e-13  dnbinom 1.6e-13
    sigma 10     max rel to exact: frmtmb 1.7e-13  stan 2.4e-13  dnbinom 8.6e-14
    sigma 1000   max rel to exact: frmtmb 2.2e-14  stan 1.7e-12  dnbinom 2.5e-14
    sigma 1e+06  max rel to exact: frmtmb 6e-12  stan 1.4e-10  dnbinom 6e-12
  d/dlog(sigma), 287 points (sigma in [1e-4, 1e3]): tape vs digamma form max rel 0.000103, tape vs Richardson 4.56e-05
    worst digamma point y 1 eta -6 sigma 0.0001: tape -2.475568939e-07 digamma -2.475824072e-07 Richardson -2.475681882e-07
  finite tape derivatives over the whole grid: 384 of 384
  rate(): 81 points, max rel vs stan 1.97e-14

dirichlet_multinomial: 126 cases
  frmtmb doubles vs stan: max rel 6.31e-07; tape vs stan: 3.68e-07
  frmtmb doubles vs exact: max abs 2.52e-07; stan vs exact: max abs 6.11e-07
  by phi, max abs error against the exact reference:
    phi      frm     stan
1 1e-04 9.00e-12 8.64e-12
2 5e-02 2.49e-13 1.14e-12
3 1e+00 5.61e-13 9.09e-13
4 8e+00 3.24e-13 9.09e-13
5 3e+02 9.09e-13 9.09e-13
6 1e+05 3.74e-10 4.92e-10
7 1e+08 2.52e-07 6.11e-07
  identical(frmtmb doubles, stan): 33 of 126
  K = 2 equals beta_binomial():
    max abs diff 2.98e-07 over 42 cases
```

Reading: negbinomial2 is within 6e-12 relative of the exact value at
every grid point; brms's Stan function is within 1.4e-10 and
`stats::dnbinom` within 2e-9 (its Poisson branch at sigma = 1e-8). The
dirichlet_multinomial half of this block is the FIRST build, which
evaluated brms's lgamma-difference formula as written: its error grew
with phi as brms's does. The final build is in the next section but
one. The worst derivative disagreement is at a
derivative of -2.48e-7 where the three values differ in the fifth
digit, which is the references' cancellation, not the tape's.

## Fit agreement with brms

Script `dev/nb2dm-brmsfit.R <case>`, logs
`dev/nb2dm-out/brmsfit-<case>.txt`, summary
`dev/nb2dm-out/brmsfit-summary.txt`. Data seed 2026 (n = 400, 20
groups), brms and frm_sample() 4 chains x 2000, seed 101; brms with
`init = 0` (see the brms defect below). Each estimate is in units of
brms's posterior sd; `diff_z_mcse` is the posterior-mean difference
over the MCSE of the difference.

```
### case nb2 (dev/nb2dm-out/brmsfit-nb2.txt)
frmtmb ML fit: 0.27 s, logLik -772.380930513 convergence 0
         par brms_mean brms_sd     ml     ml_z sample_mean sample_sd diff_z_sd
 b_Intercept    0.7523 0.04883 0.7529  0.01233      0.7533   0.04789  0.020028
         b_x    0.5588 0.04946 0.5561 -0.05501      0.5566   0.05076 -0.044705
       sigma    0.4118 0.06430 0.4025 -0.14534      0.4121   0.06657  0.004427
max |ml - brms mean| / brms sd: 0.145
max |sample - brms mean| / brms sd: 0.0447
max |sample - brms mean| / mcse: 1.79
brms rhat max: 1.002  divergences: 0
### case nb2re
max |ml - brms mean| / brms sd: 0.478
max |sample - brms mean| / brms sd: 0.0498
max |sample - brms mean| / mcse: 1
brms rhat max: 1.007  divergences: 0
### case dm
max |ml - brms mean| / brms sd: 0.0437
max |sample - brms mean| / brms sd: 0.0315
max |sample - brms mean| / mcse: 1.33
brms rhat max: 1.002  divergences: 0
### case dmre
max |ml - brms mean| / brms sd: 0.469
max |sample - brms mean| / brms sd: 0.0397
max |sample - brms mean| / mcse: 0.97
brms rhat max: 1.013  divergences: 0
```

The dm and dmre cases were rerun on the final build (the
lgamma_shift_diff() density; brms fits read from their rds): the ML
log-likelihoods are unchanged to 12 digits (-1970.99219107,
-1926.48569719) and the summaries became

```
dm:   max |ml - brms mean| / brms sd: 0.0437; |sample - brms| / sd 0.041; / mcse 1.72
dmre: max |ml - brms mean| / brms sd: 0.469;  |sample - brms| / sd 0.0325; / mcse 0.872
```

The cost, measured (`dev/nb2dm-tapecost.R`, log `tapecost.txt`: one
process, interleaved blocks of 32 fn + gr evaluations of a
`(1 | g)` fit's objective, minimum of 7 rounds, with a control arm):

```
shift plain   ctl
 1.53  0.75  1.58
ratio shift / plain: 2.04  control shift / shift: 1.033
```

(CORRECTED in fix round 1: that run was on a loaded machine.) On a
quiet machine the review measured 1.92x with a group term, where the
Laplace inner problem evaluates the density many times per outer
evaluation (`dev/nb2dm-tapecost.R` rerun, review log
`dev/nb2dm-rev-out/tapecost-lane-script.txt`: ratio 1.919, control
0.966), and 1.24x for a fixed-effects fit (`dev/nb2dm-rev-cost.R`,
n = 2000, K = 4, `phi ~ z`, 512 evaluations per block, minimum of 9
rounds, log `cost2.txt`: ratio 1.239, control 1.018). Kept: the
alternative is the noise-driven stop above the multinomial limit
measured above.

The two 0.47-sd ML gaps are the group sds (nb2re 0.5189 against a
posterior mean of 0.5741; dmre 0.5008 against 0.5524 and 0.3074
against 0.3429): an ML mode against a posterior mean of a variance
component, the ordinary shortfall (lane-rules "A shortfall in a
variance component"). frm_sample() under brms's default priors agrees
with brms's posterior to 0.05 sd and to 1.79 MCSE at most, over 20
parameters; sd ratios 0.96 to 1.08.

## Default priors against brms

Script `dev/nb2dm-priors.R`, logs `dev/nb2dm-out/priors.txt` (lane)
and `priors-before.txt` (rellib-r10). Six models: negbinomial2 y ~ x,
bf(y ~ x, sigma ~ z), y ~ x + (1 | g); dirichlet_multinomial Y |
trials ~ x, bf(..., phi ~ z), ... + (1 | g).

```
lane:    rows with a brms density that differ, all cases: 0
         (2, 2, 3, 3, 3, 5 rows where brms states a density; all equal)
base:    every case: "negbinomial2 is not a supported family" /
         "dirichlet_multinomial is not a supported family"
```

The rows brms leaves "" and frmtmb prints "(flat)", the class-wide
`sd || g` row frmtmb fills in, and frmtmb's `theta` rows are the
neighbors' behavior on rellib-r10 too (`dev/nb2dm-priors-ctl.R`, log
`priors-ctl.txt`: negbinomial and multinomial with (1 | g)).

## The bin-2 rows

Adapter (`dev/brmsport-bin2-adapters.R`): `dirichlet_multinomial` joins
`brms_adapt_catlike` and is built from its bare constructor on brms's
`ncat` categories, as `multinomial` is; the per-category map reads
`primary_dpars` (it read every dpar, which made phi a fourth `^mu`
dpar). negbinomial2 needed nothing: the adapter's
`frmtmb::brmsfamily("negbinomial2")` now builds it. The verdict files
are not edited (consolidation regenerates them).

Recorded runs of the generated log_lik and posterior_predict files
(`dev/brmsport-bin2-run.R`, one process each), before on rellib-r10
and after on the lane build; logs
`dev/brmsport-bin2-log/nb2dm-{before-,}rec-<topic>.tsv` and
`nb2dm-{before-,}run-<topic>.txt`:

```
RESULT frmtmb test-brms-suite-bin2-log_lik.R pass=84 fail=0 err=0 skip=0 blocks=12
RESULT frmtmb test-brms-suite-bin2-posterior_predict.R pass=55 fail=0 err=0 skip=0 blocks=12
rows whose `held` differs, before -> after (all rows compared: 261 log_lik, 207 posterior_predict)
assert:log_lik:288	FALSE	TRUE
assert:log_lik:510	FALSE	TRUE
setup:log_lik:287	FALSE	TRUE
setup:log_lik:508	FALSE	TRUE
setup:log_lik:509	FALSE	TRUE
value:log_lik:510	FALSE	TRUE
assert:posterior_predict:169	FALSE	TRUE
assert:posterior_predict:355	FALSE	TRUE
setup:posterior_predict:168	FALSE	TRUE
setup:posterior_predict:353	FALSE	TRUE
setup:posterior_predict:354	FALSE	TRUE
value check log_lik:510: HOLDS: 400 values (0 NA, 0 infinite), max abs diff 2.13e-14, max rel 1.97e-15
```

So the four rows flip from `cannot transfer / absent` to holding:
`log_lik:288`, `posterior_predict:169` (negbinomial2) and
`log_lik:510`, `posterior_predict:355` (dirichlet_multinomial), and no
other row moves. The value check of `log_lik:510` is against brms's
own `log_lik_dirichlet_multinomial()` on the block's prep. That
settles consolidation item 4 of the 0.72.0 round: the
dirichlet_multinomial rows hold with the family built, through the
same adapter path as `multinomial()`'s `log_lik:504`. At consolidation
`dev/rel072-bin2-moves.R` still sets `log_lik:510` and
`posterior_predict:355` to `absent`, and
`dev/brmsport-bin2-verdicts-2p.R` sets 288 and 169 to `absent`: those
lines go. (A first after-run without frmtmb.eam in the lane library
showed `log_lik:565` and `posterior_predict:420`, the wiener rows,
failing with "there is no package called 'frmtmb.eam'": the runner
reads the lane library and the user library only. With every
extension installed into the lane library they hold, as above.)

`tests.posterior_epred.R` (posterior_epred_trunc_negbinomial2,
posterior_epred_dirichlet_multinomial) and `tests.stancode.R` call the
two families too; neither file is in the bin-1 or bin-2 ledgers.

## brms 2.23.0 / Stan defect met

brms writes negbinomial2 with no sigma formula as
`neg_binomial_2_log_glm_lpmf(Y | Xc, Intercept, b, inv(sigma))`. On the
nb2 data, seed 101, default inits, one chain of four went to
b_Intercept -45.75, sigma 9.2e-24 (Rhat 1.56, 139 divergences). At
that draw the GLM function returns +2.06e13 where
`neg_binomial_2_log_lpmf` and `dnbinom` give -37774.33
(`dev/nb2dm-brms-stuck.R`, log `brms-stuck.txt`):

```
chain 2 last draw: b_Intercept -45.75 b_x 15.82 sigma 9.236e-24 | glm 2.061584302e+13 plain -37774.33341 dnbinom -37774.33341
```

A positive log density of 2e13 is an attractor the sampler cannot
leave. It is StanHeaders 2.39.1's GLM function at phi = 1.1e23, so
brms's negbinomial() (same primitive) is exposed the same way. frmtmb
does not use it. Filed in `dev/upstream-bugs.md` (Stan-1). The fit
agreement above uses `init = 0`, under which all four chains mix.

## dirichlet_multinomial at a large phi: lgamma_shift_diff()

The first build evaluated brms's formula as written, lgamma
differences. On data with no overdispersion phi runs off, and there
the differences are what the optimizer sees. `dev/nb2dm-dm-limit.R`
(12 seeds of multinomial data, n = 300, 3 categories, trials 5..40)
on that build, log `dev/nb2dm-out/dm-limit-plain.txt`:

```
 seed        ll_mn        ll_dm          excess      log_phi conv
    3 -1171.836109 -1171.836100 9.285570059e-06 17.734814472    1
    4 -1179.172736 -1179.172627 1.090267763e-04 19.057321605    1
    8 -1194.659490 -1194.659486 4.637959819e-06 17.471012351    0
    9 -1156.910518 -1156.910514 3.870506134e-06 17.110021106    0
   10 -1196.799552 -1196.798976 5.760455892e-04 21.085478975    1
   11 -1202.599270 -1202.599259 1.075331625e-05 17.741803234    0
   12 -1207.549686 -1207.549682 4.542176384e-06 17.389694860    0
```

(the five other seeds kept a finite phi, log phi 5.3 to 9.0: chance
overdispersion, where beating the multinomial is legitimate). All 7
runaway fits report a log-likelihood ABOVE the multinomial supremum
they can only approach, by up to 5.8e-4, and 3 of 7 stop with "false
convergence (8)": the student-t failure `lgamma_shift_diff()` was
written for. `dev/nb2dm-shiftdiff.R` (log `shiftdiff.txt`) measured
it against exact rising-factorial sums for shifts up to 2000 (its
documentation covers 100): worst relative error 6.3e-12 (at a = 1e-6,
where the plain difference has the same error), against 7.35e-2 for
the plain difference at a = 1e15. The density now forms both
differences with it. After (log `dev/nb2dm-out/dm-limit.txt`, same 12
seeds):

```
fits with phi run off (log phi > 12): 7 of 12
  of those, above the multinomial by more than 1e-8: 0 ; max excess -9.16464e-08
  of those, optimizer code not 0: 0
fits with finite phi: optimizer code not 0: 0
```

and `dev/nb2dm-density.R` rerun on the final build
(`dev/nb2dm-out/density.txt`):

```
dirichlet_multinomial: 126 cases
  frmtmb doubles vs stan: max rel 4.1e-07; tape vs stan: 4.1e-07
  frmtmb doubles vs exact: max abs 9.09e-12; stan vs exact: max abs 6.11e-07
  by phi, max abs error against the exact reference:
    phi      frm     stan
1 1e-04 9.09e-12 8.64e-12
2 5e-02 5.12e-13 1.14e-12
3 1e+00 4.55e-13 9.09e-13
4 8e+00 2.27e-13 9.09e-13
5 3e+02 4.55e-13 9.09e-13
6 1e+05 1.71e-12 4.92e-10
7 1e+08 1.36e-12 6.11e-07
  K = 2 equals beta_binomial():
    max abs diff 4.22e-07 over 42 cases
```

The K = 2 gap is now `beta_binomial()`'s own error at phi = 1e8
(RTMBdist's plain lgamma form); beta_binomial() was not changed. The
two tests that pin this ("the density keeps its digits as phi grows",
"on multinomial data phi runs off and the fit stays below the limit")
fail with the plain difference put back (`dev/nb2dm-naive-arm.R`, log
`naive-arm.txt`: 4 and 2 failed expectations).

## Tests

New: `tests/testthat/test-negbinomial2.R`,
`tests/testthat/test-dirichlet-multinomial.R`,
`tests/testthat/test-nb2dm-brms.R` (gated),
`extensions/frmtmb.sample/tests/testthat/test-nb2dm-draws.R`.
Changed: `test-simulate-density.R` (a negbinomial2 row in the
agreement tier, a dirichlet_multinomial count-vector test, the
coverage list, the header's count regenerated with
`dev/rel071-simcount.R`: 383 statistical assertions, log
`simcount.txt`), `test-brms-families.R` (the link-field test named
"multinomial"; it reads `multinomial_families` now).

Seen to fail before the change (logs `dev/nb2dm-out/before.txt`,
rellib-r10): test-negbinomial2.R err=11 of 11 blocks ("negbinomial2 is
not a supported family"), test-dirichlet-multinomial.R err=13 of 13
("dirichlet_multinomial is not a supported family"; corrected in fix
round 1 from "10 of 10": the file grew two large-phi blocks, 12 of 12
in the review's `rev-before.txt`, and the clamp block, 13 of 13 in
`before2.txt`),
test-simulate-density.R 1 fail 2 err, test-nb2dm-draws.R err=3 of 3.
The frmtmb.sample changes were also seen to fail on their own, with
the lane's core and rellib-r10's frmtmb.sample (a junction library,
log `t-sample-corenew-samplebase-detail.txt`): the sigma default came
out `student_t(3, 0, 2.5)` where brms has `gamma(0.4, 0.3)`;
`posterior_linpred()` returned `200 x 200` where brms's shape is
`200 x 200 x 2`; `predictive_error()` stopped on the matrix check
instead of the polytomous refusal. The two large-phi tests fail with
the plain lgamma difference (`dev/nb2dm-naive-arm.R`).

One file per process, all eight suites, lane library (every `lib:`
line: `wt-nb2dm-lib`, 414 of 414), log `dev/nb2dm-out/all.txt`, on the
build before the lgamma_shift_diff() change:

```
frmtmb 256 files pass=18445 fail=0 err=1 skip=197 warn=0
frmtmb.coupling 11 files pass=551 fail=0 err=0 skip=5 warn=0
frmtmb.eam 32 files pass=1838 fail=0 err=0 skip=5 warn=0
frmtmb.latent 11 files pass=368 fail=0 err=0 skip=2 warn=0
frmtmb.learn 15 files pass=430 fail=0 err=0 skip=13 warn=0
frmtmb.ode 11 files pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample 62 files pass=2878 fail=0 err=0 skip=8 warn=0
frmtmb.spline 16 files pass=594 fail=0 err=0 skip=1 warn=0
```

The one error was the new dirichlet_multinomial block of
test-simulate-density.R (a newparams spelling in the test, fixed).
After the density change, rerun one per process (`dm2.txt`, gated
on): test-dirichlet-multinomial 64/0/0/0 (rerun alone after its last
edit),
test-simulate-density 621/0/0/0, test-negbinomial2 65/0/0/0,
test-nb2dm-brms 20/0/0/0, test-multinomial 17/0/0/0, test-refcat
102/0/0/0, test-nb2dm-draws 29/0/0/0; test-refcat-brms (gated) 17/0/0/0.

One session per suite, three orders, final build
(`dev/release/nb2dm-sess.log`, `run-sessions.sh nb2dm-sess 8`, every
`lib:` line the lane library): 24 of 24 sessions, fail 0, err 0 in
each; frmtmb 18473 / 18575 / 18675 passes (forward / reverse /
preload), sample 2878 x 3. Every LEAK line carries `[with load of
...]`; none is untagged.

OpenBLAS with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack` and
`0.3.32 lapack`; logs `ob0326.txt`, `ob0332.txt`): test-negbinomial2,
test-dirichlet-multinomial, test-simulate-density, test-multinomial,
test-refcat, test-optstart, test-families, test-nb2dm-draws: fail 0,
err 0 under both (65, 64, 621, 17, 102, 69, 252, 29 passes). The
OpenBLAS R copies were deleted afterwards.

## R CMD check --as-cran

`dev/nb2dm-check.sh` (built and checked in `dev/nb2dm-check/`, since
deleted; logs `dev/nb2dm-out/check-frmtmb-full.txt`,
`check-sample-full.txt`), R_LIBS lane, rellib-r10, user,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH:

```
frmtmb:        Status: 1 NOTE   (HTML manual: "Skipping checking math rendering: package 'V8' unavailable")
frmtmb.sample: Status: OK
```

## Bump and floor

frmtmb minor (two new families). frmtmb.sample minor (new defaults
and draws paths); its new test file builds `dirichlet_multinomial()`
and `brmsfamily("negbinomial2")`, so its floor on frmtmb rises to the
release that carries them. No other extension changed.

## Overlaps for consolidation

R/families.R, R/compat.R (vocabulary, groups `discrete`,
`matrix_response`, `trials_families`, four rows), R/links-brms.R and
dev/famlink-gen-links.R, R/links.R, NEWS.md, man/frmtmb-families.Rd,
man/frmtmb-links.Rd, frmtmb.sample's R/sample.R and NEWS.md. With
lane cdfs: a negbinomial CDF reaches negbinomial2 with no edit; its
compat rows for negbinomial need negbinomial2 copies
(test-negbinomial2.R compares the two families' `cens()` status).

## Fix round 1 (review `dev/reviews/2026-10-20-nb2dm.md`)

1. **Small-phi Dirichlet draws (MEDIUM).** `simplex_rdirichlet()`
   (R/families-simplex.R), which dirichlet(), dirichlet2 and
   dirichlet_multinomial() share, draws `log G = log Gamma(a + 1) +
   log(U) / a` and normalizes by a softmax. `dev/nb2dm-rdir-arm.R`, log
   `dev/nb2dm-out/rdir-arm.txt`, 20000 rows, seed 5:

   ```
   phi 0.01: NaN rows of 20000, probability-scale draw 18, log-scale draw 0
   phi 0.001: NaN rows of 20000, probability-scale draw 9570, log-scale draw 0
   ```

   New test file `tests/testthat/test-simplex-rdirichlet.R`: no NaN row
   at phi 0.01 and 1e-3 for dirichlet() and dirichlet2 (and none for
   dirichlet_multinomial's sim); means, variances and the covariance of
   two categories against the analytic Dirichlet within z(1e-6) Monte
   Carlo SEs at phi 0.01, 0.5 and 20, 1e5 draws; `predict()` and
   `simulate()` on a dirichlet_multinomial fit with phi = 0.02 data
   (seed 41, n = 200; the fit's phi is below 0.1). On rellib-r10
   (`before2.txt`): pass=28 fail=22 err=2, all three blocks failing,
   the dirichlet and dirichlet2 rows with NaN (the behavioral failure,
   on families that exist there). With the old draw put back into the
   lane build (`rdir-arm.txt`): 8 fails and an error, 14 fails, and the
   predict() block errors.
   brms 2.23.0 fails the same way (`dev/nb2dm-brms-rdirichlet.R`, log
   `brms-rdirichlet.txt`): `rdirichlet()` 18 and 9570 NaN rows of
   20000 on the same seed and alpha, `posterior_predict_dirichlet_
   multinomial()` at phi 1e-3 "NA in probability vector",
   `posterior_predict_dirichlet()` 973 NaN rows of 2000. Filed as
   brms-35 in `dev/upstream-bugs.md`. NEWS: a Breaking changes bullet
   (the draws of dirichlet() and dirichlet2 change for a seed).
2. **`.eta_shape` (LOW).** Removed from `to_shape()`; the comment says
   negbinomial2 is exactly as robust as negbinomial() at the mirrored
   point. The test "the log-scale density input survives a sigma near
   the Poisson" is gone (it pinned `dnbinom_robust()`). The density
   identity test now compares with negbinomial() at shape = 1 / exp(p)
   with no `.eta_shape`, and is bitwise identical (values and Jacobian).
3. **Cost (LOW).** Corrected above: 1.92x with a group term, 1.24x with
   fixed effects only (the review's runs). No doc or NEWS text stated a
   cost.
4. **Refusal name (LOW).** `fam_negbinomial2()` sets
   `valid_y = count_y("negbinomial2")`; new test "a response that is
   not a count is refused in negbinomial2's name".
5. **Before-arm count.** Corrected in "Tests".
6. **Tests.** `expect_error(..., "trials(N) could not be evaluated on
   newdata", fixed = TRUE)`; the mocked `lcdf` returns `pnbinom()`
   probabilities, per the slot contract; new block "fitted() holds the
   Wald bounds of a count between 0 and trials" (newdata x = 3, 5, 8 at
   50 trials: category b's upper bound at x = 8 is held at the sum of
   the row's estimates, 50 to rounding, and every lower bound of b is
   above 1, which the probability clamp would cap at 1); new block in
   test-nb2dm-draws.R asserting `frmtmb.sample:::draws_count_cats`
   identical to `frmtmb:::multinomial_families`.
7. **Style.** The registry's `cumulative` alignment and the second blank
   line; `[[ ]]` in `start_ls_rules$negbinomial2`.

### Fix round 1: runs (final build)

One file per process, all eight suites (the shared simplex draw
changed), 415 of 415, every `lib:` line the lane library, log
`dev/nb2dm-out/all2.txt`, summary `all2-summary.txt`:

```
frmtmb 257 files pass=18529 fail=0 err=0 skip=197 warn=0
frmtmb.coupling 11 files pass=551 fail=0 err=0 skip=5 warn=0
frmtmb.eam 32 files pass=1838 fail=0 err=0 skip=5 warn=0
frmtmb.latent 11 files pass=368 fail=0 err=0 skip=2 warn=0
frmtmb.learn 15 files pass=430 fail=0 err=0 skip=13 warn=0
frmtmb.ode 11 files pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample 62 files pass=2879 fail=0 err=0 skip=8 warn=0
frmtmb.spline 16 files pass=594 fail=0 err=0 skip=1 warn=0
```

The lane's files there: test-simplex-rdirichlet 55/0/0/0,
test-dirichlet-multinomial 68/0/0/0, test-negbinomial2 62/0/0/0,
test-nb2dm-draws 30/0/0/0. Gated (`gated2.txt`): test-nb2dm-brms
20/0/0/0, test-simplex-brms 15/0/0/0, test-refcat-brms 17/0/0/0.

One session per suite, three orders (`dev/release/nb2dm-sess2.log`):
24 of 24, fail 0 and err 0 in every session; frmtmb 18529 / 18631 /
18731 (forward / reverse / preload), frmtmb.sample 2879 x 3; 36 LEAK
lines, 0 untagged.

OpenBLAS 0.3.26 with LAPACK (`ob0326b.txt`), the simulator files:
test-simplex-rdirichlet 55, test-dirichlet-multinomial 68,
test-negbinomial2 62, test-simulate-density 621, test-simplex-families
125, test-simplex-mvbf 26, sample test-simplex-draws 26,
test-simulators 54, test-nb2dm-draws 30 passes, fail 0 and err 0 in
each. The OpenBLAS copy was deleted.

R CMD check --as-cran (`check-frmtmb-fix1.txt`,
`check-sample-fix1.txt`): frmtmb `Status: 1 NOTE` (V8 math rendering),
frmtmb.sample `Status: OK`.

## To file

- (For consolidation's backlog, not this lane's) beta_binomial() on
  binomial data stops with "false convergence" in 4 of 8 seeds (review
  `runoff.txt`, seeds 3 to 6): the plain lgamma difference this lane
  replaced for dirichlet_multinomial. Not touched here: lane cdfs is
  editing beta_binomial().

## For consolidation

- Add "negbinomial2" to lane cdfs's `compat_cdf_discrete`, with
  `stan = "neg_binomial_2"` in its `switch()` and negbinomial2 in its
  `tmean` set (review note 1). On the review's trial merge negbinomial2
  inherits the CDF, `cens()` and `trunc()` with no edit and its fits
  equal negbinomial()'s, but `frm_compat()` reports "refused" until the
  row exists, and test-negbinomial2.R's compatibility block fails.
- frmtmb.sample's floor rises with its tests; add a "Needs frmtmb"
  line at the release.
- The bin-2 `absent` lines in `dev/rel072-bin2-moves.R` and
  `dev/brmsport-bin2-verdicts-2p.R` for the four rows go.

## Not done, and why

- frmtmb.latent `hmm()` with a dirichlet_multinomial state
  distribution: refused by its existing message (see above).

## Defects found and not fixed

- The joint-density translator of the gated brms tests
  (`stan_pars_from_fit()`, helper-brms.R) misreads uncorrelated
  per-category group terms: `multinomial()` with `(1 | g)` fails
  check C by 8.04 log units with a gradient of 2.45, and
  dirichlet_multinomial by 126.7 (seed 65, `nbdm_dm_data()`), while
  `(1 | ID | g)` passes. The model is not at fault for
  dirichlet_multinomial: frm_sample() of the (1 | g) model agrees with
  brms's posterior (case dmre above). Pre-existing for multinomial().
- RTMB-4 (`dnbinom_robust()` gradient in the size at very large sizes)
  reaches negbinomial2 as it reaches negbinomial: log(sigma) near -30
  is log(shape) near 30.

## Log

- 2026-10-20: lane started. The cdfs worktree had no changes when this
  lane started (`git status` clean), so negbinomial2 delegates to the
  negbinomial family object structurally rather than to named CDF
  functions.
