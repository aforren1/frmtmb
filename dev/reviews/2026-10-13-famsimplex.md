# Review: lane famsimplex (brms's simplex families), 2026-10-13

## Round 2 (targeted): MERGEABLE AFTER ONE NEWS FIX

This round checks fix round 1 (findings section 9) on the installed
lane build. Scripts are `dev/famsimplex-rev2-*.R`, and each log is
beside its script.

Remaining:

1. **The BREAKING note is incomplete.** The new refusal in
   `parse_one_response()` applies to every family with several
   locations, extension families included. It now refuses
   frmtmb.eam's `lba()` and `rdm()` constants on `v<k>`, and
   frmtmb.latent's `hmm()` constants on `mu<k>`. `lca()` has several
   locations too; I did not run it. On rellib-r8 each of these
   constants was silently ignored: `lba(3)` with `v3 = 0.8`, `rdm(2)`
   with `v2 = 1.4` and `hmm(2)` with `mu2 = 3` give the same logLik
   and `opt$par` as the free fit (-245.6213688723 for hmm). So the
   refusal is right for them too. But core NEWS names only
   categorical, multinomial, mixture and the simplex families, and
   frmtmb.eam's and frmtmb.latent's NEWS are unchanged. Name the
   extension families in core NEWS, and add a Breaking-changes entry
   to each of the two extensions' NEWS.

Recommended, not required: the refusal says "Write a formula for it
instead, as v3 ~ 1, which estimates it". That gives no route to the
value the user meant to hold. An offset-only formula holds a location
on this build: `bf(y ~ x, mu1 ~ 0 + offset(o))` on a mixture and
`v3 ~ 0 + offset(o)` on `lba(3)` both fit, without that location's
coefficients (`dev/famsimplex-rev2-offset.log`). One clause in the
message would give that route.

### The guard-absent side: identical to rellib-r8

`dev/famsimplex-rev2-guard.R lane|base`, then `diff` of the two logs.
The logLik is printed to 10 digits and `opt$par` to 8 digits. These
cases are identical on the two libraries:

- mixture `theta1 = 0.5`, `theta1 ~ x`, `sigma2 = 1` and
  `sigma1 = "sigma2"`;
- mixture `mu2 ~ 1`, and categorical `muc ~ 1`;
- `mu2 = "mu1"`, which gives brms's "Equating parameters of class 'mu'
  is not allowed." on both;
- `mvbf()` with `sigma = 1` in one gaussian response, and
  categorical + gaussian with `sigma = 2`;
- nonlinear `sigma = 0.3`, and a constant nlpar `b = 0.3`;
- cumulative `disc = 1`, and gaussian `sigma = 1`;
- wiener `bias = 0.5`, and `bias = 0.5, ndt = 0.2`;
- lba `A = 0.5`. The lba `ndt = 0.2` case is refused the same way on
  both, an error that predates this lane;
- rdm `A = 0.5`, and hmm `sigma1 = 0.6`.

`phi = 5` on dirichlet (450.8987887798) and `sigmac = 1` on
logistic_normal (519.0634227082) fit on the lane build. Base has no
such families, so these two have no comparison.

The only differences are the intended refusals:

- mixture `mu1`/`mu2`, categorical `mub`/`muc`, and `mvbf()`
  categorical `muc = 0`. On base the mvbf case fitted with 7
  parameters, the constant ignored.
- lba `v1`/`v3`, rdm `v2`, and hmm `mu1`/`mu2`.

The message names exactly the held locations. On
`bf(y ~ x, mu1 = 0, mu2 = 1, theta1 = 0.5)` it lists "mu1 = 0, mu2 =
1" and not theta1. No test, vignette or R source in core or the
extensions writes a constant on a location (grep), and no extension
builds `pfix` itself, so no existing suite meets the refusal.

### Spot checks (`dev/famsimplex-rev2-spot.log`)

- **Item 2.** `dharma_residuals()` refuses all three families in
  words. `residuals()` of every type (response, pearson, deviance,
  osa) and with `newdata` reaches one simplex sentence that points to
  `predict()`, and no longer points to osa or dharma.
- **Item 3.** The `method = "posterior_predict"` refusal now gives a
  true reason: one simulated composition is a row over every
  category. On dirichlet and dirichlet2 it advises `epred` or
  `dpar = "mub"` / `"mua"`. On logistic_normal it advises `dpar` only,
  and the advised `dpar` call answers.
- **Item 4.** `dev/famsimplex-rev2-ulps.R` runs the three density
  blocks of `test-simplex-families.R` with each lpdf multiplied by
  `(1 + eps)`. At eps 0 the worst values are 3.68, 1.05 and 1.83 of
  the bound 64. At eps 1e-13 all 16 expectations fail (worst 360, 160
  and 371), and at 1e-12 all 16 fail (3595, 1604 and 3722). The new
  tolerance detects an error 10 times smaller than 1e-12 and keeps a
  17-fold margin on the clean code.
- **refcat.** `dirichlet(refcat = 2)` on unnamed columns gives `mu1`
  and `mu3`. On columns a, b, c it gives brms's "must be one of a, b,
  c. Got refcat = '2'". `logistic_normal(refcat = 2)` gives `sigma1`,
  `sigma3` and `lncor(1,3)`. `refcat = 2.5` and `refcat = TRUE` become
  `"2.5"` and `"TRUE"`, as brms's `as_one_character()` makes them.
- **`?confint`** (rendered with `Rd2txt`) names `thetar_*` and
  `lncor_raw_*` as coordinates of a row-normalized Cholesky factor,
  not correlations. It also says that `summary()` (`lncor(b,c)`) and
  `hypothesis()` (`lncor__b__c`) give the natural scale.
  `?dharma_residuals` names compositions.

I did not rerun the suites. The lane's fix1, fix2 and OpenBLAS 0.3.26
records are in findings section 9. The R CMD check results there are
1 NOTE for frmtmb and OK for frmtmb.sample.

# Round 1

Reviewer, adversarial. Base d2ef6312 (0.70.0, `rellib-r8`), lane build
`wt-famsimplex-lib`, brms 2.23.0. Read-only on the lane's files. Every
script named here is `dev/famsimplex-rev-*.R` (or `.sh`) in this
worktree; its log is the `.log` beside it (local, gitignored).

## Verdict: MERGEABLE AFTER FIXES

Required, each with evidence below:

1. **A constant on a location dpar is silently dropped, and the parser
   change widens it to the first location.** `bf(y ~ x, mu1 = 0)` on a
   `mixture()` and `bf(f ~ x, mub = 0)` on a `categorical()` are refused
   on rellib-r8 and accepted on the lane, and the fit ignores the
   constant (same logLik and coefficients as without it). The new
   families inherit the older silent drop for every non-first location:
   `bf(Y ~ x, muc = 0)` on `dirichlet()` and `logistic_normal()` fits
   `muc` freely. (R1)
2. **`dharma_residuals()` dies with an internal error on a simplex fit,
   and the residuals() refusal points the user at it.** (R2)
3. **The `conditional_effects(method = "posterior_predict")` refusal
   gives a false reason on the simplex families** ("whose response is
   a category ... no further observation noise") and, on
   `logistic_normal()`, sends the user to `method = "epred"`, which is
   then refused. (R3)
4. **The two density tests use a tolerance that is not relative to the
   quantity's conditioning.** Measured margin 224 of 1000 ulps on the
   test's own seed; 1136 ulps on random draws of the same kind. (R4)
5. Small: `R/families.R:5375` is 88 columns (the findings say it was
   rewrapped to 80); `R/compat.R:839` reads `covs <-c(`. Correct the
   record's false sentence on `gamma(0.01, 0.01)` (R6).

Everything else the lane claims, I reproduced or found consistent. The
refusals have brms parity (R0), the densities agree with independent
references (R5), and the parser change leaves every previously working
fit bitwise identical (R1).

## R0. Every refusal, against brms 2.23.0 run on the lane's cached fits

`dev/famsimplex-rev-brms.R` (brms, on `dev/famsimplex-brmsfit-*.rds`)
beside `dev/famsimplex-rev-frm.R` (lane, same data and calls).

| call | brms 2.23.0 | lane |
|---|---|---|
| `residuals()` any type, any `method` | "Predictive errors are not defined for ordinal or categorical models." | refused (own sentence) |
| `predictive_error()` (draws) | same as above | refused (own sentence) |
| `pp_check(type = "stat")` and default | "'pp_check' is not implemented for this family." (`has_multicol()`, every type) | brms's sentence, every type |
| `fitted()`, `posterior_epred()`, ln | "Cannot compute expected values ... 'logistic_normal'." | brms's sentence |
| `conditional_effects()`, ln, `categorical` TRUE or FALSE | same refusal | same refusal |
| `conditional_effects(method = "posterior_predict")`, dir and ln | "Please set 'categorical' to TRUE." then "Can only use 'categorical' with method = 'posterior_epred'." | refused, wrong reason (R3) |
| `fitted(dpar = "mub")`, `posterior_epred(dpar = "mub")`, ln | answers (3 x 4, 20 x 3) | `fitted(dpar = "mub")` answers 3 x 4 |
| `conditional_effects(dpar = "mub")`, `dpar = "sigmab"` / `"phi"` | answers | answers (estimate matches; `se__` is link scale, pre-existing, R7) |
| `bayes_R2()` | refused | refused (ML) |

No refusal fires where brms answers. Shapes on newdata (rows 5, 1, 3):
`fitted()` 3 x 4 x 3 with `P(Y = a)` layers and NULL row names,
`predict()` 3 x 4 x 3 with `a, b, c`, identical to brms's dimnames on
dir, dirre and ln; dirichlet2 likewise (`dev/famsimplex-rev-misc.log`).
The categorical display's column names equal brms's; at x = -2.299,
category a: lane 0.2766 (se 0.0222), brms 0.2765 (sd 0.0227).

## R1. The parser change: blast radius and the silent constant

`dev/famsimplex-rev-parse.R lane|base`, compared bitwise by
`dev/famsimplex-rev-parsecmp.R`: categorical, categorical with
`muc ~ z`, categorical with `(1 | ID | g)`, mixture, mixture with
`mu2 ~ z`, and `mvbf()` with a categorical response give IDENTICAL
logLik, coefficients and standard errors on the two libraries (all six
TRUE). Multinomial (`dev/famsimplex-rev-mn.R`) identical as well.

Which formula wins: `bf(y ~ x, mu1 ~ z)` on a mixture gives `mu1` on
`z` and `mu2` on `x`, as brms's standata does (`X_mu1: Intercept z`,
`X_mu2: Intercept x`, `dev/famsimplex-rev-brmsparse.log`). When every
location has its own formula, the main formula's right-hand side is
unused on both sides, and an NA in its variable drops no row on either
(lane nobs 300, brms N 300).

The defect (`dev/famsimplex-rev-pfix.R`, logs `-lane` and `-base`):

```
mixture y ~ x            logLik -501.1224413371   (both)
mixture y ~ x, mu2 = 0   logLik -501.1224413371   (both: pre-existing)
mixture y ~ x, mu1 = 0   lane: -501.1224413371    base: "dpar(s) not available"
categorical, mub = 0     lane: -311.5006016534 (= no constant)  base: refused
dirichlet Y ~ x, muc = 0 456.8262864196 = dirichlet Y ~ x; muc_x estimated
logistic_normal, muc = 0 muc_x = -0.50863, estimated
```

In `parse_one_response()` the dpars loop tests `dp %in% primaries`
before `pfix`, so a location's constant is never read. The widened
`allowed` admits `pfix` on the first location for the first time.
brms 2.23.0 also ignores `mu1 = 0` on a mixture (its Stan code builds
`mu1` from `Intercept + x`), so this is brms parity, but it is a
user-written constant silently replaced by an estimate: clearly wrong,
and the tiebreaker applies. Minimal fix in this lane: widen `allowed`
for `pforms` only, and refuse `pfix` on any location of a family with
several (naming the first ones too, since the old silent drop now
reaches two new families). Pin it with a test that expects the
refusal: for `mu1 = 0` it passes on rellib-r8 and fails on the lane
build; for `muc = 0` on `categorical()` it fails on both.

## R2. `dharma_residuals()` on a simplex fit

`dev/famsimplex-rev-dharma.log`, `dev/famsimplex-rev-generics.log`:

```
dharma_residuals  ERROR [in sims[keep, , drop = FALSE]]:
                  (subscript) logical subscript too long
```

`residuals(type = "response" | "pearson")` reaches `residuals_one()`
first, whose sentence ends "On an ordinal fit, type = "osa" gives
randomized quantile residuals ..., and dharma_residuals() is the
simulation-based check". On a composition that sends the user to a
call that crashes. `dharma_residuals()` refuses `categorical`; add the
simplex families with a sentence, and let a simplex fit reach the
simplex refusal the lane wrote in `residual_values()` for every type.
The crash is new with these families: on rellib-r8 a multinomial fit
stops earlier in `fitted_point()` with a frmtmb sentence ("declares no
mean"; `dev/famsimplex-rev-dharma-mn.log`); a simplex family declares
a mean, so the call gets as far as the matrix response.

Other generics on a dirichlet and a logistic_normal fit with `(1 | g)`
(`dev/famsimplex-rev-generics.log`): `ranef`, `coef`, `VarCorr`,
`confint_varcorr`, `fitted(re_formula = NA)`, new levels,
`frm_bootstrap()`, `confint()` boot and profile, `cooks.distance`,
`influence`, `simulate(re_formula = NA)`, `emmeans(dpar = "mub")`,
`update(REML = TRUE)` all answer. Two do not, neither required here:
`marginaleffects::predictions()` on the dirichlet fit dies inside
data.table on the matrix column `Y` (a multinomial fit fails there on
base too, with marginaleffects's own message,
`dev/famsimplex-rev-basecat.log`); `drop1()` says "No linear predictor
for dpar 'mu'", as it does for `categorical()` on base.

## R3. The `method = "posterior_predict"` sentence

```
method = "predict" has no meaning on family 'dirichlet', whose response
is a category: the category probabilities conditional_effects() draws
ARE the predictive distribution, so there is no further observation
noise to add. Use method = "epred" (the default) ...
```

False for a composition: a Dirichlet or logistic-normal draw adds
noise to the mean. On `logistic_normal()` the advised `epred` is
refused in the next call. brms refuses this call too, so keep the
refusal and fix the reason for the simplex families.

## R4. Density-test tolerance

`dev/famsimplex-rev-ci.log`, `dev/famsimplex-rev-dens.log`:

```
test-simplex-families dirichlet density: worst 224 ulps against a bound of 1000
dirichlet2 density: worst 48 ulps
200 random dirichlet cases: worst relative to max(1,|ll|) 1136 ulps;
  relative to the sum of |terms| 6.674 ulps
40 random cases (K 2..6, phi to e^6): dir 3180 ulps, ln 128, dir2 230
```

The error is `lgamma()` cancellation: `lgamma(phi) - sum(lgamma(a_k))`
is large beside a small result. The test passes on its seed on
Windows, OpenBLAS 0.3.26 (lane) and 0.3.32 with LAPACK (mine, R8), but
its bound is a constant times `max(1, |ll|)`, not a ratio to the
conditioning the run measures, and macOS and glibc have their own
`lgamma()`. Scale the error by the sum of the absolute terms (7 ulps
measured) in `test-simplex-families.R` and the `log_lik()` check of
`test-simplex-draws.R`. No test of the lane sits on a rounding-decided
branch: the `refcat = NA` test's SE verdict rests on an exactly flat
direction (softmax shift), and it passed on all three arithmetics.

## R5. Densities and fits, independently of brms's Stan code

`dev/famsimplex-rev-dens.R`, seed 20261013, 40 random cases, K 2 to 6,
random reference (and `NA` for dirichlet):

- logistic_normal against `mvtnorm::dmvnorm()` on the alr scale minus
  `sum(log(y))`: worst 2.84e-14 relative.
- dirichlet and dirichlet2 against a hand `lgamma()` formula: 7.05e-13
  and 5.10e-14 relative (R4 for why dirichlet is the larger).
- Whole objective: a K = 4 `logistic_normal(refcat = "c")` fit with
  `sigmaa ~ z`, logLik 1006.9741591420 against the mvtnorm sum at its
  estimates 1006.9741591420 (6.8e-16 relative). A K = 12 unnamed
  dirichlet: logLik 9122.5713693688 against the hand sum, 2.0e-16
  (`dev/famsimplex-rev-misc.log`).
- The lncor order: on a K = 4 fit the ML correlations equal the
  residual correlations of the alr regressions to 4 digits (an identity
  when the scales carry no formula), in brms's column-by-column upper
  triangle order (b,c), (b,d), (c,d) (`dev/famsimplex-rev-vars.log`).
  `draws_lncor_block()`'s inverse uses `us_theta_cor()`'s
  `lower.tri()` order, so it holds at K >= 5 too.

Recovery under a 3 x 3 correlated `(1 | ID | g)` on K = 4, reference
LAST, 30 groups, n = 600, seeds 5001 to 5060
(`dev/famsimplex-rev-recov.R`): both families 60 of 60 fitted, 0
nonzero convergence codes, 0 warnings; every coefficient bias within
1.43 Monte Carlo SE; Wald coverage 0.90 to 0.967 at 60 replicates; sds
0.498 / 0.406 / 0.285 (dir) and 0.495 / 0.405 / 0.281 (ln) against
0.5 / 0.4 / 0.3; correlations (b,a) 0.52 and 0.53 against 0.5, (c,a)
-0.285 and -0.254 against -0.3 (2.0 MC SE, ordinary at 30 groups),
(c,b) 0.15 and 0.186 against 0.2.

## R6. Response validation and refcat

`dev/famsimplex-rev-valid.log`, `-misc.log`, `-routes.log`,
`-brmsobj.log`:

- Zeros, a value of 1, a row summing to 2, an all-NA row, an
  off-simplex row whose predictor is NA, a data.frame column, special
  characters in category names (`a.b`/`ab`/`c_d` gives `muab`,
  `mucd` as brms's `make_stan_names()`), K = 12 unnamed (`mu2` ...
  `mu12`, `refcat = "10"` correct), `update()` with renamed columns,
  `mvbf()` with a dirichlet response (logLik equals the sum of the two
  separate fits), weights against replicated rows (903.4874 both):
  all agree with brms or are correct.
- `subset()` that drops an off-simplex row: accepted, logLik equals the
  fit without the row (417.0746). `nobs()` says 300 where brms's N is
  299; gaussian does the same on base (pre-existing).
- Rounded data (4 decimals): both packages refuse (84 rows off by up to
  1e-4); renormalizing fixes it.
- The per-row tolerance refuses data brms accepts AND fits in one
  case: `logistic_normal()` with one row off by 1e-7 (brms samples; the
  log-ratio density is defined off the simplex). For `dirichlet()` brms
  accepts the data and then fails in Stan. I accept the deliberate
  divergence; the record should say the logistic-normal case is a
  refusal of a fit brms completes.
- **Parity gap, small:** `dirichlet(refcat = 2)` is refused ("with one
  string ... Got numeric"), while brms's `as_one_character()` turns it
  into `"2"`: `brms::dirichlet(refcat = 2)` passed to `frm()` works on
  unnamed columns and brms's standata accepts it. Accept a number and
  coerce it, as brms does.
- brms's `validate_prior()` accepts class `"Llncor"` as well as
  `"lncor"`; frmtmb refuses `"Llncor"`. frmtmb refuses `"Lrescor"` and
  `"L"` the same way, so this is a pre-existing convention.
- **Record correction:** findings section 5 says brms's
  `gamma(0.01, 0.01)` on `phi` "is a density set_prior() does not
  carry". It is not so: `set_prior("gamma(0.01, 0.01)", class = "phi")`
  is accepted by `frm()` (MAP phi 6.570) and by `frm_sample()`
  (`dev/famsimplex-rev-sample.log`). The missing default is the
  existing policy for every family's phi or shape, and the disclosure
  still prints "no defaults for phi" when the user gave phi a prior
  (pre-existing). Fix the sentence; the default is a follow-up for all
  families, not this lane.

## R7. frmtmb.sample, confint, and pre-existing items

- `frm_sample()` on ln (`sigmab ~ z`), ln with `lkj(4)` on lncor,
  dirichlet2, and dirichlet with `(1 | ID | g)`: all sample; the
  disclosure lists `student_t(3, 0, 2.5)` on the sigmab intercept,
  `sigmac` natural, `lkj(1)` on lncor (absent once the user's
  `lkj(4)` replaces it), as brms's prior_summary has; `lncor__b__c`,
  `sigmac` and `phi` are draws columns under brms's names; dirichlet2's
  `posterior_predict()` rows sum to 1 within 2.2e-16.
- Draws of a group block are `theta_1..3`, not brms's `sd_`/`cor_`;
  identical on base for a categorical fit
  (`dev/famsimplex-rev-sample-base.log`): pre-existing, not this lane.
- `confint()`'s `lncor_raw_1` (0.426 [0.303, 0.549] beside
  `lncor(b,c)` 0.392 [0.293, 0.484]) is the chart coordinate. It is the
  same convention as `thetar_1` for `set_rescor(TRUE)` on base, and
  `?confint` says covariance parameters are on the internal scale. Not
  a wrong answer by the package's convention, so not a fix for this
  lane; but `?confint` names `theta_*` only, and should name
  `lncor_raw_*`. Pre-existing for rescor and lncor alike:
  `confint(parm = "lncor__b__c")` is refused with "does not stand for a
  single internal one here", false at K = 3, and the advised
  `hypothesis(method = "profile")` then fails ("not linear; use method
  = 'boot'"). File it.
- `conditional_effects(dpar = "phi")`: estimate 7.364 (brms 7.312),
  `se__` 0.054 against brms's 0.383. The lane's `se__` is the link
  scale; a gaussian `sigma ~ z` fit on the lane does the same, so it is
  pre-existing.

## R8. OpenBLAS 0.3.32 with LAPACK

`dev/famsimplex-rev-ob.sh` (the 0.3.32-lapack R copy built by lane
ci701, used read-only; `OPENBLAS_NUM_THREADS=4`, gated):

```
RESULT test-simplex-families.R pass=100 fail=0 err=0 skip=0 warn=0
RESULT test-simplex-brms.R pass=15 fail=0 err=0 skip=0 warn=0
RESULT test-simplex-draws.R pass=25 fail=0 err=0 skip=0 warn=0
```

## R9. Overlap with lane famuni (read-only look at its worktree)

Textual conflicts, all trivial: the `brms_mu_link_source` and
`brms_mu_links` lists in `R/links-brms.R` and the `src` vector of
`dev/famlink-gen-links.R` (both append after `categorical`), the
`fams` vector of `compat_features_build()` in `R/compat.R`, NEWS. Take
the union of the `src` lists and rerun the generator. Neither
regeneration changes an existing family's links (both diffs only add
entries; famsimplex also adds three names to `brms_no_mixture`).
`test-simulate-density.R`: famsimplex sets the header count to 307,
and famuni adds three `sim_specs` rows that add statistical
assertions, so the count must be recomputed after the merge.

## R10. NEWS and docs

Accurate against what I ran, with two exceptions: the "Bug fixes"
bullet must say what R1's fix does to `mu1 = <constant>`, and the
`?frmtmb-families` claim that `residuals()` "is refused" should not
lead to `dharma_residuals()` until R2 is fixed. Style: the NEWS
sentences are long for ASD-STE100 ("The main formula applies to every
category, `bf(Y ~ x, muc ~ z)` gives one category its own formula, and
..." is three statements in one sentence); split them. No em dashes,
US spelling.

## Cleanup

My scratch is `dev/famsimplex-rev-*` (scripts, logs, rds) and
`dev/famsimplex-rev-ob/`; nothing large. I did not install anything.
