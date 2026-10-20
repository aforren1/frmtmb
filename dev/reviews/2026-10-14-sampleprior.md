# Review of lane sampleprior (2026-10-14)

Worktree `frmtmb-wt-sampleprior`, branch `wt-sampleprior`, uncommitted
over 47cf659d. Lane library `C:/Users/adf44/source/r/wt-sampleprior-lib`
(frmtmb and frmtmb.sample), base rellib-r9 (0.71.0, frmtmb.sample
0.19.0), brms 2.23.0. The reviewer's scripts are
`dev/sampleprior-rev-*.R`; their logs are in `dev/sampleprior-rev-log/`
(local, gitignored). No lane file was edited. No library was installed.

## Verdict

MERGEABLE AFTER FIXES.

The core of the lane holds up under an independent check of the
objective itself, not the table: on 20 models of the reviewer's choice,
the log prior that `frm_sample()` adds to the objective it hands
tmbstan equals brms's Stan `lprior`, carried to frmtmb's coordinates by
a numerical Jacobian, to at most 6.5e-11 (finding list below, "What
holds"). rellib-r9 differs on 15 of the 20. The Jacobian is applied
exactly once for every class and link tried, the seed works, and the
derived `sd_`/`cor_` columns are right.

One defect blocks the merge: `frm_sample()` now fails on every
Student-t multivariate model with `set_rescor(TRUE)`, which rellib-r9
sampled (finding 1). The table counted that row as equal because
`default_prior()` lists it; the sampler never applies it.

## Findings

### 1. HIGH: `frm_sample()` fails on every Student-t `mvbf()` with `set_rescor(TRUE)`

Evidence (`dev/sampleprior-rev-mv.R`, logs `mv-lane.txt`, `mv-r9.txt`;
data seed 20261114, 1 chain x 200, seed 1):

```
lane  mvbind student rescor            -> ERROR: A prior with no resp names no
      parameter of a multivariate model (class=nu). ... No response of this
      model has this slot
lane  mvbf student rescor              -> ERROR: (the same)
lane  fit route mvbind student rescor  -> ERROR: (the same)
lane  mvbind gaussian rescor           -> OK, draws 100 x 9 : rescor__y1__y2
r9    mvbind / mvbf / fit route        -> OK, draws 100 x 9 : nu
```

`default_prior(route = "sample")` lists `gamma(2, 0.1)` on `nu` with no
`resp`, which is why the lane's table has this row as `=`. Cause: the
new default writes the shared `nu` with `resp = ""` (sample.R,
`rs <- if (multi && !isTRUE(lp[["shared"]])) ...`), and core's
`resp_missing_refusal()` refuses class `nu` without `resp` in a
multivariate model. The rescor check in `resolve_priorlist()` refuses it
WITH `resp` too, so no spelling reaches the shared `nu`: a user's
`set_prior("gamma(2, 0.1)", class = "nu")` is refused on rellib-r9 as
well, with either spelling (`mv2-r9.txt`). That half is pre-existing;
the lane's default turns it into a failure of every such model.

Fix: in `R/priors.R`, `resp_missing_refusal()`, return `NULL` when
`!is.null(spec[["rescor_nu"]])` and the specification's `class` or
`dpar` is `"nu"`. With that patch applied in the session
(`REV_PATCH_MV=1 dev/sampleprior-rev-objective.R mv_rescor`), the
objective's log prior equals brms's to 6.5e-11, including `nu`'s
`gamma(2, 0.1)` through `logm1` and the rescor `lkj(1)`
(`objective-lane-mvpatch.txt`). Add a test that samples a Student-t
rescor model with the defaults and one with a user's
`set_prior(class = "nu")`, each seen to fail before the fix. Lane
prilabel also edits `R/priors.R`; name the overlap.

### 2. MEDIUM: `uniform()` without `lb`/`ub` is a NaN wall that stops or breaks NUTS

The new density is right as a density (finding list, item 4), but a
user who writes it as brms users do gets a sampler that cannot start or
that diverges. `dev/sampleprior-rev-uniform2.R` and `-uniform4.R`
(fit `y ~ x` on the reviewer's wide gaussian data, 2 chains x 1000,
seed 5):

```
uniform(2.2, 2.6) on b_x, formula route : ERROR the sampler returned no draws
                                          ("no more scalars to read")
uniform(2.2, 2.6) on b_x, fit route     : ERROR (the same)
uniform(2.2, 2.6), init_jitter = 0      : OK, 648 divergences of 1000
uniform(2.2, 2.6) with lb = 2.2, ub = 2.6: OK, 0 divergences
uniform(-10, 10)                        : OK, 0 divergences
```

brms 2.23.0 warns twice on the same prior ("It appears as if you have
specified a lower bounded prior on a parameter that has no natural lower
bound ... please specify argument 'lb'", and the same for `ub`;
`dev/sampleprior-rev-brmsunif.R`). The core NEWS bullet says the NaN
"rejects the value as brms's bounds do", which is not so: brms's bounds
transform the parameter, the NaN only rejects a proposal.

Fix: in `set_prior()`, a `uniform(a, b)` with no `lb`/`ub` takes
`lb = a`, `ub = b` (frm_sample() then bounds the internal parameter
through Stan's transform, as the `lb = 2.2, ub = 2.6` row shows), or at
least warns as brms does. Reword the NEWS bullet.

### 3. LOW: the `ndt` default leaves the wall in place for `shifted_lognormal()`

The same mechanism, on the lane's own default. `uniform(0, min_Y)` on a
log-linked `ndt` is NaN past `min_Y`, so NUTS meets a wall where brms
has `real<lower=0, upper=min_Y>`:

- Prior-only tape of the resolved defaults (`dev/sampleprior-rev-prioronly.R`,
  4 chains x 4000, seed 11): 4058 divergent transitions of 8000 (the
  draws are still uniform, KS p 0.086).
- On data (`dev/sampleprior-rev-shiftdiv.R`, data seed 20261116, 4 x
  2000, seed 9): 82 (formula route) and 44 (fit route) divergent of
  4000. rellib-r9 had 64 and 116 from the likelihood's own NaN, with
  every draw of `ndt` at 0.0000, so the lane is better, but not brms.
- The same default with `ub = min_Y` (`dev/sampleprior-rev-ndtub.R`,
  same data and seed, fit route): 0 divergent of 4000, ndt q05/50/95
  0.0599 0.3080 0.4649 against 0.0897 0.3192 0.4618 with the wall.

Fix: write the `ndt` default as
`natural_dpar_prior(<uniform>, "ndt", ...)` with `ub = min_Y` (core's
`internal_bound()` carries it to `log(min_Y)` on the link scale). The
frmtmb.eam `wiener()` ndt needs nothing; its link is bounded already.

### 4. LOW: an identity-link `nu` admits `0 < nu <= 1`, where brms's bound is 1

`student(link_nu = "identity")` is allowed (`dpar_links_above1`). brms
declares `real<lower=1> nu` whatever the link. The lane's default
`gamma(2, 0.1)` has a finite density below 1 and the Student-t
likelihood has no guard there (`dev/sampleprior-rev-links.R`,
`links.txt`):

```
identity nu =  0.90: default-prior nll 108.5097
identity nu =  0.50: default-prior nll 119.2109
identity nu = -0.50: default-prior nll NaN
```

The lane's bounds table read only the default links, so it did not see
this. Fix: give the `nu` default `lb = 1` when the dpar's link image
reaches below 1; `set_prior("gamma(2, 0.1)", class = "nu", lb = 1)`
reaches tmbstan as lower bound 1 on that coordinate
(`dev/sampleprior-rev-nulb.R`). The other natural classes are guarded
under an identity link by their density or the likelihood (finding 5
names the exception).

### 5. MEDIUM (pre-existing in core, now on the default path): `beta(1, 1)` and `gamma(1, r)` lose their support on the tape

rtmb-pitfalls item 22, in the densities this lane now uses as defaults.
`prior_base_logdens()` writes `(shape - 1) * log(x)`, a constant-zero
operand when the shape is 1, and the tape folds it
(`dev/sampleprior-rev-fold.R`, `fold.txt`):

```
beta(1, 1)  at -0.2, 1.2: doubles NaN, NaN, tape 0, 0
gamma(1, 1) at -0.5:      doubles NaN, tape 0.5
gamma(2, 1) at -0.5:      doubles NaN, tape NaN
beta(1, 1) on b_x, b_x = 2.0 : sampler nll 1162.417   (finite)
beta(1, 1) on b_x, b_x = -1.0: sampler nll 445.5767   (finite)
```

Under the default links (logit for `zi`, `hu`, `zoi`, `coi`, `bias`,
`quantile`; log for `bs` and the exgaussian `beta`) the link keeps the
value inside the support, and with `link_zi = "identity"` and
`wiener(link_bs = "identity")` the likelihood is NaN outside it, so the
defaults are not wrong on the cases tried. A user's `beta(1, 1)` on a
coefficient is: it is flat over the whole line on the tape.
`exponential()` has no guard at all. The lane wrote the uniform's guard
for exactly this pitfall and left the densities beside it. The lane's
own Jacobian test (`test-default-priors-dpar.R`, "a natural-scale
default is brms's density with the link Jacobian") evaluates on
doubles only.

Fix: add `nan_unless_pos(x)` to `gamma`, `inv_gamma`, `exponential`
and `lognormal`, and `nan_unless_pos(x) + nan_unless_pos(1 - x)` to
`beta`, in `prior_base_logdens()`; test each on a tape at a point
outside the support (`RTMB::MakeTape`), not on doubles.

### 6. LOW: NEWS and the `init_jitter` docs overstate what `seed` does to the caller's stream

`frmtmb.sample` NEWS: "the caller's stream is left as it was"; the
`init_jitter` doc: "the caller's random stream is left as it was".
After `frm_sample(fit, seed = 3)` the caller's `.Random.seed` has moved
(`dev/sampleprior-rev-seed.R`: `caller .Random.seed unchanged after
frm_sample(seed = 3): FALSE`), because rstan's `sampling()` draws with
`sample.int()`, as the lane's own findings and test say. What is true:
the jitter no longer draws from the caller's stream. Fix: say that, and
that rstan still advances the stream, as under `brm()`.

### 7. LOW: the fast path's two-point check misses a defect confined to a region neither point visits

`dev/sampleprior-rev-names.R` mutates `draws_cov_fast()` and asks
whether the check falls back to the slow path:

```
sd names reversed within a block               -> rejected (slow path)
cor_Intercept__x and cor_Intercept__z swapped  -> rejected (slow path)
cor as cov / var_i                             -> rejected (slow path)
cor (i, j) transposed                          -> accepted (equivalent: M is symmetric)
sd off by 1e-6 only where var > 2              -> ACCEPTED
```

The present code has no branch on the value, so this is hardening, not
a defect: on a correlated 3x3 block plus `(1 | h)` and a nested
`(1 | g:g2)` (data seed 20261118, 2 x 600, seed 21) the fast columns
equal `hyp_env_vals()` per draw to 2.2e-16 over 600 draws, all 8 names
on the fast path. Fix (cheap): in `draws_cov_values()`, also compare the
fast values with `hyp_env_vals()` at the rows holding each theta
column's minimum and maximum, and take the slow path on a mismatch.

### 8. LOW: the negbinomial agreement rests on a seed where brms's chains mix

Rerun of the lane's cases on seed 101, 4 chains x 6000 on both sides
(`dev/sampleprior-rev-agree.R`, `agree.txt`): brms's negbinomial chains
stick (ESS 4, mean log(shape) 56.3, against frmtmb's 0.5069), the
Poisson-limit failure the lane documents. The gated test asserts brms's
R-hat below 1.05 on seed 3 only, so it compares against a reference that
fails on another seed. frmtmb's density there is right (finding list,
item 1, case `negbin`, and the prior-only run). Fix: in
`test-default-priors-fit.R`, either pick a negbinomial data set on which
brms mixes on several seeds (check three), or compare frmtmb's
posterior of log(shape) with a deterministic reference (an
importance-weighted prior draw), and say which in the findings.

## What holds (each claim checked)

### Item 1: the defaults in the objective

`dev/sampleprior-rev-objective.R` captures the objective `frm_sample()`
hands tmbstan (tmbstan replaced in the session), takes its log prior as
the default tape less the `prior = "flat"` tape (both centered, so the
same coordinates), and compares differences over 6 points with brms's
Stan `lprior` written out by hand from `brms::stancode()`
(`dump*.txt`) and carried to frmtmb's coordinates by a central-difference
Jacobian (independent of frmtmb's `mu_eta`). Data seed 20261114, `x`
with mean 3 so that centering moves every intercept.

```
case            max |d frm - d brms|   lane       r9
gauss_wide      (mad(y) 9.9)           1.48e-12   1.48e-12
sigma_x         (sigma ~ x, scale 4.2) 8.88e-14   0.244
lognormal_wide  (log-scale, 4.7)       8.66e-12   8.66e-12
negbin          (inv_gamma, log)       9.39e-12   0.232
beta_phi_x      (phi ~ x, centered)    2.84e-14   0.607
student         (nu logm1, sigma)      1.25e-11   1.81
zip_zi_x        (zi ~ x, logistic)     9.77e-15   1.72
hurdle_ln       (hu beta(1,1), logit)  2.76e-11   0.55
hurdle_ln_hu1   (hu ~ 1)               7.46e-12   0.679
gev_xi1         (xi ~ 1, log1p)        5.93e-12   5.93e-12
cum_thres_gr    (thres(gr =), no centering) 2.2e-11 1.97
cum             (thresholds centered)  8.51e-12   2.33
mixture         (dirichlet(1), sigma1/2) 2.71e-11 0.428
poisson_sqrt    (sqrt mu link)         2.86e-14   2.86e-14
beta_phi_id     (link_phi identity)    2.66e-14   1.26
sigma_id_x      (link_sigma identity, sigma ~ x) 2.66e-14 0.55
shifted_ln      (ndt uniform, log)     1.19e-11   0.5
wiener          (bs, ndt scaled_logit, bias) 2.58e-11 0.285
re_sigma        ((1 + x | g), sigma ~ (1 | g): sd of sigma on mu's 9.9, lkj) 5.15e-11 5.15e-11
mv_rescor       ERROR (finding 1); 6.5e-11 with the patch   1.48
```

Also: softplus `sigma` without a formula (natural, Jacobian
`plogis`) 3.64e-14, with `sigma ~ x` (brms keeps 2.5 on a link outside
its list) 6.31e-14, identity `nu` above 1 1.15e-14
(`links.txt`); a mixture with a constant `sigma1 = 1`, whose map shifts
the `betad` positions, 2.40e-14 (`mixfixed.txt`).

### Item 2: change of variables, sampler only

A prior-only tape of exactly the resolved default entries
(`sample_resolve_priors()` then `neg_log_prior_fn()`, as
`prior_augmented_obj()` composes them, no likelihood), every other
parameter fixed by `map`, sampled by tmbstan (4 x 4000, seed 11;
`prioronly.txt`):

```
student nu ~ gamma(2,0.1)[1,]   KS D 0.0122 p 0.181  q05/50/95 3.877 16.849 46.276 target 3.729 16.858 47.495
zip zi ~ beta(1,1)              KS D 0.0149 p 0.058
shifted_ln ndt ~ U(0, min_Y)    KS D 0.0140 p 0.086  (4058 divergent, finding 3)
wiener bs ~ gamma(1,1)          KS D 0.0081 p 0.672
wiener ndt ~ U(0, min_Y)        KS D 0.0102 p 0.371
wiener bias ~ beta(1,1)         KS D 0.0107 p 0.324
beta Intercept_phi ~ t(3,0,2.5) KS D 0.0108 p 0.307  (a formula: no Jacobian)
mixture theta2 ~ dirichlet(1)   KS D 0.0136 p 0.103
```

The negbinomial `shape` gave KS p 0.000 on that run (21 divergent);
thinned by 10 over 4 x 22000 (seed 12, `prioronly2.txt`) it gives ESS
6930, KS D 0.0077 p 0.723, and the tape's density on log(shape) equals
the reviewer's to 1.8e-15 on a grid: the first run's p was
autocorrelation.

### Item 3: the posterior agreement

The lane's z already pools both chains' MCSEs:
(5.1795 - 5.2163) / sqrt(0.0166^2 + 0.0174^2) = -1.53. Rerun on seed
101 with 4 x 6000 on both sides (`agree.txt`):

```
beta      log(phi): z(mean) -0.57 | z(q05..q95) 0.16 0.59 -1.35 -1.12 -0.38 | sd ratio 0.992
student   log(nu):  z(mean)  0.55 | z(q05..q95) 0.39 0.57 -0.24 -0.25 1.25 | sd ratio 0.993
dirichlet log(phi): z(mean)  1.21 | z(q05..q95) 1.60 1.96 1.31 -0.29 0.14 | sd ratio 0.977
```

That is agreement. The negbinomial is finding 8.

### Item 4: `uniform()` on the tape

`dev/sampleprior-rev-uniform.R`, `uniform.txt`: on doubles NaN at -1.5,
-1, 2, 2.5 and -log(3) inside `uniform(-1, 2)`; on the sampler's tape
the nll is NaN at b_x = -1.5 and 2.5, and inside the uniform's part is
1.098612 = log(3) at -0.999, 0.5 and 1.999; the gradient inside equals
the defaults-only tape's, outside it is NaN; on a MAP fit's tape
(`uniform(-1, 4)`) NaN at -2 and 4.5, finite at 0.5. The constant is
the normalized one, which `er_prior_density()` needs: the lane's
`evidence-ratio.R` change adds `dunif()`/`dlnorm()` so a Savage-Dickey
ratio on a coefficient with a user's uniform or lognormal prior has a
denominator; both are normalized densities. frmtmb.sample offers no
bridge sampling.

### Item 5: seed

`dev/sampleprior-rev-seed.R` (3 chains x 400, `seed-lane.txt`,
`seed-r9.txt`): serial `seed = 3` twice identical (r9: not); seed 4
differs on all 3 chains; parallel `cores = 3` twice identical and
bitwise equal to serial (r9: differ by 11.97); with no seed,
`set.seed(9)` decides as before; formula route identical. 
`with_sample_seed()` restores an L'Ecuyer-CMRG stream and leaves no
`.Random.seed` where there was none. The call as a whole moves the
caller's stream (finding 6).

### Item 6: draw names

On the 3x3 + `(1 | h)` + `(1 | g:g2)` model: the names are
`sd_g__Intercept sd_g__x sd_g__z cor_g__Intercept__x cor_g__Intercept__z
cor_g__x__z sd_h__Intercept sd_g:g2__Intercept`, as `variables(fit)`
lists; at the ML estimate core's values equal lme4's
`VarCorr(lmer(REML = FALSE))` to the optimizer's 5e-5 (`names.txt`), so
the semantics are brms's `VarCorr()` (sd, then correlation of the
covariance). `VarCorr(ds)` means equal the column means. Mutations:
finding 7. The `theta_` columns sit where rellib-r9 has them (20 21 22
of 26, 4 5 6 under `laplace = TRUE`); the derived columns come after
the sampled ones and before `lp__` (`cols-lane.txt`, `cols-r9.txt`);
`posterior_predict()`, `log_lik()`, `summary()`, `VarCorr()`,
`hypothesis()`, `posterior_summary()`, `ranef()`, `loo()`,
`bayes_R2()` and `posterior_epred(newdata =)` run on the new matrix
(`hyp.txt`).

### Item 7: the jitter

Confirmed on 47cf659d (`extensions/frmtmb.sample/R/sample.R`, lines
587 to 597): chain 1 at the mode, chains 2 and up at
`mode + rnorm(length(mode), 0, jitter)` with `init_jitter = 0.25`, on
every coordinate of `last.par.best` (random effects included), then
clamped into the bounds. Isotropic, not SE-scaled.

### Item 8: suites

Against the lane library, gated (`FRMTMB_BRMS_FIT_TESTS=true`), one file
per process (`onefile-*.txt`):

```
RESULT test-default-priors-dpar.R: tests 7 pass 57 fail 0 error 0 skip 0 warn 0
RESULT test-draws-sd-names.R: tests 3 pass 12 fail 0 error 0 skip 0 warn 0
RESULT test-sample-seed.R: tests 3 pass 10 fail 0 error 0 skip 0 warn 0
RESULT test-loo.R: tests 31 pass 142 fail 0 error 0 skip 0 warn 0
RESULT test-sampling-ported.R: tests 34 pass 213 fail 0 error 0 skip 0 warn 0
```

frmtmb.sample's whole suite in one session (`dev/sampleprior-rev-onesession.R`):

```
lib: C:/Users/adf44/source/r/wt-sampleprior-lib/frmtmb.sample  frmtmb: C:/Users/adf44/source/r/wt-sampleprior-lib/frmtmb
ONE SESSION: files 57 tests 441 pass 3028 fail 0 error 0 skip 1 warn 0 (579 s)
```

(`onesession3.txt`, run from PowerShell.) It equals the lane's own
one-session line.

A first one-session run from Git Bash had 3 errors in `test-loo.R` and
`test-sampling-ported.R`, all "Cannot create temporary file in
C:\WINDOWS\" from brms's compile: the reviewer's shell environment, not
the lane; both files pass one per process from PowerShell (above).

### Item 9: style

No em dashes, no British spellings, no added line over 80 columns in the
diff, the new tests or the findings. Comments say why. The frmtmb.sample
NEWS marks the breaking changes in bold under "Breaking changes". Two
statements are wrong rather than badly phrased: the core NEWS bullet on
the uniform's NaN (finding 2) and the stream sentence (finding 6).

## Overlaps

`R/priors.R` (fixes 1, 2, 5; lanes prilabel and catmulti also edit it),
`extensions/frmtmb.sample/R/sample.R` (fixes 3, 4, 7),
`extensions/frmtmb.sample/R/draws-brms.R` (fix 7), both NEWS files
(fixes 2, 6), `test-default-priors-fit.R` (fix 8).

## Re-check of fix round 1

Targeted re-check of the eight findings against the fix-round build
(lane library installed 13:05, after the last source edit at 13:02).
Scripts are `dev/sampleprior-rev-*.R`; logs are
`dev/sampleprior-rev-log/fx-*.txt`.

### Verdict after fix round 1

MERGEABLE AFTER FIXES. Findings 1 to 8 are fixed. The re-check found
one new HIGH defect (finding 9). It was present in the lane's first
build too, and the first review missed it. There is one new LOW
(finding 10).

### Finding 1: fixed, and the exemption is narrow

- `sampleprior-rev-mv.R` (`fx-mv-rev.txt`): mvbind and mvbf Student-t
  rescor models, and the fit route, all sample (`nu rescor__y1__y2` in
  the draws).
- `sampleprior-rev-objective.R` (`fx-objective-rev.txt`), with no
  patch in the session: all 20 cases agree with brms's `lprior`;
  `mv_rescor` agrees to 4.34e-11.
- `sampleprior-rev-mv3.R` (`fx-mv3.txt`), each line against
  `brms::validate_prior()`:

```
student, no rescor, defaults          frmtmb captured             (brms: n/a)
student, no rescor, nu without resp   frmtmb ERROR (no resp)      brms ERROR (no parameter)
student, no rescor, nu resp = y1      frmtmb captured             brms ok (y1 user, y2 default)
gaussian rescor, nu without resp      frmtmb ERROR                brms ERROR
student rescor, nu without resp       frmtmb captured             brms ok
student rescor, sigma without resp    frmtmb ERROR                brms ERROR
```

The exemption reaches only dpar `nu` in a model with `rescor_nu`, and
each answer is brms's.

The lane's `apply_check()` sweep does not build the sampler's
objective. It calls `prior_augmented_obj()` directly. `frm_sample()`
uses `ncp_objective()` for every non-centered block, which includes
the `(1 | g)` of every base model, and the GEV and `mo()` retapes. The
sweep does check that the defaults resolve, and resolution was where
finding 1 failed. I reran the lane's table script with an apply check
that calls `frm_sample()` itself, with tmbstan captured
(`sampleprior-rev-table.R`, `fx-table.txt`). 182 of 185 objectives
build and are finite at their start. The other 3 are the
`gen_extreme_value()` scalar-xi models on the formula route, the
pre-existing defect below. On the fit route their defaults apply. The
lane's claim stands, but its sweep is weaker than it says. Wording fix
for the findings: say it builds `prior_augmented_obj()`, or call
`frm_sample()` with tmbstan captured as the reviewer's copy does.

### Finding 2: fixed

`set_prior()` warnings against brms 2.23.0's `validate_prior()` on the
same model (`sampleprior-rev-boundwarn.R`, `fx-boundwarn.txt`; frm()
adds no warning of its own on any case):

```
beta(2,2) on b coef x                 frmtmb lb+ub   brms lb+ub
lognormal(0,1) on sd                  none           none
uniform(-5,5) on Intercept            lb+ub          lb+ub
gamma(2,1) on b, lb = 0               none           none
uniform(0,10) on b, lb = 0 only       ub             ub
gamma(2,1) on Intercept dpar sigma    lb             lb
vector: gamma on b, normal on Int.    lb             lb
uniform(0,10) on sigma                ub             ub
exponential(1) on b, ub = 5 only      lb             lb
lognormal(0,1) on sd dpar sigma       none           none
gamma(2,0.1) on nu (student)          none           none
uniform(-1,1) on ar                   lb+ub          lb+ub   (fx-boundwarn2.txt)
class lb row + coef gamma row         lb             none    <- finding 10
```

brms declares `real<lower=0> delta` for equidistant thresholds and
does not warn on `lognormal(0, 1)` there (`fx-delta.txt`); frmtmb's
`prior_class_bounds()` gives `delta` its lower bound too.
`prior_class_bounds()` matches brms's `dpar_bounds()` class by class.
The defaults muffle the warning.

### Findings 3 and 4: fixed, with one exception (finding 9)

`brms_natural_bounds()` lists exactly brms's `dpar_bounds()`. The
bounds that reach tmbstan (`sampleprior-rev-natbounds.R`,
`fx-natbounds.txt`) are:

```
shifted_ln ndt, log link        upper -0.583647 = log(min_Y)
shifted_ln ndt identity link    lower 0, upper 0.55786 = min_Y
shifted_ln ndt ~ x              none
student nu logm1                none
student nu identity             lower 1
student nu ~ x, identity        none
Gamma shape identity / softplus lower 0 / none
zip zi identity / zi ~ x        [0, 1] / none
gaussian sigma identity         lower 0
beta phi identity               lower 0
mixture sigma2 identity         lower 0
wiener (bounded ndt link)       none
```

Each bound sits on the transformed coordinate. No bound is written
where a formula sits on the dpar.

### Finding 5: fixed

My `fold.R` cases on the tape now give:
- `beta(1, 1)` at -0.2 and 1.2: NaN.
- `gamma(1, 1)` at -0.5: NaN.
- A user's `beta(1, 1)` on `b_x` at 2 and -1: NaN on the sampler tape
  (`fx-fold.txt`).

On the claim "bitwise unchanged inside the support":
- The math supports it: `u / u - 1` is exactly 0 and so is its
  derivative.
- The lane's control proves the swap reaches `frm()`'s tape only.
- I added the same control on the sampler's tape
  (`sampleprior-rev-guardctl.R`, `fx-guardctl.txt`). With
  `beta(1, 1)` on `b_x` at 1.5, the guarded tape gives NaN and
  0.71.0's densities give 944.29.
- For the zinb defaults, value and gradient are bitwise equal at 20 of
  20 points inside the support.

The claim holds. NEWS nit: `exponential()` read no `log()`; it had no
guard at all, so "read their support through a `log()`" does not
describe it.

### Findings 6, 7 and 8: fixed

- **6.** NEWS and `?frm_sample` now say rstan still advances the
  session's stream.
- **7.** Region defects planted in the fast path
  (`sampleprior-rev-extreme.R`, `fx-extreme.txt`): an `sd` above its
  90th percentile and a `cor` above its 99th percentile are both
  caught. A `cor` inside its 45th to 55th percentile is missed (532
  rows off by 1e-6): extreme draws cannot cover an interior region.
  The code has no value-dependent branch, so I accept the check as the
  hardening asked for.
- **8.** I built an independent reference (`sampleprior-rev-nbref.R`,
  `fx-nbref.txt`). The priors come from brms's Stan code, not
  `default_prior()`: `student_t(3, 0, 2.7)` on the centered
  Intercept and `inv_gamma(0.4, 0.3)` on shape. Quadrature is a
  trapezoid in coordinates standardized at each log(shape)'s
  conditional mode, which is not the lane's grid. Result: mean 0.5116,
  sd 1.4362, q05/50/95 -1.010 0.218 2.819, the lane's reference to the
  printed digits.

### The pre-existing GEV defect: confirmed, not this round

A scalar-xi `gen_extreme_value()` from a formula stops with
"frm_linpred() needs a fitted model" on both the lane and rellib-r9.
The fit route samples on both, xi mean 0.108
(`sampleprior-rev-gev.R`, `fx-gev-lane.txt`, `fx-gev-r9.txt`). It is in
core's `gev_xi_chart_start()` (0.71.0's chart) and touches no prior.
File it in the backlog; it does not block this lane.

### 9. HIGH (new): under `ndt_group()` the `ndt` default caps a fraction at `min_Y`

With `ndt_group()`, frmtmb.eam's `ndt` is a fraction in (0, 1) of each
subject's fastest response, on a logit link. The default applies brms's
`uniform(0, min_Y)` to that fraction. Since fix round 1,
`brms_natural_bounds()` also hands tmbstan an upper bound `logit(min_Y)`
(the link's image top, 1, exceeds `min_Y`).

Evidence (`sampleprior-rev-ndtgroup.R`, `-ndtgroup2.R`;
`fx-ndtgroup.txt`, `fx-ndtgroup2.txt`; data seed 20261120,
`ddm_simulate()` with ndt 0.3, 4 subjects, 2 x 600, seed 4):

```
ndt link logit, image (9.4e-14, 1); default uniform(0, 0.340182823773274); upper bound -0.6624796 on ndt
ML ndt fraction 0.818; true ndt / fastest per subject 0.88 0.78 0.78 0.86
defaults               ndt fraction q05/50/95 0.2943 0.3285 0.3397  max 0.3402
uniform(0, 1) on ndt   ndt fraction q05/50/95 0.7458 0.8083 0.8664  max 0.8973
```

The default excludes the true value: the posterior is pinned against
0.34. Before fix round 1 the uniform's NaN gave the same exclusion
without the bound, so the defect has been in the lane since its first
build. brms has no `ndt_group()`, so this is about internal
consistency, not parity. lba() and rdm() document the same fraction
parameterization; only wiener() was measured.

Fix (in `brms_def_dpar_prior()` / `ndt_uniform_prior()` and
`brms_natural_bounds()`): where the `ndt` link is bounded above, take
the uniform over the link's own image, `uniform(0, linkinv(Inf))`, and
write no bound. That gives `uniform(0, min_Y)` for an ungrouped
wiener(), whose link image is (0, `min_Y`), and `uniform(0, 1)` for
the fraction. Keep `uniform(0, min_Y)` with `ub = min_Y` only for an
unbounded link (log, softplus, identity). Add a test on
`rt | dec(upper) + ndt_group(subj) ~ x` that the default's posterior
contains the ML fraction.

### 10. LOW (new): a class-level bound row does not silence a coef row's warning

`set_prior("", class = "b", lb = 0) + set_prior("gamma(2, 1)",
class = "b", coef = "x")` warns about a missing `lb`. frmtmb applies
the class row's bound to `x` (tmbstan lower 0 on `b_x` and `b_z`,
`fx-boundwarn2.txt`), and brms does not warn on this prior. The
warning runs per `set_prior()` call, before the rows are combined;
brms checks the combined frame, where `stan_base_prior()` reads the
class row's bounds. Fix: run `prior_bound_warning()` on the combined
priorlist (at resolution, or in `+`/`c()`), and read a class-level
`lb`/`ub` row as covering that class's coefficients.

### Tests on the fix-round build

Gated, one file per process, lane library, run from PowerShell:

```
RESULT test-prior-support.R: tests 4 pass 112 fail 0 error 0 skip 0 warn 0
RESULT test-default-priors-bounds.R: tests 5 pass 14 fail 0 error 0 skip 0 warn 0
RESULT test-default-priors-fit.R: tests 4 pass 15 fail 0 error 0 skip 0 warn 0
```

### Remaining fixes

1. Finding 9 (HIGH): the `ndt` default over the link's image when the
   link is bounded above, with a test under `ndt_group()`.
2. Finding 10 (LOW): the bound warning on the combined priorlist.
3. Wording: the sweep's description (finding 1) and the exponential
   NEWS nit (finding 5).

## Check of fix round 2

This is a short, targeted check of the round-2 changes, run on the lane
library: frmtmb was installed at 15:42, after the last `R/priors.R`
edit at 15:41; frmtmb.sample at 15:30, after its last edit at 14:18.
frmtmb.eam comes from rellib-r9 because the lane changed only a test
file there. Scripts are `dev/sampleprior-rev-*.R`; logs are
`dev/sampleprior-rev-log/fx2-*.txt`.

### Verdict after fix round 2

MERGEABLE. Findings 9 and 10 are fixed. The two new core functions do
what the lane says, on doubles and on the tape. The three notes below
are informational and need no change.

### Finding 9: fixed

Data from `sampleprior-rev-ndt2.R` (`fx2-ndt2.txt`): wiener data seed
20261120; lba and rdm from frmtmb.eam's simulators.

```
wiener               ndt default uniform(0, 0.340182823773274) | min rt 0.340182823773274 | bounds: none
wiener + ndt_group   ndt default uniform(0, 1)                 | bounds: none
lba                  ndt default uniform(0, 0.28632970174612)  | min rt 0.28632970174612  | bounds: none
lba + ndt_group      ndt default uniform(0, 1)                 | bounds: none
rdm                  ndt default uniform(0, 0.235316478906408) | min rt 0.235316478906408 | bounds: none
rdm + ndt_group      ndt default uniform(0, 1)                 | bounds: none
brms wiener ndt row: uniform(0, min_Y), ub min_Y (min_Y = 0.340182823773274)
wiener + ndt_group, defaults (2 x 600, seed 4): ML fraction 0.8182 |
  posterior q01/05/50/95/99 0.7117 0.7450 0.8120 0.8605 0.8731 | rhat 1.022
```

- With the default under `ndt_group()`, the posterior covers the ML
  fraction. In round 1 it was pinned at 0.29 to 0.34.
- The ungrouped wiener's default equals brms's `uniform(0, min_Y)`.
- The log-linked `shifted_lognormal()` keeps `ub = log(min_Y)`; this is
  unchanged from round 1, and the lane's test covers it.

### `support_guard()`: correct

`sampleprior-rev-guard2.R`, `fx2-guard2.txt`:

- **The tape keeps the shift.** A guard alone, taped, is 0 at v = 0.
  A folded `+ 1e-300` would give NaN there. Adding a nonzero constant
  is not one of the operations the tape folds (rtmb-pitfalls item 22).
- **Nothing changes for |v| > 1e-284.** On 2336 values with |v| from
  1e-284 to 1e300, `nan_unless_pos(v + 1e-300)` is identical to
  `nan_unless_pos(v)`: on doubles, on the tape, and for the tape's
  gradient.
- **Through `uniform(0, 1)`** (doubles and tape agree everywhere):
  - NaN at -1e-3, -1e-200, -1e-290 and -1e-299.
  - 0 at 0, 1e-310 and 1e-300, with a 0 gradient.
  - 0 at the top edge, x = 1; NaN at 1 + 2.2e-16.
  - The tape's Hessian at the edge is 0, not NaN.
- **Every guarded density, on the tape:**
  - Finite at its edge with a finite gradient: `gamma(1, 1)` and
    `exponential(2)` at 0, `beta(1, 1)` at 0 and 1, `uniform()` at both
    ends.
  - -Inf at the edge, as Stan: `gamma(2, 1)` at 0, `beta(2, 5)` at 0
    and 1.
  - NaN one step outside, every density.

Notes, no change needed:

- **Values in [-1e-300, 0) count as inside.** The guard is finite
  there, for example at -1e-301. A value that close to the edge is
  below anything a link produces, except the identity link exactly at
  the bound.
- **NaN, not -Inf, at exactly 0.** `inv_gamma()` and `lognormal()`
  give NaN at x = 0 (Inf - Inf), where Stan gives -Inf. The density
  is 0 there, both values reject the point, and a log link reaches
  x = 0 only below eta = -745.

### `natural_log_jacobian()`: correct

Checked over eta from -50 to 50 in steps of 0.25 (401 points), against
a 256-bit Rmpfr reference.

**Logit link.**
- Value, on doubles and on the tape: within 2.2e-16 relative of the
  reference; `plogis(log.p = TRUE)` sits at the same 2.2e-16.
- Tape gradient: within 1.1e-16 of 1 - 2p.
- 0.71.0's `log(mu_eta)` was -Inf at 54 of the 401 points, every
  |eta| at or above 36.75.

**frmtmb.eam's `scaled_logit`** (top = min rt = 0.3944):
- Value, on doubles and on the tape: within 2.2e-16 relative of the
  reference.
- Tape gradient: within 1.1e-16 of 1 - 2p.
- Its `linkinv()` is `top * plogis()` to rounding, not bitwise.

**Every other link is untouched.** `log`, `identity`, `softplus`,
`squareplus`, `logm1`, `log1p`, `inverse` and `sqrt` give a Jacobian
identical to `log|mu_eta|`.

### Finding 10: fixed

`sampleprior-rev-warnonce.R` and `-warnonce2.R`, `fx2-warnonce.txt`.
Each count below is lower plus upper warnings, so 2 means one warning
pair: one call warns once.

```
brms: bare set_prior(uniform on b)        0
brms: validate_prior()                    2
brms: stancode()                          2
frmtmb: bare set_prior()                  0
frmtmb: frm()                             2
frmtmb: validate_prior()                  2
frmtmb: par_template()                    2
frmtmb: update(fit, prior =)              2
frmtmb: frm_sample(formula)               2
frmtmb: frm_sample(fit)                   2
frmtmb: frm_simulate()                    2
```

brms warns when it validates the model, not at `set_prior()`. frmtmb
now does the same, once per call, with no repeat through a chain of
callers.

### Tests on the fix-round-2 build

One file per process, lane library, from PowerShell:

```
RESULT test-prior-support.R: tests 7 pass 150 fail 0 error 0 skip 0 warn 0
RESULT test-sample-ndt-default.R (frmtmb.eam, gated): tests 3 pass 13 fail 0 error 0 skip 0 warn 0
```

Under OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`; the copy in `dev/ciharden-out` was deleted after):

```
RESULT test-prior-support.R: tests 7 pass 150 fail 0 error 0 skip 0 warn 0
RESULT test-brms-priors.R: tests 12 pass 103 fail 0 error 0 skip 0 warn 0
RESULT test-prior-compat.R: tests 28 pass 199 fail 0 error 0 skip 0 warn 0
```

### Remaining

Nothing blocks the merge. File the GEV scalar-xi formula-route defect,
as the lane's "To file" does.
