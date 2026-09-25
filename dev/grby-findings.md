# Lane grby: brms's `gr(g, by = f)` and `mm(g1, g2, by = )`

Round of 2026-09-25. Branch `lane/grby`. Base: frmtmb 0.63.0.

## What brms means

brms 2.23.0 gives a by-split term one set of standard deviations and
one correlation matrix per level of the by-variable. Each level of `g`
must belong to exactly one level of `f`. The Stan code
(`dev/grby-log/brms-code.txt`, from `dev/grby-brms-code.R`, nothing
compiled) is `r_j = diag(sd_1[, Jby_1[j]]) L_1[Jby_1[j]] z_j`
(`scale_r_cor_by()`), with `Jby_1` the by-level of each level of `g`.
So the covariance of the term is block diagonal over the by-levels:
levels of `g` in one by-level share that by-level's covariance, and
levels in different by-levels are independent.

brms's refusals, read off `frame_re()` and reproduced verbatim:

- `Some levels of 'g' correspond to multiple levels of 'f'.`
- `Each grouping factor can only be associated with one 'by' variable.`
- `Grouping structure 'mm' expects 'by' to be a matrix with as many
  columns as grouping factors.`

brms's names (`dev/grby-log/brmsnames.txt`, from `rename_re()` on
`brm(empty = TRUE)`, `dev/grby-brmsnames.R`):
`sd_g__Intercept:fa`, `sd_g__x:fa`, `cor_g__Intercept:fa__x:fa`, per
by-level; `sd_g__sigma_Intercept:fa` on a distributional parameter;
`sd_mmg1g2__Intercept:cbind(f1, f2)1` for `mm(g1, g2, by = cbind(f1,
f2))`, where `cbind()` of two factors gives their integer codes. The
brief's `cor_g__Intercept__x:fa` is not brms's spelling: brms pastes the
by-level onto BOTH coefficients. `r_` names and `ranef()` columns keep
the plain coefficient (`r_g[1,Intercept]`). Priors: `get_prior()` lists
class `sd` and `cor` rows with `group = "g"` and `coef = "Intercept"`,
with no by-level, so one specification reaches every by-level.

## What changed

The design is a split, not a new covariance structure. A by-split term
becomes one ORDINARY random-effect block per by-level, over the levels
of `g` in that by-level, each with the structure the term asks for
(`R/gr-by.R` header). The objective, the Laplace machinery, `theta`,
the priors, `vcov()`, `confint()`, the sampler and the simulator read
the blocks as they read any other; the objective is untouched. The
precision over `b` keeps its sparsity: the split adds no coupling, so
the Hessian stays block diagonal by level of `g`.

- `R/gr-by.R` (new): parsing (`parse_gr_by()`), the level-to-by-level
  map with brms's refusal (`by_level_map()`), the split
  (`by_split_component()`), the one-by-variable-per-factor check, the
  merged view brms keys by grouping factor (`by_merged_blocks()`), and
  the routing of new data (`by_route_rows()`).
- `R/parse.R`: `gr()` accepts `by =`; `mm()` accepts `by =` (removed
  from the refused-argument list; `pw =` stays refused); the check
  runs beside `check_id_covstructs()`; the multi-membership help page
  gains a section.
- `R/frame.R`: the by-variable enters the model frame; the component is
  split after it is built, and the duplicate check reads the term as
  written; each block carries `by`.
- `R/brms-names.R`: `brms_re_parts()` suffixes `rnames` with the
  by-level (`Intercept:fa`) and gains `rcoef`, the plain name.
- `R/methods-fit.R`, `R/sugar.R`: `ranef()`, `coef()` and `ngrps()` read
  the merged view, one entry over every level of `g`. `VarCorr()` rows
  follow brms's `get_rnames()` order (by-level first) when every block
  of a group is by-split. `summary()` lists a group's standard
  deviations before its correlations, as `summary.brmsfit()` does; this
  also reorders groups of several non-by blocks, which were listed
  block by block before (brms lists them sd first too).
- `R/predict.R`: new-data routing for plain and multi-membership terms.
- `R/simulate-new.R`: a by-split block has no pre-brms legacy spelling.
- `R/compat.R`: feature `gr_by` (grammar) and its rules.
- `extensions/frmtmb.sample/R/draws-brms.R`: `draws_ranef_layout()`
  merges the by-levels into one column per coefficient over every
  level, and carries brms's `attr(levels, "by")`.

### New data

A level of `g` the fit saw keeps its fitted effect whatever by-value
the new row carries, and needs no by-variable at all; brms indexes the
effect by the level alone. A level the fit did not see is drawn from
the covariance of the by-level the new row names, as brms's
`get_new_rdraws()` draws it (`used_by_per_level`). Refused by name: a
new level whose by-level the fit did not estimate, a new level with a
missing by-value, a new level given two by-values in one `newdata`
(brms's message), and a by-column newdata cannot supply.

### Unused by-levels

brms keeps a by-level that no level of `g` falls in, with an `sd` its
prior alone informs. Maximum likelihood has no estimate for it, so
frmtmb leaves the by-level out and says so in a message. Both packages
drop unused factor levels from the data first (brms's
`drop_unused_levels = TRUE`), so this is reached only when the
by-expression itself creates levels.

## Validation

Script `dev/grby-validate.R`, log `dev/grby-log/validate.txt`
(seeds in the script). All at the ML optimum unless stated.

Against brms's own density. For a gaussian response the Laplace
approximation is exact, so frmtmb's `logLik()` must equal brms's
marginal density, written from its Stan code and evaluated with
`mvtnorm` on brms's own `make_standata()` (`J_1`, `Jby_1`, `Z_1_*`) at
frmtmb's estimates. This is an identity; the residual is rounding.

- `(1 | gr(g, by = f))`, 3 by-levels: frmtmb -396.4994843, brms
  density -396.4994843, difference 2.8e-13.
- `(1 + x | gr(g, by = f))`: -350.2927406 against -350.2927406,
  difference -5.7e-14.
- `(1 | mm(g1, g2, by = cbind(f1, f2)))`: -308.7559174 against
  -308.7559174, difference 1.7e-13.

Against lme4 and glmmTMB at their own optimum, each fitting one term
per by-level with indicator columns, `(0 + fa + xa | g) + (0 + fb + xb
| g) + (0 + fc + xc | g)`, which is the same model:

- gaussian `(1 | gr(g, by = f))`: frmtmb -396.4994843, lme4
  -396.4994845, glmmTMB -396.4994843; standard deviations within
  1.3e-4 relative of lme4's.
- gaussian `(1 + x | gr(g, by = f))`: frmtmb -350.2927406, lme4
  -350.2927406, glmmTMB -350.2927405; standard deviations within
  3.8e-5 relative, correlations within 7.0e-5.
- poisson `(1 + x | gr(g, by = f))`: frmtmb -555.7764683, glmmTMB
  -555.7764683; standard deviations within 1.9e-5 relative.
- `diag(1 + x | gr(g, by = f))`: -351.9323081 against glmmTMB's
  -351.9323081.
- `ar1(0 + t | gr(g, by = f))`: -129.0979522 against glmmTMB's
  -129.0979522; theta within 1.7e-6.
- `cs`, `homcs`, `toep`, `hetar1` over `gr(g, by = f)`: log-likelihood
  differences from glmmTMB 8.7e-10, 2.9e-12, 4.1e-9, 1.3e-9.
- REML `(1 | gr(g, by = f))`: -400.5261024 against lme4's REML
  -400.5261024, difference 1.2e-8; standard deviations within 2.7e-5
  relative.
- `mvbf()` with `(1 | p | gr(g, by = f))` in both responses: against
  glmmTMB on the long format with trait by by-level indicators,
  difference 2.0e-10.

glmmTMB's `homtoep` gives an NA log-likelihood on the indicator
design, and glmmTMB cannot split the num_factor() structures, so
`homdiag`, `homtoep`, `ou`, `exp`, `gau`, `mat` and
`gr(dist = "student")` are checked by factorization: with a fixed
effect and a residual sd per by-level nothing is shared, so the by-split
objective AT the two subset fits' estimates must equal the sum of their
log-likelihoods. Differences: homdiag -1.1e-13, homtoep -1.7e-13, ou
3.1e-13, exp 4.0e-13, gau 8.0e-6, mat -7.8e-8, student -1.8e-12. The gau
and mat residuals are the inner solver's on near-singular kernels: on
this 4-position design every gau and mat fit, joint and subset, reports
false convergence with a warning, so the joint fit's own optimum sits
0.0006 (gau) and 0.012 (mat) below the sum of the subset fits'. That is
where the optimizer stopped, not a different model.

Quadrature: `(1 | gr(g, by = f))`, bernoulli, `quadrature = TRUE`,
against `integrate()` over each level at frmtmb's estimates:
-155.4725458 against -155.4725458, difference 7.2e-10.

The level-to-by-level map is brms's `Jby_1` on brms's own
tests.standata.R data, for `gr(g, by = z)` and for
`mm(g, g2, by = cbind(z, z2))` (`test-gr-by.R`).

`frm_sample()` draws: `ranef()`, `ranef(summary = FALSE)`, `coef()`
and `VarCorr()` are `identical()` to brms's own methods on the same
draws through a `brm(empty = TRUE)` shim
(`extensions/frmtmb.sample/tests/testthat/test-gr-by-draws.R`).

## Speed

The objective is untouched: a by-split term is K ordinary blocks where
the term without by is one, and each block's density is the one it
always had. Nothing is added per evaluation for a model without by.

Measured (`dev/grby-timing.R`, `dev/grby-log/timing.txt`): 2000 rows,
200 levels of `g` in 2 by-levels, `(1 + x | gr(g, by = f))` against
`(1 + x | g)`, arms interleaved in one process, 8 fits per timed block,
minimum of 5 rounds. Seconds per fit: plain 0.1650, by 0.2011, control
(plain again) 0.1726. The control reads 1.046, so load moves the clock
by about 5 percent; the by-split fit costs about 1.2 times the plain
one, for twice the covariance parameters (theta 6 against 3). The
load-independent count: 8 objective and 4 gradient evaluations for the
by-split fit against 10 and 5 for the plain one. An earlier run of the
same script, before the arms were warmed up, read a control of 1.13
and is not used.

## Downstream methods (registered in `R/compat.R`, feature `gr_by`)

Verified by running each on a by-split fit
(`tests/testthat/test-gr-by.R`, `dev/grby-validate.R`, and scratch
checks recorded here): `VarCorr()`, `summary()`, `confint()` (Wald and
profile; a natural name such as `sd_g__x:fb` addresses its `theta`, on
the internal scale as for any block), `confint_varcorr()`,
`hypothesis()` with names that carry `:`, `ranef(condVar = TRUE)`,
`coef()`, `ngrps()`, `fitted()`, `predict()` (in sample, new levels,
`re_formula = NA` and `~ (1 | g)`), `simulate()`, `frm_simulate()`
with brms's names, `residuals()` including OSA, `emmeans`, `anova()`,
`drop1()`, `diagnose()`, `vcov(full = TRUE)`, `conditional_effects()`,
`frm_lp_basis()`, `update()`, REML, quadrature, `sparse_x`,
`autoscale`, `set_prior(class = "sd", group = "g")`, and
`frm_sample()` with `ranef()`, `coef()` and `VarCorr()` on the draws
identical to brms's own methods on the same draws
(`extensions/frmtmb.sample/tests/testthat/test-gr-by-draws.R`).

Refused by name: `importance` (its proposal needs every block over a
factor to carry the same levels; the existing refusal fires).

## Structures

Works, one block of the structure per by-level: `us`, `diag`,
`homdiag`, `cs`, `homcs`, `ar1`, `hetar1`, `toep`, `homtoep`, `ou`,
`exp`, `gau`, `mat`, and `gr(dist = "student")`. Refused by name:
`rr()` (its loadings span the whole factor), `equalto()` (nothing to
split), and `gr(cov = )` / `gr(prec = )`.

### Why `gr(g, by = f, cov = A)` is refused

brms accepts it and writes `scale_r_cor_by_cov()`: each level's
standard-normal draw is scaled by its by-level's covariance and THEN
correlated through the Cholesky factor of `A`,
`b = (L_A (x) I) blockdiag(Sigma_by(j)) z`. Levels in different
by-levels are then correlated through `A`, so the covariance is neither
block diagonal over the by-levels nor a Kronecker product, and the
split would silently fit a different model. It is feasible as a new
covariance structure (`u = (L_A^-1 (x) I) b` has independent levels
with by-level covariances, so the log-density is the per-level normal
of `u` plus `M log|L_A|`), with its own `vcov()`, naming and new-level
rules; that is a lane of its own, not a variant of the split.

## `mm(g1, g2, by = )`

It fell out of the same split: brms maps each POOLED level to one
by-level, read across every member column (the matrix `by`), and the
pooled levels split into one multi-membership block per by-level.
Implemented and validated against brms's density (below). A new
membership level takes the covariance of the by-level its own member
column names.

## Tests

New:

- `tests/testthat/test-gr-by.R`: the split, agreement with lme4 and
  glmmTMB (ar1), brms's density (gr and mm), brms's names, the merged
  `ranef()`/`coef()`/`ngrps()`, new-level routing with the exact
  new-level variance per by-level, every refusal, `simulate()`,
  `frm_simulate()`, brms's `Jby_1`, the registry, and the ordering of
  `VarCorr()` and `summary()`. Tolerances are relative.
- `extensions/frmtmb.sample/tests/testthat/test-gr-by-draws.R`: draws
  `ranef()`, `coef()` and `VarCorr()` identical to brms's methods.

Seen to fail on the base build (`FRMTMB_LIB=base`,
`dev/grby-log/test-gr-by-base.txt`): `RESULT frmtmb test-gr-by.R
pass=0 fail=0 err=15`, every block stopping at the base refusal "gr()
supports (x | gr(g, cov = A)) or (1 | gr(g, prec = Q))". The draws
file against this lane's frmtmb and the BASE frmtmb.sample
(`dev/grby-log/test-gr-by-draws-basesample.txt`): `pass=2 fail=2
err=1`, the ranef/coef block failing, which pins the
`draws_ranef_layout()` change on its own.

Run on this lane's build, one file per process (final code):

```
frmtmb test-gr-by.R pass=71 fail=0 err=0 skip=0
frmtmb test-methods.R pass=64 fail=0 err=0 skip=0
frmtmb test-brms-shapes.R pass=72 fail=0 err=0 skip=0
frmtmb test-brms-shapes-punch.R pass=55 fail=0 err=0 skip=0
frmtmb test-brms-shapes-punch2.R pass=22 fail=0 err=0 skip=0
frmtmb test-brms-names.R pass=138 fail=0 err=0 skip=0
frmtmb test-methods-audit.R pass=58 fail=0 err=0 skip=0
frmtmb test-sugar.R pass=45 fail=0 err=0 skip=0
frmtmb test-mm.R pass=116 fail=0 err=0 skip=0
GATED frmtmb test-brms-suite-methods.R pass=162 fail=0 err=0 skip=0
GATED frmtmb test-brms-suite-standata.R pass=87 fail=0 err=0 skip=0
frmtmb.sample test-gr-by-draws.R pass=7 fail=0 err=0 skip=0
frmtmb.sample test-brms-output.R pass=21 fail=0 err=0 skip=0
frmtmb.sample test-brms-suite-helper-copy.R pass=0 fail=0 err=0 skip=0
GATED frmtmb.sample test-brms-suite-methods.R pass=61 fail=0 err=0 skip=0
frmtmb test-compat.R pass=579 fail=0 err=0 skip=0
frmtmb test-predict-newdata.R pass=12 fail=0 err=0 skip=0
frmtmb test-re-formula-partial.R pass=39 fail=0 err=0 skip=0
frmtmb test-id-kron.R pass=60 fail=0 err=0 skip=0
frmtmb test-tre.R pass=102 fail=0 err=0 skip=0
frmtmb test-covstruct.R pass=7 fail=0 err=0 skip=0
frmtmb test-frame.R pass=22 fail=0 err=0 skip=0
frmtmb test-simulate-ergonomics.R pass=50 fail=0 err=0 skip=0
frmtmb test-naming-collisions.R pass=29 fail=0 err=0 skip=0
frmtmb test-priors-bounds-grcov.R pass=50 fail=0 err=0 skip=0
frmtmb test-aliased-grouping.R pass=33 fail=0 err=0 skip=0
frmtmb test-autoscale.R pass=46 fail=0 err=0 skip=0
frmtmb test-dry-run.R pass=23 fail=0 err=0 skip=0
frmtmb test-simulate-newdata.R pass=44 fail=0 err=0 skip=0
frmtmb test-predict-re-uncertainty.R pass=65 fail=0 err=0 skip=0
frmtmb test-glmm-gaussian.R pass=17 fail=0 err=0 skip=0
GATED frmtmb.sample test-brms-suite-helper-copy.R pass=33 fail=0 err=0 skip=0
```

## brms-suite ports

Regenerated with `/tmp/lanes/shared/brmsport-gen.R` after editing
`dev/brmsport-verdicts.tsv` (and the matching rows of
`dev/brmsport-verdicts-manual.tsv`); the generated diff touches only
`test-brms-suite-standata.R`, rows standata:705, 706, 709, 710, 712.

- standata:712, `expect_error(standata(y ~ x + (1|gr(g, by = z3)), dat),
  "Some levels of 'g' correspond to multiple levels of 'z3'")`: now
  pass.
- standata:705, 706, 709, 710 (`Nby_1`, `Jby_1`): stay "cannot
  transfer", now class `stan` rather than `absent`, since the feature
  exists and only the Stan data does not. `test-gr-by.R` asserts the
  same level-to-by-level map, and brms's `rm_wsp()` naming (`zx1` for
  level `x 1`), on brms's own data for these rows.
- Fixture 5 (`brmsfit_example5`) now uses brms's own formula,
  `count ~ Age + (1 | gr(patient, by = gender)), mu2 ~ Age`, in
  `tests/testthat/helper-brms-suite.R` and its generated copy. Gated
  `test-brms-suite-methods.R` gave pass=162 fail=0 before and after,
  and the sample package's gave pass=61 fail=0 after, so no verdict
  moved.
- `dev/brmsport-ledger.tsv` rows 461 to 465 still record the old run
  results; they are refreshed by the ledger's record run, not by hand.

## Defects found, not fixed

- An observation-level by-split term escapes the OLRE check
  (`check_re_structure()` compares a block's level count with the row
  count, and each by-level block has fewer levels than rows).
- `autoscale_plan_z()` skips by-split blocks, because a sub-block's Z
  column matches its X column only on its own by-level's rows. The fit
  is the same with and without autoscale (checked on one design);
  the by-split slope variance is simply not rescaled.
- The new-level refusal of a multi-membership by-split term names the
  unseen levels of the first by-level block it meets, not all of them.
- `variables()` lists each block's sd and cor together (sd, sd, cor per
  by-level); brms orders by class (every sd, then every cor). This is
  frmtmb's order for any model with several blocks, not specific to by.
- The fit's `ranef()` does not carry brms's `attr(levels, "by")` on its
  row names; the draws' `ranef()` does, since it is compared with
  brms's with `identical()`.
- In a group that mixes by-split and plain terms, brms's `VarCorr()`
  pastes the by-level onto the plain term's coefficient too
  (`get_rnames()` reads `by[1]` for the whole group), which names
  parameters that do not exist; frmtmb names each block's own.

## NEWS entry

* `gr(g, by = f)` fits brms's by-split group-level term: one set of
  standard deviations and correlations per level of `f`, with each
  level of `g` in exactly one level of `f`. It is fitted as one block
  per by-level, with any structure whose density is a product over the
  grouping levels (`us`, `diag`, `cs`, `ar1`, `toep`, the spatial
  ones, and `dist = "student"`). Parameters take brms's names
  (`sd_g__Intercept:fa`, `cor_g__Intercept:fa__x:fa`); `ranef()`,
  `coef()` and `ngrps()` keep one entry over every level of `g`; a new
  level in `predict(allow_new_levels = TRUE)` takes the covariance of
  its own by-level. brms's refusals are reproduced verbatim.
  `mm(g1, g2, by = cbind(f1, f2))` works the same way over the pooled
  levels. `gr(g, by = f, cov = A)` is refused by name: brms correlates
  the by-levels through `A`, which is not a by-split.
* `summary()` lists a grouping factor's standard deviations before its
  correlations when the factor has several blocks, as brms does.
