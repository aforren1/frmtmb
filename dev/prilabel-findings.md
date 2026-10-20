# Lane prilabel: group-level prior and naming parity with brms 2.23.0

Round opened 2026-10-14, base 47cf659d (frmtmb 0.71.0), branch
`wt-prilabel`, library `C:/Users/adf44/source/r/wt-prilabel-lib`; the
"before" arm is `C:/Users/adf44/source/r/rellib-r9`. Scope: the items
lane dblbar filed at 0.71.0 (`dev/test-backlog.md`, "Filed by lane
dblbar"; `dev/dblbar-findings.md`, "Not done"). Nothing is committed.

Every brms claim below was RUN on brms 2.23.0, mostly without a compile
(`default_prior()`, `stancode()`, `standata()`, `brms:::frame_re()`),
and once with a tiny fit for the names a fit carries:

- `dev/prilabel-brms-ref.R` -> `dev/prilabel-log/brms-ref.txt`: sd and
  cor rows of eleven formulas, Stan code of per-coefficient sd priors,
  Student-t terms.
- `dev/prilabel-brms-ref2.R` -> `brms-ref2.txt`: mm() with weights and
  scale, three-level nesting, `||` with a slash, `|ID|` with a slash,
  `gr(g1/g2)`, the labels brms refuses, Student-t with `|ID|`.
- `dev/prilabel-brms-ref3.R` -> `brms-ref3.txt`: per-coefficient sd
  rows and Stan code on `gr(by = )`, `mo()`, `mmc()`, a nonlinear
  parameter, a block spanning mu and sigma, a multivariate model.
- `dev/prilabel-brms-fit.R` -> `brms-fit.txt` (1 chain, 300 iterations,
  seed 1, data `dev/dblbar-data.R` seed 20261013) and
  `dev/prilabel-brms-vars.R` -> `brms-vars.txt`: `variables()`,
  `ranef()`, `coef()`, `VarCorr()`, `summary()`, `prior_summary()`,
  `ngrps()`, `re_formula` and new levels on
  `bf(y ~ x + (1 + x | g1/g2) + (1 | mm(m1, m2)), sigma ~ (1 | g))`.

The frmtmb side of the same calls is `dev/prilabel-frm-probe.R`
(`frm-probe-base.txt`, `frm-probe-lane.txt`), and the row-by-row
comparison is `dev/prilabel-sdcmp.R` (`sdcmp-base.txt`,
`sdcmp-lane.txt`).

## Item 1: a class "sd" prior per coefficient

### What brms does (`brms-ref.txt`, `brms-ref3.txt`)

- `default_prior()` lists, per prefix (resp, dpar, nlpar), the
  class-wide sd row, then per group a group row and one row per
  coefficient, sorted by `order(resp, dpar, nlpar, class, group, coef)`
  (`.default_prior()`, `R/priors.R` line 580 of the tarball). So
  `(f | g)` lists `fb`, `fc`, `Intercept` in the en_US collation, and
  `Intercept` first in the C collation testthat uses.
- `set_prior(class = "sd", group = "g", coef = "x")` writes
  `normal_lpdf(sd_1[2] | ...)` and keeps the default on `sd_1[1]`. The
  same holds on `||`, `gr(cor = FALSE)`, a factor (`fb`), `mo()`
  (`moxo`), `mmc()` (`mmcx1x2`), `gr(by = )` (one row for every
  by-level), a nonlinear parameter, a dpar block and a multivariate
  response. A group prior beside a coef prior: the coef wins on its sd.
- A coef with no dpar, on a block that spans mu and sigma
  (`(1 + x | p | g)`, `sigma ~ (1 | p | g)`), reaches mu's column alone
  (`sd_1[1]`). A group-level spec without dpar reaches the whole block.
- brms refuses a coef without group ("sd_x does not correspond to any
  model parameter") and a coef the group does not have (`sd_g__z`). It
  refuses a bound on any coef ("Prior argument 'coef' may not be
  specified when using boundaries"). frmtmb already accepts a bound on
  a class `"b"` coef (`dev/prilabel-coefbound.R`, rellib-r9: accepted,
  b_x = 0.6024 at lb = 0.6), and keeps that for class `"sd"`.

### What a per-coefficient sd prior means in frmtmb

A correlated block's `theta` is `c(log sd, correlation parameters)`
(`us`, `us_t`, `gr_cov`, `gr_prec`, `cs`, `toep`), not a Cholesky
factor, and `diag`/`diag_t` hold the log sds alone. One standard
deviation is therefore one parameter, and the density goes on
`exp(theta_k)` with the log-Jacobian, as class `"sd"` already did per
block. It is exact and does not change the correlations, which class
`"cor"` addresses; brms factors the same way (the sd vector and `L`).
In a Student-t block the parameter is the scale, which is what brms's
`sd_` is there too. A structure whose columns share one standard
deviation (`homdiag`, `ar1`, `ou`, `homcs`, `homtoep`, `exp`, `gau`,
`mat`; `dev/prilabel-sdidx.R` -> `sdidx.txt`) has none per coefficient:
a coef that names one of its columns is refused by name, and a group
prior still reaches the shared one. Smooth, gp, CAR, SPDE, `rr` and
`equalto` blocks list no coefficient rows.

### Changes

- `R/priors.R`: `sd_spec_reach()` narrows by `coef`, with an EQUAL
  prefix (brms's rule above); `block_sd_coefs()` gives brms's coef per
  standard deviation (`brms_rename()` of the column, as `get_rnames()`
  does); `sd_coef_shared()` and `sd_coef_hint()`; the resolver refuses
  a coef without group and a coef on a shared sd, and lists a group's
  coefficients when a coef matches none; `prior_table()` writes the
  coefficient rows in brms's order; the class `"sd"` branch of
  `unhonored_coef_refusal()` is removed, so the brms-table route and
  `set_prior()` accept the row; `resolve_gp_class()` reads the prefix
  only, because its coef names a sub-GP. Docs: a "Group-level standard
  deviations" section and `@param coef` in `?set_prior`, and
  `?default_prior`.
- `extensions/frmtmb.sample/R/sample.R`: a default steps aside only for
  a user specification of its own slot (`drop_superseded()`), not for
  every specification of its class. On rellib-r9 a group prior on `g`
  dropped the `student_t(3, 0, 2.5)` default from every other group,
  so `sd_h` was flat; brms keeps it (`dev/prilabel-sample-defaults.R`,
  data seed 3; theta 1, 2 = `sd_g__Intercept`, `sd_g__x`, 4 =
  `sd_h__Intercept`; `sample-defaults-base.txt`, `-lane.txt`):

```
                    base (rellib-r9)               lane
group g, coef x     ERROR (refused)                1=t  2=normal  4=t
group g             1=normal 2=normal  (4 flat)    1=normal 2=normal 4=t
class-wide          1=normal 2=normal 4=normal     same
```

  `priorlist_classes()`, `prior_class_key()` and `drop_overridden()`
  are removed. `?frm_sample` and the sampling vignette said a
  specification "takes over the classes it names ... brms's
  partial-override rule", which is not what brms does.

### Evidence

`dev/prilabel-sdcmp.R` compares, in one process, brms's and frmtmb's
sd rows (class, coef, group, resp, dpar, nlpar; `identical()`, in
order) on 20 formulas, and the set of standard deviations a
`normal(0, 0.5)` specification reaches on 25 specifications (brms's
from its Stan code through `frame_re()`, frmtmb's from
`resolve_priorlist()`). A SAME also needs a non-empty, error-free set.

```
lane (sdcmp-lane.txt):  ROWS SUMMARY 20 of 20 cases identical
                        PRIOR SUMMARY 25 of 25 specifications identical
base (sdcmp-base.txt):  ROWS SUMMARY 0 of 20 cases identical
                        PRIOR SUMMARY 0 of 25 specifications identical
```

The 25 cover every spelling in the brms list above and three
group-only specifications on the nested and mm groups (item 2).

## Item 2: nested and multi-membership group labels

### What brms does (`brms-ref.txt`, `brms-ref2.txt`, `brms-fit.txt`)

`(x | g1/g2)` is `(x | g1) + (x | g1:g2)`, outer first, and `g1/g2/g3`
adds `g1:g2:g3`; an explicit `(1 | g2:g1)` keeps `g2:g1`.
`mm(m1, m2, weights = , scale = )` is the group `mmm1m2`. brms refuses
`group = "g2:g1"` and `group = "mm(m1, m2)"` on those models, and
refuses `gr(g1/g2)` and `(1 | p | g1/g2)`. On the fit, `variables()`
lists `sd_g1__*`, `sd_g1:g2__*`, `sd_mmm1m2__*`, `sd_g__sigma_*`, then
every `cor_`; `ranef()`, `coef()`, `VarCorr()`, `summary()$random` and
`ngrps()` are keyed `g1`, `g1:g2`, `mmm1m2`, `g`; the levels of `g1:g2`
are `p1_q1`; the new-level error names `'g1:g2'` and `'p4_qNEW'`.

### Before (rellib-r9, `frm-probe-base.txt`)

`default_prior()` had the groups `g2:g1`, `g1`, `mm(m1, m2)`;
`variables()` listed `sd_g1:g2__*` before `sd_g1__*`;
`names(ranef())` and `ngrps()` were `g2:g1`, `g1`, `mm(m1, m2)`;
`VarCorr()` and `summary()` had `g1:g2` before `g1`; the levels were
`q1:p1`; `set_prior(group = "g1:g2")` and `group = "mmm1m2"` were
refused ("No random-effect SDs match"); the new-level error named
`g2:g1`.

### Changes

- `R/parse.R`: `expand_slash_groups()` rewrites every `(x | a/b)`
  (with `||`, an `|ID|` key or a covariance wrapper too) into brms's
  pieces before splitForm(), and after the twin-term check, because
  brms compares terms as written. The block's own `group_name` is then
  brms's, and every reader agrees. `slash_nested_bars()` and the
  `from_slash` flag are removed (`R/frame.R`, `R/brms-names.R`), with
  the two places that reversed reformulas's order. `mm_group_name()`
  names an mm() group by its members (`R/frame.R`, `R/gr-by.R`,
  `R/predict.R`).
- `R/methods-fit.R`: `ranef()` and `coef()` rows are brms's levels
  (`brms_levels()`, `p1_q1`), as frmtmb.sample's draws already named
  them.
- `R/confint.R` (`hyp_env_vals()`): every block's `cor_` follows every
  block's `sd_`, brms's `variables()` order (`brms-vars.txt`). This
  moved the `gr(g, by = f)` expectation of `test-gr-by.R`, whose
  interleaved order was frmtmb's, not brms's (`dev/grby-brmsnames.R`
  recorded names, not order).
- `R/predict.R`: the new-level message spells a level of `g1:g2` as
  brms does (`p4_qNEW`).

After (`frm-probe-lane.txt`): the `sd_`/`cor_` names of `variables()`
are brms's list above, in brms's order; `ranef()`, `VarCorr()`,
`ngrps()` and `summary()$random` are keyed `g1`, `g1:g2`, `mmm1m2`,
`g`; the levels are `p1_q1`; both brms group labels are accepted.

## Item 3: one Student-t latent scale per grouping factor

### What brms does (`brms-ref.txt`, `brms-ref2.txt`)

`(1 | gr(g, dist = "student")) + (0 + x | gr(g, dist = "student"))`
declares one `df_1` and `udf_1` and writes `r_1_1 = dfm_1 .* (...)` and
`r_2_1 = dfm_1 .* (...)`. A term in `sigma` on `g` reads `dfm_1` too, a
term on another factor gets `dfm_2`, and a gaussian term on `g` gets
none. `(1 + x | gr(g, ...)) + (0 + z | gr(g, ...))` shares `dfm_1`
between a correlated term and a scalar one. So the Student-t terms of
one factor are one multivariate t with a block-diagonal scale.

### Changes

Implemented, because the joint density is closed form: the quadratic
forms and log-determinants of the blocks add, and the multivariate-t
normalizer takes the summed dimension.

- `R/covstruct.R`: `student_share_sets()` marks the Student-t blocks
  of one factor (`t_share`) and refuses two `dist_nu` (brms has one df
  per factor) and two level sets; `student_t_parts()` and
  `student_shared_lpdf()`; `?frmtmb-student-re` rewritten.
- `R/objective.R`: a shared set is one density, taken at its first
  block.
- `R/predict.R` (`draw_b()`): one mixing draw per level for the set.
- `R/simulate-newdata.R`: `simulate(re_formula = )` that keeps some
  terms of a set and drops others is refused, because a dropped term
  is not independent of a kept one; dropping all of them runs.
- `R/fit.R`: `quadrature = TRUE` refuses a shared set, because the rule
  integrates one random effect at a time.
- `R/frame.R` calls `student_share_sets()`; `R/compat.R`: two rule
  texts.

Not done: brms ESTIMATES `df` under `gamma(2, 0.1)` and frmtmb fixes
`nu` (documented, `?frmtmb-student-re`, "Why nu is fixed"); an
`|ID|`-keyed Student-t term is still refused.

### Evidence (`tests/testthat/test-student-share.R`, data seed 21)

- Two scalar Student-t terms on `g` and `(1 + x || gr(g, dist =
  "student"))` are the same model: the same logLik and scales on the
  lane. On rellib-r9 the two-term fit has logLik -504.78 against
  -504.46, and log scales -0.67954, -1.35228 against -0.67939,
  -1.35772 (`dev/prilabel-tests/st2base-files/`).
- The shared density of `(1 + x | gr(g, ..., dist_nu = 4)) +
  (0 + z | gr(g, ..., dist_nu = 4))` at the fit equals the mixture
  integral over one gamma(nu/2, nu/2) precision per level by
  `integrate()` (relative tolerance 1e-9), and the taped `obj$env$f()`
  at the joint mode equals the gaussian negative log-likelihood minus
  that reference (relative tolerance 1e-9).
- `draw_b()` reproduces `z1 * s1 * w` and `z2 * s2 * w` from the same
  seed: one `w` per level for both terms.

## Item 4: the low items

- **`gr(g, id = )`** is implemented (`bar_read_id()`, `R/parse.R`):
  `(x | gr(g, id = "a"))` is `(x | a | g)`, which is what brms builds
  (`brms-ref.txt`: one cor row for the two terms). `id = NA` is no id;
  `|b|` beside `gr(id = )` and `gr(id = )` under `||` are refused.
  `mm(id = )` stays refused (an mm() block cannot share an `|ID|` key).
- **Group-level `cs()`** is not done: brms's `(cs(1) | g)` is a
  per-threshold random effect of the sequential ordinal families,
  which is a new block layout over the thresholds, not a small change.
- **`re_formula = ~ diag(x | g)`** now names the term it wraps
  (`re_formula_bars()`, `R/re-formula.R`); on rellib-r9 it named
  nothing and predicted as `re_formula = NA`.
- **A `||` term in messages**: the class `"cor"` refusals name the
  term as written, `(1 + f || g)` and `(1 + f | gr(g, cor = FALSE))`,
  through `block_written_label()` (`R/priors.R`); the block components
  carry `written_text` (`R/frame.R`). `confint_varcorr()$block` and
  the other tables keep `term_label`, which is the matching key of
  `ranef()[["..."]]` and the summary rows; not changed.

## Tests that pin the changes, seen to fail on rellib-r9

Each new file was run on the lane library and on rellib-r9 with
`dev/release/run-tests.R` through `dev/prilabel-run-par.sh` (one file
per process; logs `dev/prilabel-tests/<tier>-files/`):

```
file                         lane            rellib-r9 (tier)
test-prior-sd-coef.R         31 pass, 0/0    1 pass, 3 fail, 5 err (base2)
test-group-labels.R          32 pass, 0/0    4 pass, 12 fail, 4 err (glbase)
test-student-share.R         14 pass, 0/0    3 pass, 9 fail, 1 err (st2base)
frmtmb.sample, one new test:
test-default-priors-brms.R   36 pass, 0/0    34 pass, 1 fail, 1 err (sbase)
```

These are the round-0 files. After fix round 1 the lane counts are
`test-prior-sd-coef.R` 37, `test-group-labels.R` 43,
`test-student-share.R` 28 and `test-default-priors-brms.R` 39, all
passing (`all4`); the new tests' failures before their fixes are in
"Fix round 1".

The rellib-r9 failures are behavioral, not missing symbols, except one
error in test-student-share.R (`student_shared_lpdf` is new): the sd
rows differ from brms's tables; a per-coefficient prior is refused;
`names(ranef())` is `g2:g1`, `g1`, `mm(m1, m2)`; `group = "g1:g2"` is
refused; `gr(g, id = )` is refused; `re_formula = ~ diag(1 + x | g1)`
predicts as `NA`; the cor refusal says `1 + f | g [diag]`; the two-term
Student-t fit has a different logLik from the one-block fit; `draw_b()`
draws two mixing variables; the default on `sd_h` is dropped under a
prior on `g` (frmtmb.sample, a FAIL, ahead of the coef case's ERROR).

The absent cases: a coef on `diag()` and `cs()` blocks (each column its
own sd) is accepted and lands on that sd; a group prior on a
`homdiag()` block still reaches its one sd; `(1 | g2:g1)` keeps its
label; Student-t terms on two factors, and a gaussian term beside a
Student-t one, share nothing; one Student-t block still takes
`quadrature = TRUE` (`test-tre.R`, pass=102 fail=0); `re_formula = NA`
on a shared set runs; `id = NA` is no id.

Changed expectations, each because brms's answer differs from what
the old test pinned: `test-brms-formula-priors.R` and
`test-prior-sd-scope.R` (a group row now has its coefficient rows under
it, `(vectorized)`), `test-prior-compat.R` (the brms per-coefficient
sd row translates instead of being refused), `test-mm.R` (group
`mmg1g2`), `test-se-tier.R` (the zero sd of `(1 | Subject/a)` is
`theta_2`, Subject:a now following Subject; the test reads it from the
block), `test-gr-by.R` (sds before cors).

## Fix round 1 (review `dev/reviews/2026-10-14-prilabel.md`)

Each fix has a test that fails with the fix taken back out.
`dev/prilabel-prefix-run.R <fix>` puts the round-0 function back into
the installed lane namespace and runs one test file
(`dev/prilabel-tests/prefix/fix<k>-*.txt`):

```
fix  file (fix taken out)         result               test that fails
1    test-student-share.R         13 pass 1 fail 2 err by-term (2 tests)
2    test-prior-sd-coef.R         33 pass 1 err        own prefix only
3    test-group-labels.R          40 pass 1 fail       sd, cor, then sds
4    test-group-labels.R          37 pass 4 fail       brms's group order
5    test-group-labels.R          42 pass 1 fail       c:d_f levels
6    test-group-labels.R          40 pass 1 fail       re_formula group
7    test-default-priors-brms.R   38 pass 1 fail       class-wide note
8    test-prior-sd-coef.R         36 pass 1 fail       validate_prior()
```

(Fix 1's third failing test is the old dist_nu message, which the
round-0 mutant does not have.) With every fix in, the same files pass
(tier `r1`, then `all4`).

1. **HIGH, Student-t `gr(g, by = f)`.** Fixed. `student_share_sets()`
   now shares the latent scale per LEVEL of the factor, over the union
   of the set's levels (`t_level`, `t_n` on each block). A set whose
   blocks hold disjoint levels (a by-term alone) shares nothing and is
   not marked, so it is the 0.71.0 model, RNG stream included.
   `student_shared_lpdf()` sums each level's quadratic forms into the
   union and takes `D_j` per level; `draw_b()` draws one mixing
   variable per union level; the partial-`re_formula` refusal now
   fires only where a dropped and a kept block share a level, and
   `quadrature = TRUE` refuses only marked (overlapping) sets.
   Evidence (`test-student-share.R`, data seed 21 of
   `dev/prilabel-rev-tby.R`): the two by-term fits give logLik
   -361.6239067 and -364.5278394, the rellib-r9 values, to a relative
   1e-8, and the scalar by-term takes `quadrature = TRUE` again. For
   `(1 | gr(g, by = hg, dist = "student")) + (0 + z | gr(g, dist =
   "student"))` (all three blocks one set, as brms writes `dfm_1` into
   both terms over every level of `g`, `dev/prilabel-rev-log/misc.txt`),
   `obj$env$f()` equals a textbook bivariate-t reference at the mode
   and at 3 points perturbed by N(0, 0.3^2) (seed 7) to a relative
   1e-10, and `draw_b()` reproduces one mixing variable per level for
   the three blocks. `R/compat.R`'s `gr_by` x `student_blocks` row
   says so and is now covered by these tests.
2. **MEDIUM, `sd_coef_shared()` ignored the prefix.** Fixed: it refuses
   only a column whose (resp, dpar, nlpar) EQUALS the specification's
   (`block_col_prefix()`, the coef branch's rule of `sd_spec_reach()`).
   Tested both of the reviewer's formulas: the prior lands on theta 2
   (mu's `x`) and theta 3 (sigma's `x`), and the shared block's own
   prefix is refused.
3. **`variables()` order with a smooth.** Fixed: a smooth's `sds_` is
   written after the `cor_` flush (`sd`, `cor`, `sds`, then `sdgp`,
   `lscale`, as before). Tested on `y ~ s(x) + (1 + x | g)`.
4. **Group order.** Done, not filed: `re_brms_order()` (`R/parse.R`)
   sorts each predictor's bar terms by brms's group name with
   `order()`, as `terms_re()` does, keeping the written order within a
   group, so the blocks, `theta`, `variables()`, `ranef()`,
   `VarCorr()`, `ngrps()` and `summary()` follow. brms 2.23.0's
   `frame_re()` on three unsorted formulas (`dev/prilabel-brms-order.R`
   -> `brms-order.txt`): `g:h` before `gs`; per predictor `g`, `k` in
   mu and `h`, `k` in sigma; `(1 | gs) + (0 + x | g) + (1 | g)` gives
   `g` (x), `g` (Intercept), `gs`. Two existing expectations moved:
   `test-brms-names.R` (`g:h` before `gs`) and `test-mm.R` (`fc` before
   `mmg1g2`).
5. **Levels whose parts hold `:`.** Fixed: the frame keeps brms's name
   of each level of an interaction group (`brms_level_map()`, from the
   parts pasted with `_`, `group_values_brms()`), and `ranef()`,
   `coef()` and the new-level message read it: `c:d_f` (brms), was
   `c_d_f`. Without a map the old reading stays. The fix-5 mutant
   found a defect in that fallback (a length-1 `NA` recycled into a
   longer level vector); fixed and pinned by a unit test.
6. **`re_formula` group match.** Fixed: `re_group_key()` compares the
   group as written. `~ (1 + x | g2:g1)` on a `g1/g2` fit is refused,
   naming the fit's terms (brms keeps no term there and silently
   predicts as `NA`).
7. **frmtmb.sample, class-wide over a per-coefficient default.** The
   behavior is brms's and stays. NEWS has a second BREAKING bullet for
   it, and `frm_sample()`'s default-prior MESSAGE now says, per user
   specification, when it applies to no parameter because a more
   specific default stays there, and names the coef to write
   (`unreached_prior_notes()`). The message, not only the docs,
   because the docs are not where a user finds out that a prior did
   nothing. `?frm_sample` says it too. Tested on the reviewer's `gp(x)`
   model (seed 51), with the absent cases `coef = "gpx"` and a group
   sd prior.
8. **Docs and order.** The `?set_prior` section and the
   `?frmtmb-student-re` paragraph are rewritten in short sentences.
   `validate_prior()` and `prior_summary()` list the sd rows in brms's
   `order(class, group, resp, dpar, nlpar, coef)` (`prior_rows_sd_order()`),
   tested on the multivariate model of fit B. The other classes keep
   `default_prior()`'s order there (filed below).

## Suites (final build, after fix round 1)

One file per process, NOT_CRAN=true, the lane library first, every
extension loading the lane core (all 385 `lib:` lines name
`wt-prilabel-lib`; tier `all4`, `dev/prilabel-tests/all4.log`,
generated by awk over the RESULT lines):

```
frmtmb files=239 pass=17552 fail=0 err=0 skip=185 warn=0
frmtmb.coupling files=11 pass=542 fail=0 err=0 skip=5 warn=0
frmtmb.eam files=30 pass=1808 fail=0 err=0 skip=3 warn=0
frmtmb.latent files=10 pass=360 fail=0 err=0 skip=2 warn=0
frmtmb.learn files=15 pass=430 fail=0 err=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample files=53 pass=2704 fail=0 err=0 skip=4 warn=0
frmtmb.spline files=16 pass=594 fail=0 err=0 skip=1 warn=0
```

(Round 0's tier `all2` had frmtmb 17521 and frmtmb.sample 2701, all
passing.)

Gated (NOT_CRAN=true, FRMTMB_BRMS_FIT_TESTS=true; the 30 files that
read the gate; tier `gatedlane3` on the final build, `gatedbase` on
rellib-r9):

```
lane      frmtmb files=23 pass=2639 fail=1 err=0 skip=13
          frmtmb.sample files=7 pass=556 fail=0 err=0 skip=0
rellib-r9 test-brms-suite-priors.R pass=36 fail=0; every other file
          the same counts as the lane
```

The one failure is the ledger flip below. The 13 skips are
`test-drmtmb-agreement.R`, which reads its own gate; with
FRMTMB_DRMTMB_FIT_TESTS=true it was pass=131 fail=0 skip=0 in round 0
(`fix3.log`).

`dev/prilabel-sdcmp.R` on the final build: ROWS 20 of 20, PRIOR 25 of
25 identical to brms (`sdcmp-lane.txt`).

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
rebuilt for fix round 1, `dev/prilabel-tests/openblas2-files/`):
`test-student-share.R` 28/0, `test-prior-sd-coef.R` 37/0,
`test-group-labels.R` 43/0. Round 0 (`openblas-files/`) also had
`test-se-tier.R` 91/0 and `test-prior-sd-scope.R` 25/0.

One session, round 0 only (`dev/prilabel-onesession.R`,
`test_dir(load_package = "none")`, `dev/prilabel-log/onesession-*.txt`):
frmtmb files=229 pass=17514 fail=9 err=1, of which 5 are the
`test-compat.R` and `test-compat-register.R` failures known on 0.71.0
(lane testiso's), and one was this run's artifact (a test file edited
while the session ran); frmtmb.sample files=48 pass=2701 fail=0.

`R CMD check --as-cran` on the final build (`dev/prilabel-check.sh`,
built and checked in `dev/prilabel-check/<pkg>/`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, R_LIBS lane, rellib-r9, user
library):

```
frmtmb         Status: 1 NOTE   (HTML manual: V8 unavailable for math
                                 rendering, the expected one)
frmtmb.sample  Status: OK
```

## The ported brms-suite ledger

One row flips (do not edit `dev/brmsport-verdicts.tsv`; regenerate at
consolidation): `brms priors:74`, "default_prior returns correct prior
names for auxiliary parameters", recorded as a divergence because
frmtmb's Beta() table lacked the per-coefficient sd row
(`sd, g, Intercept`), now HOLDS (`gatedlane-files`,
`frmtmb--test-brms-suite-priors.txt`). On rellib-r9 the file is 36/0.

## Not done, and defects found but not fixed

- **brms estimates the Student-t df; frmtmb fixes nu.** Documented
  (`?frmtmb-student-re`). An `|ID|`-keyed Student-t term stays refused;
  brms shares one `dfm` there too (`brms-ref2.txt`).
- **Group-level `cs()`** (`(cs(1) | g)`): a per-threshold random effect,
  not small. Filed.
- **`mm(id = )`** stays refused: an mm() block cannot share an `|ID|`
  key in frmtmb.
- **Population-level order of `variables()`** (pre-existing): brms
  lists `b_Intercept, b_sigma_Intercept, b_x`, frmtmb
  `b_Intercept, b_x, b_sigma_Intercept` (`brms-vars.txt`,
  `frm-probe-lane.txt`).
- **`coef()` on a distributional model** (pre-existing) is keyed by
  dpar (`mu`, `sigma`) and then group; brms's is keyed by group.
- **Class `"cor"` rows of a dpar block** (pre-existing): frmtmb lists
  `cor, g, sigma` where brms lists one `cor, g` row with no dpar
  (`brms-ref.txt` against `frm-probe-lane.txt`); and brms sorts
  `cor` before `sd`, frmtmb lists `sd` first.
- **`confint_varcorr()$block` and other tables** still read
  `term_label` (`1 + f | g` for a `||` term); it is the matching key,
  so it is left.
- frmtmb.sample's draws still carry `theta_*`, not brms's `sd_` names,
  in `variables(draws)` (lane dblbar's item; not this lane's scope).
- **`ranef()` gives one entry per block** (pre-existing): two blocks on
  one factor give the names `g g`, and three shared Student-t terms on
  `g` give `g g g`; brms merges them into one `g` entry with every
  coefficient (`dev/prilabel-rev-log/order-brms.txt`).
- **`validate_prior()` and `prior_summary()` order of the classes
  other than `"sd"`** keeps `default_prior()`'s order; brms sorts
  every row by `order(class, group, resp, dpar, nlpar, coef)`. The sd
  rows are in brms's order since fix round 1.

## Files touched

Core: `R/parse.R`, `R/frame.R`, `R/priors.R`, `R/brms-names.R`,
`R/covstruct.R`, `R/objective.R`, `R/predict.R`, `R/fit.R`,
`R/simulate-newdata.R`, `R/re-formula.R`, `R/methods-fit.R`,
`R/confint.R`, `R/gr-by.R`, `R/compat.R` (two rule texts), `NEWS.md`;
`man/set_prior.Rd`, `man/default_prior.Rd`, `man/frmtmb-student-re.Rd`,
`man/frmtmb-multimembership.Rd` (roxygen). Tests: new
`test-prior-sd-coef.R`, `test-group-labels.R`, `test-student-share.R`;
changed `test-brms-formula-priors.R`, `test-prior-sd-scope.R`,
`test-prior-compat.R`, `test-mm.R`, `test-se-tier.R`, `test-gr-by.R`,
`test-brms-names.R` (fix round 1).

frmtmb.sample: `R/sample.R`, `man/frm_sample.Rd`,
`vignettes/sampling.Rmd`, `NEWS.md`,
`tests/testthat/test-default-priors-brms.R`.

Shared with other lanes: `R/priors.R` (lane sampleprior adds
`uniform`/`lognormal` to `parse_prior_dist()` and
`prior_base_logdens()`, no overlap with these hunks), `R/parse.R`,
`R/frame.R`, `R/compat.R`, `R/predict.R`, `R/fit.R`, `NEWS.md`,
`extensions/frmtmb.sample/R/sample.R` (lane sampleprior changes
`default_priors_for()` and `sample_resolve_priors()`'s tail; this lane
changes the `drop_overridden()` call two lines above that hunk and
removes `drop_overridden()`, `priorlist_classes()` and
`prior_class_key()`), `extensions/frmtmb.sample/NEWS.md`.

Dev scripts (tracked): `dev/prilabel-*.R`, `dev/prilabel-install.sh`,
`dev/prilabel-run-par.sh`, `dev/prilabel-check.sh`, and for fix round 1
`dev/prilabel-prefix-run.R` and `dev/prilabel-brms-order.R`. Logs are local
(`dev/prilabel-log/`, `dev/prilabel-tests/`, `dev/prilabel-check/`;
the OpenBLAS copy of R in `dev/ciharden-out/` was deleted after use).

## Version

Minor bumps: frmtmb (breaking labels, a new model for shared Student-t
terms) and frmtmb.sample (breaking default rule). frmtmb.sample's new
test writes a per-coefficient sd prior, so it needs the new frmtmb as
its floor.
