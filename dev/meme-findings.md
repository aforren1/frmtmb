# Lane meme: brms's `me()` prior classes `meanme`, `sdme`, `corme`

Round opened 2026-10-28. Worktree `C:\Users\adf44\source\r\frmtmb-wt-meme`,
branch `wt-meme`, base `ac4906bd` (frmtmb 0.77.0). Library
`C:/Users/adf44/source/r/wt-meme-lib`; "before" is
`C:/Users/adf44/source/r/rellib-r15`. brms 2.23.0, rstan 2.32.7.
Logs are in `dev/meme-out/` (gitignored); every number below is pasted
from them.

## 1. Measurement first

### 1.1 What frmtmb did with `me()` at 0.77.0

Script `dev/meme-before.R base`, log `dev/meme-out/before-base.txt`
(seed 11, n = 40).

- `me()` itself was complete and is brms's model (lane me, 0.64.0,
  `dev/me-findings.md`): `Xn ~ normal(Xme, noise)`, `Xme ~
  normal(meanme, sdme)`, multivariate normal with correlation `corme`
  across the terms that share a `gr =` unless `set_mecor(FALSE)`.
- Parameterization (`R/me.R`, `R/frame.R`): the latent values `Xme`
  themselves are the inner component `miss` (shared with `mi()`). The
  hyperparameters are three OUTER components: `meanme` (natural scale,
  named `meanme_mex`), `logsdme` (log of `sdme`) and `thetame` (the
  `us_chol_cor()` chart of each correlated group, unnamed:
  `thetame_1`, ...). The latent density is centered, on `Xme`, not on
  brms's standardized `zme`.
- On `frm()` the latent values are integrated out by the Laplace
  approximation beside the random effects (`random = c("b", "miss")`).
  It is exact for a gaussian response with linear `me()` terms (the
  marginal is multivariate normal), approximate otherwise
  (`?frmtmb-me`, Estimation; `dev/me-findings.md` 3.1).
- On `frm_sample()` they are sampled with every other parameter
  (`laplace = FALSE`), as brms samples `zme`. The sampler works on
  `Xme` directly (centered), brms on `zme` (non-centered); the target
  is the same posterior up to the linear map's Jacobian (section 3.1).
- Names: `bsp_mexsx`, `meanme_mex`, `sdme_mex`, `corme__mex__mez`,
  `corme_g__mexg__mezg`, `Xme_mex[i]` in `variables()` of a fit and in
  `summary()$me`. In draws: `meanme_mex`, `logsdme_mex`, `thetame_1`
  (internal, kept as for `theta_1`), derived `sdme_mex` and
  `corme__mex__mez`. DEFECT found: the grouped `corme_g__...` column was
  missing from the draws (`draws_cov_cols()` matched `corme__` only);
  `dev/meme-out/agree-lane.txt` first run, "missing columns:
  corme_g__mexg__mezg". Fixed (section 2).
- Prior classes: `set_prior(class = "meanme")`, `"sdme"`, `"corme"`
  stopped with "class = "meanme" has no faithful frmtmb spelling. ...
  has no prior slot for it". `default_prior()` listed no `meanme`,
  `sdme` or `corme` row on either route. A brms table's flat rows and
  its `lkj(1)` corme row were dropped; an edited one was refused.
- `frm_sample()` defaults: none on the hyperparameters, so flat on
  `logsdme` and on the `thetame` chart (`prior_summary(draws)` lists
  Intercept and sigma only). `sample_prior = "yes"` drew none of them;
  `"only"` refused every `me()` model naming `'meanme_1'` and "frmtmb
  has no prior slot yet".

### 1.2 What brms 2.23.0 does

Scripts `dev/meme-brms.R` (no compile), `dev/meme-brms-order.R`,
`dev/meme-brms-names.R` (one small fit each); logs
`dev/meme-out/brms.txt`, `brms-order.txt`, `brms-names.txt`.

- Stan (`stan_Xme()`): per GROUP `i` of terms sharing a `gr =`
  (numbered in order of first appearance), `vector[Mme_i] meanme_i`,
  `vector<lower=0>[Mme_i] sdme_i`, `zme` (vector per term, or a
  `Mme_i x N` matrix when correlated) and
  `cholesky_factor_corr[Mme_i] Lme_i`;
  `Xme = meanme + sdme * L * zme`, `target += std_normal_lpdf(zme)`,
  `normal_lpdf(Xn | Xme, noise)`; generated `corme_i` from `Corme_i`.
- `default_prior()` rows (`prior_Xme()`): `meanme` class row (flat),
  one `meanme` row per coef (`mex`, `mexg`, ...), `sdme` class row
  flat with `lb = 0`, one `sdme` row per coef, and for each correlated
  group `lkj(1)` on `corme` and, under `gr = g`, a `corme group = g`
  row. Cases A to H in `brms.txt`, e.g. case G (`me(x) + me(z) +
  me(xg, gr = g)`): `corme`, `meanme`, `meanme mex`, `meanme mexg`,
  `meanme mez`, `sdme lb 0`, `sdme mex`, `sdme mexg`, `sdme mez`.
- Per coefficient: `meanme` and `sdme` take `coef`; brms writes one
  statement per term once any term has its own prior
  (`normal_lpdf(meanme_1[2] | 1, 0.5)`), a vectorized one otherwise.
  `corme` takes `group` and no `coef` ("Lme_mex ~ lkj(2)" does not
  correspond to any parameter) and only `lkj` ("The only supported
  prior for correlation matrices is the 'lkj' prior"). A bound on
  `sdme`/`meanme` is a declared bound (`vector<lower=0,upper=3>`), and
  `coef` with a bound is refused ("Prior argument 'coef' may not be
  specified when using boundaries.").
- `corme` arises when two or more `me()` terms share a `gr =` (`""`
  included) and `set_mecor(FALSE)` is not set: case B has one row,
  case C (`set_mecor(FALSE)`) none, case G one (only the ungrouped
  pair is correlated; the single grouped term is not).
- `sdx`: `me(x)` with no `sdx` passes `default_prior()` and
  `stancode()`, which do not read the data, and `standata()` refuses
  "Argument 'sdx' is missing in function 'me'."
  (`dev/meme-brms-nosdx.R`, `brms-nosdx.txt`). WITHDRAWN as parity
  (review finding 4): frmtmb refused at the parse, so its
  `default_prior()` refused where brms answers. Fixed in fix round 1
  (section 7): the refusal is now at the fit. A row with `NA` in `sdx` or `x` is
  dropped ("Rows containing NAs were excluded", N 39 of 40); frmtmb
  drops it too (29 of 30 rows in `brms-nosdx.txt`). `sdx = 0` is
  refused ("Measurement error should be positive.").
- Prior statements come after every other one (after `sd_1`,
  `rescor`): `brms-order.txt`, models F2 and MV.
- Prior-draw names (`brms-names.txt`): class-wide `prior_meanme_1`,
  `prior_sdme_1`; correlation `prior_Lme_1` (brms's renaming of `L_` to
  `cor_` does not reach `Lme`); a coefficient's own prior splits the
  vector into `prior_meanme_1__1`, `__2`, which `rename_prior()` turns
  into `prior_meanme_1_meanme_mex`, `prior_meanme_1_meanme_mez` (it
  replaces the trailing `_<k>` with the whole new name). With brms's
  defaults only `Lme_1` is drawn (fit f2: `Intercept sigma meanme_1
  sdme_1 Lme_1` with proper priors; defaults: `Intercept sigma
  Lme_1`).
- `"only"`: refused while a statement is flat, naming the Stan
  parameter: `'bsp'`, then `'meanme_1'`, then `'sdme_1'`, then `'Lme_1'`
  for a flat corme. brms does NOT refuse a vector with only some terms
  priored (`meanme coef = mez` alone compiles); frmtmb refuses that, as
  it already does for class `b` (documented in `?frm_sample`). Under
  `"only"` brms keeps `normal_lpdf(Xn | Xme, noise)` and the latent
  density outside `if (!prior_only)`.
- brms's `prior_summary()`/`validate_prior()` show the internal class
  `Lme` (`lkj_corr_cholesky(2)`) instead of `corme`, as they show `L`
  for `cor`. frmtmb keeps `corme` in its tables, as it keeps `cor`.

### 1.3 brms's documentation: `me()` is soft deprecated

`?brms::me` (dev/brms-suite/brms/R/formula-sp.R, lines 6 and 19): "(Soft
deprecated) ... `me` terms are soft deprecated in favor of the more
general and consistent `mi` terms." brms 2.23.0 gives no runtime
warning (`me()` at line 48 has none). Implication: brms keeps `me()`
working and documented, so parity is still owed to code written for
brms, but brms will not grow it. frmtmb should match it and go no
further (no `me()` features beyond brms's); new measurement-error work
belongs to `mi()`. `?frmtmb-me` already points to `mi()`.

### 1.4 Earlier decisions in `dev/`

- `dev/me-findings.md` section 6: "Priors on meanme, sdme, corme are
  refused, like `simo`: frmtmb has no slot for them ... it is the next
  step if wanted."
- `dev/test-backlog.md`, lines 1202-1204 (lane me) and 4349-4355
  ("Not done, with the reason", 0.77.0: "Decide whether to add them"),
  `dev/priordraws-findings.md` section 14 item 5 ("To file: brms's
  `meanme` and `sdme` as prior classes in core"), `dev/round-handoff.md`
  1375-1376, `dev/round-20261028.md` line 1096 (user question).
- `dev/core-boundary.md`: no entry for `me()`; it is in core since
  0.64.0 and the prior classes sit with the rest of the prior machinery
  in core (`R/priors.R`), so nothing moves.
- `SPEC.md` line 468 lists `me()` as shipped. `dev/feature-gaps.md`
  lines 68 and 423 predate `me()` (they say it was done as `mi(sdx)`);
  stale, not edited.

## 2. What changed

Core (`R/`):

- `R/priors.R`: `meanme`, `sdme`, `corme` added to
  `frmtmb_prior_classes` and `brms_direct_prior_classes` (own lines);
  their refusals removed from `brms_prior_class_refusal()`, and the
  translator's drop of brms's `lkj(1)` corme row removed.
  `set_prior_one()`: `corme` is an LKJ class (lkj only, no lb/ub);
  `resp`, `dpar`, `nlpar` refused on the three classes.
  `unhonored_coef_refusal()`: `group` refused on `meanme`/`sdme`,
  `coef` on `corme`. `prior_class_bounds()`: `sdme` has its own
  `lb = 0` (no unbounded-prior warning for `exponential()`).
  `prior_table()`: brms's rows (section 1.2). `resolve_priorlist()`:
  new `resolve_me_class()`: `meanme` on component `meanme` (identity),
  `sdme` on `logsdme` with scale `"sd"` (density on sdme plus the log
  Jacobian, as class `sd`), `corme` an `lkj_dist(kind = "chol")` on the
  group's `thetame` positions (the chart `rescor` uses), with `claim()`
  so a group's row replaces the class row; bounds go to `lower`/`upper`
  under the template names (`logsdme_mex` takes `log(lb)`). Refusals
  name what the model has (coefs, `gr =` variables, "set_mecor(FALSE)
  fixes it at zero", "this model has no me() term"). Docs of
  `set_prior()` updated.
- `R/me.R`: `?frmtmb-me` gains a Priors section; examples.
- `R/compat.R`: rule `me() x prior` from "conditional" to "works".

One prior slot per parameter: the hyperparameters live in their own
components (`meanme`, `logsdme`, `thetame`), which are not among
`theta_components`, so class `"theta"` cannot reach them and
`cov_class_coefs()` needs no rows for them: each parameter has exactly
one slot, its class (plus its coef or group row, which replaces the
class row by `prior_specificity_order()`). `prior = "flat"` is brms's
`""`: `sdme` stays flat on sdme, `corme` becomes flat on the Cholesky
factor (`lkj_dist(flat = TRUE)`). The route mark, `validate_prior()`
and `prior_summary()` read the same table and resolver, unchanged.

frmtmb.sample (`extensions/frmtmb.sample/R/`):

- `sample.R` `default_priors_for()`: brms's defaults, `(flat)` with
  `lb = 0` on `sdme` held as a flat density (the log Jacobian, as for
  `delta`), `lkj(1)` on `corme` when a group is correlated. `meanme`
  needs no entry. The `"only"` refusal of every `me()` model removed;
  the generic flat-statement check now names `'meanme_1'`/`'sdme_1'`
  as brms does. `announce_default_priors()` prints a flat default as
  `(flat), lb = 0` (it printed `flat()`, for `delta` too). Docs.
- `prior-draws.R`: `me_prior_stmts()` adds brms's statements after
  every other one, with brms's names (section 1.2).
- `draws-brms.R`: the derived-column pattern takes `corme_`, so the
  grouped `corme_<gr>__...` is in the draws (bug fix).

BREAKING: the sampled posterior of every `me()` model changes (flat on
`sdme` and `lkj(1)` on `corme` instead of flat on their internal
charts), and on `frm()` a brms table's `corme lkj(1)` row now applies.

## 3. Evidence

### 3.1 L1: brms's Stan `log_prob`, several points, value and gradient

`tests/testthat/helper-me-lp.R`, gated test "check C: row 23b" in
`tests/testthat/test-brms-likelihood.R`; numbers from `dev/meme-l1.R`,
log `dev/meme-out/l1.txt`, data seed 23, points seed 23, n = 40.
Every class has a proper prior on both sides (b normal(0,1), Intercept
and sigma student_t(3,0,2.5), meanme normal(0,2), sdme exponential(1),
plus in case 2 `meanme coef mezg normal(1,0.5)`, `sdme coef mexg
student_t(3,0,1)`, `corme group g lkj(2)`, in case 3 `sd
student_t(3,0,2.5)`), so `frm()`'s MAP objective is the posterior in
frmtmb's coordinates and brms's `log_prob(adjust_transform = TRUE)` the
posterior in Stan's. Identity checked:
`-f(theta) = log_prob(U(theta)) + log|det dU/dtheta| + C`, with `U`
built by `stan_pars_from_fit()` and its Jacobian by central differences
(independent of both densities); gradient by the chain rule,
`-grad f = J' grad log_prob + grad log|det J|`, the last term from a
closed form checked against the numeric log-determinant at every
point. Four points per case (the mode and three at sd 0.1 around it):

    case                     spread(C)/scale  max grad rel  |Lnum-Lana|
    one me()                 3.51e-12         7.59e-11      1.65e-11
    two me(gr = g), corme    6.37e-12         4.26e-10      5.25e-11
    me() + (1 | h)           4.09e-12         1.76e-10      2e-11

`C` is -0.693147180527 (case 1) and -1.38629436080 (cases 2 and 3):
brms's `- student_t_lccdf(0 | ...)` normalization, `log(2)` per
half-Student-t (sigma; sdme mexg; sd), which frmtmb leaves out as for
class `sd`. The test asserts `k = -C / log(2)` is an integer to 1e-8.

The MAP objective's prior term is also checked against a closed form
written from the densities (`test-me-priors.R`, "frm() adds brms's
densities"): normal, half-t on `exp(logsdme)` with `+ logsdme`,
exponential, and the LKJ(2) density of `rho = t / sqrt(1 + t^2)` with
its Jacobian, at the MAP and ML modes: relative error below 1e-10.

### 3.2 Prior tables equal to brms's `default_prior()`

`test-me-priors.R` ("the me() rows are brms's own, read from brms"):
five shapes (one, two, set_mecor(FALSE), two grouped, mixed), the
`class|coef|group|lb` of every me() row identical to brms 2.23.0's live
`default_prior()`. Sampling route: `test-me-prior-draws.R` ("the
sampling route lists brms's default rows"): `prior` strings of the me()
rows identical to brms's (`""` read as `(flat)`) on the mixed shape,
`lkj(1)` on corme, `lb = 0` on the sdme class row.

WITHDRAWN (review finding 1, a defect, not a convention): on the
sampling route the `corme group = g` row printed the inherited
`lkj(1)`, where brms prints `""`, as frmtmb's `cor` and `sd` group rows
did. Passed back, that printed value overrode an edited class row.
Fixed in fix round 1 (section 7): every group row prints `(flat)`.

### 3.3 Fit agreement, `frm_sample()` against brms

`dev/meme-agree.R` (`brms`, `lane`, `base` arms), comparison
`dev/meme-agree-cmp.R`, log `dev/meme-out/agree-cmp.txt`. Data seed
2026: 10 levels of `g` with 6 rows each; model `y ~ me(xg, sg, gr = g) +
me(zg, sg, gr = g) + me(x, sx)`; brms's default priors on both sides;
4 chains, 4000 iterations, 1000 warmup, seed 1 (brms: 10 divergent
transitions; lane 0; base 20). z = (frmtmb - brms) /
sqrt(MCSE^2 + MCSE^2) of the mean, sd, 5 and 95 percent quantiles:

    lane: max |z| 2.52 over 14 variables x 4 statistics (56)
          (meanme_mex sd -2.52, sdme_mexg sd 2.09; all others < 2.2)
    base: max |z| 13.63 over 13 variables x 4 (corme_g missing)
          sdme_mexg mean -12.64, sdme_mezg -13.63, sdme_mex -5.23,
          Xme_mexg[1] -9.85, b_Intercept -9.39, bsp_mexgsggrEQg 8.00

The base arm is the 0.77.0 sampler: flat on `log(sdme)` and on the
corme chart. 56 statistics at the null give a largest |z| near 2.5, so
2.52 is ordinary.

### 3.4 Guard absent: models without `me()`

`dev/meme-guard.R` (base and lane arms), logs
`dev/meme-out/guard-base.txt`, `guard-lane.txt`: four shapes without
`me()` (`(1 + x | g)` with b, sd and cor priors; a multivariate
rescor model; `ar(cov = TRUE)`; `y | mi(sy) ~ x`): both routes' prior
tables, every estimate and logLik at 17 digits, the objective and its
gradient at a shifted point, `prior_summary()`, and a seeded
`frm_sample(sample_prior = "yes")` run (digest of the draws rounded to
1e-10, the column names and the prior draw names). `diff`: identical,
281 lines each; draws digest `caaa455998804d753502b7284af73f67` in
both.

### 3.5 Seen to fail on rellib-r15

`dev/meme-run1.R <pkg> <file> base`, logs `dev/meme-out/base-*.txt`:

    test-me-priors.R       pass=7  fail=10 error=7   (lane 64/0/0)
    test-me.R              pass=67 fail=1  error=1   (lane 71/0/0)
    test-me-prior-draws.R  pass=2  fail=8  error=2   (lane 24/0/0)
    test-prior-draws.R     pass=82 fail=0  error=1   (lane 87/0/0)

Behavioral failures (not only missing symbols): the me() rows absent
from `default_prior()` on both routes; on the sampling route the
spread of `-log prior + log|J|` 1.17 (flat on the chart, not on the
correlation) and its constant off `log(pi^2/2)` by a ratio of 2.73;
`th$dist$flat` FALSE under `prior = "flat"`; the grouped `corme`
missing from the derived draws; `logLik(f1) == logLik(f0)` (brms's
corme row dropped); `"only"` refused with "frmtmb has no prior slot
yet".

## 4. Test runs

All with the lane library first (`lib:` lines name
`C:/Users/adf44/source/r/wt-meme-lib` for frmtmb and for every
extension, each installed from this worktree).

One file per process, `dev/release/run-par.sh`, job lists
`dev/meme-jobs-all.txt` and `dev/meme-jobs-gated.txt`, tallies by
`dev/meme-tally.sh`, logs `dev/release/meme-plain.log`,
`dev/release/meme-gated.log`:

    plain (NOT_CRAN=true): MEME-PLAIN ran 464 of 464
    frmtmb          files=286 pass=20431 fail=0 err=0 skip=206 warn=0
    frmtmb.coupling files=11  pass=551   fail=0 err=0 skip=5   warn=0
    frmtmb.eam      files=32  pass=1838  fail=0 err=0 skip=5   warn=0
    frmtmb.latent   files=12  pass=681   fail=0 err=0 skip=2   warn=0
    frmtmb.learn    files=15  pass=430   fail=0 err=0 skip=13  warn=0
    frmtmb.ode      files=11  pass=549   fail=0 err=0 skip=1   warn=0
    frmtmb.sample   files=74  pass=3571  fail=0 err=0 skip=15  warn=0
    frmtmb.spatial  files=8   pass=478   fail=0 err=0 skip=6   warn=0
    frmtmb.survival files=15  pass=552   fail=0 err=0 skip=0   warn=0

    gated (FRMTMB_BRMS_FIT_TESTS=true), the 49 files that skipped:
    MEME-GATED ran 49 of 49
    frmtmb files=33 pass=2787 fail=0 err=0 skip=15; frmtmb.sample
    files=7 pass=731 skip=1; eam 13/skip 3; learn 71/skip 2; spatial
    193/skip 0; coupling, latent, ode 0 pass (scale tier only).
    Every remaining skip is another gate: 15 FRMTMB_SCALE_TESTS, 13
    FRMTMB_DRMTMB_FIT_TESTS, 1 FRMTMB_FUZZ.
    test-brms-likelihood.R pass=494 (row 23b included),
    test-brms-priors.R pass=102.

The lane's files in the plain tier: test-me-priors.R 64, test-me.R 71,
test-me-prior-draws.R 24, test-prior-draws.R 87,
test-default-priors-brms.R 40, all fail 0.

One session, three orders, `dev/release/run-sessions.sh meme-sess 6`,
log `dev/release/meme-sess.log`: MEME-SESS ran 27 of 27, fail 0 and
err 0 in every session; frmtmb files=276/286 (the ten
`test-brms-suite-*.R` files skip at file level and leave no row, as in
earlier records, e.g. 269/279 in `dev/acprior-findings.md`),
frmtmb.sample 69/74, every extension n/n. 46 LEAK lines, all tagged
"[with load of ...]"; no untagged LEAK line.

OpenBLAS 0.3.26 with LAPACK (`dev/ciharden-openblas.sh 0.3.26 lapack`,
logs `dev/meme-out/ob-*.txt`): test-me-priors.R 64, test-me.R 71,
test-me-prior-draws.R 24, test-prior-draws.R 87,
test-prior-only-objective.R 10, all fail 0. The gated L1 row was not
run under OpenBLAS.

`update()` composes (`dev/meme-update.R`, `dev/meme-out/update.txt`):
dropping `me(z, sz)` drops the `corme` prior with brms's message,
dropping both terms drops both priors; `anova()` refuses a fit whose
`me()` prior adds a density, as for any prior.

R CMD check --as-cran (`dev/meme-check.sh`, `dev/meme-check/`):
frmtmb "Status: 1 NOTE" (the V8 math-rendering NOTE of the HTML
manual; tests 24 min, OK); frmtmb.sample "Status: 1 NOTE" (the same
V8 NOTE; tests 429 s, OK). The `.Rcheck` trees and tarballs were
deleted afterward; `check.txt` stays in each directory.

## 4b. Release needs

- frmtmb: minor bump (new prior classes, BREAKING table behavior).
- frmtmb.sample: minor bump; floor frmtmb at the new core version,
  because `default_priors_for()` writes `class = "sdme"` and
  `"corme"`, which core 0.77.0's `set_prior()` refuses by name.
- No other extension changes; none needs a new floor.

## 5. Not done, and why

- WITHDRAWN (review finding 2): `coef` with a bound was accepted on
  `meanme`/`sdme`, and its prior draws ignored the bound. Fixed in fix
  round 1 (section 7) for every class on `frm_sample()`.
- `"only"` with a vector priored on some terms only: brms samples it
  (flat on the rest), frmtmb refuses, as for class `b` (0.77.0 rule).
- brms's `Lme` internal class in its `validate_prior()` and
  `prior_summary()`: frmtmb shows `corme`, as it shows `cor` for `L`.
- Ledger: no row of `dev/brmsport-verdicts.tsv` flips. The only brms
  suite test of these classes is `tests.stancode.R:1483` ("noise-free
  terms appear in the Stan code"), which reads Stan code and is not in
  the ported ledger; `brmsfit-methods:423/425` stay "cannot transfer"
  for their own reasons (fixture, stanfit).
- `dev/feature-gaps.md` lines 68 and 423 are stale about `me()`; not
  edited (historical).

## 6. Overlaps with other lanes

`R/priors.R` (class vectors on their own lines, `prior_table()`,
`resolve_priorlist()`), `R/me.R`, `R/compat.R`, `NEWS.md`, `man/`;
frmtmb.sample `R/sample.R` (`default_priors_for()`,
`announce_default_priors()`, the `"only"` block, docs), `R/prior-draws.R`
(`prior_plan()` end, new `me_prior_stmts()`), `R/draws-brms.R`
(`draws_cov_cols()` pattern), `NEWS.md`, `man/frm_sample.Rd`. rescorpp
touches `R/sample.R` and the draws readers; the edits here are in
different functions.

## 7. Fix round 1 (review `dev/reviews/2026-10-28-meme.md`)

### 7.1 brms, measured again

`dev/meme-fix-brms-grp.R`, log `dev/meme-out/fix-brms-grp.txt`:

- brms 2.23.0 prints `""` on EVERY group row: `cor group g`, `corme
  group g`, `sd group g` and its coef rows, and a level of `thres(gr =
  th)` (`Intercept group a`, `delta group a lb 0`).
- `coef` with `lb`/`ub`: `set_prior()` accepts it, and `stancode()`
  (so `brm()`) refuses "Prior argument 'coef' may not be specified when
  using boundaries." on `b coef w ub 0.1`, `sd group g coef x ub 2`,
  `meanme coef mexg lb 0`, `sdme coef mexg ub 0.5` and `sdme coef mexg`
  bounds only. A class-level bound (`b lb -1`, `sdme ub 3`) and a
  group's (`sd group g ub 2`) are taken. The check is
  `brms/R/priors.R` line 1330, on any row with a coef and a bound.
- `me(x)` without `sdx`: `default_prior()` gives `b`, `b mex`,
  `Intercept`, `meanme`, `meanme mex`, `sdme lb 0`, `sdme mex`,
  `sigma`, and `validate_prior()` answers too.

### 7.2 What changed

- `R/priors.R` `prior_table()`: a group row (a `group`, no `coef`)
  takes no registered default, so it prints `(flat)`, brms's `""`.
  This is generic: `sd`, `cor`, `corme`, grouped `Intercept`
  thresholds and `delta` group rows. Passed back, such a row is
  `set_prior("", group = ...)`, which applies no density and so
  inherits the class row (the resolver keeps the class entry; no
  default targets a group, so `keep_default_bounds()` does not make it
  flat). BREAKING.
- `?set_prior` Hard bounds: frm() keeps a bound on one coefficient,
  `frm_sample()` refuses it.
- `R/me.R`, `R/frame.R`: the missing-`sdx` refusal moved from
  `parse_me_call()` to `me_build_frame(check_noise =)`, called with
  `check_response`, so the prior tables answer and `frm()` refuses with
  the same words.
- frmtmb.sample `R/sample.R`: `refuse_coef_bounds()`, called on the
  stacked user prior (the call's and the fit's own) in
  `sample_resolve_priors()`, refuses `coef` with `lb`/`ub` with brms's
  sentence, for every class. Not refused: class `"theta"` (frmtmb's
  own, no brms row), frmtmb's own `class = "Intercept", nlpar =`
  spelling (stored with a coef), and a row restating `lb = 0` on `sd`,
  `sdgp`, `lscale` or `sdme` (the gp() coef rows of a sampling-route
  table carry it). `?frm_sample` (`prior`) says both rules. BREAKING.

Decision to flag: the refusal is on `frm_sample()` only. `frm()` has no
`lower`/`upper` of its own and a coefficient's bound is its only box on
one coefficient; 14 test files across core and frmtmb.sample use it
(grep `coef = ..., lb =`). brms has no maximum-likelihood route, so it
offers no parity there. The wrong prior draws (finding 2) were on the
sampling route, which is now closed for every class.

### 7.3 Evidence, seen to fail on rellib-r15

Per-block runner `dev/meme-run1-blocks.R`, logs
`dev/meme-out/f2-*-{lane,base}.txt`:

    test-prior-table-groups.R  lane pass=23 fail=0  error=0
                               base pass=2  fail=6  error=3
    test-me-prior-draws.R      lane pass=40 fail=0  error=0
                               base pass=3  fail=13 error=4
    test-me-priors.R           lane pass=67 fail=0  error=0
                               base pass=7  fail=10 error=8

Behavioral failures on base among the new blocks:
- "without me(), the sd and cor group rows round-trip too": fail 3 (the
  group rows print `student_t`/`lkj(1)`; the resolved cor entry has
  eta 1, not 4; the sd entries are `t`, not `exponential`).
- "the group rows print brms's (flat)": fail 2 (`sd||g`, `cor||g`).
- "a bound on one coefficient is refused": fail 2 (b and sd accepted).
- "the announcement prints delta, sdme and corme": fail 3 (`delta`
  read `flat()`; no `sdme` or `corme` line).
The blocks that need the new classes on a `me()` model error on base
with the old class refusal (the weak form); the no-`me()` blocks above
are their behavioral counterparts.

Lane, the round trip of the review's case: `corme`, `cor` class rows
edited to `lkj(4)`, `sd` to `exponential(3)`: resolved `thetame` eta 4,
the block's two sds `exponential` and its correlation eta 4; prior
draws (2 chains of 500, seed 4) KS against Beta(4, 4) and
exponential(3) above 0.01, and against Beta(1, 1) below 0.01
(test-prior-table-groups.R). A class-level `sdme ub = 0.95`: posterior
and prior draws at most 0.95, KS against the truncated exponential(1)
above 0.01 (test-me-prior-draws.R).

### 7.4 Filed, not fixed

Repros in `dev/meme-filed.R`, log `dev/meme-out/filed.txt` (lane
library, seed 5):

1. `update(<draws>, prior = )` replaces the first call's prior instead
   of merging it, as brms's `update.brmsfit()` does
   (`dev/brms-suite/brms/R/update.R` lines 144-155). After
   `frm_sample(prior = meanme normal(0, 1) + sdme exponential(1))` then
   `update(prior = corme lkj(4))`, `prior_summary()` lists `sdme ""`
   and `corme lkj(4)`; `meanme normal(0, 1)` and `sdme exponential(1)`
   are gone.
2. `variables()` orders the `me()` hyperparameters by term, brms by
   group: on `me(xg, gr = g) + me(z) + me(wg, gr = g)` frmtmb gives
   `meanme_mexg meanme_mez meanme_mewg`, brms `meanme_mexg meanme_mewg
   meanme_mez`. Same set of names.
3. `default_prior()` gives `sigma` no `lb = 0` (`y ~ x`: `(flat) sigma
   lb NA`); brms lists `lb = 0`.
4. `mi()` latent values are `miss_1`, `miss_2` in draws; brms names
   them `Ymi_x[1]`, `Ymi_x[2]`.

Not filed (fixed): the `b` and `sd` coefficient bounds of finding 2.

### 7.5 Guard absent, after the fix

`dev/meme-guard.R lane`, log `dev/meme-out/guard-lane-fix1.txt`,
against `guard-base.txt`: the only differing lines are the four
sampling-route group rows (`cor group g`, `sd group g`, and `sd group g`
of each response of the multivariate model), now `(flat)`. Estimates,
logLik, objective and gradient, `prior_summary()`, and the seeded
`frm_sample(sample_prior = "yes")` draws digest
(`caaa455998804d753502b7284af73f67`) are identical.

### 7.6 Test runs after the fix

Whole suite with gates on (NOT_CRAN, FRMTMB_BRMS_FIT_TESTS), one file
per process, `dev/release/run-par.sh meme-fix1-gated`, log
`dev/release/meme-fix1-gated.log`: MEME-FIX1-GATED ran 464 of 464.

    frmtmb          files=286 pass=22886 fail=1 err=0 skip=15
    frmtmb.coupling files=11  pass=551   fail=0 err=0 skip=5
    frmtmb.eam      files=32  pass=1844  fail=0 err=0 skip=3
    frmtmb.latent   files=12  pass=681   fail=0 err=0 skip=2
    frmtmb.learn    files=15  pass=501   fail=0 err=0 skip=2
    frmtmb.ode      files=11  pass=549   fail=0 err=0 skip=1
    frmtmb.sample   files=74  pass=3937  fail=2 err=0 skip=1
    frmtmb.spatial  files=8   pass=561   fail=0 err=0 skip=0
    frmtmb.survival files=15  pass=552   fail=0 err=0 skip=0

Every skip is another gate (15 scale, 13 drmTMB, 1 fuzz). The three
failures were two tests that pinned the old group-row print
(`test-get-prior-route.R` "route = 'sample' reports the registered
defaults", `test-prior-route.R` sd and cor rows). Both now assert the
class row's default and the group row's `(flat)`; rerun
(`dev/meme-out/fix1-rerun-*.txt`): test-get-prior-route.R pass=82
fail=0, test-prior-route.R pass=11 fail=0; on rellib-r15
(`fix1-base-*.txt`) 81/1 and 9/2. New file test-prior-table-groups.R,
gated runner: pass=23 fail=0 (`fix1-gated-ptg.txt`).

One session, three orders, frmtmb and frmtmb.sample
(`dev/release/run-sessions.sh meme-fix1-sess 6 frmtmb frmtmb.sample`,
log `dev/release/meme-fix1-sess.log`): MEME-FIX1-SESS ran 6 of 6.

    frmtmb        forward files=276/286 pass=20435 fail=0 err=0 skip=206
    frmtmb        reverse files=276/286 pass=20585 fail=0 err=0 skip=206
    frmtmb        preload files=276/286 pass=20678 fail=0 err=0 skip=208
    frmtmb.sample forward files=70/75   pass=3612  fail=0 err=0 skip=15
    frmtmb.sample reverse files=70/75   pass=3612  fail=0 err=0 skip=15
    frmtmb.sample preload files=70/75   pass=3612  fail=0 err=0 skip=15

34 LEAK lines, all tagged "[with load of ...]"; no untagged one.

R CMD check --as-cran after the fix (`dev/meme-check.sh`):
frmtmb "Status: 1 NOTE" (V8 math rendering; tests 22 min, OK);
frmtmb.sample "Status: 1 NOTE" (the same; tests 425 s, OK). The
`.Rcheck` trees and tarballs were deleted; `check.txt` stays.

OpenBLAS was not rerun in this round: the fix touches the prior table,
a refusal and the announcement, no numerical code.

## 8. Fix round 2 (short check, S1 and S2)

### 8.1 S1: the nonlinear-intercept spelling

`refuse_coef_bounds(fit, pl)` (frmtmb.sample `R/sample.R`) now reads
the stacked user prior against the fit. A specification written
`set_prior(class = "Intercept", nlpar = "a")` (stored by core as class
`"b"`, `coef = "Intercept"`, with `$written` set) is returned as
brms's class row `class = "b", nlpar = "a"` when `a`'s predictor is an
intercept alone, density and bound included; a class row written both
ways keeps the later one. With a bound and a predictor that has more
than an intercept, it is refused in brms's words, as brms's own
`class = "b", coef = "Intercept", nlpar = "a", lb = 0` is. Without a
bound on such a predictor it stays a coefficient row, as before. Core
is unchanged in this round.

Test `extensions/frmtmb.sample/tests/testthat/test-prior-nlpar-intercept.R`
(data seed 2, sampler seed 3, 2 chains of 300 draws), per block
(`dev/meme-run1-blocks.R`, logs `dev/meme-out/f3-pni-{lane,base}.txt`):

    lane       pass=17 fail=0
    rellib-r15 pass=8  fail=9

Behavioral failures on base: the table's class row has `lb NA` and
`(flat)`; the prior draws reach -1.7 under `lb = 0` and -1.9 under
`lb = 0.4`; KS against the truncated normal(1, 1) p below 1e-17 for
both;
no prior draw is named `b_a`; the slope case is not refused (2). On
the lane: resolved bounds and densities identical to brms's spelling
(the guard-absent case), prior draws at least 0 and 0.4, posterior
draws at least 0 and 0.4, KS p above 0.01.

`frm()` with the same spelling, on `a ~ 1` and on `a ~ 1 + w`
(`dev/meme-fix2-frm.R`, `dev/meme-out/fix2-frm-{base,lane}.txt`):
logLik and every estimate printed to 15 digits are identical on
rellib-r15 and the lane.

### 8.2 S2: NEWS

frmtmb.sample's refusal bullet now says "for the brms classes" and
names the three cases not refused (class `theta`, a restated `lb = 0`,
the nonlinear-intercept spelling on `a ~ 1`); a new bullet records the
prior draws of that spelling. The group-row bullet names the
`thres(gr = )` `Intercept` and `delta` rows. The split
"Kolmogorov-Smirnov" in frmtmb.sample's NEWS is joined; core's NEWS
had it on one line already ("with a / Kolmogorov-Smirnov").

### 8.3 Runs after fix round 2

Gates on, one file per process (`dev/release/run-par.sh
meme-fix2-gated`, job list `dev/meme-jobs-fix2.txt`: all of
frmtmb.sample, and every core and extension file that names
`frm_sample` or frmtmb.sample), log `dev/release/meme-fix2-gated.log`:
MEME-FIX2-GATED ran 117 of 117.

    frmtmb          files=36 pass=5172 fail=0 err=0 skip=0
    frmtmb.eam      files=3  pass=136  fail=0 err=0 skip=0
    frmtmb.ode      files=1  pass=15   fail=0 err=0 skip=0
    frmtmb.sample   files=76 pass=3981 fail=0 err=0 skip=1 (scale tier)
    frmtmb.spatial  files=1  pass=25   fail=0 err=0 skip=0

test-prior-nlpar-intercept.R 17, test-me-prior-draws.R 40,
test-prior-table-groups.R 23, all fail 0.

frmtmb.sample in one session, three orders (`run-sessions.sh
meme-fix2-sess 3 frmtmb.sample`): forward, reverse and preload each
files=71/76 pass=3629 fail=0 err=0 skip=15. Every LEAK line is tagged
"[with load of ...]"; no untagged one.

R CMD check --as-cran of frmtmb.sample: "Status: 1 NOTE" (V8 math
rendering; tests OK). Core did not change in this round.
