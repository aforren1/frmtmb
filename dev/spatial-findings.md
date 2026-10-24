# Lane spatial: the covariance-structure seam, guardrail 1, frmtmb.spatial

Base 19b162b9 (frmtmb 0.73.1), worktree `frmtmb-wt-spatial`, branch
`wt-spatial`, library `C:/Users/adf44/source/r/wt-spatial-lib` in front
of `rellib-r11`. Nothing is committed. Logs are under
`dev/spatial-log/` (gitignored); every number below was pasted from
them. The parts stand alone: `dev/spatial-out/AB.diff` is the whole
change after parts A and B, before any of part C (54 files), and part C
is recorded in its own section with its own patches.

## Part A. The covariance-structure seam

### What changed

`frmtmb_register_covstruct()` (new, exported, `R/covstruct-seam.R`,
`?frmtmb_register_covstruct`) registers a covariance structure from a
package's `.onLoad()`. It follows `frmtmb_register_aterm()`: it checks
every argument before it commits anything, refuses a name core's own
structures or formula grammar hold (`covstruct_reserved_names()`),
replaces an entry the same caller registered before, and adds the
structure to the compatibility vocabulary as kind `"covstruct"` in the
same step, so a structure `frm()` accepts is always one `frm_compat()`
can describe. A refusal leaves the registry and the vocabulary as it
found them (tested).

The registry is now an environment, `frmtmb_covstruct_registry$reg`,
filled at build time from `covstruct_core` (the built-in entries) and
extended at run time. Every consumer reads it through `cs_entry(name)`,
`cs_entries()` and `cs_names()` (R/covstruct.R); no code reads
`covstruct_registry` any more. `cs_entry()` of a name nothing
registered refuses, and names the extension when core's home table
lists one.

ONE entry shape, used by the built-in entries and a registered one
alike (`cs_check_entry()`, and a test that every built-in passes it):

| slot | read by |
|---|---|
| `npar(blk)`, `start(blk)`, `sd_idx(blk)` | frame assembly, priors, autoscale, SE boundary reading, confint, frm_simulate, frmtmb.sample's `block_sd_idx()` |
| `nll(b, theta, blk)` (the log density; public argument `lpdf`) | objective, importance |
| `vcov(theta, blk)` | VarCorr, confint, new levels, conditional_effects, simulate |
| `natural(theta, blk)` | confint_varcorr, summary, `log_sd_theta_index()`, and `cs_is_field()`: VarCorr, ngrps, brms names, priors, hypothesis |
| `draw(theta, blk)` | `draw_b()`, which simulate() and frm_simulate() use |
| `levels_coupled` | predict new levels (`lp_extra_var()`, `old_level_pick_add()`), conditional_effects, sandwich units, importance, simulate(re_formula =) |
| `expand` (apply, jacobian, condsd) on blocks with `expands = TRUE` | `expand_b()`, the delta-method jacobian, ranef() condSD, the batch and bound predicates |
| `term` (args, parse, build, levels_of, levels_complete, usage) | parser, frame, newdata, simulate(newdata =) |
| `from_natural`, `cor_spec(blk)`, `chol_sd`/`chol_L`/`chol_A` | frm_simulate(newparams =), LKJ priors, frmtmb.sample non-centering |

Every list of structure names that stood for a property became the
property: `c("smooth", "gp", "hsgp", "car", "spde")` (6 sites in core,
one in frmtmb.sample that I did not change, see below) is
`cs_is_field()`; `c("gr_cov", "gr_prec", "car", "spde")` (6 sites) is
`cs_levels_coupled()`; the spde-only reads are `cs_levels_of()` and
`cs_levels_complete()`; `block_is_esicar()` is `block_expands()`; the
`car`/`spde`/`gp`/`hsgp`/`rr` branches for `npar` and `start` in
frame.R are gone. The field specials are parsed generically
(`parse_field_term()`, `field_specials()`): the spec slots
`carterms`/`spdeterms` are one `fieldterms` list, which frmtmb.latent's
`hmm_dpar_is_constant()` reads (changed there, one line).

### The design decision

Two seam shapes had real trade-offs:

1. ONE function, with the formula special as a `term =` slot of the
   structure (chosen). A field and its density cannot drift apart, the
   frame validates what `build()` returns (`field_term_build()`), and
   the registrant learns one call. Cost: the function has 13 arguments.
2. TWO functions, `frmtmb_register_covstruct()` for the density and a
   `frmtmb_register_special()` that maps a predictor special onto a
   structure. More general (a special could build a block of an
   existing structure, or a non-covariance term), but every pair needs a
   consistency check the first shape gets by construction, and no
   current or planned special needs the generality.

The second decision: entry functions take the BLOCK, not its dimension
(`npar(blk)` where it was `npar(dim)`), because a field's parameter
count depends on its fixed data (car's type) and an rr block's on its
rank. `by_dim()` keeps the built-in one-liners as they were. Cost: 9
test files called `$npar(9L)` and similar and now pass `list(dim = 9L)`
(listed below).

### Ported and left

Ported onto the shared shape: every built-in structure's `npar`,
`start`, `sd_idx` and `cor_spec` (block argument); `smooth`, `gp`,
`hsgp`, `car` and `spde` report through `natural`; `gr_cov`,
`gr_prec`, `car` and `spde` draw through `draw` and declare
`levels_coupled`; `car` and `spde` are fields through `term`, and
esicar's centering is car's `expand`; `rr`, `gp`, `hsgp`, `car` and
`spde` count and start their parameters in their entries.

Left special-cased in core, with the reason:

- `rr`: its `b` lives in a factor space of rank `d`, so `b` and the
  coefficient vector have different lengths and the map depends on
  theta. `expand` as published is a fixed linear map of equal length;
  generalizing it for rr is possible but not cheap.
- `gp()` and `hsgp()` as specials: their parser, the kriging of new
  positions and the HSGP basis are their own code, and neither is a
  block over a grouping factor.
- Student-t latents (`us_t`, `diag_t`): the shared mixing scale across
  blocks is a cross-block construction.
- `equalto`: still named in `imp_crosslevel()` and in the
  non-centering switch.
- Name lists that are allow-lists, not properties, stay: quadrature
  (`us`, `diag`, `homdiag`), importance (`imp_covstructs`), autoscale
  (`us`, `diag`, `us_t`, `diag_t`) and `check_re_structure()`'s
  `structured` set. Each refuses or skips an unknown structure, which
  is the right default for a registered one.
- frmtmb.sample's `draws_cov_fast()` keeps its own copy of the field
  list. It validates its map against `hyp_env_vals()`, which now uses
  `cs_is_field()`, and falls back when they disagree, so a registered
  field is slower there, not wrong. A reader export would let it ask
  core; not done, to keep frmtmb.sample unchanged.

### Evidence

`tests/testthat/test-covstruct-seam.R` (43 expectations, all pass;
`dev/spatial-log`, one process): two structures registered through the
exported function only, each in a test that calls
`local_frmtmb_registries()`, and a later test that asserts both names
are gone from the registry and the vocabulary.

- `xcs(0 + t | g)`, a toy exchangeable block with a FIXED correlation of
  0.5: `logLik()` equals the dense gaussian marginal evaluated at the
  fit's estimates to 1e-9 relative (an identity, so rounding only), and
  the dense likelihood maximized by `optim()` is not higher than the
  fit's by more than 1e-7 relative. VarCorr, `log_sd_theta_index()`,
  `confint_varcorr()` (the fixed correlation is reported at 0.5 with
  no width), an unseen level, `simulate()` and `frm_compat_features()`
  read the entry.
- `toyfld(A, gr = g)`, a field term over a fixed correlation matrix:
  the same model as core's `(1 | gr(g, cov = A))`, which is code that
  does not go through the field path. logLik agrees to 1e-7 relative,
  the coefficients to 1e-3 standard errors, the sd to 1e-3 relative;
  `ngrps()` is `NULL`; an unseen location is predicted at the
  population level (`fitted()`), as `gr(cov = A)` predicts it; a bar
  spelling and a misspelled argument are refused by name.
- Registration refusals (a core name, a grammar name, a non-syntactic
  name, a name the vocabulary holds as another kind, a negative
  `npar`, a `vcov` of the wrong arity, a `term` without `"gr"`) leave
  `cs_names()` and `frm_compat_features()` unchanged; the absent case,
  a correct registration and its replacement, takes.

Seen to fail on rellib-r11: the file cannot load there
(`object 'frmtmb_covstruct_registry' not found`), the weak form. The
behavioral failure is the wrapper defect below.

Suites after parts A and B (`dev/spatial-log/ab-all.log`, one file per
process, NOT_CRAN, each log's `lib:` line shows frmtmb 0.73.1 from the
lane library):

    frmtmb files 264 pass 18911 fail 0 err 0 skip 202
    frmtmb.coupling files 11 pass 551 fail 0 err 0 skip 5
    frmtmb.eam files 32 pass 1838 fail 0 err 0 skip 5
    frmtmb.latent files 11 pass 368 fail 0 err 0 skip 2
    frmtmb.learn files 15 pass 430 fail 0 err 0 skip 13
    frmtmb.ode files 11 pass 549 fail 0 err 0 skip 1
    frmtmb.sample files 64 pass 2971 fail 0 err 7 skip 9
    frmtmb.spline files 16 pass 594 fail 0 err 0 skip 1

The 7 frmtmb.sample errors were its tests reading
`frmtmb:::covstruct_registry`; after changing them to
`frmtmb:::cs_entry()` (`ab-sample2.log`): test-loo.R pass=134 fail=0
err=0 skip=2, test-reparam.R pass=379 fail=0 err=0 skip=0.

Tests rewritten for the entry shape (no expectation dropped):
test-sparsear1.R, test-covstruct-floor.R, test-lkj.R, test-tre.R,
test-student-share.R, test-simulate-ergonomics.R, test-compat.R,
test-importance.R (`imp_crosslevel()` is a function),
test-predict-re-uncertainty.R (`block_expands()`), helper-brms.R,
test-brms-likelihood.R, test-car-spde.R and test-review-v28.R
(`aux` for `aux_car`); helper-registry.R saves and restores
`frmtmb_covstruct_registry`, and test-registry-helper.R's probe
compares it.

## Part B. Guardrail 1 for specials and addition terms

### What changed

`R/brms-vocabulary.R` (new) holds brms 2.23.0's vocabulary as fixed
data (50 families, 15 addition terms, 19 specials) and the tables that
complete `extension_family_homes`:

- `extension_special_homes`: formula specials an extension supplies.
  Empty after part B; part C adds `car` and `spde`.
  `check_extension_specials()` runs on every top-level formula term
  before anything else reads it, on CALL positions
  (`called_functions()`), so a column named `car` is data.
- `extension_aterm_homes`: `dec = "frmtmb.eam"`.
  `check_unsupported_aterm()` runs before the general "not supported"
  list.
- `brms_family_refusals` (`custom`: use `custom_family()`),
  `brms_aterm_refusals` (`bhaz`: the basis is set in
  `cox(df =, degree =, intercept =)`), `brms_special_refusals` (`sar`,
  `fcor`, the existing refusal, now read from the table).

The refusal text is the Wiener family's: "The dec() addition term is in
the extension package frmtmb.eam, which is not loaded. Install it and
run library(frmtmb.eam), then use dec()".

### The guardrail-1 test

`tests/testthat/test-guardrail1.R` (24 expectations at the end of the
lane round, 33 after fix round 1; all pass). Every brms
name is exercised by a real call (`brmsfamily()`, the `family =
"<name>"` path, and `frm(..., dry_run = "spec")` for each term and
special) in a CHILD process holding core alone, and classified from the
refusal it gets: core, an extension the refusal names, a refusal by
name, or a generic catch-all. It asserts no generic outcome, that all
three good outcomes occur, that the set refused by name equals the set
the three refusal tables declare (so a probe that its own misspelling
refuses cannot pass for a refusal), and that the classifier calls a
generic message generic. Then, for each extension a refusal names that
is installed, a child with it loaded runs the same calls and none may
be refused: the absent-condition case. It is generic: a move that adds
a home needs no edit to the test, and an extension not installed skips
its arm. A further test checks the vocabulary against an installed brms
(`.family_*` less the `.family_info()` helper, the `resp_*` exports,
the specials among the exports).

Classification, `dev/spatial-guard-report.R`, logs
`guard-before.txt` (rellib-r11) and `guard-after-B.txt`, names that
are not "core":

    before (0.73.0)                         after part B
    family  custom   generic                 refused
    family  wiener   ext:frmtmb.eam          ext:frmtmb.eam
    aterm   bhaz     generic                 refused
    aterm   dec      generic                 ext:frmtmb.eam
    special sar      refused                 refused
    special fcor     refused                 refused

So the test fails on rellib-r11 on three names (custom by both
spellings, bhaz, dec), behaviorally.

### Found and fixed on the way: a wrapper that names no structure

reformulas finds a bar inside ANY call, so a wrapper no structure
answers to was dropped and the term fitted with the default `us`
block, with no message. On rellib-r11 (`dev/spatial-wrapper-probe.R`,
`dev/spatial-log/wrapper-probe.txt`, seed 1, 50 rows, 10 groups):

    y ~ x + foo(1 | g)      -> FIT covstruct = us logLik = -76.03093204
    y ~ x + us_t(1 | g)     -> FIT covstruct = us logLik = -76.03093204
    y ~ x + (1 | g)         -> FIT covstruct = us logLik = -76.03093204

This is a silent gap of exactly the kind guardrail 1 is about: a bar
structure from an unloaded extension would have fitted as `us`. The
comment in parse.R saying `us_t(x | g)` "lands in the not-supported
message" was wrong. `check_bar_wrapper()` now refuses such a call by
name (and points `us_t` at `gr(g, dist = "student")`); `diag(1 | g)`
and calls that hold no bar are untouched (test in
test-covstruct-seam.R).

## Part C. car() and spde() into frmtmb.spatial

### Measured before the move (guardrail 3)

All on the base commit 19b162b9 or on rellib-r11 (0.73.0), ungated
(NOT_CRAN only).

Lines in R/ (`dev/spatial-measure.R`, `measure-loc.txt`): the
top-level objects that are car's or spde's, with the comment and
roxygen lines directly above each, plus the code lines elsewhere in R/
that name them.

| | own objects | lines naming it elsewhere in R/ (17 files) |
|---|---|---|
| car | 446 | 66 |
| spde | 157 | 49 |

Most of the "elsewhere" lines are shared: lists such as
`c("gr_cov", "gr_prec", "car", "spde")` name both.

Test time (`dev/spatial-run-tests.R`, the testthat `real` time of each
`test_that()` block, `dev/spatial-log/before-core/`): the whole core
suite is 2742 blocks, 1711.23 s (2059.87 s of per-file wall time
including package loads, 262 files, 18846 passes, 202 skips).

| | blocks | seconds | share of the suite |
|---|---|---|---|
| car | 21 in test-car-spde.R (5.50 s) and 7 elsewhere (1.83 s) | 7.33 | 0.43% |
| spde | 2 in test-car-spde.R (0.25 s) and 5 elsewhere (1.94 s) | 2.19 | 0.13% |
| both | registry block (1.92 s), call-valued gr (0.25 s), data2 (0.47 s, shared with gr(cov =)) | 2.64 | 0.15% |

The "elsewhere" blocks: car in test-predict-re-uncertainty.R (0.60),
test-re-formula-partial.R (0.08), test-review-v28.R (0.63 over four
blocks), test-brms-likelihood.R (0.52 ungated; two more blocks are
gated); spde in test-review-v28.R (0.86 over four blocks) and
test-tmb-examples.R (1.08).

Churn over the last five rounds, `git log 4f5ea39f..19b162b9` (0.68.1
to 0.73.1, the 0.69 to 0.73 rounds, 37 commits): ONE commit changed a
line of R/ naming car or spde (ed824ffb, which added the two to the
list `sd_coef_unaddressable`); two changed test-car-spde.R (071184ee,
a tolerance for the least-squares starts; 7eb86c43, the boundary
message). Both are stable.

Dependencies only they need: none in Imports (Matrix, methods and RTMB
are core's). fmesher, in core's Suggests, is read only by the spde
example in test-tmb-examples.R and by `vignettes/frmtmb.Rmd`.

Ported brms rows that depend on them (`dev/spatial-brmsrows.R`,
`brmsrows.txt`, a row depends on car when its line or its test block
calls `car(` or `cor_car(`): 0 of 599 in `dev/brmsport-verdicts.tsv`
and 0 in `dev/brmsport-bin2-verdicts.tsv`. brms's car() tests are in
its stancode and standata files, whose car blocks were never ported.
spde() is not a brms special.

Reading of the numbers: both are small and stable, need no dependency
of their own, and no ported brms row reads them. By guardrail 3 they
are weak candidates, spde the weaker. They move anyway because the
brief says to, and the move is what proves the seam.

### The move

`extensions/frmtmb.spatial` (new, version 0.1.0): `R/car.R`,
`R/spde.R`, `R/zzz.R`, `R/frmtmb.spatial-package.R`; help pages
`?frmtmb-car` (alias `car`) and `?frmtmb-spde` (alias `spde`); a
vignette, `vignette("spatial", package = "frmtmb.spatial")`, that holds
core's "Spatial Gaussian Markov random fields" section; DESCRIPTION,
NEWS.md, README.md, `_pkgdown.yml`, `tests/testthat.R` and seven test
files. The code moved from core's R/covstruct.R unchanged in
arithmetic. The only rewrites are `car(con_sd =)`'s evaluation
(`car_con_sd_arg()`, the same checks and words as core's internal
`eval_spec_arg()`), a local `%||%`, and each entry written as the
arguments of `frmtmb_register_covstruct()`. `.onLoad()` registers
`car`, then `spde`, then their compat rows.

Guardrail 2, grep of the package's R/ (pasted):

    $ grep -rnE "frmtmb:::|getFromNamespace|asNamespace|:::" extensions/frmtmb.spatial/R/
    extensions/frmtmb.spatial/R/car.R:313:#' brms's scaling factor for BYM2 (brms:::.car_scale): the geometric mean
    extensions/frmtmb.spatial/R/frmtmb.spatial-package.R:8:#' `frmtmb:::`, no `getFromNamespace()`, no `asNamespace("frmtmb")`.

Both hits are comments. The NAMESPACE imports `compat_rule_builder`,
`frm_message`, `frm_stop`, `frmtmb_register_compat` and
`frmtmb_register_covstruct` from frmtmb, all exported. New core exports
for the move: none beyond part A's `frmtmb_register_covstruct()`. The
package's TESTS use `frmtmb:::` (draw_b, rr_jacobians, the batch
predicates in test-car-postfit.R, and the copied helper-brms.R), as
other extensions' tests do.

Core after the move: `car` and `spde` are listed in
`extension_special_homes` (R/brms-vocabulary.R), so without the
package `frm()` refuses `car()`, `spde()`, `(1 | car(...))` and
`(1 | spde(...))` with "The car() special is in the extension package
frmtmb.spatial, which is not loaded. Install it and run
library(frmtmb.spatial), then use car()". Their compat rows left core's
R/compat.R; the groups `crosslevel_blocks` and `wide_blocks` lost them,
and the `gmrf` group, which held nothing else, is gone. fmesher left
core's Suggests (no core file reads it now) for frmtmb.spatial's.
`brms_whole_term_specials` keeps `car`, which is a fact of brms's
grammar; `check_extension_specials()` runs before it. Core's vignettes
(frmtmb.Rmd, brms-migration.Rmd, inputs.Rmd) and the `data2` help say
the two terms are frmtmb.spatial's; the frmtmb.Rmd section is now a
pointer, with no code to guard.

car and spde are separable: `dev/spatial-out/Ccar.diff` is the tree
after the car move alone (87 files against the base) and
`dev/spatial-out/C.diff` after both (92 files); `tree-Ccar/` and
`tree-C/` beside them hold the changed files at each stage. The spde
stage depends on the car stage only for the package skeleton and for
`%||%`, which is in car.R.

### Fits are unchanged

`dev/spatial-fitcompare.R` fits the moved tests' models (the four car
types on the 4 x 4 lattice of test-car-spde.R, seed 1; icar under
poisson; esicar and icar on a graph of two 2 x 2 components and a
singleton, seed 3; spde on the 1-D chain, seed 5; spde with a Weibull
censored response on an fmesher mesh, seed 126) on rellib-r11 and on
the lane build with frmtmb.spatial loaded, each arm in one process. It
compares logLik, every estimate, every standard error
(`vcov(full = TRUE)`), `fitted()` at the data, `fitted()` on newdata
with an unseen location (`allow_new_levels = TRUE`), the population
prediction there (`re_formula = NA`), the `ranef()` conditional SDs,
`confint_varcorr()` and `simulate(nsim = 3, seed = 7)`
(`dev/spatial-log/fitcompare-C.txt`):

    before: C:/Users/adf44/source/r/rellib-r11 0.73.0 | after: C:/Users/adf44/source/r/wt-spatial-lib 0.73.1 frmtmb.spatial
    gauss_escar    logLik=identical est=identical se=identical fitted=identical fitted_new=identical pop_new=identical condsd=identical vc=identical sim=identical
    gauss_esicar   logLik=identical est=identical se=identical fitted=identical fitted_new=identical pop_new=identical condsd=identical vc=identical sim=identical
    gauss_icar     logLik=identical est=identical se=identical fitted=identical fitted_new=identical pop_new=identical condsd=identical vc=identical sim=identical
    gauss_bym2     logLik=identical est=identical se=identical fitted=identical fitted_new=identical pop_new=identical condsd=identical vc=identical sim=identical
    pois_icar      logLik=identical est=identical se=identical fitted=identical fitted_new=identical pop_new=identical condsd=identical vc=identical sim=identical
    disc_esicar    logLik=identical est=identical se=identical fitted=identical fitted_new=identical pop_new=identical condsd=identical vc=identical sim=identical
    disc_icar      logLik=identical est=identical se=identical fitted=identical fitted_new=identical pop_new=identical condsd=identical vc=identical sim=identical
    spde_chain     logLik=identical est=identical se=identical fitted=identical fitted_new=identical pop_new=identical condsd=identical vc=identical sim=identical
    spde_weibull   logLik=identical est=identical se=identical fitted=identical fitted_new=identical pop_new=identical condsd=identical vc=identical sim=identical

"identical" is `identical()` on the values, so bit for bit. The same
holds after parts A and B alone (`fits-AB.rds`) and after the car stage
(`fitcompare-Ccar.txt`).

frmtmb.sample on a car fit and an spde fit, through the seam
(`dev/spatial-sample-check.R`, `frm_sample(chains = 1, iter = 300,
warmup = 150, seed = 11)`, `sample-compare-C.txt`):

    car draws 150x11 identical: TRUE colnames identical: TRUE msgs identical: TRUE
    spde draws 150x14 identical: TRUE colnames identical: TRUE msgs identical: TRUE

Both stay centered, with frmtmb.sample's reason. frmtmb.sample's R/ is
unchanged; its `ncp_reasons` table still names the two.

### frm_compat() with the package loaded

`dev/spatial-compat-dump.R` and `dev/spatial-compat-compare.R`, every
pair, rellib-r11 with the seven extensions against the lane build with
the eight (`compat-compare-c2.txt`):

    pairs 11350 -> 11350 ; features 164 -> 164
    features identical (name, key, kind): TRUE
    pairs only before: 0  only after: 0  status differs: 0  note differs: 16

Every pair keeps its status. The 16 notes that differ are all ties
between two blanket rules, `"<a> x *"` against `"<b> x *"`, which the
resolver breaks by position (the later rule wins). Core's car and spde
blankets came before core's blankets for `equalto`, `gr_cov`, `rr`,
`smooth`, `gp`, `hsgp` and `mmc()`, and before frmtmb.sample's for
`frm_sample`. A contributed rule is appended, so car's and spde's
blankets now win those ties: 8 pairs each. With frmtmb.spatial loaded
before frmtmb.sample the `frm_sample` pairs match again (14 differ,
`compat-compare-c2-first.txt`). The pair `car x spde` matches, because
the package registers spde after car, as core's file had it. To match
every note, a contributor would have to lose a tie, which
`frmtmb_register_compat()` cannot express. I did not change the
resolver. This needs the user.

### Tests, accounted for

Counts from the logs: base tests on rellib-r11 (`before-core.log`
ungated, `before-gated.log` with FRMTMB_BRMS_FIT_TESTS) against the
lane build (`c-all.log` ungated, `c1-gated.log` and `c2-gated.log`
gated).

| file (base) | before | after, in core | after, in frmtmb.spatial |
|---|---|---|---|
| test-car-spde.R | 138 | file removed | test-car.R 118 + 3 (its registry rows); test-spde.R 13 + 3 + 1 (gr_prec's note, kept with the block) |
| test-review-v28.R | 74 | 20 | test-car-v28.R 25, test-spde.R 29 |
| test-data2.R | 31 | 29 | test-car-postfit.R 1, test-spde.R 1 |
| test-predict-re-uncertainty.R | 78 | 68 | test-car-postfit.R 10 |
| test-re-formula-partial.R | 44 | 41 | test-car-postfit.R 3 |
| test-tmb-examples.R | 28 | 27 | test-spde.R 1 |
| test-brms-likelihood.R (gated) | 543 | 476 | test-brms-likelihood-car.R 67 |
| frmtmb.sample test-reparam.R | 379 | 377 | test-car.R, a block of 3 (the 2 moved, plus `covstruct_has_chol()`) |

Every moved `test_that()` block is verbatim: dev/spatial-move-tests.R
copies each with the comment lines above it and prefixes unexported
frmtmb calls with `frmtmb:::`. The blocks that tested car and spde
together ("car() and spde() accept a call-valued gr", "every structural
special resolves from data2", "the spatial GMRF blocks are declared in
the registry") were split by term. Dropped: none. test-review-v28.R
opens with a file-level `skip_on_cran()`, which its moved blocks do not
carry in their new files, so they now run on CRAN too.

Changed counts that are not moves: test-compat.R 772 to 769, one fewer
group checked (the removed `gmrf` group) and two fewer rule patterns
(core's car and spde rules; the contributed rules are checked at
registration and by test-compat-extensions.R, whose list now includes
frmtmb.spatial). New: test-covstruct-seam.R 43, test-guardrail1.R 24,
frmtmb.spatial's test-guardrail.R 11 and test-bracket-access.R 1 (the
counts of the lane round; "Fix round 1" below has the counts after it).

frmtmb.spatial gated: 290 passes over 7 files. Ungated: 227, with the
two brms Stan rows skipped, which are the 2 skips core's
test-brms-likelihood.R lost (50, from 52).

### Plumbing

- `.github/workflows/check-frmtmb-spatial.yaml`, modeled on
  check-frmtmb-spline.yaml (no siblings).
- `.github/workflows/test-coverage.yaml`: a `frmtmb.spatial` row in the
  matrix (`siblings: ''`).
- `.github/workflows/pkgdown.yaml` and `dev/release/build-docs.R`:
  frmtmb.spatial in the list of extensions installed and built. Core's
  `_pkgdown.yml`: a navbar entry and the `frmtmb_register_covstruct`
  topic. The package's own `_pkgdown.yml`.
- `dev/release/run-sessions.sh` and `dev/release/run-session.R`:
  frmtmb.spatial in the default package list and in the list every
  session needs installed and preloads.
- `tests/testthat/test-compat-extensions.R`: frmtmb.spatial among the
  extensions the child loads.
- The floor script `dev/rel073-floor.sh` is the 0.73.0 release's and
  names frmtmb.sample alone; a release floor run for frmtmb.spatial
  needs its own arm. Measured here: frmtmb.spatial installed over
  rellib-r11 (frmtmb 0.73.0), with its floor lowered, stops at lazy
  loading with `object 'frmtmb_register_covstruct' is not exported by
  'namespace:frmtmb'` (`dev/spatial-log/floor.txt`). The base 0.73.1
  has no such export either. So the floor is the core release that
  carries this lane's seam. DESCRIPTION says `frmtmb (>= 0.73.1)` only
  because a lane does not bump core's version; it must become that
  release's version at consolidation.

## Suites, sessions, OpenBLAS and checks (final tree)

One file per process, NOT_CRAN, every log's `lib:` line showing frmtmb
0.73.1 from the lane library (`dev/spatial-log/c-all.log`, summary
`c-all-summary.txt`):

    frmtmb files 263 pass 18700 fail 0 err 0 skip 200 warn 0
    frmtmb.coupling files 11 pass 551 fail 0 err 0 skip 5 warn 0
    frmtmb.eam files 32 pass 1838 fail 0 err 0 skip 5 warn 0
    frmtmb.latent files 11 pass 368 fail 0 err 0 skip 2 warn 0
    frmtmb.learn files 15 pass 430 fail 0 err 0 skip 13 warn 0
    frmtmb.ode files 11 pass 549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample files 64 pass 3020 fail 0 err 0 skip 9 warn 0
    frmtmb.spline files 16 pass 594 fail 0 err 0 skip 1 warn 0
    frmtmb.spatial files 7 pass 227 fail 0 err 0 skip 2 warn 1

The one escaped warning was the spde toy model of frmtmb.spatial's
test-guardrail.R (an identity mesh, which leaves theta flat); the test
now uses the 1-D chain mesh and runs clean (`c4-gated.log`: pass=11
warn=0). frmtmb.latent was installed from this worktree (hmm.R reads
`fieldterms`); the other extensions are rellib-r11's builds.

One session per suite, three orders, NOT_CRAN
(`dev/release/run-sessions.sh spatial-sess 8 <all nine>`,
`dev/release/spatial-sess.log`): fail 0 and err 0 in all 27 sessions.
36 LEAK lines, every one tagged `[with load of ...]` (an extension's
registrations, ode's `rxTempDir`, and `GDAL_DATA`, which sf sets when
fmesher loads it in test-spde.R); 0 untagged.

OpenBLAS 0.3.26 with LAPACK, 4 threads (`dev/ciharden-openblas.sh
0.3.26 lapack`, `dev/ciharden-scan.sh`, every gate on), on 21 core
files that read the entries this lane changed (covstruct, the seam,
guardrail 1, gp, smooths, gr(by =), |ID| Kronecker, LKJ, Student-t,
simulate, the data2 and review-v28 remainders,
predict-re-uncertainty, confint) and frmtmb.spatial's 7: 28 of 28
files, every count equal to the reference-BLAS arm, fail 0
(`dev/ciharden-log/spatial-ob.sum` against `spatial-ref.sum`). The R
copy was deleted after.

R CMD check --as-cran (built and checked in `dev/spatial-check/`):

    frmtmb.spatial  Status: 1 NOTE   (NEWS.md has only a development
                                      heading, which R's news reader
                                      does not count as an entry)
    frmtmb.latent   Status: 1 NOTE   (V8 math rendering, expected)
    frmtmb.sample   Status: OK
    frmtmb          Status: 1 NOTE   (V8 math rendering, expected)

The first frmtmb.spatial check had 1 WARNING and 3 more NOTEs, all
fixed: a link to `frmtmb::predict()` (no such topic), `_pkgdown.yml`
at the top level (now in `.Rbuildignore`, copied from frmtmb.spline),
RTMB and stats in Imports with nothing imported, and codetools reading
the `build` slot's resolver, then called `data()`, as `utils::data()`.
The resolver is now named `lookup` in the seam's documentation and in
every caller.

## Stage patches

- `dev/spatial-out/AB-final.diff`: parts A and B alone, at their final
  state (59 files against 19b162b9 after fix round 1, written by
  `dev/spatial-ab-diff.sh`). It is `AB.diff`, the snapshot
  the A and B suites ran on, plus the A and B fixes made later
  (`dev/spatial-stage-fix.R`: the `lookup` name, one wrapped test line,
  frmtmb.latent's NEWS bullet). It applies to the base and every R file
  of the result parses.
- `dev/spatial-out/Ccar.diff`: the snapshot after the car move alone,
  on which the full suites ran (`ccar-all.log`; one frmtmb.sample
  error, `test-reparam.R` asking `ncp_eligible()` of the unregistered
  `car`, fixed in the final tree).
- `dev/spatial-out/C.diff`: the final tree; `dev/spatial-snapshot.sh C`
  regenerates it.

## Fix round 1 (review dev/reviews/2026-10-23-spatial.md)

The fix scripts are `dev/spatial-fix-*.R`; each takes a tree and was
run on the worktree and on the A and B stage tree, so that
`AB-final.diff` carries the A and B fixes and nothing of part C.

- A1. `frmtmb_register_covstruct()` runs `cs_check_entry()` on the
  entry it builds. The frame runs `cs_check_block()` once per block:
  `start` length and finiteness, `sd_idx` range, `vcov` finite,
  dim x dim and named by the block's columns (unless `natural`), and
  for a registered structure `lpdf` finite at zero and larger there
  than at 1e3 (a negated density is refused by name). Seven refusals
  and the absent case (built-ins and two toys) in
  test-covstruct-seam.R.
- A2. The OLRE check in `R/fit.R` reads `cs_levels_coupled()`. Tests:
  frmtmb.spatial test-car.R (6 x 6 lattice, one row per area, escar,
  icar, bym2: no warning) and test-covstruct-seam.R (a coupled toy: no
  warning; `(1 | g)` with one row per level still warns). NEWS bug fix.
- A3. Each block records `covstruct_package`. The helper is
  `require_structure_packages(x)` in `R/extension-home.R`, which is
  splitsurv's file verbatim (`require_extension(pkg, what)`,
  `require_model_packages()`) plus this one caller-side function. It is
  called by `build_objective()`, `sdr_of()` (summary, vcov, Wald
  confint), print, ranef, VarCorr, predict, fitted, frm_linpred,
  simulate, update and confint_varcorr; `cs_entry()`'s fallback uses
  `require_extension()` too. Class `frmtmb_extension_missing`. Tests:
  frmtmb.spatial test-guardrail.R, a child process with saved car and
  spde fits, 11 methods x 2 fits, in two arms: frmtmb.spatial
  installed but unloadable (a stub in front of every library) and not
  installed (a library without it; skipped under R CMD check, which
  puts frmtmb's dependencies in the library that holds the package).
  Core: a toy from the uninstalled package `frmtmbnotinstalled` and a
  session-registered toy.
- A4. Reserved: formula functions of base, stats, splines, utils,
  methods and mgcv (`poly`, `ns`, `s`, `t2`, `te`, `offset`, ...),
  the brms vocabulary (`gr`, `mm`, `mo`, `me`, `mi`, `cs`, `gp`,
  families and addition terms), core's family registry, every
  extension's home names and the registered addition terms.
- A5. Replacement only by the registering package, or with
  `replace = TRUE`; documented under `replace` and `package`; tested.
- A6. `frmtmb_register_covstruct()` writes the structure's
  `importance` and `quadrature` rows of `frm_compat()` (refused), keyed
  by the structure so a replacement replaces them. The notes are
  constants in `R/compat.R` that the groups `crosslevel_blocks` and
  `wide_blocks` share, so frmtmb.spatial's `zzz.R` no longer writes
  those four rows. A first form (kind rules in core) won no pair and
  failed the contract test; it is gone.
- A7. `bar_wrapper_head()` walks the whole call, so `foo(bar(1 | g))`
  is refused; `re_formula =` gets its own wording. Tests with the
  absent cases.
- A8. `frmtmb_local_registries()` is exported (`R/registries.R`, with
  docs and a pkgdown entry); the test helper wraps it. The comment at
  R/covstruct.R:1592 is fixed. The 21/24 in this file now reads 24
  (lane round) and 33 (after this round).
- B1. `(1 | gr(g))` is `(1 | g)`: same logLik, VarCorr and variables
  (test-guardrail1.R). The probes take every special in its plainest
  brms spelling and in the others brms documents. NEWS corrected.
- B2. The bhaz probe uses `gaussian()`. The pointer to `cox()` is
  `bhaz_family_pointer` in `R/brms-vocabulary.R`, the one place to
  change it. The cox wording is unchanged.
- B3. The bhaz refusal says `bhaz(gr = )` has no frmtmb equivalent.
- B4. The extension arm calls `skip()` naming the absent packages.
- B5. The bare-string arm calls `frm(family = "<name>")`.
- C2. Contributed rules resolve in the order of the owning package's
  name. With frmtmb.spatial loaded last and first: 11350 pairs, 164
  features, status differs 0, notes differ 16 against the base in both
  orders (the 16 accepted), 0 between the two orders
  (`compat-compare-f1b-*.txt`).
- C4. `dev/spatial-fitcompare.R` prints a `lib:` line per package:
  `lib: frmtmb.spatial C:/Users/adf44/source/r/wt-spatial-lib`; all 9
  models identical in every quantity (`fitcompare-F1.txt`).
- T1. test-importance.R 226 to 224: the block "every covariance
  structure importance refuses reads refused" loops over the registry,
  which no longer holds car and spde (27 to 25 there); the two moved to
  frmtmb.spatial test-car.R with quadrature (4 expectations). Nothing
  dropped. test-bracket-access.R 33 to 36: one more extension's sources
  in "no extension file reads ..." (+1) and the new "every extension
  in the tree runs the guard on itself" (+2). test-dirichlet-
  multinomial.R 66 to 67: the base log was 0.73.0 on rellib-r11, where
  one expectation failed (the fitted() clamp fixed in 0.73.1); not a
  lane change.
- X1. test-setprior.R, test-v07.R, test-v11.R and test-v12.R are CRLF
  again (`dev/spatial-crlf.R`); no changed file differs from the base
  only by line endings.
- S1. `dev/spatial-ab-diff.sh` writes `AB-final.diff` from the stage
  tree against `git archive 19b162b9`, with `a/<path>` headers for new
  files too; `git apply --check` on the base passes.

Suites after the round, one file per process, NOT_CRAN and
FRMTMB_BRMS_FIT_TESTS (`dev/spatial-log/f1c-all.log`, 430 of 430):

    frmtmb files 263 pass 21211 fail 0 err 0 skip 14 warn 0
    frmtmb.coupling files 11 pass 551 fail 0 err 0 skip 5 warn 0
    frmtmb.eam files 32 pass 1844 fail 0 err 0 skip 3 warn 0
    frmtmb.latent files 11 pass 368 fail 0 err 0 skip 2 warn 0
    frmtmb.learn files 15 pass 501 fail 0 err 0 skip 2 warn 0
    frmtmb.ode files 11 pass 549 fail 0 err 0 skip 1 warn 0
    frmtmb.sample files 64 pass 3290 fail 0 err 0 skip 1 warn 0
    frmtmb.spatial files 7 pass 390 fail 0 err 0 skip 0 warn 0
    frmtmb.spline files 16 pass 594 fail 0 err 0 skip 1 warn 0

frmtmb.spatial test-guardrail.R is 103 there (57 before its second,
unloadable arm). The A and B stage tree, installed alone (no
frmtmb.spatial in the lane library), ran core's 264 files: pass 21485,
fail 0, err 0, skip 14, warn 0 (`ab2-core.log`).

One session per suite, three orders, NOT_CRAN
(`dev/release/spatial-f1.log`): frmtmb, frmtmb.spatial, frmtmb.latent
and frmtmb.sample, 12 of 12, fail 0 and err 0 in every session; 27
LEAK lines, all tagged `[with load of ...]`, 0 untagged.

R CMD check --as-cran: frmtmb `Status: 1 NOTE` (V8 math rendering,
expected); frmtmb.spatial `Status: 1 NOTE` (NEWS.md development
heading). frmtmb.latent's R code did not change this round, so it was
not checked again.

`AB-final.diff` now has 59 files and `C.diff` 100.

## Fix round 2 (review section "Short check of fix round 1")

`dev/spatial-fix-r2.R` and `dev/spatial-fix-r2b.R`, run on the
worktree and the stage tree.

- F1. `cs_check_block()` no longer compares `lpdf` at `b = 1000` with
  `b = 0`; it keeps the finiteness check at zero. No curvature test
  replaces it: a bimodal density with modes at -1 and 1 has a minimum
  at zero, and an improper one is flat along its null space, so no test
  of the shape at zero tells a correct density from a negated one; the
  author's own test does. The `lpdf` parameter's documentation, the
  internal comment and NEWS say so. test-covstruct-seam.R: the
  sign-flip refusal became a refusal of an `lpdf` that is `-Inf` at
  zero, and a new block fits three correct structures: the reviewer's
  `rwone` (intrinsic RW1, `Q1 = crossprod(diff(diag(30)))`, seed 7,
  `y ~ 0 + x + rwone(Q1, gr = t)`), a mean of 600 (logLik equal to
  `(1 | g)`'s, relative 1e-6) and a two-mode mixture fitted to data
  drawn from it (sd within a factor of 2 of 0.3; on the iid data of the
  other tests it does not converge, which is the data, not the check).
- F2. The base/stats/splines/utils/methods/mgcv lookup is gone. The
  refused formula functions are a fixed list, `covstruct_formula_terms`
  in `R/covstruct-seam.R` (poly, polym, ns, bs, s, t2, te, ti, offset,
  I, C, cbind, Surv, strata, cluster, the transformations log ... scale,
  and the recodings factor ... relevel), plus the brms vocabulary, the
  family registry and the other extensions' names, as before. Tests:
  `toeplitz`, `kernel`, `matern` and `rw1` register; `poly`, `ns`,
  `log` and `s` are refused. In a fresh process `library(frmtmb)` and
  then `library(frmtmb.spatial)` leave mgcv unloaded
  (`loadedNamespaces()`), now a test in frmtmb.spatial's
  test-guardrail.R.

Runs, one file per process, NOT_CRAN and the brms gate
(`dev/spatial-log/r2.log`): core test-covstruct-seam.R 131,
test-guardrail1.R 33, test-compat.R 769, test-compat-extensions.R 2,
test-compat-register.R 103, test-registry-helper.R 7,
test-message-uniqueness.R 6, test-bracket-access.R 36; frmtmb.spatial
test-guardrail.R 104, test-car.R 132, test-spde.R 48,
test-brms-likelihood-car.R 67, test-car-postfit.R 14, test-car-v28.R
25, test-bracket-access.R 1. All fail 0, err 0, warn 0. The stage tree
installed alone (`ab3.log`): test-covstruct-seam.R 131, test-compat.R
772, test-guardrail1.R 30, test-car-spde.R 138 and the other five as
in core, all fail 0.

## What I did not do, and what needs the user

- frm_compat() notes: 16 pairs that tie between two blanket rules now
  show car's or spde's note (statuses all unchanged). A resolver change
  could restore them; I left the resolver alone.
- frmtmb.sample's `draws_cov_fast()` keeps its own list of field
  structures; a registered field is correct there (the map falls back)
  but slower. A core reader export would close it.
- `rr`'s factor-space expansion, `gp()`/`hsgp()` and the Student-t
  sharing are not on the published entry shape (see Part A).
- Floors and bumps: frmtmb minor (a new export, new features, a
  breaking move). frmtmb.spatial is new at 0.1.0, with its floor the
  frmtmb release that carries the seam. frmtmb.latent needs the same
  floor (its `hmm_dpar_is_constant()` reads `fieldterms`); a patch bump
  on its own, more with the mvnlatent lane's move. frmtmb.sample: only
  its tests changed (two read `frmtmb:::cs_entry()`, one list lost car
  and spde), so no bump is needed and its floor does not move.
- The brmsport ledger was not regenerated; no row depends on car or
  spde (Part C, measured), so no row flips.
