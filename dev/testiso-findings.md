# Lane testiso: test isolation in one session

Base 47cf659d (frmtmb 0.71.0 plus the test-parse.R fix), 2026-10-14.
Before arm: `rellib-r9` with the base tests read from the main checkout.
After arm: this worktree's tests, frmtmb and frmtmb.spline installed
into `C:/Users/adf44/source/r/wt-testiso-lib` (no package code is
changed after fix round 2, so both are the same build as `rellib-r9`),
the other extensions from `rellib-r9`.

## Summary

- test-compat.R and test-compat-register.R pass with no extension, with
  each extension, and with all seven preloaded, and leave every
  registry as they found it. Every changed assertion was seen to fail
  under a mutant of the behavior it guards, with and without
  extensions loaded.
- What the narrowed compat tests no longer check in an arbitrary
  session is checked with all seven loaded, in a child process
  (test-compat-extensions.R, fix round 1). A contributed method row
  must name an exported function or, in core's own spelling
  `<function>_<argument>`, a formal argument of one; frmtmb.spline's
  `frm_curve_contrast` (the contrast argument of `frm_curve()`)
  passes. Fix round 1 renamed that row instead; fix round 2 reverted
  the rename, so frmtmb.spline is unchanged and needs no bump.
- One-session result: all 24 sessions (8 suites; forward, reverse,
  preload) read `fail=0 err=0` after fix round 1, with no LEAK line a
  test made itself, and core's three orders again after fix round 2
  (section "Fix round 2"). Base: core failed 9 forward, 8 plus 1 error
  reverse, and 19
  plus 1 error preload, and 25 LEAK lines across the three orders came
  from tests. The seven extension suites passed in every order on both.
- The audit found five tests or helpers that left a registry changed
  (test-parse.R, test-cens-trunc.R, helper-rl.R, frmtmb.ode
  test-compat.R, frmtmb.sample test-compat-preflight.R), and 39
  `data()` calls in 20 files that wrote to the global environment.
  Each now saves and restores, never deletes.
- `dev/release/run-session.R` and `dev/release/run-sessions.sh` run a
  suite in one session three ways and print a LEAK line for each
  registry, removed notice, option, environment variable, search-path
  entry or global a file leaves changed. The release gate is in
  `dev/round-handoff.md` ("Consolidation") and the rule in
  `dev/lane-rules.md` ("Added after the 0.71.0 round").
- CI: coverage ran core against whatever extensions the restored cache
  held (pkgdown seeds it). By the user's decision,
  `.github/workflows/test-coverage.yaml` now measures core with every
  frmtmb package removed, and each extension in its own matrix job,
  with core and the siblings it names installed from the checkout, and
  `FRMTMB_SAMPLER_GATES: "false"` as the other CI jobs of frmtmb.sample
  set it. Only a push tests the workflow.
- No package code changed. No NEWS entry and no version bump: the
  changes are tests, helpers, dev and CI.

## Fix round 2 (re-check appended to the review)

- **R1. The frmtmb.spline rename is reverted** (`R/zzz.R` and
  `NEWS.md` rewritten from `git show HEAD:`; `git status` shows no
  frmtmb.spline file changed, and frmtmb.spline was reinstalled into
  wt-testiso-lib from that tree). The key `frm_curve(contrast =)`
  broke `R/compat.R`'s rule that a key carries no `()`, and core spells
  its own argument-setting rows `<function>_<setting>`.
  `compat_contract_violations()` now accepts a contributed method row
  `m` that is an exported function, or that splits at some underscore
  into `f` and `a` with `f` a function a loaded frmtmb package exports
  and `a` a formal argument of `f` other than `...`. Every split is
  tried, because function names hold underscores. The comment at
  test-compat.R:363 now says that frm_curve_contrast is the contrast
  argument of frm_curve() and where contributed rows are checked.
  Mutants (`dev/testiso-contract-mutants.R`, all seven loaded, spline
  0.10.2 with its `frm_curve_contrast` row; log
  `dev/testiso-log/contract-mutants2.txt`):

<!-- BEGIN GENERATED: dev/testiso-contract-mutants.R (round 2) -->

```
CONTRACT none | spline 0.10.2 wt-testiso-lib | frm_curve_contrast row: TRUE | violations: 0 
CONTRACT typo | spline 0.10.2 wt-testiso-lib | frm_curve_contrast row: TRUE | violations: 2 
   unresolved rule side: hmmm 
   rule wins no pair: gaussian x hmmm 
CONTRACT method | spline 0.10.2 wt-testiso-lib | frm_curve_contrast row: TRUE | violations: 1 
   method row names no exported function or argument: zz_contrib_method 
CONTRACT misarg | spline 0.10.2 wt-testiso-lib | frm_curve_contrast row: TRUE | violations: 1 
   method row names no exported function or argument: frm_curve_contrsat 
CONTRACT dots | spline 0.10.2 wt-testiso-lib | frm_curve_contrast row: TRUE | violations: 1 
   method row names no exported function or argument: frm_sample_... 
CONTRACT dead | spline 0.10.2 wt-testiso-lib | frm_curve_contrast row: TRUE | violations: 2 
   rule wins no pair: royston_parmar x cens() 
   rule wins no pair: royston_parmar x cens() 
```

<!-- END GENERATED -->

  `misarg` is the misspelled argument row the review asked for
  (`frm_curve_contrsat`); `dots` (`frm_sample_...`) closes the gap the
  re-check named, since `frm_sample()` has `...` among its formals. The
  control passes with `frm_curve_contrast` present.
- **R2. `FRMTMB_SAMPLER_GATES: "false"`** is in the `extension-coverage`
  job's `env:`, with the reason, as `check-frmtmb-sample.yaml:59` and
  `pkgcheck.yaml:28` set it. The check-frmtmb-* workflows set three
  variables in all: `NOT_CRAN` (already set), `FRMTMB_SAMPLER_GATES`
  (now carried over) and `_R_CHECK_DONTTEST_EXAMPLES_` (sample's job;
  read by R CMD check's examples, which covr does not run, so not
  carried). The extensions' tests read no other variable a workflow
  sets (`Sys.getenv()` calls in `extensions/*/tests`: the FRMTMB_SCALE_*,
  brmsport, BRMS_FIT and LEARN_LP_LOG gates, which no CI job sets). No
  other CI step was added. YAML validation (`dev/testiso-yaml.R`, which
  now also asserts the matrix job's NOT_CRAN and FRMTMB_SAMPLER_GATES;
  log `dev/testiso-log/yaml2.txt`):

<!-- BEGIN GENERATED: dev/testiso-yaml.R (round 2) -->

```
jobs: test-coverage extension-coverage 
test-coverage / Remove every frmtmb package from the library : parses 
test-coverage / Measure coverage : parses 
test-coverage cache-version: coverage-core-1 
extension-coverage / Install core and the siblings from this checkout : parses 
extension-coverage / Measure coverage : parses 
extension-coverage cache-version: coverage-${{ matrix.pkg }}-1 
matrix: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline 
frmtmb.coupling | names:  | siblings:  | same set: TRUE | order ok: TRUE 
frmtmb.eam | names: frmtmb.sample | siblings: frmtmb.sample | same set: TRUE | order ok: TRUE 
frmtmb.latent | names:  | siblings:  | same set: TRUE | order ok: TRUE 
frmtmb.learn | names: frmtmb.eam | siblings: frmtmb.eam | same set: TRUE | order ok: TRUE 
frmtmb.ode | names:  | siblings:  | same set: TRUE | order ok: TRUE 
frmtmb.sample | names: frmtmb.eam frmtmb.latent frmtmb.ode | siblings: frmtmb.latent frmtmb.ode frmtmb.eam | same set: TRUE | order ok: TRUE 
frmtmb.spline | names:  | siblings:  | same set: TRUE | order ok: TRUE 
extension-coverage env: GITHUB_PAT=${{ secrets.GITHUB_TOKEN }} NOT_CRAN=true CODECOV_TOKEN=${{ secrets.CODECOV_TOKEN }} PKG=${{ matrix.pkg }} SIBLINGS=${{ matrix.siblings }} FRMTMB_SAMPLER_GATES=false 
```

<!-- END GENERATED -->

- **Text.** round-handoff.md's consolidation sentence says "three
  ways"; test-compat.R:363 as above.
- **Optional, done: frmtmb.sample in the helper's proof.**
  test-registry-helper.R now loads frmtmb.learn and frmtmb.sample
  inside the scope in one child and outside any scope in another, and
  compares the registries and frmtmb.sample's active bindings (count,
  and that each resolves to a function). The first run FAILED: the
  scoped child's features came out in a different order from the plain
  one (`frm_sample` first), because the helper reran `.onLoad()` in an
  order derived from imports, and frmtmb.sample imports neither learn
  nor eam. A contributed rule's precedence is its position, so that
  was a real inconsistency, not a cosmetic one. The helper now records
  the order the frmtmb packages load in with an onLoad hook per
  package not yet loaded (and puts the previous hooks back), and
  reruns in that order. Log of the failing run:
  `dev/release/testiso-r3.log` (`test-registry-helper.R pass=6
  fail=1`); after the fix, `testiso-r3b.log`.

Changed files, one per process (`dev/release/run-par.sh`,
`FRMTMB_RELLIB` = wt-testiso-lib):

<!-- BEGIN GENERATED: dev/release/testiso-r3.log, testiso-r3b.log -->

```
RESULT test-compat-extensions.R pass=2 fail=0 err=0 skip=0 warn=0
RESULT test-registry-helper.R pass=6 fail=1 err=0 skip=0 warn=0
RESULT test-compat.R pass=797 fail=0 err=0 skip=0 warn=0
RESULT test-surface.R pass=47 fail=0 err=0 skip=0 warn=0
TESTISO-R3 ran 4 of 4
--- after the load-order fix (testiso-r3b):
RESULT test-registry-helper.R pass=7 fail=0 err=0 skip=0 warn=0
RESULT test-compat-register.R pass=103 fail=0 err=0 skip=0 warn=0
RESULT test-parse.R pass=44 fail=0 err=0 skip=0 warn=0
RESULT test-cens-trunc.R pass=69 fail=0 err=0 skip=0 warn=0
RESULT test-rl-example.R pass=96 fail=0 err=0 skip=2 warn=0
TESTISO-R3B ran 5 of 5
```

<!-- END GENERATED -->

Core in one session, three orders, final tree
(`run-sessions.sh testiso-after3 3 frmtmb`; log
`dev/release/testiso-after3.log`, lines without a load tag):

<!-- BEGIN GENERATED: dev/release/testiso-after3.log -->

```
SESSION frmtmb forward files=228/238 pass=17449 fail=0 err=0 skip=185 warn=0
SESSION frmtmb reverse files=228/238 pass=17547 fail=0 err=0 skip=185 warn=0
SESSION frmtmb preload files=228/238 pass=17643 fail=0 err=0 skip=186 warn=0
TESTISO-AFTER3 ran 3 of 3
```

<!-- END GENERATED -->

The seven extension suites were not rerun in this round: no extension
file differs from round 1's tree except frmtmb.spline, which is back
to its base, and spline's suite passed all three orders on the base
build (before table) and test-surface.R passes on the reverted build
(above).

## Fix round 1 (review `dev/reviews/2026-10-14-testiso.md`)

1. **Contributed rows are checked with every extension loaded.** New
   `tests/testthat/test-compat-extensions.R` loads all seven in a
   child Rscript (this session is left unchanged) and runs
   `compat_contract_violations()` (new
   `tests/testthat/helper-compat-contract.R`): no unresolved rule side,
   no rule that wins zero pairs (no excuse), and every contributed
   `method` row names a function a loaded frmtmb package exports, or an
   argument of one spelled `f(arg =)` whose `arg` is in `formals(f)`.
   It skips unless all seven are installed (so under R CMD check and on
   the core coverage job). It flagged frmtmb.spline's
   `frm_curve_contrast`, which named no function. Fix round 1 renamed
   the spline feature to `"frm_curve(contrast =)"` (with a NEWS
   bullet); **fix round 2 reverted that rename** and made the check
   accept core's `<function>_<argument>` spelling instead (section
   "Fix round 2"), so the paragraphs below record round 1's state.
   Round 1's evidence (`dev/testiso-contract-mutants.R` as it was then,
   log `dev/testiso-log/contract-mutants.txt`; all seven loaded):

<!-- BEGIN GENERATED: dev/testiso-contract-mutants.R -->

```
CONTRACT none lane spline 0.10.2 C:/Users/adf44/source/r/wt-testiso-lib | violations: 0 
CONTRACT typo lane spline 0.10.2 C:/Users/adf44/source/r/wt-testiso-lib | violations: 2 
   unresolved rule side: hmmm 
   rule wins no pair: gaussian x hmmm 
CONTRACT method lane spline 0.10.2 C:/Users/adf44/source/r/wt-testiso-lib | violations: 1 
   method row names no exported function or argument: zz_contrib_method 
CONTRACT badarg lane spline 0.10.2 C:/Users/adf44/source/r/wt-testiso-lib | violations: 1 
   method row names no exported function or argument: frm_curve(no_such_arg =) 
CONTRACT dead lane spline 0.10.2 C:/Users/adf44/source/r/wt-testiso-lib | violations: 2 
   rule wins no pair: royston_parmar x cens() 
   rule wins no pair: royston_parmar x cens() 
CONTRACT none base spline 0.10.2 C:/Users/adf44/source/r/rellib-r9 | violations: 1 
   method row names no exported function or argument: frm_curve_contrast 
```

<!-- END GENERATED -->

   `typo` and `method` are the reviewer's mutants; `badarg` and `dead`
   are mine. The last line is the base build (rellib-r9 spline 0.10.2):
   the check names `frm_curve_contrast`. The test itself, one per
   process: lane `pass=2 fail=0`, base (rellib-r9 first) `pass=1
   fail=1` with the violation `method row names no exported function
   or argument: frm_curve_contrast` (`dev/release/testiso-ce-lane.log`,
   `testiso-ce-base.log`).
2. **Global and search-path leaks are never tagged**, and there is a
   third order, `preload` (alphabetical, all seven extensions loaded
   first), in `run-session.R` and `run-sessions.sh` (24 sessions). The
   release gate in `dev/round-handoff.md` and the rule in
   `dev/lane-rules.md` now name forward, reverse and preload. On the
   base tests, preload shows the frmtmb.sample preflight leak the first
   version could not see (`test-compat-preflight.R compat features:
   12 -> 13 entries`, untagged) and 19 failures plus 1 error in core
   (tables below).
3. **The runner reports removed notice tags** (`LEAK <file> notices:
   removed ...`); an added tag is not a leak.
4. **`local_frmtmb_registries()` keeps an extension that loads inside
   its scope registered.** After the restore it runs the `.onLoad()`
   of every `frmtmb.*` namespace that appeared, a package after the
   frmtmb packages it imports. This keeps the registries consistent
   with the loaded namespaces, which a `stop()` would not; every
   extension's `.onLoad()` only registers, apart from frmtmb.sample's
   `frm_install_generics()`, which is idempotent. helper-rl.R's
   whole-run scope now ends with every extension the run loaded
   registered again, so a `devtools::test()` session keeps them. New
   `tests/testthat/test-registry-helper.R`: a child loads frmtmb.learn
   (which loads frmtmb.eam) inside the scope and another outside any
   scope, and the registries come out identical and `y | dec(d)`
   parses.
5. **A renamed registry fails loudly.** The helper and the runner read
   each registry with `get(nm, envir = ns, inherits = FALSE)` and
   `stopifnot(is.environment())`. test-registry-helper.R asserts that
   the namespace's top-level environments are exactly the helper's six
   registries (`frmtmb_test_registries`) plus three caches
   (`frmtmb_test_caches`), and that an absent name errors.
6. **CI diagnosis corrected** (section "CI" below): pkgdown seeds the
   cache and coverage carries it, not the check-frmtmb-* jobs; an
   exact-key hit freezes the set; the R CMD check mask is the internal
   `R_check_suggests_only`, not an environment variable.
7. **test-confint-anova.R**: the comment and the findings now say
   `update()` alone re-evaluates the call.
8. **run-sessions.sh**: the handoff runs tier `rel<ver>`, not
   `rel<ver>-sessions`, and a non-number N is refused (`bash
   dev/release/run-sessions.sh x frmtmb` prints "N must be a number,
   got: frmtmb" and exits 1).

Changed files in this round, one per process with the final tree
(`dev/release/run-par.sh`, `FRMTMB_RELLIB` = wt-testiso-lib; log
`dev/release/testiso-r2.log`; test-compat-extensions.R rerun after a
fix to its child's output as `testiso-ce-lane`):

<!-- BEGIN GENERATED: dev/release/testiso-r2.log, testiso-ce-lane.log -->

```
RESULT test-compat-extensions.R pass=1 fail=1 err=0 skip=0 warn=0
RESULT test-registry-helper.R pass=5 fail=0 err=0 skip=0 warn=0
RESULT test-compat.R pass=797 fail=0 err=0 skip=0 warn=0
RESULT test-compat-register.R pass=103 fail=0 err=0 skip=0 warn=0
RESULT test-parse.R pass=44 fail=0 err=0 skip=0 warn=0
RESULT test-cens-trunc.R pass=69 fail=0 err=0 skip=0 warn=0
RESULT test-confint-anova.R pass=35 fail=0 err=0 skip=0 warn=0
RESULT test-rl-example.R pass=96 fail=0 err=0 skip=2 warn=0
RESULT test-surface.R pass=47 fail=0 err=0 skip=0 warn=0
TESTISO-R2 ran 9 of 9
--- testiso-ce-lane (after the child-output fix):
RESULT test-compat-extensions.R pass=2 fail=0 err=0 skip=0 warn=0
```

<!-- END GENERATED -->

## One-session results, before and after

`NOT_CRAN=true`, no gate, every extension installed
(`dev/release/run-sessions.sh`, summarized by `dev/testiso-counts.R`).
"own LEAK" counts LEAK lines that no namespace load explains.

### Before (rellib-r9, base tests from the main checkout)

The forward and reverse sessions ran the first version of the runner. The preload sessions were added in fix round 1 and ran the final runner on the same base tests.

<!-- BEGIN GENERATED: dev/testiso-counts.R testiso-before -->

From `dev/release/testiso-before-sessions/` (24 sessions).

| package | order | files | pass | fail | err | skip | warn | own LEAK | frmtmb from |
|---|---|---|---|---|---|---|---|---|---|
| frmtmb | forward | 226/236 | 17435 | 9 | 0 | 185 | 0 | 5 | rellib-r9 |
| frmtmb | preload | 226/236 | 17647 | 19 | 1 | 186 | 0 | 5 | rellib-r9 |
| frmtmb | reverse | 226/236 | 17680 | 8 | 1 | 185 | 0 | 5 | rellib-r9 |
| frmtmb.coupling | forward | 11/11 | 542 | 0 | 0 | 5 | 0 | 0 | rellib-r9 |
| frmtmb.coupling | preload | 11/11 | 542 | 0 | 0 | 5 | 0 | 0 | rellib-r9 |
| frmtmb.coupling | reverse | 11/11 | 542 | 0 | 0 | 5 | 0 | 0 | rellib-r9 |
| frmtmb.eam | forward | 30/30 | 1808 | 0 | 0 | 3 | 0 | 0 | rellib-r9 |
| frmtmb.eam | preload | 30/30 | 1808 | 0 | 0 | 3 | 0 | 0 | rellib-r9 |
| frmtmb.eam | reverse | 30/30 | 1808 | 0 | 0 | 3 | 0 | 0 | rellib-r9 |
| frmtmb.latent | forward | 10/10 | 360 | 0 | 0 | 2 | 0 | 0 | rellib-r9 |
| frmtmb.latent | preload | 10/10 | 360 | 0 | 0 | 2 | 0 | 0 | rellib-r9 |
| frmtmb.latent | reverse | 10/10 | 360 | 0 | 0 | 2 | 0 | 0 | rellib-r9 |
| frmtmb.learn | forward | 15/15 | 430 | 0 | 0 | 13 | 0 | 0 | rellib-r9 |
| frmtmb.learn | preload | 15/15 | 430 | 0 | 0 | 13 | 0 | 0 | rellib-r9 |
| frmtmb.learn | reverse | 15/15 | 430 | 0 | 0 | 13 | 0 | 0 | rellib-r9 |
| frmtmb.ode | forward | 11/11 | 549 | 0 | 0 | 1 | 0 | 2 | rellib-r9 |
| frmtmb.ode | preload | 11/11 | 549 | 0 | 0 | 1 | 0 | 2 | rellib-r9 |
| frmtmb.ode | reverse | 11/11 | 549 | 0 | 0 | 1 | 0 | 2 | rellib-r9 |
| frmtmb.sample | forward | 48/53 | 2697 | 0 | 0 | 4 | 0 | 1 | rellib-r9 |
| frmtmb.sample | preload | 48/53 | 2697 | 0 | 0 | 4 | 0 | 2 | rellib-r9 |
| frmtmb.sample | reverse | 48/53 | 2697 | 0 | 0 | 4 | 0 | 1 | rellib-r9 |
| frmtmb.spline | forward | 16/16 | 594 | 0 | 0 | 1 | 0 | 0 | rellib-r9 |
| frmtmb.spline | preload | 16/16 | 594 | 0 | 0 | 1 | 0 | 0 | rellib-r9 |
| frmtmb.spline | reverse | 16/16 | 594 | 0 | 0 | 1 | 0 | 0 | rellib-r9 |

Failing and erroring tests (expectation count):

- frmtmb--forward: FAIL test-compat-register.R [a declared but unresolved rule side is reported, not dropped] x4
- frmtmb--forward: FAIL test-compat-register.R [expects = accepts a name another package already supplies] x1
- frmtmb--forward: FAIL test-compat.R [declared post-fit methods are callable] x1
- frmtmb--forward: FAIL test-compat.R [every rule wins at least one pair] x1
- frmtmb--forward: FAIL test-compat.R [rule patterns name real kinds, groups, and features] x2
- frmtmb--preload: ERROR test-compat-register.R [an unknown spelling with no near miss is refused without a guess] x1
- frmtmb--preload: FAIL test-compat-register.R [a declared but unresolved rule side is reported, not dropped] x4
- frmtmb--preload: FAIL test-compat-register.R [a refused registration leaves the vocabulary as it found it] x1
- frmtmb--preload: FAIL test-compat-register.R [a registered addition term becomes a compat feature] x1
- frmtmb--preload: FAIL test-compat-register.R [a rules = builder that throws says which argument failed] x1
- frmtmb--preload: FAIL test-compat-register.R [an unknown spelling with no near miss is refused without a guess] x1
- frmtmb--preload: FAIL test-compat-register.R [expects = admits another package's feature without adding it] x1
- frmtmb--preload: FAIL test-compat-register.R [expects = refuses a name this very call supplies] x1
- frmtmb--preload: FAIL test-compat.R [declared post-fit methods are callable] x8
- frmtmb--preload: FAIL test-compat.R [declared specials still parse] x1
- frmtmb--reverse: ERROR test-compat-register.R [an unknown spelling with no near miss is refused without a guess] x1
- frmtmb--reverse: FAIL test-compat-register.R [a refused registration leaves the vocabulary as it found it] x1
- frmtmb--reverse: FAIL test-compat-register.R [a registered addition term becomes a compat feature] x1
- frmtmb--reverse: FAIL test-compat-register.R [a rules = builder that throws says which argument failed] x1
- frmtmb--reverse: FAIL test-compat-register.R [an unknown spelling with no near miss is refused without a guess] x1
- frmtmb--reverse: FAIL test-compat-register.R [expects = refuses a name this very call supplies] x1
- frmtmb--reverse: FAIL test-compat.R [declared addition terms are the ones the parser accepts] x3

LEAK lines a test made itself (not a namespace load):

- frmtmb--forward: test-aliased-grouping.R globals: +sleepstudy -
- frmtmb--forward: test-brms-agreement.R globals: +inhaler +kidney -
- frmtmb--forward: test-glmm-binomial.R globals: +cbpp -
- frmtmb--forward: test-multiple-methods.R globals: +nhanes -
- frmtmb--forward: test-parse.R compat features: 4 -> 7 entries
- frmtmb--preload: test-aliased-grouping.R globals: +sleepstudy -
- frmtmb--preload: test-brms-agreement.R globals: +inhaler +kidney -
- frmtmb--preload: test-glmm-binomial.R globals: +cbpp -
- frmtmb--preload: test-multiple-methods.R globals: +nhanes -
- frmtmb--preload: test-parse.R compat features: 12 -> 15 entries
- frmtmb--reverse: test-v19.R globals: +sleepstudy -
- frmtmb--reverse: test-v14.R globals: +cbpp -
- frmtmb--reverse: test-parse.R compat features: 1 -> 4 entries
- frmtmb--reverse: test-multiple-methods.R globals: +nhanes -
- frmtmb--reverse: test-brms-agreement.R globals: +inhaler +kidney -
- frmtmb.ode--forward: test-compat.R compat rules: 1 -> 2 entries
- frmtmb.ode--forward: test-compat.R frame checks: 1 -> 2 entries
- frmtmb.ode--preload: test-compat.R compat rules: 9 -> 10 entries
- frmtmb.ode--preload: test-compat.R frame checks: 3 -> 4 entries
- frmtmb.ode--reverse: test-compat.R compat rules: 1 -> 2 entries
- frmtmb.ode--reverse: test-compat.R frame checks: 1 -> 2 entries
- frmtmb.sample--forward: test-parallel-chains.R globals: +sleepstudy -
- frmtmb.sample--preload: test-compat-preflight.R compat features: 12 -> 13 entries
- frmtmb.sample--preload: test-parallel-chains.R globals: +sleepstudy -
- frmtmb.sample--reverse: test-parallel-chains.R globals: +sleepstudy -

<!-- END GENERATED -->

### After fix round 1 (this worktree, final tree)

<!-- BEGIN GENERATED: dev/testiso-counts.R testiso-after2 -->

From `dev/release/testiso-after2-sessions/` (24 sessions).

| package | order | files | pass | fail | err | skip | warn | own LEAK | frmtmb from |
|---|---|---|---|---|---|---|---|---|---|
| frmtmb | forward | 228/238 | 17447 | 0 | 0 | 185 | 0 | 0 | wt-testiso-lib |
| frmtmb | preload | 228/238 | 17641 | 0 | 0 | 186 | 0 | 0 | wt-testiso-lib |
| frmtmb | reverse | 228/238 | 17545 | 0 | 0 | 185 | 0 | 0 | wt-testiso-lib |
| frmtmb.coupling | forward | 11/11 | 542 | 0 | 0 | 5 | 0 | 0 | wt-testiso-lib |
| frmtmb.coupling | preload | 11/11 | 542 | 0 | 0 | 5 | 0 | 0 | wt-testiso-lib |
| frmtmb.coupling | reverse | 11/11 | 542 | 0 | 0 | 5 | 0 | 0 | wt-testiso-lib |
| frmtmb.eam | forward | 30/30 | 1808 | 0 | 0 | 3 | 0 | 0 | wt-testiso-lib |
| frmtmb.eam | preload | 30/30 | 1808 | 0 | 0 | 3 | 0 | 0 | wt-testiso-lib |
| frmtmb.eam | reverse | 30/30 | 1808 | 0 | 0 | 3 | 0 | 0 | wt-testiso-lib |
| frmtmb.latent | forward | 10/10 | 360 | 0 | 0 | 2 | 0 | 0 | wt-testiso-lib |
| frmtmb.latent | preload | 10/10 | 360 | 0 | 0 | 2 | 0 | 0 | wt-testiso-lib |
| frmtmb.latent | reverse | 10/10 | 360 | 0 | 0 | 2 | 0 | 0 | wt-testiso-lib |
| frmtmb.learn | forward | 15/15 | 430 | 0 | 0 | 13 | 0 | 0 | wt-testiso-lib |
| frmtmb.learn | preload | 15/15 | 430 | 0 | 0 | 13 | 0 | 0 | wt-testiso-lib |
| frmtmb.learn | reverse | 15/15 | 430 | 0 | 0 | 13 | 0 | 0 | wt-testiso-lib |
| frmtmb.ode | forward | 11/11 | 549 | 0 | 0 | 1 | 0 | 0 | wt-testiso-lib |
| frmtmb.ode | preload | 11/11 | 549 | 0 | 0 | 1 | 0 | 0 | wt-testiso-lib |
| frmtmb.ode | reverse | 11/11 | 549 | 0 | 0 | 1 | 0 | 0 | wt-testiso-lib |
| frmtmb.sample | forward | 48/53 | 2697 | 0 | 0 | 4 | 0 | 0 | wt-testiso-lib |
| frmtmb.sample | preload | 48/53 | 2697 | 0 | 0 | 4 | 0 | 0 | wt-testiso-lib |
| frmtmb.sample | reverse | 48/53 | 2697 | 0 | 0 | 4 | 0 | 0 | wt-testiso-lib |
| frmtmb.spline | forward | 16/16 | 594 | 0 | 0 | 1 | 0 | 0 | wt-testiso-lib |
| frmtmb.spline | preload | 16/16 | 594 | 0 | 0 | 1 | 0 | 0 | wt-testiso-lib |
| frmtmb.spline | reverse | 16/16 | 594 | 0 | 0 | 1 | 0 | 0 | wt-testiso-lib |

Failing and erroring tests (expectation count):

- none

LEAK lines a test made itself (not a namespace load):

- none

<!-- END GENERATED -->

The first after run (16 sessions, `testiso-after`, before fix round 1) also read fail=0 err=0 everywhere with no own LEAK line.

### The two files with each extension preloaded

`dev/testiso-compatload.sh` (one fresh session per preload set; logs `dev/testiso-log/compatload-before/` and `compatload-after5/`):

Before:

<!-- BEGIN GENERATED: dev/testiso-compatload.sh before -->

| preloaded | test-compat.R | test-compat-register.R | registries unchanged |
|---|---|---|---|
| all | pass 1018 fail 9 err 0 | pass 91 fail 10 err 1 | TRUE |
| frmtmb.coupling | pass 800 fail 2 err 0 | pass 103 fail 0 err 0 | FALSE |
| frmtmb.eam | pass 900 fail 0 err 0 | pass 96 fail 5 err 1 | FALSE |
| frmtmb.latent | pass 799 fail 0 err 0 | pass 98 fail 5 err 0 | FALSE |
| frmtmb.learn | pass 1008 fail 2 err 0 | pass 96 fail 5 err 1 | TRUE |
| frmtmb.ode | pass 796 fail 3 err 0 | pass 98 fail 5 err 0 | FALSE |
| frmtmb.sample | pass 797 fail 4 err 0 | pass 98 fail 5 err 0 | FALSE |
| frmtmb.spline | pass 801 fail 3 err 0 | pass 103 fail 0 err 0 | FALSE |
| none | pass 797 fail 0 err 0 | pass 103 fail 0 err 0 | FALSE |

<!-- END GENERATED -->

After:

<!-- BEGIN GENERATED: dev/testiso-compatload.sh after5 -->

| preloaded | test-compat.R | test-compat-register.R | registries unchanged |
|---|---|---|---|
| all | pass 993 fail 0 err 0 | pass 103 fail 0 err 0 | TRUE |
| frmtmb.coupling | pass 797 fail 0 err 0 | pass 103 fail 0 err 0 | TRUE |
| frmtmb.eam | pass 895 fail 0 err 0 | pass 103 fail 0 err 0 | TRUE |
| frmtmb.latent | pass 797 fail 0 err 0 | pass 103 fail 0 err 0 | TRUE |
| frmtmb.learn | pass 993 fail 0 err 0 | pass 103 fail 0 err 0 | TRUE |
| frmtmb.ode | pass 797 fail 0 err 0 | pass 103 fail 0 err 0 | TRUE |
| frmtmb.sample | pass 797 fail 0 err 0 | pass 103 fail 0 err 0 | TRUE |
| frmtmb.spline | pass 797 fail 0 err 0 | pass 103 fail 0 err 0 | TRUE |
| none | pass 797 fail 0 err 0 | pass 103 fail 0 err 0 | TRUE |

<!-- END GENERATED -->

## What was wrong, and the fix

### test-compat-register.R

The tests registered `wiener`, `dec()` and `hmm`, which frmtmb.eam and
frmtmb.latent register when they load, and asserted that the session
held no forward reference (`attr(frm_compat(), "unresolved")` NULL),
which is false when frmtmb.sample (expects hmm, lca) or frmtmb.ode
(expects frm_sample) is loaded. The fix:

- test-only names no extension registers: `tfam`, `tterm()`, `thmm`,
  and `xyzzy()` for "an unknown spelling with no near miss";
- the unresolved assertions are relative to `un0`, the session's own
  unresolved set at the start of the test: after the registration it
  is `sort(c(un0, "thmm"))`, after the feature is supplied it is `un0`
  again, and the `frmtmb_compat` class is present exactly when `un0`
  is non-empty;
- `local_registries` is now the shared `local_frmtmb_registries()`
  (`tests/testthat/helper-registry.R`), which also restores the
  vocabulary cache, the frame checks, the prior defaults and the
  notice state.

### test-compat.R

Four tests assumed a core-only vocabulary:

- "every rule wins at least one pair": a contributed rule naming a
  feature another package owns is dead until that package loads
  (frmtmb.sample's `hmm x frm_sample` without frmtmb.latent). Such a
  rule is excused only when it is contributed (row index past
  `compat_core_rule_count()`) AND one of its sides is in
  `compat_unresolved_sides()`. A core rule is never excused, and
  neither is a contributed rule dead for another reason.
- "rule patterns name real kinds, groups, and features": now checks
  core's rules only. A contributed side is refused at registration
  unless it is in the vocabulary or in `expects =`
  (test-compat-register.R pins that), so the check on contributed
  sides was a load-order test, not a check of the table.
- "declared specials still parse": pins core's eight specials exactly,
  over the non-contributed rows; frmtmb.ode adds `frm_ode()` and
  `frm_lincmt()`.
- "declared post-fit methods are callable": core's rows only, on the
  search path as before. A contributed method is in a namespace an
  earlier file may load without attaching, and one row names an
  argument setting: frmtmb.spline's `frm_curve_contrast` is the
  `contrast` argument of `frm_curve()`, not a function. The
  contributed rows are checked in test-compat-extensions.R (fix rounds
  1 and 2).

`compat_contributed_names()` and `compat_core_rule_count()` are test
helpers at the top of the file. The second relies on contributed rules
being appended after core's, which `R/compat.R` documents ("contributor
seam").

The pass count of test-compat.R is unchanged with no extension (797)
and falls with extensions loaded (1018 to 993 with all seven), because
the pattern loop and the methods loop no longer assert on contributed
rows.

### Mutation checks

`dev/testiso-mutants.R` replaces one function in the frmtmb namespace,
or registers a planted rule set, and runs one file; `dev/testiso-mutants.sh`
runs each mutant with no preload, with frmtmb.sample, frmtmb.eam,
frmtmb.ode, and with all seven. Log: `dev/testiso-log/mutants-r2.txt`.

<!-- BEGIN GENERATED: dev/testiso-mutants.sh r2 -->

mutant: none preloaded: none file: test-compat.R pass 797 fail 0 err 0 
mutant: none preloaded: frmtmb.sample file: test-compat.R pass 797 fail 0 err 0 
mutant: none preloaded: frmtmb.eam file: test-compat.R pass 895 fail 0 err 0 
mutant: none preloaded: frmtmb.ode file: test-compat.R pass 797 fail 0 err 0 
mutant: none preloaded: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline file: test-compat.R pass 993 fail 0 err 0 
mutant: none_reg preloaded: none file: test-compat-register.R pass 103 fail 0 err 0 
mutant: none_reg preloaded: frmtmb.sample file: test-compat-register.R pass 103 fail 0 err 0 
mutant: none_reg preloaded: frmtmb.eam file: test-compat-register.R pass 103 fail 0 err 0 
mutant: none_reg preloaded: frmtmb.ode file: test-compat-register.R pass 103 fail 0 err 0 
mutant: none_reg preloaded: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline file: test-compat-register.R pass 103 fail 0 err 0 
mutant: contrib_dead preloaded: none file: test-compat.R pass 793 fail 4 err 0 
  FAILS [no rule is stated twice for the same unordered pattern pair]
  FAILS [every rule wins at least one pair]
  FAILS [no pair is decided by file order alone]
  FAILS [a more specific rule beats a broader one]
mutant: contrib_dead preloaded: frmtmb.sample file: test-compat.R pass 793 fail 4 err 0 
  FAILS [no rule is stated twice for the same unordered pattern pair]
  FAILS [every rule wins at least one pair]
  FAILS [no pair is decided by file order alone]
  FAILS [a more specific rule beats a broader one]
mutant: contrib_dead preloaded: frmtmb.eam file: test-compat.R pass 891 fail 4 err 0 
  FAILS [no rule is stated twice for the same unordered pattern pair]
  FAILS [every rule wins at least one pair]
  FAILS [no pair is decided by file order alone]
  FAILS [a more specific rule beats a broader one]
mutant: contrib_dead preloaded: frmtmb.ode file: test-compat.R pass 793 fail 4 err 0 
  FAILS [no rule is stated twice for the same unordered pattern pair]
  FAILS [every rule wins at least one pair]
  FAILS [no pair is decided by file order alone]
  FAILS [a more specific rule beats a broader one]
mutant: contrib_dead preloaded: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline file: test-compat.R pass 989 fail 4 err 0 
  FAILS [no rule is stated twice for the same unordered pattern pair]
  FAILS [every rule wins at least one pair]
  FAILS [no pair is decided by file order alone]
  FAILS [a more specific rule beats a broader one]
mutant: core_unresolved preloaded: none file: test-compat.R pass 796 fail 2 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
mutant: core_unresolved preloaded: frmtmb.sample file: test-compat.R pass 796 fail 2 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
mutant: core_unresolved preloaded: frmtmb.eam file: test-compat.R pass 894 fail 2 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
mutant: core_unresolved preloaded: frmtmb.ode file: test-compat.R pass 796 fail 2 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
mutant: core_unresolved preloaded: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline file: test-compat.R pass 992 fail 2 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
mutant: core_special preloaded: none file: test-compat.R pass 794 fail 3 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
  FAILS [declared specials still parse]
mutant: core_special preloaded: frmtmb.sample file: test-compat.R pass 794 fail 3 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
  FAILS [declared specials still parse]
mutant: core_special preloaded: frmtmb.eam file: test-compat.R pass 892 fail 3 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
  FAILS [declared specials still parse]
mutant: core_special preloaded: frmtmb.ode file: test-compat.R pass 794 fail 3 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
  FAILS [declared specials still parse]
mutant: core_special preloaded: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline file: test-compat.R pass 990 fail 3 err 0 
  FAILS [every rule wins at least one pair]
  FAILS [rule patterns name real kinds, groups, and features]
  FAILS [declared specials still parse]
mutant: core_method preloaded: none file: test-compat.R pass 797 fail 1 err 0 
  FAILS [declared post-fit methods are callable]
mutant: core_method preloaded: frmtmb.sample file: test-compat.R pass 797 fail 1 err 0 
  FAILS [declared post-fit methods are callable]
mutant: core_method preloaded: frmtmb.eam file: test-compat.R pass 895 fail 1 err 0 
  FAILS [declared post-fit methods are callable]
mutant: core_method preloaded: frmtmb.ode file: test-compat.R pass 797 fail 1 err 0 
  FAILS [declared post-fit methods are callable]
mutant: core_method preloaded: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline file: test-compat.R pass 993 fail 1 err 0 
  FAILS [declared post-fit methods are callable]
mutant: unresolved_drop preloaded: none file: test-compat-register.R pass 99 fail 4 err 0 
  FAILS [a declared but unresolved rule side is reported, not dropped]
mutant: unresolved_drop preloaded: frmtmb.sample file: test-compat-register.R pass 99 fail 4 err 0 
  FAILS [a declared but unresolved rule side is reported, not dropped]
mutant: unresolved_drop preloaded: frmtmb.eam file: test-compat-register.R pass 99 fail 4 err 0 
  FAILS [a declared but unresolved rule side is reported, not dropped]
mutant: unresolved_drop preloaded: frmtmb.ode file: test-compat-register.R pass 99 fail 4 err 0 
  FAILS [a declared but unresolved rule side is reported, not dropped]
mutant: unresolved_drop preloaded: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline file: test-compat-register.R pass 99 fail 4 err 0 
  FAILS [a declared but unresolved rule side is reported, not dropped]
mutant: commit_early preloaded: none file: test-compat-register.R pass 101 fail 2 err 0 
  FAILS [a refused registration leaves the vocabulary as it found it]
mutant: commit_early preloaded: frmtmb.sample file: test-compat-register.R pass 101 fail 2 err 0 
  FAILS [a refused registration leaves the vocabulary as it found it]
mutant: commit_early preloaded: frmtmb.eam file: test-compat-register.R pass 101 fail 2 err 0 
  FAILS [a refused registration leaves the vocabulary as it found it]
mutant: commit_early preloaded: frmtmb.ode file: test-compat-register.R pass 101 fail 2 err 0 
  FAILS [a refused registration leaves the vocabulary as it found it]
mutant: commit_early preloaded: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline file: test-compat-register.R pass 101 fail 2 err 0 
  FAILS [a refused registration leaves the vocabulary as it found it]
mutant: aterm_nofeature preloaded: none file: test-compat-register.R pass 74 fail 19 err 3 
  FAILS [a registered addition term becomes a compat feature]
  FAILS [one display name carries one kind]
  FAILS [a refused frmtmb_register_aterm() registers nothing at all]
  FAILS [the refusal reads in the term registrant's own terms]
  FAILS [a kind beginning with a vowel takes the right article]
  FAILS [four registrations cost one vocabulary build]
  FAILS [the cached vocabulary is the vocabulary]
  FAILS [the cache is validated against the registry, not a flag]
  FAILS [a refused registration leaves the cache as it found it]
mutant: aterm_nofeature preloaded: frmtmb.sample file: test-compat-register.R pass 74 fail 19 err 3 
  FAILS [a registered addition term becomes a compat feature]
  FAILS [one display name carries one kind]
  FAILS [a refused frmtmb_register_aterm() registers nothing at all]
  FAILS [the refusal reads in the term registrant's own terms]
  FAILS [a kind beginning with a vowel takes the right article]
  FAILS [four registrations cost one vocabulary build]
  FAILS [the cached vocabulary is the vocabulary]
  FAILS [the cache is validated against the registry, not a flag]
  FAILS [a refused registration leaves the cache as it found it]
mutant: aterm_nofeature preloaded: frmtmb.eam file: test-compat-register.R pass 74 fail 19 err 3 
  FAILS [a registered addition term becomes a compat feature]
  FAILS [one display name carries one kind]
  FAILS [a refused frmtmb_register_aterm() registers nothing at all]
  FAILS [the refusal reads in the term registrant's own terms]
  FAILS [a kind beginning with a vowel takes the right article]
  FAILS [four registrations cost one vocabulary build]
  FAILS [the cached vocabulary is the vocabulary]
  FAILS [the cache is validated against the registry, not a flag]
  FAILS [a refused registration leaves the cache as it found it]
mutant: aterm_nofeature preloaded: frmtmb.ode file: test-compat-register.R pass 74 fail 19 err 3 
  FAILS [a registered addition term becomes a compat feature]
  FAILS [one display name carries one kind]
  FAILS [a refused frmtmb_register_aterm() registers nothing at all]
  FAILS [the refusal reads in the term registrant's own terms]
  FAILS [a kind beginning with a vowel takes the right article]
  FAILS [four registrations cost one vocabulary build]
  FAILS [the cached vocabulary is the vocabulary]
  FAILS [the cache is validated against the registry, not a flag]
  FAILS [a refused registration leaves the cache as it found it]
mutant: aterm_nofeature preloaded: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline file: test-compat-register.R pass 74 fail 19 err 3 
  FAILS [a registered addition term becomes a compat feature]
  FAILS [one display name carries one kind]
  FAILS [a refused frmtmb_register_aterm() registers nothing at all]
  FAILS [the refusal reads in the term registrant's own terms]
  FAILS [a kind beginning with a vowel takes the right article]
  FAILS [four registrations cost one vocabulary build]
  FAILS [the cached vocabulary is the vocabulary]
  FAILS [the cache is validated against the registry, not a flag]
  FAILS [a refused registration leaves the cache as it found it]

<!-- END GENERATED -->

### The other leaks the audit found

Each was measured by `dev/testiso-filecheck.R`, which loads all seven
extensions first (so the file's own loads register nothing new) and
compares the registries before and after one test file
(`dev/testiso-log/filecheck.txt`):

<!-- BEGIN GENERATED: dev/testiso-filecheck.R; base = main checkout, after = worktree; the last two lines with TESTISO_SEED_NOTICE=1, base then after -->

```
base: test-compat-preflight.R pass 86 fail 0 err 0 | changed: features 
  features: 12 -> 13
base: test-compat.R pass 15 fail 0 err 0 | changed: rules checks 
  rules: 9 -> 10
  checks: 3 -> 4
base: test-parse.R pass 44 fail 0 err 0 | changed: features 
  features: 12 -> 15
base: test-cens-trunc.R pass 69 fail 0 err 0 | changed: nothing 
base: test-compat-register.R pass 91 fail 10 err 1 | changed: nothing 
base: test-compat.R pass 1018 fail 9 err 0 | changed: nothing 
after: test-compat-preflight.R pass 86 fail 0 err 0 | changed: nothing 
after: test-compat.R pass 15 fail 0 err 0 | changed: nothing 
after: test-parse.R pass 44 fail 0 err 0 | changed: nothing 
after: test-cens-trunc.R pass 69 fail 0 err 0 | changed: nothing 
after: test-compat-register.R pass 103 fail 0 err 0 | changed: nothing 
after: test-compat.R pass 993 fail 0 err 0 | changed: nothing 
test-cens-trunc.R pass 69 fail 0 err 0 | changed: notices 
  notices: 1 -> 0
test-cens-trunc.R pass 69 fail 0 err 0 | changed: nothing 
```

<!-- END GENERATED -->

- test-parse.R: `on.exit(frmtmb_aterm_registry$reg[["side"]] <- NULL)`
  deleted by name and left the `side()`, `pair()` and `tag()` compat
  features, because `frmtmb_register_aterm()` also writes a feature.
  In the reverse-order before run this failed test-compat.R's
  "declared addition terms are the ones the parser accepts" three
  times. Now `local_frmtmb_registries()`.
- test-cens-trunc.R `arm_cens_notice()`: cleared the notice state at
  both ends, so the notices earlier files had given were lost. Now it
  saves and restores (seeded check: 1 notice before, 0 after on base).
- helper-rl.R: sourcing `inst/rl/rw-delta.R` registers `reward()` for
  the whole run; the registries are now restored at the run's end
  through `teardown_env()` when `testthat::is_testing()`.
- frmtmb.ode test-compat.R: re-running `.onLoad()` appended a second
  compat rule and a SECOND FRAME CHECK, which then ran twice on every
  ODE fit for the rest of the session. Now saves and restores both
  (base R `on.exit()`, because frmtmb.ode does not suggest withr).
- frmtmb.sample test-compat-preflight.R: registered a bare `frm_ode`
  structure for the rest of the session. Now restored with
  `withr::defer()`. The one-session runner could not see this one:
  the same file loads frmtmb.ode, so the feature count change was
  attributed to that load (a limitation of the LEAK rule, below).
- `data(x, package = )` without `envir =` writes to the global
  environment. 39 calls in 19 core files and one frmtmb.sample file now
  pass `envir = environment()`. test-confint-anova.R needed the data at
  file level instead, because `update()` evaluates the stored call,
  which names `sleepstudy`, from where it is called
  (`anova(refit = TRUE)` refits from the stored frame and does not; the
  first version of this sentence named both, review finding 7). With
  the data inside the file's `local()` block the test's own `update()`
  calls did not find it: 3 errors in `testiso-data` below, fixed and
  rerun as `testiso-ca`.

Every changed file was run one per process (`dev/release/run-par.sh`).
`testiso-data` ran after the `data()` and compat edits; the test-parse,
test-cens-trunc, helper-rl and extension edits came after it and are in
`testiso-misc`. The one-session after arm ran the final tree.

<!-- BEGIN GENERATED: dev/release/run-par.sh; logs dev/release/testiso-*.log -->

`testiso-data` (gates `NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS`):

```
RESULT test-parallel-chains.R pass=4 fail=0 err=0 skip=0 warn=0
RESULT test-aliased-grouping.R pass=33 fail=0 err=0 skip=0 warn=0
RESULT test-boot.R pass=49 fail=0 err=0 skip=0 warn=0
RESULT test-brms-agreement.R pass=191 fail=0 err=0 skip=0 warn=0
RESULT test-brms-likelihood.R pass=543 fail=0 err=0 skip=0 warn=0
RESULT test-brms-priors.R pass=103 fail=0 err=0 skip=0 warn=0
RESULT test-ce-bands.R pass=169 fail=0 err=0 skip=0 warn=0
RESULT test-compat-register.R pass=103 fail=0 err=0 skip=0 warn=0
RESULT test-compat.R pass=797 fail=0 err=0 skip=0 warn=0
RESULT test-confint-anova.R pass=23 fail=0 err=3 skip=0 warn=0
RESULT test-glmm-binomial.R pass=10 fail=0 err=0 skip=0 warn=0
RESULT test-glmm-gaussian.R pass=16 fail=0 err=0 skip=0 warn=0
RESULT test-interop.R pass=45 fail=0 err=0 skip=0 warn=0
RESULT test-lkj.R pass=212 fail=0 err=0 skip=0 warn=0
RESULT test-methods.R pass=65 fail=0 err=0 skip=0 warn=0
RESULT test-multiple-methods.R pass=113 fail=0 err=0 skip=0 warn=0
RESULT test-open-issues.R pass=46 fail=0 err=0 skip=0 warn=0
RESULT test-sugar.R pass=45 fail=0 err=0 skip=0 warn=0
RESULT test-v07.R pass=25 fail=0 err=0 skip=0 warn=0
RESULT test-v11.R pass=22 fail=0 err=0 skip=0 warn=0
RESULT test-v14.R pass=81 fail=0 err=0 skip=0 warn=0
RESULT test-v19.R pass=23 fail=0 err=0 skip=0 warn=0
TESTISO-DATA ran 22 of 22
```

`testiso-ca` (gates `NOT_CRAN`):

```
RESULT test-confint-anova.R pass=35 fail=0 err=0 skip=0 warn=0
TESTISO-CA ran 1 of 1
```

`testiso-misc` (gates `NOT_CRAN`):

```
RESULT test-cens-trunc.R pass=69 fail=0 err=0 skip=0 warn=0
RESULT test-parse.R pass=44 fail=0 err=0 skip=0 warn=0
RESULT test-rl-example.R pass=96 fail=0 err=0 skip=2 warn=0
RESULT test-compat.R pass=15 fail=0 err=0 skip=0 warn=0
RESULT test-compat-preflight.R pass=86 fail=0 err=0 skip=0 warn=0
TESTISO-MISC ran 5 of 5
```


<!-- END GENERATED -->

## The one-session runner

`dev/release/run-session.R <package> [forward|reverse|preload] [test dir]`
runs one suite through testthat's own `test_files()` with the file list
testthat would use (`find_test_scripts()`), reversed for `reverse`, so
helpers and setup files run once per session as under covr and
`R CMD check`. `preload` is alphabetical order with all seven
extensions loaded before the first file. It requires all eight
packages installed and prints which library each came from. After each
file it compares frmtmb's five registries (compat features and rules,
addition terms, frame checks, prior defaults), the notice tags (a
REMOVED tag is a leak; a fit adds one legitimately), `options()`, the
environment variables, the search path and the global environment, and
prints a LEAK line for a difference. It reads each registry with
`get(inherits = FALSE)` and stops if one is not an environment, so a
renamed registry cannot compare NULL with NULL. A namespace load is
printed as LOAD and cannot be undone, so a registry difference in a
file that loaded an frmtmb namespace, or an option or variable
difference in a file that loaded any namespace, is marked
`[with load of ...]`. A global or a search-path difference is never
marked: a load creates no global and attaches nothing (review finding
2; the first version marked them). In `preload` no file loads an frmtmb
namespace, so every registry LEAK there is the test's own; that closes
the limitation the first version had (the sample preflight case
above). `dev/testiso-filecheck.R` is the instrument for one file.

`dev/release/run-sessions.sh <tier> [N] [package ...]` runs 24
sessions (three orders per suite), N at a time (default 8; a non-number
N is refused), like `run-par.sh`, and writes
`dev/release/<tier>.log` with each SESSION line, its failures and its
LEAK lines. `TEST_ROOT` reads the tests from another checkout (the
before arm used the main checkout). A full run takes about 21 minutes
of wall time here: core 1214 s forward and 1249 s reverse, frmtmb.eam
585 s, frmtmb.sample 410 s, the rest under 2 minutes each (sums of the
TIME lines of the before arm).

The runner sets only what the caller sets: these runs used
`NOT_CRAN=true` and no gate, as the coverage job does. 10 core files
and 5 frmtmb.sample files skip at file level under that setting and
have no result row, which is the `files=226/236` and `48/53`.

## CI: the coverage job ran against stale extensions

Corrected in fix round 1 from the review (`dev/reviews/2026-10-14-testiso.md`,
finding 6), which read the run logs and `gh cache list`. The first
version of this section said the `check-frmtmb-*` jobs seeded the
cache; that was wrong.

Run 38563139754 (0.71.0, the failing one) restored its R library from
the cache key `Ubuntu 24.04.5 LTS-R version 4.6.1 (2026-06-24)-x86_64-1-664bc0b8...`
through the restore-key prefix, not its own key (`...-858c7898...`).
That library held `frmtmb 0.69.0`, `frmtmb.coupling 0.6.0`,
`frmtmb.eam 0.11.2`, `frmtmb.latent 0.6.0`, `frmtmb.learn 0.7.2`,
`frmtmb.ode 0.7.2`, `frmtmb.sample 0.17.0` and `frmtmb.spline 0.10.1`,
all "local", installed 2026-10-13 (the session-info block of the
setup-r-dependencies step). The checkout had eam 0.12.0, latent 0.6.1,
sample 0.19.0 and spline 0.10.2.

The mechanism:

1. `setup-r-dependencies` keys the cache on the OS, the R version, a
   cache version and a hash of the resolved dependencies, and restores
   by the prefix up to the cache version. A job whose exact key misses
   restores the newest cache with that prefix, from any workflow. The
   prefix is the same for every ubuntu-24.04 R release job.
2. `pkgdown.yaml` installs all seven extensions from the checkout into
   that cached library (its "Install the seven extensions" step), so
   its saved cache holds them. The `check-frmtmb-*` jobs install core
   and only the siblings their package needs, and `check-r-package`
   installs the checked package into a temporary library; only
   pkgdown's library holds coupling, learn and spline.
3. Coverage carries the set forward. Per the review: the cache 38563139754
   restored (`...-664bc0b8...`, created 03:24:41Z) was saved by the
   previous coverage run, 38502130389, which had restored
   `...-4e42e1a1...`, saved by coverage run 38457152845; run
   38579297702 restored `...-fd87d7bf...`, which matches pkgdown run
   38563139900. So pkgdown seeds the library and coverage keeps it.
4. The root DESCRIPTION names no extension (not in Suggests, no
   Remotes), so pak neither installs nor updates them. covr installs
   core from the checkout into its own path, so core is current; the
   extensions are whatever the restored library holds.
5. Run 38579297702 saved coverage's EXACT key (`...-858c7898...`). The
   next run with the same dependency lock hits that key exactly and
   restores 47cf659d's extensions until the lock changes. So a separate
   `cache-version` alone does not end the staleness; only removing the
   extensions, or installing them from the checkout, does.

31 core test files name an extension (most of them frmtmb.sample). In
run 38563139754 they ran against frmtmb.sample 0.17.0, where every
other suite used 0.19.0. The 0.71.0 failure itself did not depend on
the versions: any frmtmb.eam registers `dec()`.

Why run 38579297702 (47cf659d) passed, `[ FAIL 0 | WARN 0 | SKIP 206 |
PASS 17379 ]`, although test-compat.R fails locally in alphabetical
order (before arm, 9 failures): locally test-arg-refusal.R loads
frmtmb.sample before test-compat.R sorts, through a source scan that
needs the package's `dev/` tree. On CI that test skips ("dev/ not
reachable from the test dir"), so nothing loads an extension before
test-families-extreme.R, which sorts after both compat files. The
compat tests were order-dependent on CI too and passed by sort order.

The ubuntu-24.04 (release) job of `R-CMD-check.yaml` on the same commit
(run 38563140196) restored the same cache and listed the same eight
stale packages, and passed. It runs `rcmdcheck(args = c("--no-manual",
"--as-cran"))`, and under `--as-cran` `tools:::.check_packages` sets
its internal flag `R_check_suggests_only` (the review read this in R
4.6.1's source), so the tests run with a library that holds only the
packages the DESCRIPTION names; frmtmb.eam cannot load there. covr has
no such mask.

### The change made (decided by the user)

Coverage measures core ALONE, and each extension has its own coverage
run, because they are separate packages. `.github/workflows/test-coverage.yaml`:

- **Core job** (`test-coverage`): its own `cache-version:
  coverage-core-1`, and a step after `setup-r-dependencies` that
  removes every installed `frmtmb*` package (and stops if one is left).
  covr installs core from the checkout itself. The tests that need an
  extension skip, as under R CMD check. Upload flag `frmtmb`.
- **Extension jobs** (`extension-coverage`, a matrix over the seven,
  `fail-fast: false`): `extra-repositories` for RTMBode and
  `libgsl-dev`, copied from `pkgdown.yaml`; `cache-version:
  coverage-<pkg>-1`; one step that removes every `frmtmb*` package from
  the restored library, installs core from the checkout, then the
  siblings the package names (each with its Depends and Imports), then
  the package's own Imports and Suggests, and prints the frmtmb
  versions installed; `covr::package_coverage("extensions/<pkg>")`;
  the cobertura file names rewritten from `R/` to
  `extensions/<pkg>/R/`, because core has files of the same names;
  upload with flag `<pkg>` when `CODECOV_TOKEN` exists, else the
  coverage is printed in the log. Both jobs keep the 90-minute timeout
  and the "Show testthat output" step. Since fix round 2 the matrix job
  also sets `FRMTMB_SAMPLER_GATES: "false"`, as the other CI jobs of
  frmtmb.sample do.
- Siblings, from each DESCRIPTION (`dev/testiso-yaml.R`, which also
  parses the YAML with `yaml::read_yaml()` and every Rscript step with
  `parse()`; log `dev/testiso-log/yaml.txt`): eam gets frmtmb.sample
  (Suggests); learn gets frmtmb.eam (Imports); sample gets
  frmtmb.latent, frmtmb.ode and frmtmb.eam (Suggests); coupling,
  latent, ode and spline name none. For every entry the set equals the
  frmtmb packages the DESCRIPTION names, and each sibling's own frmtmb
  imports come before it.

<!-- BEGIN GENERATED: dev/testiso-yaml.R -->

```
jobs: test-coverage extension-coverage 
test-coverage / Remove every frmtmb package from the library : parses 
test-coverage / Measure coverage : parses 
test-coverage cache-version: coverage-core-1 
extension-coverage / Install core and the siblings from this checkout : parses 
extension-coverage / Measure coverage : parses 
extension-coverage cache-version: coverage-${{ matrix.pkg }}-1 
matrix: frmtmb.coupling frmtmb.eam frmtmb.latent frmtmb.learn frmtmb.ode frmtmb.sample frmtmb.spline 
frmtmb.coupling | names:  | siblings:  | same set: TRUE | order ok: TRUE 
frmtmb.eam | names: frmtmb.sample | siblings: frmtmb.sample | same set: TRUE | order ok: TRUE 
frmtmb.latent | names:  | siblings:  | same set: TRUE | order ok: TRUE 
frmtmb.learn | names: frmtmb.eam | siblings: frmtmb.eam | same set: TRUE | order ok: TRUE 
frmtmb.ode | names:  | siblings:  | same set: TRUE | order ok: TRUE 
frmtmb.sample | names: frmtmb.eam frmtmb.latent frmtmb.ode | siblings: frmtmb.latent frmtmb.ode frmtmb.eam | same set: TRUE | order ok: TRUE 
frmtmb.spline | names:  | siblings:  | same set: TRUE | order ok: TRUE 
```

<!-- END GENERATED -->

**Only a push tests this workflow.** It is parsed and its R steps
parse, but no part of it has run on a runner; the user pushes.

**For the user, from the review:** before this change, coverage was the
only CI job that ran core's suite in one session with extensions
present, which is the only CI place where this lane's class of defect
(a test that depends on what an earlier file loaded) could show. Core
alone no longer has that run on CI. The local release gate covers it:
`dev/release/run-sessions.sh` runs every suite forward, reverse and
with all seven extensions preloaded (`dev/round-handoff.md`,
"Consolidation"). The extension jobs run each extension's suite in one
session with its named siblings only.

## R CMD check --as-cran

`dev/testiso-check.sh` (built with vignettes in `dev/testiso-check/`,
`R_LIBS` = wt-testiso-lib, rellib-r9, user library;
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`; logs
`dev/testiso-check/check-*.txt`):

- frmtmb: Status: 1 NOTE (HTML manual, "package 'V8' unavailable", the
  known environmental NOTE). Tests `[ FAIL 0 | WARN 0 | SKIP 208 |
  PASS 17375 ]`, 21 minutes.
- frmtmb.ode: Status: 1 NOTE (the same V8 NOTE). Tests `[ FAIL 0 |
  WARN 0 | SKIP 2 | PASS 545 ]`.
- frmtmb.sample: Status: OK. Tests `[ FAIL 0 | WARN 0 | SKIP 12 |
  PASS 2683 ]`.

After fix round 1 (`bash dev/testiso-check.sh frmtmb.spline frmtmb`;
the ode and sample checks above are from before the round, whose
changes to them were none):

- frmtmb.spline: Status: 1 NOTE (the V8 NOTE). Tests `[ FAIL 0 | WARN 0
  | SKIP 3 | PASS 586 ]`. This checked round 1's rename, which fix
  round 2 reverted; frmtmb.spline is again identical to its base, so
  it was not checked again.
- frmtmb: Status: 1 NOTE (the V8 NOTE). Tests `[ FAIL 0 | WARN 0 |
  SKIP 210 | PASS 17377 ]`, two skips more than before the round. The
  two new tests that need the extensions (test-compat-extensions.R and
  the child test of test-registry-helper.R) skip when R CMD check
  hides them, which accounts for two; the per-test skip list was not
  kept, so that attribution is from the code, not the log.

## Not done, and found but not fixed

- Core's own argument-setting rows (`confint_profile` and the others)
  are checked by the map in test-compat.R, not by the
  `<function>_<argument>` rule the contract applies to contributed
  rows: core's settings name a value (`profile`, `osa`), not a formal
  argument, so the rule does not fit them.
- The reviewer's R3 note (no CI job runs test-compat-extensions.R; a
  step in pkgdown.yaml could) and the other CI suggestions are for the
  user to decide; this lane added no CI step beyond R2.
- The reviewer's option of one registry list in `R/`, so that the
  frmtmb.ode and frmtmb.sample tests stop carrying their own
  save-and-restore code, was not taken: it adds package code for tests
  only, and the two copies are a few lines each.
- **A reload duplicates registrations.** `frmtmb_register_frame_check()`
  and `frmtmb_register_prior_defaults()` append without a key, so a
  `pkgload::load_all()` of an extension (which re-runs `.onLoad()`)
  makes each of its frame checks run twice per fit for the rest of the
  session. Compat rules duplicate too (harmless, as the ode test
  notes). A user who reloads in development pays only time, so it is
  filed, not fixed: a fix needs a key per registrant (package and
  function name) and is package behavior outside this lane.
- **test-brms-suite-emmeans.R** calls `library(emmeans)` at file level
  in the gated tier, which attaches emmeans for every later file of a
  gated one-session run. The file is generated by `dev/brmsport-gen.R`,
  which belongs to the brmsport ledger; not edited here.
- **Once-per-session state in extensions**: frmtmb.ode's `ode_warned`
  (a warning raised once per key per session) has no test seam. No
  test asserts that warning today; a test that did would depend on
  file order.
- **RNG state** is shared across files in one session. Both orders pass
  on every suite, so no test depends on it today; nothing guards it.
- The one-session runs used `NOT_CRAN=true` only. The gated tiers
  (`FRMTMB_BRMS_FIT_TESTS`, fuzz, scale) were not run in one session;
  they are not what covr or R CMD check runs.

## Files

Fix round 1 added: `tests/testthat/helper-compat-contract.R`,
`test-compat-extensions.R` and `test-registry-helper.R` (new);
`extensions/frmtmb.spline/R/zzz.R` and `NEWS.md` (reverted in fix
round 2: no frmtmb.spline file differs from the base now);
`.github/workflows/test-coverage.yaml`; and the lane scripts
`dev/testiso-contract-mutants.R`, `dev/testiso-yaml.R`,
`dev/testiso-extdeps.R`, `dev/testiso-results.R`,
`dev/testiso-splice.R` and `dev/testiso-paste.R` (the last three write
generated blocks into this file).

Tests: `tests/testthat/helper-registry.R` (new), `helper-rl.R`,
`test-compat.R`, `test-compat-register.R`, `test-parse.R`,
`test-cens-trunc.R`, `test-confint-anova.R`, and the `data()` calls in
`test-aliased-grouping.R`, `test-boot.R`, `test-brms-agreement.R`,
`test-brms-likelihood.R`, `test-brms-priors.R`, `test-ce-bands.R`,
`test-glmm-binomial.R`, `test-glmm-gaussian.R`, `test-interop.R`,
`test-lkj.R`, `test-methods.R`, `test-multiple-methods.R`,
`test-open-issues.R`, `test-sugar.R`, `test-v07.R`, `test-v11.R`,
`test-v14.R`, `test-v19.R`; `extensions/frmtmb.ode/tests/testthat/test-compat.R`;
`extensions/frmtmb.sample/tests/testthat/test-compat-preflight.R` and
`test-parallel-chains.R`.

Release harness: `dev/release/run-session.R`, `dev/release/run-sessions.sh`
(new); `dev/round-handoff.md`, `dev/lane-rules.md`.

Lane scripts: `dev/testiso-compatload.R` and `.sh` (the coordinator's
`compatload.R`, copied and extended), `dev/testiso-mutants.R` and `.sh`,
`dev/testiso-filecheck.R`, `dev/testiso-counts.R`, `dev/testiso-check.sh`,
and `dev/testiso-envs.R` (lists every top-level environment in the eight
namespaces, which is how the registry list was checked as complete:
frmtmb has five registries, four caches and the notice state;
frmtmb.eam a cache, frmtmb.ode two once-per-session stores,
frmtmb.sample one call-state store).
