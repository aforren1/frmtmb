# Lane codesize: where core's size is, and what would reduce it

Round of 2026-11-16, base b255350e (frmtmb 0.83.0 plus CI files). A
measurement lane: no package code changed. `dev/codesize-findings.md` is
built from `dev/codesize-findings-src.md` by `dev/codesize-assemble.R`,
which pastes each script's generated block from its log. Logs and dumps
(`dev/codesize-*.log`, `dev/codesize-out/`, `dev/codesize-cov/`) are local
and gitignored; the blocks below are the record.

## Scripts

| Script | Answers | Output |
|---|---|---|
| `dev/codesize-parse.R` | parses core `R/`, `tests/`, `extensions/`, `inst/`; holds the subsystem map | `dev/codesize-out/parsed.rds` |
| `dev/codesize-inventory.R` | Q1 | `files.tsv`, `subsystems.tsv`, `functions.tsv`, `rd.tsv` |
| `dev/codesize-comments.R`, `dev/codesize-roxtags.R` | Q1, Q7: what the comments hold | logs |
| `dev/codesize-deadcode.R` | Q2 (static), Q4 (single use) | `callgraph.tsv`, `callgraph-novalue.tsv` |
| `dev/codesize-covr-env.R`, `dev/codesize-covsum.R` | Q2 (coverage) | `dev/codesize-cov/cov.rds`, `results.tsv`; `coverage-*.tsv` |
| `dev/codesize-dup.R` | Q3 | `dup-files.tsv`, `dup-regions.tsv` |
| `dev/codesize-guard.R` | Q4 (guard against modeling) | log |
| `dev/codesize-growth.sh`, `dev/codesize-growth.R` | Q5 | `growth.tsv`, `releases.txt` |
| `dev/codesize-candidates.R` | Q6 | log |

Run order: parse, inventory, then the others in any order; covsum and
candidates after the coverage run. One R process each; the coverage run
was the second process while the static scripts ran.

The subsystem map is in `dev/codesize-parse.R` (`sub_map`). The nine
subsystems of the brief are there; the rest of core is in six more groups
(fit engine and starts, inference extras, curve inference, model terms,
post-fit methods, infrastructure) so that every file has one home.
`fit-end.R` is counted under the checks.

## Summary

1. Size: of 118,169 lines, 64,884 are code, 39,065 roxygen (15,980 on
   exported objects, 23,085 on internal ones) and 10,989 other comments.
   The internal roxygen makes no Rd page (all `@noRd`; 5 of 133 Rd
   pages are internal, one the package page); it is prose. The largest
   subsystems are the families (16,705 lines), the fit engine and starts
   (14,962), prediction and draws (13,651) and parse and frame (11,389).
2. Dead code: almost none. One function is unreachable (9 code lines),
   8 more are reached only from tests (79). The suite covers 94.4% of
   relevant lines; about 25 internal functions (about 300 code lines)
   have no test once covr's registry artifact is set aside.
3. Duplication: low. 5.3% of code lines lie in a repeated 60-token
   window, a third of those in one data table (`compat_hand_rules_tbl()`);
   3.7% elsewhere. Conditions go through `frm_stop()`; what costs lines
   is message text (10.2% of code lines are inside a condition call).
4. Fragmentation: 990 internal functions have one call site, 254 of
   them under 10 lines (1,528 code lines). Guard and check code is about
   22.5% of function code (14,202 lines).
5. Growth: 46,171 to 118,169 lines from 0.60.0 to 0.83.0. Code,
   roxygen and comments grew in step. The fit-end check files grew from
   44 lines to 8,372; the fit engine x4.3. Drivers: brms features, and
   rules for how a fit ends at a limit, each with its message, proof
   conditions and BREAKING note.
6. Moves: the rule allows at most `importance.R` (1,878, after a seam)
   and, weakly, `spectral.R` (654) out of `R/`; the clearest domain
   material (cognitive-modeling and RL examples, tests and vignettes)
   is outside `R/`.
7. Reduction: mechanical and low-risk work removes 10% to 15% (ranked
   list, Q7). A large cut needs a decision on the fit-end judgement
   layer.

## Q1. Where the size is

A line is "code" if it holds any token that is not a comment. A roxygen
block belongs to the next top-level expression; it is "exported" when
that object is exported or a registered S3 method, or the block is a
topic page, and "internal" otherwise (including every `@noRd` block).
Byte counts on disk include CRLF line ends (5269 KB); the git blobs are
5154 KB (Q5 block).

<!-- BLOCK codesize-inventory -->

The internal roxygen makes no Rd page: 2,163 of the 2,184 internal
functions have a roxygen block, and every one of those blocks has
`@noRd`. Only 5 of the 133 Rd pages are internal (`@keywords internal`),
and one of them is the package page. What those blocks hold:

<!-- BLOCK codesize-roxtags -->

<!-- BLOCK codesize-comments -->

So the internal roxygen is prose ("why" text, derivations, cases), not
parameter lists: 31 `@param` lines in all. History and provenance
markers are rare (2.6% of comment lines).

## Q2. Dead code

### Static call graph

<!-- BLOCK codesize-deadcode -->

Reading the static result: with every symbol counted as a possible call
(the conservative choice, so that `do.call(f)`, `lapply(x, f)`,
`match.fun(f)` and registry entries are not missed), almost nothing is
unreachable. One function is dead (`notify_once()` in `R/zzz.R`, 9 code
lines). Eight (79 code lines) are reached only from tests or an
extension's tests (`ncp_block_chol()` from
`extensions/frmtmb.sample/tests/testthat/test-reparam.R`). No function is
reached only through a string. The 38 functions never called by name are
all reached as values; I checked 7 of them by grep (`ce_plan_part`,
`loglik_rows_same`, `layer_eta`, `gp_start`, `smooth_group_block_ids`,
`subset_rows`, `hyp_parse`): each is passed to `lapply()`, `vapply()`,
`do.call()` or a registry slot. Dynamic lookups by a pasted name
(`get(paste0(...))`, `match.fun(paste0(...))`) do not occur in `R/`.
Limits: a symbol edge also matches a local variable of the same name, so
the static graph can only under-count dead code. Dead code at the level
of branches is what the coverage run measures.

### Coverage

`dev/codesize-covr-env.R` ran all 349 core test files in one R process
under `covr::environment_coverage()` (5766 s wall, NOT_CRAN=true, gated
brms fit tests off). It fits the brief's "about 2 hours", so I did not
fall back to `file_coverage()` samples. Covered: 94.4% of 37,281
relevant lines (covr counts the first line of each expression, so a
multi-line call is one line); 2,100 lines no test runs.

Two test results differ from a one-file run. I take both to be
artifacts of one traced session, not package defects, but did not
verify it: `test-no-superassign.R` failed 1 (it scans the live function
bodies of the namespace, which covr has rewritten), and
`test-bracket-access.R` had 4 errors in the session.

A second artifact limits the "functions no test enters" list. A
function stored as a value in a registry list (`family_registry`, the
covstruct slots) is a separate copy in the lazy-load database, and covr
traces only the namespace binding. So `fam_inverse_gaussian()` reads as
never run although `test-brms-families.R` fits `inverse.gaussian`.
The block separates the 12 internal functions that are reached only as
values. Among "the rest", code that runs at load or build time
(`.onLoad()`, `by_dim()`, `spatial_entry()`) and helpers called only
from registry-held closures (`cs_chol_cor()`, `us_chol_scaled()`,
`homogeneous_sd()`) are the same artifact. What is left is about 25
internal functions and 300 code lines with no test: chart helpers for
priors and the GEV shape, the `loo_*` refusals, `quad_breakdown_message()`,
`sep_column_cert()`, `predictor_away()`, the truncated-mean helpers of
the extreme-value families, and some warnings. The 91 exported functions
that no test enters are mostly small S3 methods (print, chart, format)
and the `frmtmb_multiple` wrappers.

<!-- BLOCK codesize-covsum -->

## Q3. Duplication

<!-- BLOCK codesize-dup -->

Reading it: 1,063 of the 3,467 duplicated lines are rows of one data
table, `compat_hand_rules_tbl()` in `R/compat.R` (1,081 code lines of
`r("feature", "feature", "status", "note")` calls), which a token window
sees as repeats. Outside `R/compat.R`, 2,404 code lines (3.7%) lie in a
repeated 60-token window. Whole-function clones are small: 91 code lines
beyond one copy (strict), 327 (loose), and most of the loose ones are
exported family constructors and `frmtmb_multiple` or `loo()` method
wrappers, which mirror brms's API one function per name. The repeated
blocks worth a merge: the zero-inflated and hurdle constructors in
`R/families.R` (three copies each of 15 to 20 lines, at 2579/2637/2738,
2680/2805/3331 and 3121/3166/3216), `se_ridge_curvature()` and
`se_ridge_ritz()` (16 lines), two covstruct registry blocks
(`R/covstruct.R` 724/769 and 1122/1173), `build_objective()` and
`frm_linpred()` against `R/predict.R` (17 to 23 lines), and the shared
operator table of `R/links-brms.R` and `R/generic-owners.R` (42 to 51
lines).

Condition construction is already central: 1,378 code lines call
`frm_stop()`, 7 call base `stop()`. What costs lines is the message
text: 6,603 code lines (10.2%) are inside a condition call. Branches on
a family's name are rare (29 lines compare a family name, 4 switch on
it): families reach the core through slots.

## Q4. Single-use helpers and guard code

The single-use counts are in the deadcode block above ("Single-use
helpers"): 990 of the 2,161 live internal functions have one call site,
and 254 of those are under 10 code lines (1,528 code lines, 1,650 roxygen
lines). A small single-use helper carries about as many roxygen lines as
code lines (median 6 code, 5 roxygen).

<!-- BLOCK codesize-guard -->

Guard and check logic, by file and name, is 22.5% of function code
(14,202 lines). It overlaps the 6,603 lines inside condition calls (Q3).
The modeling core is the rest, and its largest functions are the
assemblers (`assemble_frame()`, 1,781 code lines in one function;
`resolve_priorlist()`, 921; `fit_assembled()`, 525).

## Q5. Growth

<!-- BLOCK codesize-growth -->

What drove it. Code, roxygen and comments grew in step (x2.67, x2.45,
x2.40), so the growth is not a documentation habit alone. The subsystem
that grew fastest is the fit-end judgement layer: the check files went
from 44 lines to 8,372, and guard-named functions in other files (Q4)
add about 8,000 more. The NEWS of the biggest releases says the same:
59% of the 1,020 bullets of 0.61.0 to 0.83.0 name a check, warning,
refusal, limit, boundary, standard error or start (a keyword proxy,
`dev/codesize-growth.R`). In the sample I read:

- 0.83.0 (+6,185 lines; se-check.R +742, predict-brms.R +643,
  gp-start.R +513): its first section, "Convergence and standard
  errors", is four BREAKING bullets, each a new rule of the SE check for
  a way a fit ends at a limit (coefficients that run off together, a
  zero-inflation probability that separates, NaN SEs for rows a limit
  direction moves, a gradient along a level direction). 50 of its 98
  bullets are BREAKING.
- 0.82.0 (+4,784; se-check.R +591, fit.R +516, predict-brms.R +477):
  the "limit" verdict and its four helpers, the quadrature settle rule,
  the prior-conflict escape, the draw law per side.
- 0.79.0 (+7,610; autocor-latent.R +1,440, predict.R +772,
  par-template.R +745): a feature release, brms's latent residuals
  (`ar()`, `ma()`, `cosy()`, `unstr()` on non-gaussian families), plus
  the defects the equivalence survey measured.
- 0.80.0 (+5,892; quadrature-nested.R +992, predict.R +744,
  par-template.R +653): a feature (nested quadrature) plus its guards.
- 0.71.0 (+5,174; families-extreme.R +970, families-simplex.R +711,
  multiple.R +677): new brms families.

So two things drove the growth at about the same rate: brms features,
which add a family, term or method and its tests, and edge-case
judgement at the end of a fit, which adds a rule, its proof conditions,
its message and its BREAKING note per case found. The second has no
counterpart in brms, which samples and does not have to judge a
maximum. Each NEWS bullet costs 60 to 70 lines of `R/` in the large
releases.

## Q6. Candidates to move out of core

Under `dev/core-boundary.md`: question 1 (extensions build on it) keeps
the seams, curve inference (`frmtmb.survival` uses it), the compat
table and the brms post-fit surface in core. Question 2 (most users meet
it) keeps the families, the downstream glue (emmeans, marginaleffects,
insight) and the post-fit methods. The candidates under question 3 (a
model class or workflow for one domain), measured under guardrail 3:

<!-- BLOCK codesize-candidates -->

(Test seconds are summed over the core test files that name a
candidate's functions, timed under covr in one process: an order of
magnitude only. An exported function that tests reach through another
name is not counted, so a 0 means "no test file names it", not "no test
runs it".)

Verdicts under the rule:

| Candidate | Lines in `R/` (all / code) | To | What in core calls it | Verdict |
|---|---|---|---|---|
| Domain material outside `R/`: `inst/bcm/` (4 files, 1,284 lines), `inst/rl/` (414), `tests/testthat/test-bcm-*.R` + `helper-bcm.R` (3,679), `helper-rl.R` (151) and `test-rl-example.R`, vignettes `bayesian-cognitive-modeling.Rmd` (522), `reinforcement-learning.Rmd` (560), `habit.Rmd` (387), data `habit_prep` | 0 in `R/` (`R/data.R` documents `habit_prep`) | `frmtmb.learn` (RL, habit), `frmtmb.eam` or `frmtmb.latent` (the cognitive-modeling ports) | nothing in `R/`; core tests and vignettes only | strongest under question 3: a domain's workflow with no core caller. It does not shrink `R/`. |
| `spectral.R`: `frm_periodogram()`, `whittle()`, `frm_series_draw()`, and `vignettes/spectral.Rmd` (476) | 654 / 254 | none fits; a new time-series extension | nothing (0 callers) | question 3 says extension; guardrail 3 says weak (small, no change since 0.78.0, no dependency of its own). Leave it. |
| `importance.R`: the importance-sampling correction of the Laplace approximation | 1,878 / 1,136 | `frmtmb.sample` (the second inference engine) | 9 functions: `fit_assembled()`, `quad_reanchored_fn()` (fit.R), the nested-quadrature path, `confint.frmtmb_fit()`, `print.summary.frmtmb_fit()` | a seam is needed first (a fit-end hook that replaces the objective and its log-likelihood report). Medium; the largest single `R/` saving among the moves. |
| `multiple.R`: `frm_multiple()` and 34 methods | 1,819 / 1,131 | a new extension | `posterior_summary.frmtmb_multiple()` only | a workflow, but not one domain's (missing data is in core's list). Keep. |
| `interop.R`, `insight.R` | 3,097 / 1,679 | none | 11 callers | post-fit glue most users meet. Keep. |
| curve inference | 2,098 / 945 | none | `frmtmb.survival` builds on it | question 1: core. Keep. |
| extreme-value families | 1,292 / 623 | none | the family registry, priors, CDFs | general families: core. Keep. |

Moving everything that the rule allows takes at most `importance.R` and
`spectral.R` out of `R/`: 2,532 lines, 2.1% of `R/`. Moves do not
answer the size question.

## Q7. Ranked reduction opportunities

Ranked by saving against risk. Line savings are of `R/` (118,169
lines). "Parse-identical" means: parse every file of `R/` before and
after with `keep.source = FALSE` and compare the deparsed top-level
expressions with `identical()`; then no behavior can change and no test
run is needed beyond one smoke file.

| # | Opportunity | Est. saving | Risk | How a lane proves behavior unchanged |
|---|---|---|---|---|
| 1 | Write the 2,435 `@noRd` blocks as plain `#` comments: drop the 2,438 `@noRd` tag lines and the 3,653 empty `#'` separator lines (the 16,336 prose lines stay). Needs a user decision: house style puts internal docs in roxygen. | 2,438 to 6,091 (2.1% to 5.2%) | none to behavior | parse-identical; `roxygenise()` leaves `man/` and `NAMESPACE` byte-identical |
| 2 | Remove the dead and test-only functions: `notify_once()` and the 8 reached only from tests (keep `frmtmb_compat_validate()` if its test is wanted; move `ncp_block_chol()`'s use into frmtmb.sample's test helper) | about 160 (88 code, 62 roxygen, comments) | none | identical results of the core suite and the frmtmb.sample suite, one file per process |
| 3 | Move the data table `compat_hand_rules_tbl()` (`R/compat.R` 1069 to 2360: 1,081 code, 173 comment, 30 roxygen lines; 528 `r()` rows) to a TSV in `inst/` read at load | about 1,250 out of `R/` (content moves, not removed) | low | `identical()` of `frm_compat()`'s full table with core alone and with all eight extensions loaded; `test-compat*.R`, `test-compat-extensions.R`, `test-guardrail1.R` |
| 4 | Merge the repeated blocks of Q3: the zero-inflated and hurdle constructors in `R/families.R` (three sets of three), `se_ridge_curvature()`/`se_ridge_ritz()`, the covstruct registry pairs, the operator table shared by `R/links-brms.R` and `R/generic-owners.R`, the `loo()`/`waic()`/`LOO()` and `frmtmb_multiple` wrappers, `prior_*()` constructors | 500 to 1,200 (2,404 code lines lie in repeated windows outside compat.R; about half can go) | low to medium | identical suite results; the equivalence survey pairs on the families touched (`dev/equiv-measure.R`); the ported brms ledger rows of those families unchanged |
| 5 | Edit internal commentary to "why, not how": 16,336 prose lines in `@noRd` blocks and 10,989 other comments; only 1,287 comment lines carry history markers, so most of it is explanation | about 5,000 at a 20% cut | none to behavior; a cost in review time and in knowledge if done carelessly | parse-identical |
| 6 | Shorten condition messages: 6,603 code lines (10.2%) are inside condition calls, about 4 lines per call; move the long explanations to a help topic per condition class and keep one sentence plus the class | 2,000 to 3,000 | medium: tests that match message text by regular expression break and must be read one by one | identical condition classes over the suite (record class per expectation); `test-conditions*.R`; equivalence survey verdicts unchanged |
| 7 | Inline the 254 single-use helpers under 10 code lines (1,528 code, 1,650 roxygen lines) where the call site reads as well | 1,000 to 1,500 | low to behavior; a readability judgement | identical suite results |
| 8 | Act on the coverage gaps: about 25 internal functions (about 300 code lines) that no test runs; test each or delete it | a few hundred | medium: a guard that no test reaches may still be reachable | identical suite results plus a new test for each one kept; equivalence survey; ported ledger |
| 9 | Moves under the core boundary (Q6): `importance.R` to `frmtmb.sample` after a seam; the domain material outside `R/` to `frmtmb.learn` and `frmtmb.eam` | 1,878 in `R/` (importance), about 7,000 package lines outside `R/` | medium: a new floor, cross-package defects, release cost (`dev/core-boundary.md` "Costs") | `test-guardrail1.R`; the ported ledger and the release sessions on the ecosystem (`dev/release/run-sessions.sh`); the importance tests moved and passing in the extension |
| 10 | Consolidate the fit-end judgement layer: 14,202 code lines classified as guard or check (22.5%), the subsystem that grew from 44 to 8,372 lines in its own files | unknown; the only lever that reaches tens of percent | high: these rules exist because of silent wrong answers; needs a user decision on which verdict classes and escapes to keep | equivalence survey verdicts unchanged (`dev/equiv-data/results.tsv`), identical warning and message classes over every suite, SE tests under OpenBLAS (`dev/ciharden-openblas.sh`), ported ledger |

Items 1 to 7 together are about 12,000 to 18,000 lines (10% to 15%)
with no change to what a fit computes. None of them brings `R/` near
brms's 44,729 lines: the code alone (64,884 lines) is 2.2 times brms's
29,772, and the part of it that brms has no counterpart for, the
judgement at the end of a maximum-likelihood fit, is item 10.

## What I did not do

- No package code changed, so no `R CMD check`, no NEWS bullet, no
  backlog item closed.
- Coverage ran core's suite in ONE R process with
  `covr::environment_coverage()` on a keep.source build of the worktree
  installed in `C:/Users/adf44/source/r/wt-codesize-lib`, not
  `covr::package_coverage()`, because package_coverage's parent plus
  child would have used both of the lane's R processes and left none for
  the static scripts. Running all files in one session is what CI's
  coverage job does, so a test that leaks state across files can differ
  from a one-file run. The extensions in rellib-r22 were on the library
  path, unlike CI, so tests that use an installed extension ran.
- The guard classification (Q4) and the NEWS keyword split (Q5) are
  name and keyword rules, stated in their blocks. They estimate; they do
  not prove.
- Duplication used one window size (60 tokens, about 8 to 12 lines).
  A smaller window finds more, mostly idiom.
