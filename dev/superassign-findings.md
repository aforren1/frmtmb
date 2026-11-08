# Lane superassign: no `<<-` in the package code

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1). Worktree
`frmtmb-wt-superassign`, library `wt-superassign-lib`. Status: DONE.

## The finding, and how pkgcheck reads it

pkgcheck 0.3.4.5 (the version in the CI log of run 42116483431, commit
03004e4d) reports "Package uses global assignment operator ('<<-')".
The check is `pkgchk_left_assign()` in pkgcheck's
`R/checks-left-assign.R` (read with `gh api` from
ropensci-review-tools/pkgcheck). It does not use pkgstats. It parses
every `*.R`, `*.r`, `*.q`, `*.s` file under `R/`, takes
`utils::getParseData()`, and tallies the tokens `LEFT_ASSIGN` and
`EQ_ASSIGN` by their text. The check fails when the `<<-` tally is
above 0 after two exemptions: a `<<-` inside a `setRefClass()`
definition, and `<<- memoise` inside `.onLoad()`. frmtmb uses neither.
The same function also scores "inconsistent assignment operators"
(both `<-` and `=` used as assignment).

`dev/superassign-leftassign.R` repeats that tally without pkgcheck.
It also counts `->>` (`RIGHT_ASSIGN`), which pkgcheck does not.

Before (log `dev/superassign-leftassign-before.log`, unedited tree):

    frmtmb           TOTAL <<- 76 in 21 files; ->> 0; '<-' 18811; '=' 0
                     pkgcheck global_assign: FAIL
    frmtmb.ode       TOTAL <<- 1 in 1 files      global_assign: FAIL
    frmtmb.sample    TOTAL <<- 6 in 4 files      global_assign: FAIL
    frmtmb.survival  TOTAL <<- 1 in 1 files      global_assign: FAIL
    coupling, eam, latent, learn, spatial: 0

The brief's "74 uses" counts lines: two lines held two `<<-` each
(`frame.R` `frame_raw_vars()`, `simulate-new.R` `nat_slots()`), so 74
lines carry 76 tokens.

After (log `dev/superassign-leftassign-after.log`):

    every package: TOTAL <<- 0; ->> 0; '=' 0
    pkgcheck global_assign: pass; left_assign consistency: pass

A second instrument agrees: goodpractice's lint for the operator,
`lintr::undesirable_operator_linter()` restricted to `<<-` and `->>`
(`dev/superassign-lintr.R`), counts 76 in core R/, 1 in ode, 6 in
sample, 1 in survival on the base sources (main checkout, log
`dev/superassign-lintr-base.log`) and 0 everywhere on the lane's
(`dev/superassign-lintr.log`).

## Other pkgcheck findings, read cheaply

From the summary of the last completed pkgcheck CI run (42116483431,
commit 03004e4d; log saved as `dev/superassign-pkgcheck-ci.log`):

- "The following function has no documented return value:
  [frmtmb-gp]" and "These functions do not have examples:
  [frmtmb-gp]". The base commit 3cd2a0c8 ("add return value & example
  to gp() doc") addresses both; its own pkgcheck run was in progress
  when this lane started.
- "Package fails continuous integration checks": pkgcheck scores the
  latest R-CMD-check and test-coverage runs, which were red.
- Info lines, not crosses: goodpractice's lintr flags `<<-` and `:::`
  together under one linter (1246 lines), and several style linters
  (`expect_true(a && b)`, `x %in% 'value'`, `c()` with one argument).
  This lane removes the `<<-` part only.
- Left-assign consistency passes: 0 uses of `=` as assignment.

## What replaced each `<<-`

The rule used, site by site:

- A closure that must keep state between calls (a condition handler, an
  optimizer callback, an active binding, a memo, a builder) writes to an
  environment created for it, `x <- new.env(parent = emptyenv())`, with
  `x$field <- value`. The environment is created in the function that
  owns the state, beside the closure, and named for what it holds. No
  site uses `assign(..., envir = parent.frame())` or
  `environment()`-of-the-caller.
- A handler whose only job was to report what happened returns it:
  `tryCatch(..., error = function(e) e)` and an `inherits()` test.
- A recursive builder whose callers are all in one place passes the
  accumulator through (`parts <- build(nm, parts)`, `out <- put(out,
  full, i)`).

The functions changed, by file. Every one is listed so the consolidator
can resolve hunks.

| File | Function | Sites | Replacement |
|---|---|---|---|
| R/autoscale.R | autoscale_choose (inner `capture`) | 2 | env `got` |
| R/brms-shapes.R | fit_fd_se (memo `geta`) | 1 | `geta` is an env |
| R/ce-levels.R | ce_level_plan | 1 | env `tally` |
| R/compat.R | compat_rule_builder | 1 | env `acc` |
| R/compat.R | compat_hand_rules_tbl | 1 | env `acc` |
| R/confint.R | par_alias_index | 2 | env `got` |
| R/confint.R | varcorr_trans_rows | 1 | env `acc` |
| R/confint.R | hyp_env_vals | 1 | env `acc` |
| R/curve-cov.R | curve_catch_span | 1 | env `got` |
| R/fit.R | nlminb_trial_fn | 6 | env `st` |
| R/fit.R | newton_polish | 3 | env `rp` (M, Mi, decided) |
| R/fit.R | mo_search | 1 | env `tried` |
| R/frame.R | resolve_group_calls | 1 | env `acc` |
| R/frame.R | lookup_structural | 1 | handler returns `list(err = e)` |
| R/frame.R | frame_raw_vars | 2 | env `seen` |
| R/frame.R | assemble_frame (`add_part`, aterm codes) | 3 | env `rhs_acc`, env `codes_acc` |
| R/generic-owners.R | frm_bind_generic | 2 | env `memo`; doc sentence updated |
| R/gp-by.R | gp_brms_natural | 4 | env `acc` |
| R/hazard-containers.R | frm_hazard_reads | 1 | env `acc` |
| R/importance.R | build_importance_objective | 1 | env `shift` |
| R/me.R | me_newdata_extra | 1 | `obs` is an env |
| R/me.R | me_node_value | 1 | env `memo` |
| R/methods-fit.R | VarCorr.frmtmb_fit (`take`) | 1 | env `cur` |
| R/par-template.R | nl_restart_round | 1 | env `h1` |
| R/par-template.R | nl_fit_off_saddle | 1 | env `hq` |
| R/parse.R | mo_bar_rewrite | 1 | env `acc` |
| R/parse.R | expand_double_verts | 1 | env `queue` |
| R/parse.R | parse_linpred | 2 | env `track` |
| R/parse.R | rhs_split_ac | 1 | env `found` |
| R/predict.R | lp_basis_nl | 2 | `parts` threaded through `build()` |
| R/priors.R | as_priorlist | 1 | env `acc` |
| R/priors.R | priorlist_from_rows | 1 | handler returns `e` |
| R/priors.R | prior_table | 1 | env `acc` |
| R/priors.R | cov_class_coefs | 1 | env `acc` |
| R/priors.R | resolve_priorlist | 19 | env `st` (assigned, lower, upper) |
| R/ps.R | ps_extract | 1 | env `found` |
| R/re-formula.R | re_formula_bars | 1 | env `found` |
| R/simulate-new.R | nat_slots | 3 | env `acc` |
| ext ode R/ode.R | frm_ode (inner `run_solve`) | 1 | env `got` |
| ext sample R/loo.R | log_lik.frmtmb_draws | 1 | env `seen` |
| ext sample R/methods-draws.R | draws_laplace_watch | 1 | env `seen` |
| ext sample R/methods-draws.R | draws_unseen_level | 1 | handler returns `list(what =)` |
| ext sample R/prior-draws.R | prior_plan | 1 | env `acc` |
| ext sample R/sample.R | default_priors_for | 1 | env `acc` |
| ext sample R/sample.R | tmbstan_build_broken | 1 | env `memo` in its local() |
| ext survival R/cox-family.R | cox_bhaz_args | 1 | `out` threaded through `put()` |

Also changed, because a test reads the state the last row moved:
`extensions/frmtmb.sample/tests/testthat/helper-sampling.R`,
`with_tmbstan_hpp()`, which resets the memo. It wrote `cached` in the
closure's enclosure; it now writes `cached` in `env[["memo"]]`.

Mechanical edits used `dev/superassign-rename.R`, which rewrites SYMBOL
tokens by parse data in a line range (never strings, comments, argument
names or `x$name`), and every hunk was read after it. No added line is
over 80 columns. Line endings stay CRLF, as checked out
(`dev/superassign-eol.R`).

`resolve_priorlist()` is the large hunk: 61 references to `assigned`,
`lower` and `upper` became `st$...`, because the class resolvers and the
main loop both write them and read each other's writes (`claim()`).

## The guard

`tests/testthat/test-no-superassign.R`, a test and not a CI grep. Why:
a grep in a workflow would also hit `<<-` in a comment or a string,
would run only where that workflow runs, and would touch the workflow
files the ci lane owns this round. The test runs in every suite run,
under `R CMD check` and under covr.

Two scanners, as test-bracket-access.R has two:

- `superassign_ns_scan()` walks every function of the installed
  namespace, bodies and default arguments, inner functions included,
  for a `<<-` call. It runs everywhere, including `R CMD check` on a
  tarball. A `->>` parses to a `<<-` call, so it sees both spellings.
  A function whose top-level environment is another package's (an
  active binding resolved to its owner's generic) is not judged.
- `superassign_src_scan()` applies pkgcheck's own token rule to the
  sources of core and of every extension. It skips where there is no
  source tree, as the bracket backstop does.

A control test shows each scanner a hit and the non-hits beside it (a
`<<-` in a string, in a comment, `acc$n <- 1`), including a hit that
sits only in an inner function's default argument.

Seen to fail on the base build. `dev/superassign-guard-base.R`
evaluates the test file's scanner definitions and runs them on
rellib-r18 and on the base sources (main checkout, 3cd2a0c8, R/ and
extensions/ unmodified there). Log `dev/superassign-guard-base.log`:

    lib: C:/Users/adf44/source/r/rellib-r18 version: 0.79.1
    NS functions with <<-: 38      (the 38 core functions in the table)
    SRC <<- tokens in .../frmtmb/R: 76 in 21 files
    SRC frmtmb.ode: 1, frmtmb.sample: 6, frmtmb.survival: 1, others 0

The test file itself, run against the base build with this worktree's
sources (so only the namespace scan can see the defect), logs
`dev/superassign-guard-before.log`:

    lib: frmtmb 0.79.1 C:/Users/adf44/source/r/rellib-r18
    DETAIL FAIL [no frmtmb function uses `<<-`] ... [1]
      "VarCorr.frmtmb_fit" [2] "as_priorlist" ...
    RESULT test-no-superassign.R pass=13 fail=1 err=0 skip=0 warn=0

and against the lane build, `dev/superassign-guard-after.log`:

    RESULT test-no-superassign.R pass=14 fail=0 err=0 skip=0 warn=0

The guard's positive condition was also seen false before it was made
right: its first spelling returned NULL where `character(0)` was
expected, so the namespace assertion and the control failed on the
lane build too (2 failures); `as.character()` on each scanner's result
fixed both.

## Behavior

Instrument: `dev/superassign-trace.R` traces, on exit, every changed
function that returns data (and `run_optimizer()`, which reads the state
`nlminb_trial_fn()` keeps) while one test file runs, hashes a canonical
form of each return value (functions and environments as placeholders,
formulas and calls deparsed, doubles to 9 significant digits for the
cross-process ulp noise of rtmb-pitfalls item 21), and writes the
(function, hash) sequence. `dev/superassign-trace-run.sh` runs each
file under the base build (rellib-r18) and the lane build, one process
at a time, `NOT_CRAN=true`, gates off. `dev/superassign-trace-compare.R`
compares the two sequences call by call. Job list
`dev/superassign-trace.jobs` (31 files: 26 core, 4 frmtmb.sample, 1
frmtmb.survival). Log `dev/superassign-trace-compare.log`:

    FILES 31, with any difference 0
    every file: same call count in both arms, 0 differing calls, and the
    same pass/fail/err/skip counts in both arms

Calls per traced function, lane arm, summed over the 31 files:

    as_priorlist 1204, assemble_frame 827, ce_level_plan 34,
    compat_hand_rules_tbl 307, cov_class_coefs 109, cox_bhaz_args 66,
    curve_catch_span 94, default_priors_for 44, draws_unseen_level 43,
    expand_double_verts 944, fit_fd_se 9, frame_raw_vars 791,
    frm_compat 146, frm_hazard_reads 5, gp_brms_natural 56,
    hyp_env_vals 4443, lookup_structural 49, lp_basis_nl 32,
    me_newdata_extra 9180, mo_search 1496, nat_slots 25,
    newton_polish 1494, nl_fit_off_saddle 30, par_alias_index 32,
    prior_plan 29, prior_table 83, priorlist_from_rows 2, ps_extract 42,
    re_formula_bars 45, resolve_group_calls 464, resolve_priorlist 291,
    run_optimizer 1883, VarCorr.frmtmb_fit 9, varcorr_trans_rows 18

Traced but not reached by these 31 files: autoscale_choose,
mo_bar_rewrite, nl_restart_round, rhs_split_ac. A second batch
(`dev/superassign-trace2.jobs`: test-mo-factor-re, test-nl-family,
test-predfix, test-grad-verdict, test-latres) reaches all four, with
the compare rerun over all 36 files:

    FILES 36, with any difference 0
    autoscale_choose 8, mo_bar_rewrite 9, nl_restart_round 38,
    rhs_split_ac 7
    (test-latres 4796 calls, test-nl-family 4008, test-predfix 616,
    test-mo-factor-re 259, test-grad-verdict 199; 0 differing in each)

So all 33 changed core functions that return data (the other 5 return
closures, below), and 3 of the 9 changed extension functions, returned the same value on
the base and lane builds at every call these 36 files make.

Not traced, because what they return is a closure or a flag: the
closures of nlminb_trial_fn (judged through run_optimizer, whose
result carries `nonfinite_trials`), compat_rule_builder (judged
through frm_compat), build_importance_objective (judged through the
fits of test-importance.R), me_node_value, frm_bind_generic,
draws_laplace_watch, tmbstan_build_broken, log_lik.frmtmb_draws and
frm_ode. These are covered by the full suites below.

## Timings

Script `dev/superassign-bench.R`, seed 20261106, one process, arms
interleaved in random order, each block grown past 1.2 s of
`system.time()`, minimum per-call time over 7 blocks. Arms: `base`
(the 3cd2a0c8 source of the function, from the main checkout, evaluated
in the lane namespace and byte-compiled), `control` (the lane's source
built the same way), `lane` (the installed function). Base against
control is the before/after comparison; lane against control is the
same code twice and measures the noise.

Sites timed: the ones that run more than once per fit.
`nlminb_trial_fn()$fn` runs once per objective evaluation; the
active binding of `frm_bind_generic()` once per call of a borrowed
generic; `compat_hand_rules_tbl()` builds about 650 rows one at a time
(cached once per session; timed uncached); `prior_table()` builds the
default prior table row by row.

First run (`dev/superassign-bench.log`) had the base and control arms
interpreted and the lane arm byte-compiled; it is not a comparison and
is kept only because it is where the construction problem showed.
Runs 2 and 3 (`dev/superassign-bench-1.log`, `-2.log`), s per call:

    site                     run  base      control   lane      lane/ctl
    nlminb_trial_fn()$fn      2   1.612e-06 1.869e-06 1.860e-06 1.00
                              3   4.616e-06 5.379e-06 2.861e-06 0.53
    active binding access     2   2.708e-06 2.766e-06 2.804e-06 1.01
                              3   2.975e-06 4.578e-06 4.005e-06 0.87
    compat_hand_rules_tbl()   2   1.812e-01 1.975e-01 1.800e-01 0.91
                              3   4.225e-01 5.225e-01 5.275e-01 1.01
    prior_table()             2   3.828e-02 3.141e-02 4.141e-02 1.32
                              3   5.312e-02 7.500e-02 5.188e-02 0.69

What this shows: no difference resolvable above the noise. The
lane/control column is the same code timed twice, and it spans 0.53 to
1.32 on this machine, which twelve lanes are loading; the
base/control ratios (0.86 to 1.38) sit inside that span. Run 3 is about
2.5 times slower than run 2 on the same work. The magnitudes bound the
cost: the optimizer wrapper and the binding cost 2 to 5 us per call
either way; the compat table is built once per session; the default
prior table of a model with two blocks and a sigma formula is tens of ms
with or without the change. Identity checks in the same script: compat
table 647 rows, `identical()` base and lane TRUE; prior table 15 rows,
TRUE.

Not done: a load-independent count (allocations or AD nodes) for these
sites. None of the rewritten sites is on the tape; the importance
objective reads `shift$offs` once per taping, not per evaluation. The
one-row-at-a-time appends (compat, prior table, priors) were kept as
appends into an environment rather than preallocated: their row count
is not known in advance, each runs once per session or per prior
table, and the run above shows no cost to remove.

## Tests

### One file per process, every suite

`dev/superassign-run-par.sh` (a copy of dev/release/run-par.sh pointed
at `dev/superassign-run-tests.R`, which is dev/release/run-tests.R with
rellib-r18 in place of rellib-r17), job list `dev/superassign-full.jobs`
(500 files), two processes at a time, `NOT_CRAN=true`,
`FRMTMB_BRMS_FIT_TESTS=true`, `FRMTMB_RELLIB=wt-superassign-lib`. Log
`dev/superassign-full.log`, per-file logs `dev/superassign-suite-full/`.
Every one of the 500 `lib:` lines names wt-superassign-lib for the
package and for frmtmb.

    FULL ran 500 of 500; NO RESULT lines: 0
    files=500 pass=35196 fail=0 err=0 skip=29 warn=0
    frmtmb           files=309 pass=25575 fail=0 err=0 skip=15 warn=0
    frmtmb.coupling  files=11  pass=551   fail=0 err=0 skip=5  warn=0
    frmtmb.eam       files=33  pass=1854  fail=0 err=0 skip=3  warn=0
    frmtmb.latent    files=12  pass=683   fail=0 err=0 skip=2  warn=0
    frmtmb.learn     files=16  pass=545   fail=0 err=0 skip=2  warn=0
    frmtmb.ode       files=11  pass=549   fail=0 err=0 skip=1  warn=0
    frmtmb.sample    files=85  pass=4326  fail=0 err=0 skip=1  warn=0
    frmtmb.spatial   files=8   pass=561   fail=0 err=0 skip=0  warn=0
    frmtmb.survival  files=15  pass=552   fail=0 err=0 skip=0  warn=0

Every skip is a tier gate outside this run: the scale tier
(`FRMTMB_SCALE_TESTS`: test-curve-scale 1, and test-scale in coupling
5, eam 3, latent 2, learn 2, ode 1, sample 1), the drmTMB agreement
tier (`FRMTMB_DRMTMB_FIT_TESTS`: test-drmtmb-agreement 13) and the
grammar fuzz tier (`FRMTMB_FUZZ`: test-fuzz 1). The last two reach
the parser and the fit path, so they were run separately with their
gates on (`dev/superassign-extra.jobs`, log `dev/superassign-extra.log`,
lane build on both `lib:` lines):

    RESULT test-drmtmb-agreement.R pass=131 fail=0 err=0 skip=0 warn=0
    RESULT test-fuzz.R pass=2 fail=0 err=0 skip=0 warn=0

The scale tier was not run: it measures cost, and the timings above are
the cost evidence for this change.

### One session per suite, three orders

`dev/superassign-run-sessions.sh` (copy of dev/release/run-sessions.sh
pointed at `dev/superassign-run-session.R`, rellib-r18 behind the lane
library), one session at a time, same environment as the full run.
Per-session logs `dev/superassign-suite-sess/`. Every `lib:` line names
wt-superassign-lib for all nine packages. Every LEAK line carries
`[with load of ...]`; untagged LEAK lines: 0 in every session.

    package          each of forward, preload, reverse
    frmtmb.coupling  11/11 pass=551  fail=0 err=0 skip=5 warn=0
    frmtmb.eam       33/33 pass=1854 fail=0 err=0 skip=3 warn=0
    frmtmb.latent    12/12 pass=683  fail=0 err=0 skip=2 warn=0
    frmtmb.learn     16/16 pass=545  fail=0 err=0 skip=2 warn=0
    frmtmb.ode       11/11 pass=549  fail=0 err=0 skip=1 warn=0
    frmtmb.sample    85/85 pass=4326 fail=0 err=0 skip=1 warn=0
    frmtmb.spatial   8/8   pass=561  fail=0 err=0 skip=0 warn=0
    frmtmb.survival  15/15 pass=552  fail=0 err=0 skip=0 warn=0
    frmtmb  forward  309/309 pass=25575 fail=0 err=0 skip=15 warn=0
    frmtmb  reverse  309/309 pass=25725 fail=0 err=0 skip=15 warn=0
    frmtmb  preload  309/309 pass=25818 fail=0 err=0 skip=17 warn=0
    (core pass counts differ by order because some tests assert more
    when an extension is loaded; preload adds 2 skips of "when it is
    not loaded" tests; LEAK lines 48, 49 and 38, all tagged)

### OpenBLAS

No arithmetic changed, but `nlminb_trial_fn()` and `newton_polish()`
sit in the optimizer, so the identity check was repeated under R with
OpenBLAS 0.3.26 for BLAS and LAPACK (`dev/ciharden-openblas.sh 0.3.26
lapack`, the ubuntu-24.04 runner's configuration): the trace of
`dev/superassign-trace-run-ob.sh` over 7 optimizer-heavy files
(`dev/superassign-traceob.jobs`), base and lane both under OpenBLAS.
Log `dev/superassign-traceob-compare.log`:

    FILES 7, with any difference 0
    test-autoscale 214 calls, test-importance 1220, test-mo-terms 109,
    test-nl-family 3985, test-optrap 476, test-predfix 617,
    test-se-reparam 508; same pass counts in both arms, 0 failures

The call counts differ from the reference-BLAS run on two files
(importance 1220 against 1219, nl-family 3985 against 4008): the BLAS
takes the optimizer down a different path, the same path in both arms.
0.3.32 was not run: no rounding-decided branch is near this change.
`dev/ciharden-out/` was deleted afterwards.

The gated brms tests compiled the Stan programs missing from the
seeded cache (240 programs at the start, it grows as they compile);
that, not the change, is why this run took about 5 hours.

### R CMD check --as-cran

`dev/superassign-check.sh`, once per changed package, built with
vignettes and checked in `dev/superassign-check/`,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX on PATH,
`R_LIBS` = lane library, rellib-r18, user library. Logs
`dev/superassign-check/check-<pkg>.log` (the .Rcheck trees and
tarballs were deleted afterwards).

    frmtmb           Status: 1 NOTE   (tests 38 min, OK)
    frmtmb.ode       Status: 1 NOTE   (tests 92 s, OK)
    frmtmb.sample    Status: 1 NOTE   (tests 11 min, OK)
    frmtmb.survival  Status: 1 NOTE   (tests 25 s, OK)

Each NOTE is "checking HTML version of manual ... Skipping checking
math rendering: package 'V8' unavailable", the machine's known one.

## Bump and floor

Patch for all four packages (frmtmb, frmtmb.ode, frmtmb.sample,
frmtmb.survival): no export, argument, default or result changes. No
new floor: no extension calls anything this lane added, and each
extension's change stands on its own code. The other five extensions
are unchanged.

NEWS: a bullet under a new `# <pkg> (development version)` heading in
NEWS.md of frmtmb, frmtmb.ode, frmtmb.sample and frmtmb.survival.

## Files touched

Code: the 21 core files and 6 extension files in the table above.
Tests: `tests/testthat/test-no-superassign.R` (new),
`extensions/frmtmb.sample/tests/testthat/helper-sampling.R`.
NEWS.md of the four packages. Lane tooling (new, under dev/):
superassign-leftassign.R, -lintr.R, -rename.R, -eol.R, -guard-base.R,
-install.R, -trace.R, -trace-run.sh, -trace-run-ob.sh,
-trace-compare.R, -bench.R, -check.sh, -run-tests.R, -run-par.sh,
-run-session.R, -run-sessions.sh, job lists -full.jobs, -trace.jobs,
-trace2.jobs, -traceob.jobs, -extra.jobs, and this file.

roxygenise changed only `Config/roxygen2/version` (8.1.0 to 8.1.1) in
all nine DESCRIPTION files; those were restored to the base content,
because the bump is not this lane's and every lane would make it.
No Rd or NAMESPACE file changed.

## What I did not do, and why

- No CI grep. The guard is a test (see "The guard"); the workflow
  files belong to the ci lane this round.
- Extensions do not get a namespace-scan test of their own. pkgcheck
  runs on core alone; core's source scan covers every extension's R/
  wherever the monorepo is checked out (local runs). Under an
  extension's own `R CMD check` nothing scans it. Adding the same
  test to eight suites, or exporting the scanner as was done for
  `frm_hazard_reads()`, is a choice for the user.
- The one-row-at-a-time appends were not preallocated (see Timings).
- Other goodpractice lints in the pkgcheck log (`:::`, expectation
  style) were not touched: not this lane's.
- The scale tier was not run, and OpenBLAS 0.3.32 was not run.
- A load-independent cost count was not made; the clock could not
  resolve a difference on a machine twelve lanes were loading.

## Defects found, not fixed

- `frame_raw_vars()` collects `off_only` (names read only inside
  `offset()`) and never reads it. Kept, because removing dead code is
  not this lane's change; the consolidator may drop it.
- dev/release/run-tests.R and run-session.R still name rellib-r17;
  this lane used copies naming rellib-r18 rather than edit shared
  harness files.

## Review follow-ups (dev/reviews/2026-11-06-superassign.md)

1. `superassign_ns_scan()` now also walks every list and environment
   held in the namespace, recursively, and each function's enclosure;
   each environment is walked once (keyed by `format(env)`), so a cycle
   ends. The control test adds a `<<-` in a nested registry list
   (`reg$fcor$parse`) and in a self-referring environment (`box$fn`).
   Against the reviewer's mutant build (rev-superassign-lib), log
   `dev/superassign-guard-mutant.log`: `pass=15 fail=1`; the namespace
   scan now reports `autocor_core$fcor$parse`,
   `frmtmb_autocor_registry$reg$fcor$parse`, `rev_mut_list$f`,
   `rev_mut_env$f`, `.__T__show:methods$revMutS4`, `fam_gaussian`,
   `family_registry$gaussian`, `rev_mut_backtick`, `rev_mut_rassign`.
   The dynamic spellings (bquote, eval(parse()), assign()) stay unseen,
   as the review accepted.
2. The core NEWS bullet now says the extension check holds only when
   the tests run from the source repository.
3. helper-sampling.R: the orphan "# Prepending a" line is reflowed.
4. Test header: "in a comment-free line of code or not" is replaced by
   "any `<<-` written as code in R/ fails it, inside a function or at
   top level; one in a comment or a string is not a token".

Final run on the lane library (`dev/superassign-guard-final.log`):
`RESULT test-no-superassign.R pass=16 fail=0 err=0 skip=0 warn=0`.
