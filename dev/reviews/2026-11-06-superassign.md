# Review of lane superassign, 2026-11-06

Adversarial review. Worktree `C:\Users\adf44\source\r\frmtmb-wt-superassign`, base
3cd2a0c8 (frmtmb 0.79.1). Lane build `wt-superassign-lib`, before-build `rellib-r18`.
Base sources are the main checkout (`C:\Users\adf44\source\r\frmtmb`, 3cd2a0c8; only
`.github/workflows/test-coverage.yaml` is modified there). Scripts and logs are
`dev/review-superassign-*`. The mutant build is in `C:/Users/adf44/source/r/rev-superassign-lib`,
the mutant tree and the trial merges in `C:/Users/adf44/source/r/rev-superassign-src`
(`m1/`, `nosrc/`, `merge/`, `merge2/`). No package, test or lane file was edited. No git
operation beyond `diff`, `status`, `log` and `show`, with `GIT_OPTIONAL_LOCKS=0`. At most
three R processes at once. `dev/ciharden-openblas.sh` was not run by this review (see the
note at the end).

## Verdict: SHIP

No blocker and no major. Every rewritten site is equivalent to the code it replaces, and
the evidence below finds no difference at full precision. All the remaining items are minor.
They are a blind spot in the guard, two documentation details, and one merge conflict
with the readers lane that has an obvious resolution.

## What was checked, with the evidence

### 1. Semantic equivalence, site by site

I read every hunk of the 21 core files and the 6 extension files.

- **Where each environment is created.** Each new environment is created in the same
  frame as the variable it replaces, so its lifetime and sharing do not change.
  `tmbstan_build_broken` memo: once per package load, inside the same `local()`, as
  `cached` was. `draws_laplace_watch`: once per call, as `first` was. The active
  binding: once per `frm_bind_generic()` call, so there is one memo per binding.
  `compat_rule_builder`: one per builder. `build_importance_objective`: one per call.
  `nlminb_trial_fn`: one per call. `newton_polish`: one `rp` per call, with M, Mi and
  decided in it.
- **Dangling references.** `dev/review-superassign-globals.R` takes each changed
  function from the base and lane sources and compares its free variables
  (`codetools::findGlobals`). A local that is read after its definition moved into an
  environment shows up as a new free variable. Log `dev/review-superassign-globals.log`:
  `RESULT changed functions 46, with new free variables 0`.
- **Shadowing.** `dev/review-superassign-shadow.R` looks in each changed function for
  any other binding of the new environment's name, such as a nested formal, a `for`
  variable or a second assignment. A nested `st <- ...` would make `st$assigned` read
  the wrong object. Result: 0. The only other definitions are the intended per-call
  ones inside `capture()` (autoscale) and `run_solve()` (ode).
- **Snapshots.** Some sites copy the environment's field back into the old local
  (`rows <- acc$rows`, `held <- hq$w`, `bad <- acc$bad`, `out <- acc$out`,
  `rhs_comb <- rhs_acc$comb`). I read each one to confirm that no writer runs after the
  copy. In `nl_fit_off_saddle`, `hold` is installed once (line 806) and `release()`
  reads `hq$w` itself, not the copy. In `assemble_frame`, the last `add_part()` is at
  line 2989 and the copy is at 2993. `varcorr_trans_rows`, `nat_slots`, `as_priorlist`
  and `prior_plan` also copy after their last writer.
- **Handlers that now return the condition.** All three were already `tryCatch`
  handlers that unwind (`lookup_structural`, `priorlist_from_rows`,
  `draws_unseen_level`), so control flow does not change. `lookup_structural`
  initializes `e2 <- NULL` before the `tryCatch`, as before. `draws_unseen_level` uses
  `isTRUE(known)` / `is.list(known)`, and nothing below them reads `what`. No
  `withCallingHandlers` site was converted to `tryCatch`. Those sites keep
  `invokeRestart("muffle...")` and write to an environment.
- **List to environment.** `geta` in `fit_fd_se` and `obs` in `me_newdata_extra` are now
  environments. One edge differs: `env[[""]]` raises an error, but `list()[[""]]`
  returns NULL. The keys are `coef_block_key()` (a dpar name or `resp_dpar`) and
  `as.character(k)`, and neither can be empty, so this is noted and not filed.
- **Recursive accumulators.** In `lp_basis_nl`, the order of `parts` is post-order, as
  before, and an empty `nl_pars` leaves `list()`. In `cox_bhaz_args`, `put()` returns
  `out` and all three call sites assign it back. `names(out)` in `step1` reads the
  updated `out`, as it read the `<<-`-updated one before.
- **`resolve_priorlist`.** The free-variable check confirms that no reference to a
  stale `assigned`, `lower` or `upper` is left (all 61 references now use `st$`). The
  P8 probe in section 2 exercises the lb, ub and assigned interplay.

### 2. Behavioral probes of the closures the trace did not cover

`dev/review-superassign-probe.R` and `-probe2.R` run once against the base build
(rellib-r18) and once against the lane build. They save every result with doubles as
`%.17g` strings, and `-probe-cmp.R` compares the two. The `lib:` lines of all four
logs name the intended library. Log `dev/review-superassign-probe-cmp.log`:

| probe | what it exercises | base vs lane |
|---|---|---|
| P1 | `nlminb_trial_fn`: trials 5, 3, NaN, Inf, NaN, 2, NA, 1, NaN; a second tracker alive and interleaved; tracker A2 whose fn calls tracker B2 (nested) | identical (count 3, tally 3, best 1; B count 0; B2 count 3) |
| P2 | Gamma identity-link fits (seed 31) that meet nonfinite trials: 3, 1 and 5 mapped trials from three starts | identical at 17 digits (objective 91.041090378540929, ...) |
| P3 | `importance = 400L` poisson mixed fit (seed 7), twice in one process | identical; rep 1 = rep 2 = 164.99688439870019 |
| P4 | `frm_ode()`: failing solve, clean rerun, failing solve in the other group | identical values, warnings and failing groups; only `frm_ode_failures()$when` (a timestamp) differs |
| P5 | `frm_bind_generic()`: two bindings of `ndraws` to posterior; owner loaded, unloaded, reloaded; bound twice | identical (fallback, owner, fallback, fresh owner; a stale instance is never served) |
| P6 | two `compat_rule_builder()`s interleaved | identical (2 rows and 1 row) |
| P7 | two `me_node_value()` closures alive; calls to `what()` counted | identical (4 calls for 7 reads; the memos are separate) |
| P8 | `resolve_prior_input()` on 6 prior lists with 13, 3, 4, 3, 1 and 2 rows: class, group and coef bounds, bounds-only rows, order overrides, an lb = 0 sd, a refusal | identical entries, lower, upper and messages |
| P9 | mixed `(1 + x \| g)`, `b1 * exp(b2 * x)`, `mo()`, `sigma ~ x` fits (seed 5) | identical at 17 digits, including `mo_search` |
| nested | a custom optimizer that fits a second model with `frm()` inside the outer fit's optimizer | identical at 17 digits |
| allfit | `frm_allfit()` | 951 of 954 numbers identical; the 3 that differ are the `seconds` column |
| nlsaddle | `b1 * exp(b2 * x)` from its default start (`nl_fit_off_saddle`, restart rounds) | identical |

On the question of whether the trace's 9-significant-digit hash was too coarse to see a
change in the optimizer path: every fit above was compared at 17 digits, including the
three fits with nonfinite trials, and none differs. The trace's resolution was not what
let any change through.

`log_lik.frmtmb_draws` (`once`) and `draws_laplace_watch` were not probed directly,
because each needs draws, and so a Stan build. Each creates its state once per call, as
before. `test-loo.R` and `test-loo-predict.R` pass (section 6).

### 3. The guard

`dev/review-superassign-guard-mutant.log` runs the lane's test file from a mutant tree
(`rev-superassign-src/m1`) against a mutant build of that tree (`rev-superassign-lib`,
`lib:` line confirmed). The mutants are: a `<<-` in `fam_gaussian()`'s lpdf closure; a
`<<-` in `autocor_core$fcor$parse` (a top-level registry list); a `<<-` in
frmtmb.survival's `cox_bhaz_args()`; and `R/zz-review-mutants.R`, which holds the
`` `<<-`(x, 1) `` call form, `->>`, `bquote`, `eval(parse(text =))`,
`assign(envir = parent.env(environment()))`, a function in a top-level list, a function
in an environment, and an S4 `show` method. Result: `pass=11 fail=3`.

| mutant | namespace scan | source scan | pkgcheck's rule |
|---|---|---|---|
| lpdf closure in `fam_gaussian` | found | found (families.R:1658) | counted |
| `autocor_core` registry list | **missed** | found (autocor-seam.R:1072) | counted |
| extension `cox_bhaz_args` | not scanned | found (cox-family.R:335) | counted (survival) |
| `` `<<-`(x, 1) `` | found | **missed** | **not counted** |
| `1 ->> x` | found | found | not counted (RIGHT_ASSIGN) |
| function in a top-level list | **missed** | found | counted |
| function in an environment | **missed** | found | counted |
| S4 method | **missed** | found | counted |
| `bquote`, `eval(parse())`, `assign(envir =)` | missed | missed | not counted |

The two scanners together catch every static spelling. Only the dynamic spellings get
through, and pkgcheck cannot see those either. That is acceptable, and I report it only.
`assign(..., envir = parent.frame())` is the same hidden coupling. The guard's advice
text forbids it, but no rule checks for it. Whether the round should also forbid it is
for the user to decide.

`dev/review-superassign-guard-nosrc.log`: when the test file runs from a tree with no
`../../R`, the namespace test runs and both source tests skip with "package sources are
not available (installed-package run)". This is the same reason
`test-bracket-access.R` gives. `.github/scripts/skip-gate.sh` only fails on
`{pkg} cannot be loaded` and `{pkg} is not installed` skips, so these skips pass it.

### 4. pkgcheck itself

pkgcheck is not installed. `dev/review-superassign-pkgcheck-upstream.R` is upstream's
`R/checks-left-assign.R` (ropensci-review-tools/pkgcheck main, DESCRIPTION 0.3.4.005),
fetched with `gh api` and sourced verbatim. `dev/review-superassign-pkgcheck.R` calls
its own `pkgchk_left_assign()`, `output_pkgchk_global_assign()` and
`output_pkgchk_left_assign()`. Log `dev/review-superassign-pkgcheck.log`:

    base    frmtmb           <<- tokens  76 in 21 files; global_assign FAIL; left_assign pass (<- 18811, = 0)
    base    frmtmb.ode       <<- tokens   1 in  1 files; global_assign FAIL
    base    frmtmb.sample    <<- tokens   6 in  4 files; global_assign FAIL
    base    frmtmb.survival  <<- tokens   1 in  1 files; global_assign FAIL
    lane    every package    <<- tokens   0; global_assign pass; left_assign pass
    mutant  frmtmb           <<- tokens   5 in 3 files; global_assign FAIL

The lane's counts (76 in 21 files; 1, 6, 1) are confirmed by pkgcheck's own function.
On the mutant, pkgcheck's rule does not count the `` `<<-`() `` call form or `->>`.
The lane's source scan counts `->>` and the namespace scan counts both, so the guard is
stricter than pkgcheck.

### 5. Merge hazards

Every other worktree's diff was read on 2026-11-06. The lanes are still active, so this
is a snapshot. `dev/review-superassign-merge.R` maps each hunk to its top-level
definition in the base file. Log `dev/review-superassign-merge.log`:

- `<<-` or `->>` added by any other lane in `R/` or `extensions/`: **0**.
- Functions that superassign changes and another lane also changes: 7.

| function | other lanes | textual |
|---|---|---|
| `compat_hand_rules_tbl` | fitedge, readers, loglik, nestquad, osa, mfxse | apart (108 to 404 lines); they add `r(...)` rows through the unchanged `r()` closure |
| `assemble_frame` | bigx | apart (651 lines) |
| `me_node_value` | readers | **conflict** |
| `nl_restart_round` | loglik | apart (5 lines) |
| `nl_fit_off_saddle` | fitedge, nlstart, loglik | apart (29 to 78 lines) |
| `lp_basis_nl` | readers | near (one unchanged line between) |
| frmtmb.sample `NEWS.md` | readers | conflict (top-of-file heading, trivial) |

The other candidates the brief names are not touched by any other lane:
`newton_polish`, `nlminb_trial_fn`, `mo_search`, `VarCorr.frmtmb_fit` (loglik, nestquad
and latperf edit `methods-fit.R` elsewhere), `draws_unseen_level` (readers edits
`methods-draws.R` elsewhere) and `nl_restart_round` apart from loglik.

Trial merges. `dev/review-superassign-trialmerge.sh` runs `git merge-file`, pairwise, on
LF-normalized copies of every shared R file, and `dev/review-superassign-mergecheck.R`
reads the results. 28 merges: 27 are clean and 1 has conflicts (readers `R/me.R`, 2
hunks, both in `me_node_value`). Every clean merge parses, has 0 `<<-`, and has **0**
free variables that neither parent has. A merged function reads nothing that one side
moved. Readers' new `lp_basis_nl` code reads `parts` after the build loop
(`nl_basis_extra(object, parts, ...)`) and still works with the threaded accumulator.
Its merged body has `build(nm, parts)` and `return(parts)`.

`dev/review-superassign-trialmerge2.sh` merges the other lanes in sequence, with
superassign last and then first. In both orders superassign adds 0 conflicts in every
file except `R/me.R`, where it adds 2. The other conflicts it shows are between other
lanes and not this lane's (`R/fit.R`: nestquad 1, arstep 1; `R/par-template.R`:
loglik 1).

### 6. Tests

`dev/review-superassign-tests.sh`, one file per process, two at a time, `NOT_CRAN=true`,
`FRMTMB_BRMS_FIT_TESTS=true`. Every `lib:` line names wt-superassign-lib for the package
and for frmtmb (11 of 11). Log `dev/review-superassign-tests.log`:

    test-no-superassign.R          pass=14  fail=0 err=0 skip=0 warn=0
    test-priors-autocor-classes.R  pass=182 fail=0 err=0 skip=0 warn=0
    test-priors-bounds-grcov.R     pass=53  fail=0 err=0 skip=0 warn=0
    test-compat.R                  pass=759 fail=0 err=0 skip=0 warn=0
    test-nl.R                      pass=92  fail=0 err=0 skip=0 warn=0
    test-me-newdata.R              pass=178 fail=0 err=0 skip=0 warn=0
    test-importance.R              pass=344 fail=0 err=0 skip=0 warn=0
    sample test-loo.R              pass=146 fail=0 err=0 skip=0 warn=0
    sample test-loo-predict.R      pass=153 fail=0 err=0 skip=0 warn=0
    sample test-sample-direct.R    pass=135 fail=0 err=0 skip=0 warn=0
    sample test-sample-seed.R      pass=10  fail=0 err=0 skip=0 warn=0

(`test-priors*.R` matches two files in this tree. `test-importance*.R` matches one.)

### 7. NEWS

There are four development sections, each under its own `# <pkg> (development
version)` heading (frmtmb, frmtmb.ode, frmtmb.sample, frmtmb.survival), and no version
numbers. No lane name and no `dev/` path appears in any of the four. Each section is
patch-level and states that the results do not change. See MINOR 2 for one overclaim.

## MINOR 1. The namespace scan does not reach the 209 function literals held in top-level registries

`superassign_ns_scan()` walks `ls(ns)` and judges only objects that are functions. A
function literal inside a top-level list or environment is never walked. Examples are
`autocor_core`, `covstruct_core`, `frmtmb_links`, `nl_stan_funs` and `start_ls_rules`.
`dev/review-superassign-literals.log` counts **209** such literals in core `R/`
(covstruct_core 109 including its sub-entries, frmtmb_links 60, nl_stan_funs 20,
start_ls_rules 17, autocor_core 3). `dev/review-superassign-nsreach.log`: 526 function
values are held in lists or environments in 9 objects of the installed namespace. Some
of these duplicate top-level functions, such as the S3 table and `family_registry`.

The mutant shows the effect: a `<<-` in `autocor_core$fcor$parse` is missed by the
namespace scan. In the monorepo the source scan catches it. Under `R CMD check` on a
tarball the source scan skips, so there only pkgcheck stands between the code and a
release. The test header says the namespace scan "runs wherever the package does",
which is true, but it reaches less than the header implies.

Fix (cheap): in `superassign_ns_scan()`, also descend into non-namespace lists and
environments held in the namespace, deduplicated by environment identity, as
`dev/review-superassign-nsreach.R` does. Add the registry case to the control test.

## MINOR 2. The core NEWS bullet overstates what guards the extensions

"The test `test-no-superassign.R` fails if `<<-` comes back, in core or in an
extension." The extension check is a source scan, so it runs only where core's tests run
from the monorepo checkout. Under an extension's own `R CMD check`, nothing scans the
extension. The findings say this ("What I did not do"), but the NEWS bullet does not.
Fix: "in core, and in an extension when the tests run from the source repository", or
drop the clause.

## MINOR 3. Merge conflict with readers in `me_node_value`

Readers rewrites `me_node_value()` (multi-direction `z`) and removes its memo. That
removes the `<<-` too. Resolution: take readers' side for both conflict hunks. The
result has no `<<-` and needs no environment. Whether dropping the memo costs anything
is a question for the readers review, not this one.

## MINOR 4. A broken comment line in `helper-sampling.R`

The rewrapped comment leaves "# Prepending a" alone on its own line, followed by
"# synthesized library was tried first ...". Rewrap the paragraph.

## MINOR 5. The test header's description of pkgcheck's rule is unclear

"so a `<<-` anywhere in R/ fails it, in a comment-free line of code or not". pkgcheck
counts parse tokens, so a `<<-` in a comment or a string does not fail it. The control
test shows exactly this. Suggested wording: "so any `<<-` in the code of R/ fails it; a
`<<-` in a comment or a string does not".

## Recommended merge order

Superassign touches more files than any other lane, but it is mechanical, and trial
merges show that its hunks merge in either order. Merge it **last**, after readers.
Resolve `R/me.R` by taking readers' `me_node_value()`. Then run
`tests/testthat/test-no-superassign.R` and `dev/review-superassign-mergecheck.R` (or the
free-variable check of `dev/review-superassign-globals.R`) on the consolidated tree.
Those two runs catch any `<<-` that another lane adds after this snapshot (the lanes are
still running), and any reference to a variable this lane moved. If a hunk does conflict
later, re-run the lane's rewrite on the merged function (`dev/superassign-rename.R`)
instead of hand-merging both sides. Between the other lanes, the trial merge shows
conflicts that are not this lane's: `R/fit.R` (nestquad, arstep) and `R/par-template.R`
(loglik).

## Not done, and why

- No direct probe of `log_lik.frmtmb_draws`'s `once` handler or of two
  `draws_laplace_watch` closures alive at once, because both need Stan draws. Each
  creates its state once per call, as before (section 1), and test-loo*.R pass.
- `R CMD check` and the full suites were not rerun. The lane's logs give them.
- `dev/ciharden-openblas.sh` was not run. The lane ran it successfully; the readers lane
  had the same script denied by the permission system.
