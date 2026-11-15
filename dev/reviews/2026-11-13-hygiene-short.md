# Short check of lane hygiene's fix round 1 (round of 2026-11-13, base ec8178d3)

Checks only the fixes for F1 to F5 and the nits of
`dev/reviews/2026-11-13-hygiene.md`. Scripts and logs are
`dev/reviews/hygiene-short-*`. No lane file was edited. Git: `show`,
`hash-object`, `rev-list` on a bare clone; the only git writes went to a
scratch repository and a bare clone in the session scratchpad.

## Verdict: FIX

F2, F3, F4, F5 and the nits are fixed. F1 is fixed for the review's
case only. The key still merges two closures that one factory made for
different values when the value is not in the closure's own
environment, or when the value is an environment or an external
pointer. When that happens, the second registration replaces the first
and its refusal is lost with no message. The exported documentation and
NEWS say that this cannot happen.

## 1. F1, `registration_key()` (R/registries.R)

The installed lane build has the source's code (`INSTALLED
registration_key == lane source: TRUE`, `register_once` TRUE). Script:
`hygiene-short-regkey.R <mine|base>`. Logs:
`hygiene-short-regkey-{mine,base}.log`. Each pair registers two checks.
The first check refuses the fit and the second does nothing.

    case                                          base: added, fit   mine: added, fit
    factory inside a factory, outer value differs      2, refused         1, fitted
      same, registered in reverse order                2, refused         1, refused
    local() inside a factory (one factory)              2, refused         1, fitted
    two closures on two empty environments,
      first environment set to refuse after             2, refused         1, fitted
    two closures on two different environments          2, refused         2, refused
    two closures on two external pointers               2, refused         1, fitted
      (identical(p1, p2) FALSE; serialize() equal TRUE)

(b) The nested-factory failure is a merge, not a duplicate. So is the
plain `local()` idiom: `function(mode) local(function(spec, frame)
...)` is one factory, and the closure's own environment is the empty
`local()` frame. `?frmtmb_register_frame_check` ("Two closures that one
factory made for different values are two checks, and both run") and
the NEWS bullet both say the opposite. The note in the `@noRd` roxygen
about parents is not visible to an extension author.

(c) `serialize()` writes an environment by its contents and an external
pointer without its address. So two distinct environments with equal
contents, or two distinct pointers, give equal keys. The key is unique
only when serialization fails, and these values do not make it fail.

(d) The same closure, plain and byte-compiled: 1 (correct). A captured
helper, plain and byte-compiled: 2. The same code with and without
srcref: 1. A captured helper with and without srcref: 2. The two cases
that give 2 are duplicates, which are harmless.

(e) `pkgload::load_all()` of frmtmb.ode twice (from a copy in the
scratchpad): checks 3, 3. `pkgload::unload()` + `library()`: 3.
`unloadNamespace()` + `library()`: 3. Base: 4, 5, 6, 7. frmtmb.sample
`load_all()` twice: providers 2, 2. `unload()` + `library()`: 2.
`unloadNamespace()` + `library()`: 2. Base: 2, 3, 4, 5. The one extra
provider in mine comes from `load_all()` after `library()`: the
provider captures `libname`, which differs between the two. This is a
bounded duplicate, not a merge. The review's reproducer gives the
lane's numbers again (`hygiene-short-regonce-{base,mine}.log`):

    base  FACTORY checks 2, refused; providers 2; RELOAD ode 4; sample 2
    mine  FACTORY checks 2, refused; providers 2; RELOAD ode 3; sample 1
          LOAD eam+ode+sample 0.28 s, +4.5 MB; COST 0.02 s, +0.0 MB

(a) Cost, with one check that captured a 49.6 MB data frame:

    mine  register it                    0.140 s, peak +50 MB
          register it again              0.530 s, peak +100 MB
          any later registration         0.120 s, peak +50 MB each
          10 small closures after it     1.300 s, peak +150 MB
          library(eam, ode, sample)      0.720 s, peak +212 MB
    base  the same five                  0.000 s (x4), 0.290 s, +94 MB
    a fit with it registered: 0.01 s (mine), not affected

For each registration, the code computes the key of every entry again,
so the cost is (entries) x (captured size). This does not happen today:
eam and ode register namespace functions, and sample's provider
captures two strings. LOW. Keep the key with the entry if this ever
matters.

Fix (recommended): go up from `environment(fn)` to the first named
environment, and serialize the bindings of each frame on the way. Make
the key unique when the `refhook` of `serialize()` sees a reference
object (an unnamed environment or an external pointer). A duplicate is
harmless; a merge is not. A prototype (`hygiene-short-altkey.R`, log
`hygiene-short-altkey.log`, not the lane's code) gives the wanted
answer in all 9 cases. It keeps apart: nested, `local()`, two
environments, two pointers, and one factory with two values. It merges:
one factory value twice, a nested factory with the same values, an
`.onLoad()` closure made twice, and a namespace function twice. Add the
`local()` and nested pairs to `test-register-once.R` as cases that must
stay two.

## 2. F2, `.github/scripts/coverage-plan.py`

Script: `hygiene-short-plan.py`. Log: `hygiene-short-plan.log`. The
scratch repository is in the scratchpad. The lane's script ran as the
workflow runs it (EVENT=push, BEFORE, HEAD):

    git mv R/y.R extensions/frmtmb.eam/R/y.R
      git default: extensions/frmtmb.eam/R/y.R; script: 2 files, core true, all eight
    tests/testthat/_snaps/print.md only            core true, all eight
    extensions/frmtmb.ode/tests/testthat/_snaps/o.md only
                                                   core false, ode sample
    README.md only                                 core false, none
    diff.renames=copies in the repository config, git mv R/z.R dev/attic/z.R
      git default: dev/attic/z.R; script: 2 files, core true, all eight

The 60-commit dry run (bare clone of main, tip ec8178d3, uses the
script's own `changed_files()`): `jobs run 398 of 540; all nine: 43;
parents up to 13`. Confirmed.

## 3. F3, the texts

`vignettes/inputs.Rmd:132-141`, NEWS.md ("Documentation") and the
header of `test-weights-constant.R` now apply the statement to a
maximum-likelihood fit with no random effect and no prior, and name
the terms that the weight does not multiply. The REML formula is
correct. The weighted gaussian log-likelihood is c times the per-row
log density. Integrating the fixed effects adds `-1/2 log|c X'X /
sigma^2|` = `p log sigma` + const, so the restricted log-likelihood is
`-(c n - p)/2 log sigma^2 - c RSS / (2 sigma^2)`. Its stationary point
is `sigma^2 = c RSS / (c n - p)`. At c = 2, n = 80 and p = 3,
`sqrt(2 x 77 / 157)` = 0.9904, which is the review's measured ratio.
ASCII only. No new line is longer than 80 columns in the changed hunks.

## 4. F4, F5, nits

- F4: `dev/hygiene-findings.md:520-532` says that the bin-1 suite is not
  in the plain tier, and quotes the gated counts of the review. Correct.
- F5: `dev/hygiene-findings.md:381-396` says 398 of 540 and 43, with the
  method. These match section 2.
- The 3164 CLOSED line (`dev/test-backlog.md:3185-3188`) says that the
  RNG half was declined, and why.
- `dev/release/checklist.md` step 3 says "once core's `DESCRIPTION`
  carries the release version". It no longer refers to a step that does
  not exist.
- The backlog has 16 lane-hygiene CLOSED lines (205 CLOSED lines in
  total, against 189 at ec8178d3). This matches the 16 closed rows of
  the findings' item table. The review's "17" was a miscount, so
  nothing to fix.
- After F1 is fixed again, the 3154 CLOSED line and NEWS must describe
  the final key.

## 5. The frmtmb.latent forward session

`hygiene-short-latent.sh` uses the lane library and
`R_PARALLEL_PORT=12417`. Log: `hygiene-short-latent.txt`.
`SESSION frmtmb.latent forward files=13/13 pass=697 fail=0 err=0
skip=2 warn=0`, `EXIT 0 after 116 s`, no LEAK line. The session
completes.

The lane's changes do not explain the AD-overloads leak:

- RTMB attaches "AD-overloads" in `attachADoverloads()`.
  `detachADoverloads()` removes it only outside an AD context
  (`hygiene-short-adover.log`). So the leak means that a tape was cut
  off partway.
- frmtmb.latent registers no frame check and no prior default. Only a
  comment in `zzz.R` names the seam. Its tests do not call
  `local_frmtmb_registries()`. In the dead run, test-hmm.R loaded no
  frmtmb.* namespace. `register_once()` and the restore of the hooks
  therefore did not run in that file. The `compat.R` rows do not touch
  RTMB taping.
- The dead run differs in test-hmm.R itself: 11.9 s against 48.2 s now.
  hmmTMB, rstan, tmbstan and ggplot2 never loaded, and the hmmTMB
  blocks printed no convergence lines. Both hmmTMB blocks were skipped
  or failed before the reference fit, which is consistent with a
  process that was short of resources on a busy machine. The runner
  prints DETAIL lines only at the end, so the cause of the dead run is
  lost. That is a runner gap, not a lane defect. To keep the cause next
  time, print DETAIL lines for each file.

## 6. Test files on the lane library

`hygiene-short-tests.sh`, `hygiene-short-tests.log`. One file for each
process. Every `lib:` line names wt-hygiene-lib:

    test-register-once.R     pass=14  fail=0 err=0 skip=0 warn=0
    test-registry-helper.R   pass=10  fail=0 err=0 skip=0 warn=0
    test-structure.R         pass=128 fail=0 err=0 skip=0 warn=0
    test-compat-register.R   pass=103 fail=0 err=0 skip=0 warn=0
    test-weights-constant.R  pass=8   fail=0 err=0 skip=0 warn=0
    test-codemeta.R          pass=4   fail=0 err=0 skip=0 warn=0

## Re-check after fix round 2 (F1 only)

### Verdict: FIX (one line)

The walk stops at ANY named environment, but only namespaces, the
global environment, base and empty are safe places to stop. A closure
whose own environment or parent has a name (`attr(e, "name")`, an
`attach()`ed list) is keyed by code alone, and two different checks
merge. The refusal is then lost. This case is contrived, but it breaks
the rule "merge only what is provably the same". It also differs from
the exported documentation, which says "up to the first namespace,
global or base environment". Everything else holds.

The installed code is the lane source: a deparse comparison of
`registration_key`, `register_once` and `registration_named_env` gives
TRUE for each. The Rd file is newer than R/. Two R processes at most.

### 1. The earlier cases on the final build

Log: `hygiene-short-regkey-mine-r2.log`, made with
`hygiene-short-regkey.R mine`.

    nested factory, either order          2, refused
    local() in one factory                2, refused
    two empty environments, then refuse   2, refused
    two environments, different contents  2, refused
    two external pointers                 2, refused
    same closure plain / compiled         1;  captured helper 2
    same code srcref / none               1;  captured helper 2
    load_all(ode) x2, unload()+library(), unloadNamespace()+library()
                                          checks 3, 3, 3, 3
    load_all(sample) x2, unload()+library(), unloadNamespace()+library()
                                          providers 2, 2, 2, 2 (one is the
                                          documented libname extra)

### 2. Cost with a 49.6 MB capture

    register it                   0.110 s, peak +50 MB
    register it again             0.110 s, peak +51 MB
    a small closure after it      0.000 s, +0 MB
    10 small closures after it    0.000 s, +1 MB
    library(eam, ode, sample)     0.160 s, +95 MB (base, round 1: 0.290 s, +94 MB)

Later registrations now cost nothing, because each key is stored on its
entry.

### 3. Can it still merge two different checks?

Log: `hygiene-short-regkey2.log`, made with `hygiene-short-regkey2.R`.
In each pair the first check refuses. "added 1" means the two checks
merged.

    environment deep in a nested list, 2 envs, first turned to refuse  2, refused
    the same nested list object twice                                   2 (duplicate)
    S4 slot TRUE vs FALSE                                               2, refused
    S4 equal slots twice                                                1 (same)
    reference class, 2 objects, first turned to refuse                  2, refused
    R6-style environment object, same                                   2, refused
    value two unnamed frames up differs                                 2, refused
    equal chain values, first chain's frame turned to refuse after      1, fitted
    own frame `state <- FALSE`, first set TRUE via assign() after       1, fitted
    own env named "cfg", values differ                                  1, fitted   <- merge
    own envs named "cfgA" vs "cfgB", values differ                      1, fitted   <- merge
    parent env named "cfg", values differ                               1, fitted   <- merge
    -0 vs 0 (bytes differ)                                              2, refused
    NA_real_ vs NaN (bytes differ)                                      2
    NaN vs -(0/0) (bytes equal, identical() TRUE)                       1 (same value)
    1:3 vs c(1L, 2L, 3L); UTF-8 vs latin1 of one string; 1L vs 1        2 each (duplicates)

- Named environments are a merge, and it is harmful. This is the FIX.
  The walk stops when `nzchar(environmentName(env))`, and it compares
  neither the values in that environment nor its name or identity.
  Fix: stop only at a namespace, the global environment, base or
  empty, and make the key unique (a fresh environment) when the walk
  meets any other named environment. Do not walk on from it, because
  above an attached environment are the search path's package
  environments. Add the "cfg" pair to `test-register-once.R`. No
  extension registers such a closure, so nothing else moves.
- Two frames with equal values that are mutated from outside after
  registration can merge. This follows from the rule: a reload must
  merge on equal frame values, because `.onLoad()` makes new frames.
  The exported doc states the rule ("hold the same values"). It is
  harmful only if code assigns into a registered closure's environment
  later. Accept it. Optionally, the doc can say "at registration".
- `serialize()` writes doubles bit for bit, so -0 and 0, and NA and
  NaN, stay apart. The only equal-bytes pair found is two NaNs that
  `identical()` also calls equal, so that merge is harmless. ALTREP
  against plain vectors, encodings, types and a nested list that holds
  an environment give duplicates only.

### 4. Docs and NEWS

`?frmtmb_register_frame_check` (Rd regenerated), NEWS, the comments in
`R/structure.R` and `R/sampling-api.R`, and the 3154 CLOSED line
describe the rule. Two statements are not exact:

- "up to the first namespace, global or base environment": the code
  stops at any named environment. The FIX above makes the code match
  this text.
- "A reload ... registers exactly that: ... a closure its `.onLoad()`
  makes from the same arguments": this is not exact for two idioms of
  that kind (`ONLOAD` lines). `check <- function(...) ...;
  frmtmb_register_frame_check(check)` puts the closure in its own frame.
  A closure next to a formula carries an environment. Each of these is
  added again on every reload (added 2). An anonymous closure merges
  (added 1). None of the eight extensions uses these idioms. A
  duplicate is the harmless direction, but the doc should say "unless
  that frame also holds a function, a formula or another value that
  carries an environment".

### 5. test-register-once.R

Logs: `hygiene-short-ro-r2-{mine,base}.log`. One process each, and each
`lib:` line names the library under test.

    wt-hygiene-lib  pass=19 fail=0 err=0 skip=0 warn=0
    rellib-r21      pass=11 fail=8 err=0

Each failure on rellib-r21 is a length one higher than expected,
because base appends:

- the reload of `mk()` (2 vs 3, 3 vs 5)
- the factory value registered again (2 vs 3)
- the nested pair registered again (2 vs 3, 4 vs 5)
- the three provider counts (1 vs 2, 2 vs 3, 4 vs 5)

These are the right reasons. Nit: `expect_identical(strip(fns[[n0 +
1L]]), second)` also passes on base, where that entry is the FIRST
`mk()` closure. It proves the position relative to `other`, not that
the new closure replaced the old one. The length checks cover that.
