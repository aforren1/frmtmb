# Lane cdfname: the CDF slot `lcdf` renamed to `cdf` (round of 2026-11-16)

Item 3088 of `dev/test-backlog.md` (R8, low), from
`dev/hygiene-findings.md` "3045: the `lcdf` slot" and its question 1.
User decision 2026-11-14: rename, and refuse the old name with an error
that names `cdf`; no alias, no deprecation period.

Base: b255350e (frmtmb 0.83.0). Library: `C:/Users/adf44/source/r/wt-cdfname-lib`
(core, frmtmb.eam, frmtmb.survival, frmtmb.coupling from this worktree;
the other five extensions read from rellib-r22, see "Extensions" below).

## 1. Readers and writers of the slot, before the rename

Counted with `grep -rIc lcdf` on the worktree at b255350e
(`dev/cdfname-before-counts.txt`, generated). Corrected in fix round 1
from the review's count: outside every `dev/` directory, 49 files and
180 lines; 47 files and 171 lines without the NEWS files. The three
ported-suite scripts in `dev/` add one line each (that file's total,
52 files and 183 lines, includes them). `extensions/frmtmb.eam/dev/dev-findings.md`
has 5 more, historical; left as they are.

These include names that only contain the string and are not the slot
(left alone, see section 2): Stan's and brms's `poisson_lcdf()`,
`gamma_lcdf()`, `asym_laplace_lcdf()` in comments and compat text; TMB's
`oneStepPredict()` columns `nlcdf.lower` and `nlcdf.upper` in
`R/predict.R`; frmtmb.eam's `ddm_rt_lcdf2()` and `ddm_rt_lcdf_b()`.

The rename script changed 146 lines in 48 files (its run output,
`Rscript dev/cdfname-rename.R`, "RENAME done: 146 line(s)"):

| Where | Files | Lines | What |
|---|---|---|---|
| core `R/` | 9 | 57 | `frmtmb_family()`'s argument and slot, `fam_lcdf()`, the mixture's slot, the core families that fill it (families.R, families-extreme.R), the readers (objective.R censoring and truncation, predict.R OSA and `pit`, families-cdf.R, frame.R's refusal of `cens()`/`trunc()`), roxygen in ad-env.R, parse.R, compat.R |
| core `tests/` | 7 | 14 | test-brms-lik-cdfs, test-cdfs, test-cens-lccdf, test-cens-trunc, test-cens-windows, test-families-extreme, test-negbinomial2 |
| core `man/` | 4 | 13 | frmtmb_family, frmtmb_ad_overload, frmtmb_register_aterm, frmtmb_register_compat |
| core `vignettes/`, `inst/bcm/`, `SPEC.md` | 3 | 4 | |
| frmtmb.eam | 17 | 50 | every family that fills the slot (wiener, rdm, gddm, lba), the ndt seam (ddm-shared.R), the contaminant wrapper, compat rows in zzz.R, docs, tests |
| frmtmb.survival | 6 | 6 | cox(), royston_parmar(), two tests, one vignette |
| frmtmb.coupling | 1 | 1 | one compat row's text |
| `dev/` (ported-suite tooling) | 3 | 3 | brmsport-bin2-adapters.R (comment), brmsport-bin2-survey.R and brmsport-bin2-verdicts-2p.R (reason text) |

No registry holds the slot name: the compat registry holds text only,
and no `frmtmb_register_*()` takes a family. The latent, learn, ode,
sample and spatial extensions do not name the slot.

## 2. The rename

`dev/cdfname-rename.R` does it. Run from the package root:

    Rscript dev/cdfname-rename.R            # edit
    Rscript dev/cdfname-rename.R --dry-run  # list only

then `roxygen2::roxygenise()` on core and on each extension it touched.
It is idempotent (a second run on this tree: "RENAME done: 0 line(s)").
For the consolidation: run it on the merged tree after squant's family is
in; it renames that family's `lcdf =` with the rest. It keeps each file's
line endings (CRLF stays CRLF).

Rules, in order, per line:
1. `an lcdf` / ``an `lcdf` `` / `an \code{lcdf}` become `a cdf` and so on.
2. In package code (`*/R/*.R`, non-comment lines) `$lcdf` becomes
   `[["cdf"]]`: `$` on a `frmtmb_family` refuses a prefix, so on a
   hand-built family with no `cdf` element `fam$cdf` would be an error
   naming `cdf_log`. One hit: the mixture's component check.
3. `fam_lcdf` becomes `fam_cdf`; frmtmb.eam's `ddm_cont_lcdf` becomes
   `ddm_cont_cdf`.
4. `lcdf` with no letter, digit, `_` or `.` before it and no letter or
   digit after it becomes `cdf`. This takes `lcdf`, `$lcdf`, `"lcdf"`,
   `\code{lcdf}`, and frmtmb.eam's local `lcdf_rows`, `lcdf_edge`,
   `lcdf_upper` (all three return probabilities, so the `l` misled there
   too).

Not touched: NEWS files (history); `tests/testthat/test-cdf-slot-name.R`
and any line holding `old_cdf_slot` (the refusal must spell the old
name); prefixed names (`poisson_lcdf`, `nlcdf.lower`, `ddm_rt_lcdf2`,
`ddm_rt_lcdf_b`: Stan, TMB, and eam helpers that do return a log CDF);
`dev/` scripts other than the three the ledger runs (records of past runs
against the builds of their day; about 60 of them read `$lcdf` and will
not run on the new core as written); generated files
(`dev/brmsport-bin2-groups.tsv`, `dev/equiv-data/results.tsv`).

Hand edits after the script: the `@param cdf` text ("Despite the `l` in
its name" removed); `fam_cdf()`'s comment; the SPEC.md column alignment.

Acceptance grep, `grep -rIn '\blcdf\b'` outside `dev/` and `.git/`, after
the change: `R/families.R` 1 line (`old_cdf_slot <- "lcdf"`),
`tests/testthat/test-cdf-slot-name.R` 12 lines, and the NEWS files
(core 8, eam 7, survival 2). Prefixed names as listed above remain by
design.

## 3. The refusal, and where it sits

`old_cdf_slot`, `check_cdf_slot_name()` and `refuse_family_dots()` in
`R/families.R`. Error class `frmtmb_cdf_slot_renamed`; every message
names `cdf`.

Routes, each tested in `tests/testthat/test-cdf-slot-name.R`, and the
function that raised the error on each, from
`dev/cdfname-route-probe.R` (log `dev/cdfname-route-probe.log`):

| Route | Raised in |
|---|---|
| `frmtmb_family(lcdf = )` (and `custom_family()`, its alias) | `refuse_family_dots()`: `frmtmb_family()` gained `...` for this only; any other argument there is "unused argument" |
| `frm(family = fam)`, `bf(family = )`, `bf() + fam`, a list of families for `mvbf()` | `as_frmtmb_family()` |
| a `mixture()` component | `as_frmtmb_family()`, through `mixture()`'s `lapply()` |
| a `bf()` whose family was edited after `as_frmtmb_family()` saw it | `assemble_frame()`, which checks every response's family first: every fit, refit, `update()`, `predict(newdata)`, `simulate(newdata)`, prior and template path builds its frame there (7 callers) |
| a `family_finalize()` that returns a family with `lcdf` | `assemble_frame()`, right after the call |
| a read of `fam$lcdf` | `$.frmtmb_family` (else it returns NULL, which every reader takes as "no CDF") |

Registration: no `frmtmb_register_*()` takes a family object
(`autocor`, `covstruct`, `aterm`, `compat`, `frame_check`,
`prior_defaults`), and extension families reach core only as family
objects through the routes above, so registration is not a route.
Not covered: a `[["lcdf"]]` read by code outside core (a plain list
read; it returns NULL).

Evidence that the test pins the change: on rellib-r22 (0.83.0),
`dev/cdfname-run-test.R frmtmb tests/testthat/test-cdf-slot-name.R base`
gives `pass=2 fail=11 error=3` (log `dev/cdfname-t-slot-base.log`): the
old build accepts `lcdf` on every route and has no `cdf` slot
(`$cdf` is a prefix error naming `cdf_log`). On this lane's build:
`pass=19 fail=0 error=0 skip=0 warn=0` (`dev/cdfname-t-slot-new.log`).
The guard's condition absent: the same family with `cdf =` fits a
right-censored model and matches core `poisson()`'s log-likelihood
(last test of the file).

## 4. `lccdf` and `cdf_log`

- `lccdf`: returns `log(1 - F(q))`. The name is Stan's suffix for the
  log complementary CDF (`normal_lccdf`) and is accurate. Leave it.
- `cdf_log`: returns `list(lF = log F, lS = log S)`, both tails. "cdf on
  the log scale" is half of what it returns. Mildly misleading, but now
  that `cdf` is the probability scale, `cdf_log` reads as its log-scale
  sibling, which is how the docs present it. Recommendation: leave it.
  If renamed, `log_tails` says what it returns.

## 5. Extensions

- frmtmb.eam, frmtmb.survival: fill the slot; need the new core (the
  old core would read no CDF from them, and they cannot construct a
  family on the old core: `frmtmb_family(cdf = )` partial-matches
  `cdf_log` there). Floor: `frmtmb (>= <this round's version>)`.
- frmtmb.coupling: one compat row's text names the slot. No code change;
  floor unchanged.
- latent, learn, ode, sample, spatial: do not name the slot and build
  their families at call time (no family object is stored at install),
  so the rellib-r22 builds run unchanged on this core.

## 6. Compat rows changed (text only, no verdict change)

- frmtmb.coupling: `cross_wishart` / `cens()`.
- frmtmb.eam (`extensions/frmtmb.eam/R/zzz.R`, ten rows): `wiener` /
  `trunc()`; `gddm` / `cens()`, `trunc()`; `lba` / `cens()`, `trunc()`;
  `rdm` / `ndt_group()`, `cens()`, `trunc()`; `wiener_gng` / `cens()`,
  `trunc()`.
- core: the `inverse.gaussian` / `trunc()` row of `R/compat.R` (one
  word), and the `frmtmb_register_compat()` example.

The frame refusal of `cens()`/`trunc()` for a family with no CDF now
says "through the cdf argument of frmtmb_family()". The equivalence
record (`dev/equiv-data/results.tsv`, pairs p010, p011, p157) holds the
old wording; the verdicts ("refusal correct") do not depend on it. Not
remeasured. The ported-suite reason text in two generators changes; no
verdict row flips.

## 7. Tests

### One file per process, every suite

`dev/release/run-par.sh cdfname-all dev/cdfname-all.jobs 2` with
`NOT_CRAN=true FRMTMB_BRMS_FIT_TESTS=true
FRMTMB_RELLIB=<lane lib>` (log `dev/release/cdfname-all.log`, per-file
logs in `dev/release/cdfname-all-files/`): "CDFNAME-ALL ran 550 of 550".
Totals summed from the RESULT lines:

| Package | Files | Pass | Fail | Err | Skip | Warn |
|---|---|---|---|---|---|---|
| frmtmb | 350 | 28775 | 0 | 0 | 15 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1857 | 0 | 0 | 3 | 0 |
| frmtmb.latent | 13 | 700 | 0 | 0 | 2 | 0 |
| frmtmb.learn | 16 | 545 | 0 | 0 | 2 | 0 |
| frmtmb.ode | 11 | 549 | 0 | 0 | 1 | 0 |
| frmtmb.sample | 92 | 4651 | 0 | 0 | 1 | 0 |
| frmtmb.spatial | 8 | 634 | 0 | 0 | 0 | 0 |
| frmtmb.survival | 16 | 575 | 0 | 0 | 0 | 0 |

Every `lib:` line names frmtmb from the lane library (550 of 550);
eam, survival and coupling also load from it, the other five from
rellib-r22. The 29 skips are other tiers' gates only: 15
`FRMTMB_SCALE_TESTS`, 13 `FRMTMB_DRMTMB_FIT_TESTS`, 1 `FRMTMB_FUZZ`.

### One session per suite (alphabetical)

`bash dev/cdfname-sessions.sh` (`dev/release/run-session.R <pkg>
forward`, `NOT_CRAN=true`, lane library first; log
`dev/cdfname-sessions.log`, per-suite logs in `dev/cdfname-sessions/`):

    frmtmb          files=337/350 pass=25665 fail=0 err=0 skip=243 warn=0
    frmtmb.coupling files=11/11   pass=551   fail=0 err=0 skip=5   warn=0
    frmtmb.eam      files=33/33   pass=1851  fail=0 err=0 skip=5   warn=0
    frmtmb.latent   files=13/13   pass=700   fail=0 err=0 skip=2   warn=0
    frmtmb.learn    files=16/16   pass=474   fail=0 err=0 skip=13  warn=0
    frmtmb.ode      files=11/11   pass=549   fail=0 err=0 skip=1   warn=0
    frmtmb.sample   files=87/92   pass=4286  fail=0 err=0 skip=18  warn=0
    frmtmb.spatial  files=8/8     pass=551   fail=0 err=0 skip=6   warn=0
    frmtmb.survival files=16/16   pass=575   fail=0 err=0 skip=0   warn=0

The 13 core and 5 sample files without a row are the brms-gated files
(`test-brms-suite-*`, `*-brms.R`), skipped at file level without
`FRMTMB_BRMS_FIT_TESTS`; they ran in the gated per-file run above.
Every LEAK line carries a tag: "[with load of ...]" (an extension's
`.onLoad()` registrations, or a third-party package's environment
variable). No untagged LEAK line.

### OpenBLAS 0.3.26 with LAPACK

`bash dev/ciharden-openblas.sh 0.3.26 lapack` (Rblas: 74 symbols to
OpenBLAS; Rlapack: 616), `OPENBLAS_NUM_THREADS=4`, one file per process
through `dev/cdfname-run-test.R` (logs `dev/cdfname-ob-*.log`): the new
test file and the renamed CDF readers' tests.

    test-cdf-slot-name.R pass=19  fail=0 error=0 skip=0 warn=0
    test-cdfs.R          pass=234 fail=0 error=0 skip=0 warn=0
    test-cens-lccdf.R    pass=31  fail=0 error=0 skip=0 warn=0
    test-cens-trunc.R    pass=100 fail=0 error=0 skip=0 warn=0
    test-cens-windows.R  pass=38  fail=0 error=0 skip=0 warn=0
    test-negbinomial2.R  pass=69  fail=0 error=0 skip=0 warn=0

`dev/ciharden-out/` deleted afterwards. The change renames code and adds
a refusal; it changes no arithmetic, so 0.3.32 was not run.

### R CMD check --as-cran

`bash dev/cdfname-check.sh core frmtmb.eam frmtmb.survival
frmtmb.coupling` (built with vignettes in `dev/cdfname-check/`,
`NOT_CRAN=true`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, R_LIBS = lane
lib, rellib-r22, user lib; log `dev/cdfname-check.log`):

    frmtmb           Status: 1 NOTE   (tests 54 min, OK)
    frmtmb.eam       Status: 1 NOTE   (tests 11 min, OK)
    frmtmb.survival  Status: 1 NOTE   (tests 28 s, OK)
    frmtmb.coupling  Status: OK

Each NOTE is the known "Skipping checking math rendering: package 'V8'
unavailable" of the HTML manual. `dev/cdfname-check/` deleted
afterwards.

The rename script, rerun on the finished tree: "RENAME dry run: 0
line(s) in 1019 candidate file(s)".

## 8. NEWS, bump, floors

- Core: BREAKING bullet under `# frmtmb (development version)`.
  Minor bump (an exported argument renamed).
- frmtmb.eam, frmtmb.survival: BREAKING bullet; patch or minor bump at
  the consolidation's choice (the API of the extension does not change,
  only what it needs). Their DESCRIPTION floors MUST rise to this
  round's frmtmb version, set by the consolidation (not edited here:
  the version does not exist yet, and a floor above 0.83.0 would not
  install against the base). Why: their families now call
  `frmtmb_family(cdf = )`, and on 0.83.0 `cdf =` partial-matches the
  formal `cdf_log`, so the CDF is installed as the log pair. The base
  run of the new test shows what follows: a fit stops with
  "`[[.default`(lp, "base")`: subscript out of bounds". Loud, but the
  floor is what should stop it. eam today: `frmtmb (>= 0.82.0)`;
  survival: `frmtmb (>= 0.75.0)`.
- frmtmb.coupling: one bullet; no floor change.

## 9. Backlog

Item 3088 (R8) closed with one line under it in `dev/test-backlog.md`.
No new items filed.

## 10. Not done

- The ~60 historical `dev/` scripts that read `$lcdf` are not renamed.
- The eam helpers `ddm_rt_lcdf2()` / `ddm_rt_lcdf_b()` keep their names
  (they return log CDFs).

## 11. Fix round 1 (review `dev/reviews/2026-11-17-cdfname.md`, verdict FIX)

### Item 1 (MEDIUM): every reader of the slot refuses an old-shaped family

Choice: one helper at the readers, not one post-fit entry point. A
saved fit reaches the slot through several methods that build no frame,
and a new method would have to remember an entry check; a reader cannot
miss the check if the slot can only be read through it.
`fam_cdf_slot(fam)` (`R/families.R`) is `check_cdf_slot_name(fam)[["cdf"]]`.
Every reader of the slot in core now goes through it:

- `fam_cdf()` (so `osa_log_pair()`, `sim_cdf_pair()`, and the
  objective's censoring and truncation terms);
- `osa_has_cdf()` (`R/predict.R`), the route choice of
  `residuals(type = "osa")`;
- `sim_response()`'s `invert` test (`R/families.R`), before the
  truncated draws;
- the `cens()`/`trunc()` refusal in `assemble_frame()` (`R/frame.R`);
- `mixture()`'s component test.

The only other read is frmtmb.eam's `ddm_ndt_install()`, on a family it
has just built.

Post-fit methods checked, by reading what they call and by
`dev/reviews/cdfname-review-rds-b.R` on the reviewer's real 0.83.0 RDS
fits (log `dev/cdfname-fix1-rds-b.log`):

- `residuals(type = "osa")`: refuses on uncensored weibull, censored
  poisson and the censored custom family. On censored weibull it gives
  the 0.83.0 values (-0.41853283764 0.01319703446 -0.47996173141, as in
  the base log). That route never asks whether the family has a CDF:
  `cdf_default` short-circuits on the censoring block. It scores
  through the fit's stored tape.
- `update()` / refit / `predict(newdata)` / `simulate(newdata)`: refuse
  in `assemble_frame()`.
- `simulate()`, `posterior_predict()`-style draws and
  `conditional_effects(method = "posterior_predict")` go through
  `sim_response()`. Under `trunc()` they refuse (new test). Without
  `trunc()` they read no CDF and give the 0.83.0 draws.
- `logLik()`, `fitted()`, `predict()`, `summary()`, `vcov()` and
  `confint()` do not read the slot. The review measured `logLik`,
  `fitted`, `predict` and `simulate` equal to 0.83.0.

Test: `tests/testthat/test-cdf-slot-name.R` builds the old-shaped fit
in the test with `old_shape()`, which renames the `cdf` element of
every `frmtmb_family` inside the fit. No RDS file is used.
- "an old-shaped fit is refused where a method reads the CDF": an
  uncensored weibull fit, seed 3088, n = 120. The unmodified fit's
  OSA residuals are all finite (the case where the guard's condition is
  absent); the old-shaped one refuses.
- "an old-shaped truncated fit is refused where draws invert": poisson
  with `trunc(ub = 8)`, seed 3089; `simulate()` refuses.

Seen failing first on the round-0 lane build, before this round was
installed (`dev/cdfname-fix1-t-before.log`): pass=21 fail=2 error=2.
- The OSA test errored, unclassed: "failed in oneStepPredict(): $
  operator is invalid for atomic vectors", after 120 `integrate()`
  "non-finite function value" messages.
- The truncated `simulate()` returned without an error.
- The two item 2 and item 4 tests failed too.

After this round: pass=27 fail=0 error=0 skip=0 warn=0
(`dev/cdfname-fix1-t-after.log`), with no `integrate()` message.

### Item 2 (LOW): `$cdf` on a family without the element

- `$.frmtmb_family` returns NULL for a name that is a formal of
  `frmtmb_family()` (every slot it defines) when the object does not
  hold it. Any other prefix is still refused (`f$cdf_l` names
  `cdf_log`, tested).
- frmtmb.eam's `ddm_ndt_install()` now writes `fam[nm] <-
  list(raw[[nm]])`, so it keeps a NULL element.
- Tests: core "an absent `cdf` element reads as NULL, not as
  `cdf_log`", and eam `test-rdm-gng.R` "go/no-go's absent cdf reads as
  NULL, also after the ndt seam".
- The review's eam probe rerun (`dev/cdfname-fix1-eam.log`):
  "wiener_gng has cdf element: TRUE", "after ndt install has cdf
  element: TRUE", "finalized $cdf value NULL". Before the fix
  (`dev/reviews/cdfname-review-eam.log`) both were FALSE and `$cdf` was
  the `cdf_log` error.

### Items 3 to 5

- Item 3: the counts in section 1 are corrected.
- Item 4: `refuse_family_dots()` now counts unnamed extras ("unused
  arguments: `x`, 2 unnamed"). Tested.
- Item 5: the floors are in section 8.

NEWS: core bullet extended (old saved fits); new core bullet for `$`;
eam bullet extended.

### Reruns

Per file, the five suites the change reaches, with `NOT_CRAN=true` and
`FRMTMB_BRMS_FIT_TESTS=true` (`dev/release/run-par.sh cdfname-fix1
dev/cdfname-fix1.jobs 2`; log `dev/release/cdfname-fix1.log`): "CDFNAME-FIX1
ran 502 of 502".

| Package | Files | Pass | Fail | Err | Skip | Warn |
|---|---|---|---|---|---|---|
| frmtmb | 350 | 28783 | 0 | 0 | 15 | 0 |
| frmtmb.coupling | 11 | 551 | 0 | 0 | 5 | 0 |
| frmtmb.eam | 33 | 1861 | 0 | 0 | 3 | 0 |
| frmtmb.sample | 92 | 4651 | 0 | 0 | 1 | 0 |
| frmtmb.survival | 16 | 575 | 0 | 0 | 0 | 0 |

Every `lib:` line names frmtmb from the lane library. Skips are other
tiers' gates: 13 drmTMB, 10 scale, 1 fuzz.

One session per suite, alphabetical (`OUT=dev/cdfname-sessions-fix1
EXTS="eam survival coupling sample" bash dev/cdfname-sessions.sh`; log
`dev/cdfname-sessions-fix1.log`): frmtmb 337/350 pass=25673, coupling
11/11 pass=551, eam 33/33 pass=1855, sample 87/92 pass=4286, survival
16/16 pass=575; all fail=0 err=0 warn=0. The missing files are the
brms-gated ones, as before. The LEAK lines are the same tagged set as
round 0. There is no untagged LEAK line.

OpenBLAS 0.3.26 with LAPACK (Rblas 74 symbols, Rlapack 616;
`OPENBLAS_NUM_THREADS=4`; logs `dev/cdfname-ob-fix1-*.log`):
test-cdf-slot-name.R pass=27, eam test-rdm-gng.R pass=225,
test-osa-cdf.R pass=71, test-osa-inference.R pass=34,
test-cens-trunc.R pass=100; all fail=0 error=0 skip=0 warn=0.
`dev/ciharden-out/` deleted.

R CMD check --as-cran on the two packages changed this round
(`bash dev/cdfname-check.sh core frmtmb.eam`; log
`dev/cdfname-check.log`, round 0's kept as `dev/cdfname-check-r0.log`):
frmtmb "Status: 1 NOTE", frmtmb.eam "Status: 1 NOTE". Each NOTE is the
V8 math-rendering one. frmtmb.survival and frmtmb.coupling did not
change this round; round 0 gave them 1 NOTE and OK.
`dev/cdfname-check/` deleted.

Rename script rerun: "RENAME dry run: 0 line(s)". Acceptance grep is
unchanged: `\blcdf\b` outside `dev/` is in `R/families.R` (1, the
`old_cdf_slot` line), the refusal test (13) and the NEWS files.
