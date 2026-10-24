# Review of lane mvnlatent (`mixture_mvn()` into frmtmb.latent), 2026-10-23

Adversarial review. Worktree `C:\Users\adf44\source\r\frmtmb-wt-mvnlatent`,
base 19b162b9 (frmtmb 0.73.1), all changes uncommitted. Lane build
`C:/Users/adf44/source/r/wt-mvnlatent-lib`, read only; base build
`rellib-r11`, read only. No package file was edited. The reviewer's
scripts are `dev/review-mvnlatent-*.R` and `.sh`; their logs are in
`dev/review-mvnlatent-out/` (gitignored).

## Verdict

**FIX.** The move is complete and I could not break it. The code moved
without a change, and the fits are bitwise equal: I ran them again. No
core code line names the family. The compat table holds the same rows.
Every moved expectation is accounted for. I found no blocker. One major
finding is work for the consolidation, not for the lane: the floor in
frmtmb.latent's DESCRIPTION is known to be false. The other findings
are about docs and the record: the new export's help page claims
glmmTMB's parameterization, which is false at d >= 4; the seam has no
inverse, so the moved test forks it; core-boundary.md names a
prerequisite that the lane did not meet (correctly); and the record
leaves out some numbers.

On the numbers, I agree with the lane: this is a weak candidate under
guardrail 3. The move takes 316 lines out of core in net, removes 0.3%
of the core suite's time, and drops one Suggests (mclust). Its cost is
a release in lockstep (finding 5). I would keep it, because the rule's
question 3 decides it and `us_sigma()` is a useful export on its own.
DROP-MOVE is a defensible choice for the user.

## First, the instrument

`dev/review-mvnlatent-libcheck.R` deparses every function in the
worktree's `R/` and compares it with the installed namespace in the
lane library (log `libcheck.txt`):

```
frmtmb 0.73.1 C:/Users/adf44/source/r/wt-mvnlatent-lib : source fns 1978  differing 0  missing 0
frmtmb.latent 0.6.2 C:/Users/adf44/source/r/wt-mvnlatent-lib : source fns 67  differing 0  missing 0
frmtmb.sample 0.21.0 C:/Users/adf44/source/r/wt-mvnlatent-lib : source fns 285  differing 0  missing 0
```

So the lane library is the worktree source.

## What I verified and found correct

- **Moved code is identical.** The non-comment lines of base
  `R/families.R` 7229-7586 and of `extensions/frmtmb.latent/R/mixture-mvn.R`
  (up to the compat section) differ only in two blank lines. The nine
  rule texts are identical (sorted diff of the `r(...)` lines).
- **Fits.** I ran `dev/mvnlatent-fitcompare.R` again, once on
  rellib-r11 and once on the lane library, in parallel processes:
  `frmtmb 0.73.0 .../rellib-r11 ; mixture_mvn from frmtmb` and
  `frmtmb 0.73.1 .../wt-mvnlatent-lib ; mixture_mvn from frmtmb.latent`,
  59 fits each. `identical(A, B)` on the whole saved list is TRUE. Mine
  are also identical to the lane's `fits-base1.rds` and `fits-lane1.rds`.
- **Compat sets.** I wrote an independent comparison
  (`dev/review-mvnlatent-compat.R`; unordered pair key plus status and
  note). It uses six dumps: base and lane, each with core only, with
  core plus latent, and with all seven extensions:
  all 11350 = 11350 rows, 0 only-before, 0 only-after, vocabulary
  164 = 164; latent 7735 = 7735, 0 / 0, 134 = 134; core with the
  mixture_mvn rows dropped from before 7339 = 7339, 0 / 0, 131 = 131.
  The columns are the same (`feature_a,kind_a,feature_b,kind_b,status,note`).
- **Guardrail 3.** Recounted from `git show 19b162b9:R/families.R`
  lines 7229-7586: code 235, comment 117, blank 6 (the lane counts 3,
  because it leaves out the 3 separators), 358 in all. R/ has 86204
  lines in 70 files, and `mixture_mvn` has 38 mentions in 9 files.
  Churn over `4f5ea39f..19b162b9` (37 commits): `git log -L` on the
  block finds 0 commits, and on the test file 0. `-G` finds 3 commits
  (7b2bf010, afb1cb27, 83708f40), and each changed only a comment or a
  roxygen line. mclust: only `test-mvn-mixture.R`, the srr sentence in
  `helper-reference.R`, and prose. `grep -c mixture_mvn
  dev/brmsport-verdicts.tsv` gives 0. The time share, recomputed from
  `before.log`: 262 files, a sum of 2045.6 s, mvn file 6.3 s, 0.308%.
- **No core branch names the family.** `grep -rn -i "mvn\|sigmaraw"
  R/` on the lane tree finds only `mvn_draw_cov` (a generic
  multivariate normal draw) and comment or roxygen prose. Base had no
  code line outside the family and compat.R either (`git grep`). The
  core reads of `fam[["mix"]]` (`K`, `log_pi`, `comp_lpdf` with an
  optional fifth `extra` argument, `comp_dpars`, `comp_sim`, `ord`,
  `ref`) are generic.
- **Guardrail 2.** frmtmb.latent's `R/`, `NAMESPACE`, `DESCRIPTION` and
  `vignettes/` have no `frmtmb:::`, `getFromNamespace` or
  `asNamespace`. A codetools scan of the four moved functions found
  that every frmtmb export they use is imported:
  `compat_rule_builder frm_stop frmtmb_family frmtmb_structure
  mixture_multimodal_refusals mixture_posterior structure_supports_all
  us_sigma`, with none used but not imported. So the family works with
  frmtmb.latent loaded and not attached (next item).
- **A case the lane did not test** (`dev/review-mvnlatent-cases.R`, log
  `cases.txt`). I attached frmtmb only and called
  `frmtmb.latent::mixture_mvn()`. The namespace was loaded and not
  attached. The `family =` and `bf() +` spellings give identical
  logLik. `predict()`, `fitted()`, `mixture_probs()`, `pp_mixture()`,
  `simulate()`, `conditional_effects(dpar = "theta1")`, `summary()`,
  `confint()`, `frm_allfit()`, `anova()`, `update()` (logLik
  identical), `update(formula. = . ~ 1)` and `update(newdata =)` all
  work. The REML refusal names the family, and
  `conditional_effects(method = "predict")` refuses by name. I also
  saved a fit made with `library(frmtmb.latent)` and an unqualified
  `mixture_mvn()` call, and read it back in a session that attaches
  frmtmb only (`dev/review-mvnlatent-reload.R`, log `reload.txt`).
  `readRDS()` loads frmtmb.latent, and `predict()`, `simulate()`,
  `mixture_probs()`, `frm_compat()`, `update(fit)` and
  `update(newdata =)` all work.
- **The moved test's start helper.** For all 8 models,
  `mvn_start_vector(fit, st)` is identical to
  `unlist(frmtmb:::make_start(fit$frame, st))`. The template order
  matches the order of `obj$par`.
- **`us_sigma()` on the tape.** I compared the Jacobian with central
  differences (h = 1e-6, seeds 1 and 3). The maximum relative error is
  4.77e-11 at d = 1 and 6.45e-11 at d = 3, with `floored` both FALSE
  and TRUE. A wrong length is refused while taping ("a 2 x 2 matrix
  takes 3 parameters, not 4"). The two core callers (`gr_cov`,
  `gr_prec`) handle d = 1 before they call, and they pass
  `blk[["dim"]]`-sized segments. So neither behavior change reaches a
  core caller, and the fits above confirm it.
- **Tests, lane library first** (`dev/review-mvnlatent-par.sh`; every
  `lib:` line shows frmtmb, frmtmb.latent and frmtmb.sample from
  `wt-mvnlatent-lib`):
  `test-core-boundary.R` pass=10, `test-us-sigma.R` 18,
  `test-compat.R` 771, sample `test-simulators.R` 54, and every
  frmtmb.latent file: bracket-access 1, conditions-census 11,
  conditions 6, hmm-aterms 2, hmm-starts 90, hmm 119, lca 105,
  message-uniqueness 4, mixture-mvn 308, prior-dpar 13, scale 0
  (2 skips, the scale gate), structure-latent 20. Fail, error and warn
  are 0 everywhere. In one session (`dev/release/run-session.R`):
  frmtmb.latent forward, reverse and preload, each 12/12 files,
  pass=679, fail 0, err 0, skip 2, no LEAK line.
  frmtmb.sample reverse: 59/64 files, pass=3022, fail 0, err 0,
  skip 9. Every LEAK line is tagged `[with load of ...]`.
  `test-simulators.R` leaves "compat features: 1 -> 4" and "compat
  rules: 1 -> 4" (the lane reports the features line only).
- **Test accounting.** The moved file differs from
  `git show 19b162b9:tests/testthat/test-mvn-mixture.R` in exactly two
  changed lines (250 and 305), plus the added header and helpers
  (3-43) and the added blocks (423-540). The per-file counts I ran
  match `count.txt` (308, 13, 10, 18, 771). 308 = 277 + 3 (REML and
  profile) + 7 (structure) + 6 (unpinned seams) + 5 ("this package's")
  + 10 (compat rows). The core `test-structure.R` loop body is the
  same 7 expectations that the latent block asserts.
- **Spellings on both builds** (`dev/review-mvnlatent-strings.R`, log
  `strings.txt`, four arms). I checked `family = "mixture_mvn"`,
  `"Mixture_MVN"`, `brmsfamily("mixture_mvn")` and
  `brmsfamily("mixture_mvn", K = 2, D = 2)`. Each gives "is not a
  supported family" on base and on lane, with and without the
  extension. `family = mixture_mvn` (the function, not a call) gives
  the constructor's own "needs K >= 2" refusal wherever the function is
  visible. The lane's table in section 3 is correct.
- **0.6.2 over the new core.** Confirmed (`dev/review-mvnlatent-oldlatent.R`,
  log `oldlatent.txt`): ".onLoad failed in loadNamespace() for
  'frmtmb.latent' ... the rule `hmm x mixture_mvn` names 'mixture_mvn',
  and the registry has no such feature".
- **Docs.** `pkgdown::check_pkgdown()` gives "No problems found" for
  core and for frmtmb.latent, with pandoc on PATH. `_pkgdown.yml`
  drops `mixture_mvn` and adds `us_sigma`. The latent site has no
  reference index, so it lists every topic. `?mixture` names the
  package and `library(frmtmb.latent)`. Core NEWS is under the
  development heading, marked BREAKING, names the package and its
  install line, and has no lane name, `dev/` path or version number.
  The README moves the family to the frmtmb.latent row and paragraph.
- **Merges with the other lanes.** I ran trial `git merge-file` merges
  (CR stripped) against `frmtmb-wt-spatial`. They are clean for
  `R/covstruct.R`, `NAMESPACE`, `R/families.R`, `R/compat.R`,
  `_pkgdown.yml`, `R/predict.R`, `R/structure.R`, the simulate files,
  `R/priors.R`, `DESCRIPTION`, `test-core-boundary.R`, `README.md` and
  `man/frmtmb-extension-api.Rd`. Only `NEWS.md` conflicts, because both
  lanes add a top section. Against `frmtmb-wt-splitsurv`, one hunk of
  `R/compat.R` conflicts: two adjacent deletions near the
  `cens() x cox` rule. To resolve it, keep the cox rule where
  splitsurv leaves it and delete the two mixture_mvn rules.

## Findings

### 1. minor: `?us_sigma` claims glmmTMB's parameterization, which is false at d >= 4

The new help page says "in glmmTMB's parameterization". It then states
that the correlation entries fill the strict lower triangle of `L` in
column-major order, and that statement is correct for frmtmb. glmmTMB
fills the UPPER triangle of `R` in column-major order and takes
`crossprod(R)`, which is `L`'s lower triangle in ROW-major order. The
two orders agree for d <= 3 and differ from d = 4 up
(`dev/review-mvnlatent-glmmtmb.R`, log `glmmtmb.txt`, glmmTMB 1.1.15.2,
seeds 3 and 4):

```
d = 3  equal to glmmTMB get_cor: TRUE
d = 4  equal to glmmTMB get_cor: FALSE
```

A third party who feeds glmmTMB's `theta` into the new export gets a
different matrix, with no message, at d >= 4. The internal function
had the same order before, so it is not a regression. But the claim is
new and public. `test-us-sigma.R` pins the order only at d = 2, where
there is one entry and so no order to pin.

Fix: drop "in glmmTMB's parameterization" from `?us_sigma`, or say
"the same map as glmmTMB's for d <= 3; from d = 4 the entries are in
column-major order of `L`'s lower triangle, where glmmTMB uses
row-major order". Add an expectation at d = 4 with distinct nonzero
entries that pins the documented order against a hand-built `L`.

### 2. minor: the seam has no inverse, and the moved test forks it

The reason given for the export is that "a copy would fork the
parameterization". But the moved test needs the inverse map
(covariance to `theta`), which is not exported, so it carries a copy,
`mvn_theta_cor()` (`test-mixture-mvn.R` lines 36-41). `mclust_extras()`
uses that copy to start the fit at mclust's solution. A third party who
starts a free-covariance family from a sample covariance needs the same
inverse. So does the spatial lane's `frmtmb_register_covstruct()` for a
`us`-like structure: its `from_natural =` needs the inverse, and its
`chol =` needs the factor `us_chol_L()`. Both are core internals.

`us_sigma()` does not conflict with the spatial seam. The name is not
in the spatial lane's NAMESPACE, and the covstruct.R merge is clean.
But it belongs to that seam: the two are the "AD-safe pieces of an
unstructured block" that a registrant uses.

Fix (at consolidation, or file it as a follow-up): export the inverse
(for example `us_theta(S)`, which returns the `theta` whose
`us_sigma()` is `S`) on the same help page. Link `?us_sigma` and
`?frmtmb_register_covstruct` to each other. Then the moved test can use
the export in place of its copy.

### 3. minor: the lane did not meet core-boundary.md's stated prerequisite (correctly), and the user must amend the document

`dev/core-boundary.md`, "Proposed moves", names "an
`extension_family_homes` entry" as the prerequisite for this move. The
lane added none. Its facts are correct (strings.txt): no string
spelling reached the family on rellib-r11, so an entry would not mend
a gap. An entry would ADD spellings. `family_from_name()` passes named
arguments through, so with an entry `brmsfamily("mixture_mvn", K = 2,
D = 2)` would start to work, and `family = "mixture_mvn"` would give
the K/D refusal. The lane says only the second ("which only errors"),
so its argument is incomplete. That is still no reason to add the
entry. I agree with no stub and no entry, for reasons 1-3 of the
findings. The only way that a user without the extension can meet the
family is the call `mixture_mvn(...)`, and that fails in R before core
runs: "could not find function". No core mechanism can intercept it
short of a stub that masks.

Fix: in the record, correct the argument to say that an entry would
create `brmsfamily("mixture_mvn", K =, D =)`. Ask the user to change
the core-boundary.md row to "none; no string spelling reached it".

### 4. minor: the record leaves out some numbers

- The net change to core `R/` is +101 / -417 = **-316 lines**
  (`git diff --numstat 19b162b9 -- R/`), not the 355 moved.
  `us_sigma()` and its docs add 56 lines to `covstruct.R`, and the new
  `?mixture` paragraph and the prose add 20 to `families.R`. Guardrail
  3 asks for the size of what moves, but the user decides on what core
  loses.
- The time share comes from one loaded run. The after-arm sum for the
  same 262 files is 1678.9 s, against 2045.6 s before (`after.log`,
  `before.log`). That swing of 367 s is 58 times the 6.3 s that the
  move removes. "About 0.3%" is right as an order of magnitude. Say
  that it is a single loaded measurement.
- `count.txt` uses a `test-core-boundary.R` line that was "replaced"
  in `after.log` after a separate rerun. The rerun is legitimate, but
  the record should name the log of the rerun, not an edited line.
- The core preload session shows skip=204, against 202 forward. One of
  the two is the new guardrail block (finding 7). The record does not
  explain them.

Fix: add these sentences to section 1 and section 5 of the findings.

### 5. major (consolidation): frmtmb.latent's floor is false in the tree, and the release must be in lockstep

frmtmb.latent's `DESCRIPTION` still says `frmtmb (>= 0.72.0)`. The lane
measured that the package does not install over 0.73.0, because
`us_sigma` is not exported there (`floor-install.txt`). The lane did
not change the floor, because lane rules forbid version numbers. In the
other direction, 0.6.2 does not load over the new core (confirmed
above). Then a user who upgrades core alone loses `hmm()` and `lca()`
too.

My judgment on the handling: a lockstep release with a raised floor is
enough, because frmtmb has no existing users. Core should NOT accept
rules that name an absent feature. The refusal at registration is
deliberate: it turns a misspelled rule side into an error and not a
row that is dropped without a word. A rule that names a feature that
core removed looks the same as a misspelling. Do not weaken that for
one transition with no users. The `expects =` route that the lane
names is the right tool for a future release: from then on,
frmtmb.latent can declare `expects =` for any feature that another
package supplies.

Fix (consolidation): raise frmtmb.latent's Depends floor to the core
version that exports `us_sigma()`, and say it in latent's NEWS (the
"Needs" line is already there with no number). Add a frmtmb.latent arm
to the floor script that replaces `dev/rel073-floor.sh`. Bump core
minor and latent minor, as the lane says.

### 6. minor: codemeta.json still lists mclust as a dependency of core

`codemeta.json` lines 292-300 list `mclust` in core's software
requirements, and core's DESCRIPTION no longer has it.

Fix: regenerate `codemeta.json`, or remove the entry at the release.

### 7. nit: the guardrail-1 block skips whenever frmtmb.latent is loaded

`test-core-boundary.R`, "mixture_mvn() is not in core, and no string
spelling reaches it", calls `skip_if(isNamespaceLoaded("frmtmb.latent"))`
before both string refusals. So it skips in the core preload session
(`mvnlatent-sess-sessions/frmtmb--preload.txt` line 898). But the two
refusals do not depend on the extension (strings.txt, "lane ... latent
attached yes"). Only `expect_false(exists("mixture_mvn"))` depends on
it, and that depends on ATTACHMENT, not on loading.

Fix: run the string refusals with no skip. Guard only the `exists()`
line, with `skip_if("package:frmtmb.latent" %in% search())`.

### 8. nit: the tape test checks the value only, and "floored as the random-effect blocks do" is not exact at d = 1

`test-us-sigma.R`, "us_sigma() runs on the tape", compares `tp(th)`
with `f(th)` at one point. A tape that gave wrong derivatives would
pass. My finite-difference check (item above) passes, so this is a gap
in the test, not a defect in the function.

At d = 1, `floored = TRUE` gives `(exp(x) + 1e-150)^2`. The 1 x 1 path
of the random-effect blocks uses `sd2_floored()`,
`exp(2x) + 1e-300` (`R/covstruct.R` 789 and 811). The two agree at
underflow, so the difference has no numerical effect. But the help
page claims "as the random-effect blocks of frmtmb do".

Fix: compare `tp$jacobian(th)` with a central difference, as a ratio to
the Jacobian's scale. Either make the d = 1 floored path use
`sd2_floored()`, or reword the help page.

### 9. nit: the refusal for a non-whole `d` prints a value that hides why

`us_sigma(c(0, 0, 0), 2 + 1e-12)` gives "`d` must be a positive whole
number, not numeric "2"" (cases.txt). The printed value looks whole.

Fix: print `d` with 15 significant digits in this message, or say "not
a whole number: 2.000000000001".

### 10. nit: the fam$mix component interface is named, not specified (pre-existing)

The lane says that the family reaches core "only through generic slots
(`fam$mix`, ...)". That is true. But `?frmtmb-extension-api` names the
"`fam$mix` component interface" and does not say what is in it. Core
reads the fields `K`, `log_pi(dpars)`,
`comp_lpdf(y, dpars, aterms, k[, extra])`, `comp_dpars(dpars, k)`,
`comp_sim`, `ord` and `ref`. A third party cannot write a fourth such
family from the help pages. The gap is not this lane's, because
`lca()` relied on it already. But the move makes the slot an
inter-package contract.

Fix: file it in `dev/test-backlog.md` or the next interface round.
Document the fields under `mixture_posterior()` in
`?frmtmb-extension-api`.

### 11. nit: smaller points

- The latent `test-prior-dpar.R` block that moved there reaches three
  core internals through the existing `pd_reached()`
  (`frmtmb:::prior_design`, `resolve_priorlist`, `par_template_names`).
  The helper was there before, but the extension test now depends on
  one more core internal. Note it in the out-of-tree inventory.
- The core `test-unpinned-seams.R` rewrite drops the softmax check on
  `conditional_effects(dpar = "theta1")`. The check is kept in latent's
  copy. Core's own `mixture()` covers the `dpar_response` seam
  elsewhere. Not a loss, but the core block now asserts only
  `is.finite()` on `dpar = "mu"`.
- The latent vignette says "The family moved here from frmtmb, and its
  fits did not change in the move." That is history, not tutorial
  content (Diataxis), and it will age badly. Its
  `stopifnot(abs(... ) < 1e-3)` is an absolute tolerance. A vignette is
  not a test, so this is a nit.
- NEWS (core and latent) says the fits are "identical to the last bit".
  That was measured on this machine with the reference BLAS. It is
  true, but say "on the same machine".
- `test-mixture-mvn.R` line 465, a `test_that()` title, is 93
  columns. The four long lines in `R/mixture-mvn.R` are the rule
  strings that moved, and their length follows `R/compat.R`.
- Pre-existing, not this lane's: `predict(newdata =)` refuses this
  family with "its draws are structured (a hidden state sequence, a
  group-level latent class, or a correlated residual)"
  (`R/predict-brms.R` 615). None of the three applies to a multivariate
  normal mixture.

## Files the reviewer wrote

`dev/reviews/2026-10-23-mvnlatent.md` (this file),
`dev/review-mvnlatent-run1.R`, `dev/review-mvnlatent-par.sh`,
`dev/review-mvnlatent-cases.R`, `dev/review-mvnlatent-reload.R`,
`dev/review-mvnlatent-compat.R`, `dev/review-mvnlatent-strings.R`,
`dev/review-mvnlatent-oldlatent.R`, `dev/review-mvnlatent-glmmtmb.R`,
`dev/review-mvnlatent-libcheck.R`; logs in `dev/review-mvnlatent-out/`
and `dev/review-mvnlatent-out-tests.txt`.

## Short check of fix round 1

Verdict: **SHIP**, with one nit and one item for the consolidation
(F5, unchanged). Script `dev/review-mvnlatent-fix1.R`; logs in
`dev/review-mvnlatent-out/fix1/` and the per-file
`dev/review-mvnlatent-out/<pkg>__<file>.txt`.

- **Instrument.** The lane library is the worktree source again:
  frmtmb has 1979 source functions (one more, `us_theta`) and 0
  differ; latent has 67 and 0 differ; sample has 285 and 0 differ.
- **Fits.** `dev/mvnlatent-fitcompare.R` on the rebuilt lane library:
  59 fits, and `identical()` with my earlier base save is TRUE.
- **Tests, lane library first.** Every `lib:` line shows
  `wt-mvnlatent-lib`.
  - The three that were asked for: `test-us-sigma.R` 52 pass (8
    blocks), `test-core-boundary.R` 10, latent `test-mixture-mvn.R`
    308.
  - Other files the fix can reach: core `test-compat.R` 771,
    `test-covstruct.R` 7, `test-covstruct-floor.R` 940,
    `test-data2.R` 31, `test-id-kron.R` 61, `test-gr-by.R` 71,
    `test-structure.R` 128, `test-prior-location-dpar.R` 55,
    `test-unpinned-seams.R` 30, sample `test-simulators.R` 54, and all
    12 latent files (679 in all).
  - Fail, error and warn are 0 everywhere.
- **F1 fixed.** `?us_sigma` no longer claims glmmTMB's
  parameterization in general. It says that the map is glmmTMB's up to
  d = 3 and that the order differs from d = 4. My d = 4 check against a
  hand-built `L`, filled column by column, is TRUE.
- **F2 fixed.** `us_theta()` is exported on the same help page, and
  the moved test uses it in place of its copy.
  - Round trips, seed 99, d = 1 to 6, both directions: the largest
    relative error is 5.7e-16. A matrix with rho = 0.999 round-trips
    to 2.2e-16.
  - It refuses, with a message that names `us_theta()`: a non-square
    matrix, a character matrix, a vector, NA, Inf, an asymmetric
    matrix, a zero or negative diagonal, a 0 x 0 matrix, and the zero
    matrix.
  - It accepts: a symmetric matrix at scale 1e-300, an asymmetry of
    1e-14 (rounding), and an integer matrix.
- **F4 fixed.** Recomputed: `git diff --numstat 19b162b9 -- R/` gives
  +182 / -417, net -235, as section 1 now says.
- **F6 fixed.** `codemeta.json` has 0 mentions of mclust and still
  parses.
- **F7 fixed.** The string refusals run with no skip. Only
  `exists()` is skipped, and only when `package:frmtmb.latent` is on
  the search path.
- **F8 fixed.** The tape test compares the AD Jacobian with central
  differences, as a ratio. `us_sigma(-400, 1, floored = TRUE)` is
  identical to `sd2_floored(-400)`.
- **F9 fixed.** The message says "not 2.000000000001". For
  `d = "a"` and `d = c(2, 3)` it falls back to `arg_desc()`.
- **F10 fixed.** `?frmtmb-extension-api` has "The mixture component
  interface", with `K`, `log_pi()` and `comp_lpdf()` (the fifth
  argument only with `extra_pars`). The `fit_extras()` that it names is
  exported. The fields it says only `mixture()` reads match my earlier
  grep of the core reads.
- **F11 fixed.** The test title is shorter, the vignette history
  sentence is gone, the vignette assertion is a ratio (1e-7), and both
  NEWS files say "on the same machine". The two items left as they
  were (`pd_reached()`, and the pre-existing `predict(newdata =)`
  message) are as I asked.

Remaining:

- **nit:** for a matrix with a positive diagonal that is not positive
  definite (`matrix(c(1, 2, 2, 1), 2)`, or the singular `matrix(1, 2,
  2)`), `us_theta()` refuses with the internal message "The requested
  correlation matrix is not positive definite". The message does not
  name `us_theta()`, and it says "correlation" when the user passed a
  covariance. The test pins only "not positive definite". Fix: catch
  the error in `us_theta()` and say "us_theta(): `S` must be positive
  definite".
- **nit:** the edits for "on the same machine" left `NEWS.md` line 13
  at 87 columns and `extensions/frmtmb.latent/NEWS.md` line 18 at 89.
  Reflow both.
- **F5 (consolidation, unchanged).** Raise frmtmb.latent's floor and
  add a latent arm to the floor script.
