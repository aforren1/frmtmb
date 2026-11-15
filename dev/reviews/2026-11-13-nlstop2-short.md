# Short check of lane nlstop2's fix round (round of 2026-11-13, base ec8178d3)

Scope: the fixes for F1 to F6 of `dev/reviews/2026-11-13-nlstop2.md`
only; the lane is not reviewed again. Scripts
`dev/reviews/nlstop2-short-list.R` (one test file, every expectation
listed with its test and source line) and `dev/reviews/nlstop2-short-run.sh`;
logs `dev/reviews/nlstop2-short-out/`. File compares are `git show
ec8178d3:<path>` against the tree with `diff --strip-trailing-cr` (the
tree has CRLF line ends) and `git hash-object`.

Verdict: **FIX, record only.** F1 to F5 are done and hold: the revert
is clean, the explin test asserts the same set on every BLAS, the
backlog is right, and the check reads "Status: 2 NOTEs" (environment).
F6 needs three edits to `dev/nlstop2-findings.md` (section 4). With
those edits made, SHIP. No code change and no rerun needed.

## 1. F5: the huber() route, reverted

- `R/fit.R`: hash-object b68a5bdb..., equal to ec8178d3's. `man/frm.Rd`:
  36acd839..., equal.
- `huber_quad_route` / `quad_route`: no match in `R/`, `tests/`, `man/`,
  `NEWS.md`, `extensions/`, `dev/test-backlog.md`.
- `NEWS.md`: one BREAKING bullet in the development section (the pole
  restart). The huber() item is a `## Documentation` bullet, not BREAKING.
- `R/fit-restart.R` against ec8178d3 (289 lines, now 576), all 12 hunks:
  - base 187 to 194: roxygen of `nl_limit_message()`, the Newton step (limit check).
  - base 248 to 281 (10 hunks, new 249 to 300): `nl_limit_message()` body, `val()`, `f(w, t)` with
    the Newton step across the move, the quadratic gate on `w`, `best`
    (limit check).
  - after base 288, new 308 to 575: new `nl_denoms()`, `nl_pole_walls()`, `nl_pole_pressed()`,
    `nl_pole_escape()`, `nl_pole_message()`, `nl_pole_warning()`,
    `nl_refit_held()` (pole restart). No "huber" in the file.
- `R/families.R` against ec8178d3: one hunk, 9324 to 9346, `huber()`
  `@details`. `man/frmtmb-families.Rd`: the same text, one hunk (244 to
  260).
- `?huber` text against the logs:
  - "100 rows in 20 groups, heavy-tailed noise": `dev/nlstop2-huberkink.R`
    (20 groups of 5, t3 noise).
  - "11 of 20 ... up to 0.188 and 0.272 ... 6 and 5 gave no warning":
    recounted from `dev/nlstop2-out/huberkink-before.txt` (gap > 1e-4):
    scale 1 short 11 silent 6 max 0.188; scale 1000 short 11 silent 5 max
    0.272.
  - "On the same 40 fits it reached the maximum every time":
    `huberquad-before.txt`, gap at most 1.81e-8 and 9.41e-8, code 0 on all
    40. "The maximum" is the best of 8 jittered starts (nit, below).
  - "7 to 9 ... 39 to 49": `cost-huber.txt` 6.82, 38.88; `cost-huber2.txt`
    8.57, 48.63; 2000 rows are 100 groups of 20 (`dev/nlstop2-cost.R`).
  - "sigma ~ (1 | g) ... 6 of 10 data sets; huber(link = "log") ... 4 of
    5": `dev/reviews/nlstop2-review-out/huber2-after.txt`, 6 of 10 ERROR
    "could not marginalize" at scale 1 and 6 of 10 at scale 1e3; log link 4
    of 5.
  - predict(), fitted(), residuals(), frm_linpred(se.fit = TRUE),
    avg_slopes() warn; vcov_cluster() refuses: `huber2-after.txt` rows 43
    to 54 (`frmtmb_modes_conditional_se`, `frmtmb_warning`,
    `frmtmb_interop_se_part`; "vcov_cluster() cannot use a fit made with
    quadrature = TRUE").
  - It says where quadrature fails and recommends it only for a scalar
    term in the mean with the identity link. Plain, no overclaim.
- Nits, not blocking: `?huber` "6 of 10 data sets" is 6 of 10 at each of
  two scales; NEWS "it reached the maximum on every fit of a scalar (1 |
  g) term" reads as a general claim; "on all 40 measured fits" would be
  exact.

F5: done.

## 2. F1: the explin test, and 58 against 59

The test "a ridge that bends is read after a step across it" now asserts,
with no branch, a warning (line 245) and
`expect_false(is.null(nl_limit_message(got$value)))` (line 246).

Runs (`nlstop2-short-run.sh a`, three processes at once):

| log | build | result |
|---|---|---|
| `fit-restart-ref.txt` | lane, reference BLAS | pass 59, fail 0, err 0, skip 0, warn 0 |
| `fit-restart-ob0332.txt` | lane, OpenBLAS 0.3.32 + LAPACK, `OPENBLAS_CORETYPE=Haswell` (`dev/ciharden-out/Rob0.3.32-lapack`, the lane's copy) | pass 58, fail 0, err 0, skip 0, warn 0 |
| `fit-restart-before.txt` | rellib-r21, reference BLAS | pass 35, fail 17, err 1 |

Expectation lists compared by (class, source line, count): the only
difference is ONE `expectation_success` at line 203, present on reference
BLAS, absent under 0.3.32. Lines 245 and 246 (explin) are in both lists.

Line 203 is inside `if (has_class(got$conds, "frmtmb_nl_limit"))` in the
test "a body on a ridge to its linear limit says so" (bernoulli mm 55).
That branch is at ec8178d3 lines 195 to 200, unchanged by the lane. Its
comment names the two endings (at the linear limit with the limit
warning, or 0.0065 below with the standard-error warning; 17 and 3 of 20
ulp draws; OpenBLAS takes the second). Both endings get the same
unconditional assertions (lines 197, 198: a warning, and
`nl_limit_message()` non-NULL); the branch only adds the value check that
exists on the at-limit ending. That is a conditional that pins one of two
valid endings with a justification, so it passes the rule. The lane's
first-pass 69 against 71 was this branch plus the explin code-0 branch the
fix removed, not the explin branch alone as the findings say.

F1: done.

## 3. F2 to F4: backlog

`dev/test-backlog.md`, "Filed in the round of 2026-11-13" / "From lane
nlstop2" (line 4797, directly above "## Filed at the 0.82.0 release"):
- F2: "LOW. Cumulative `a * xp^b` seed 22 ... 1.03 short", -408.3803
  against -407.3511 (`accsum-after2.txt` line 59: gap 1.029229). The 5079
  CLOSED line's "filed again" is now true.
- F3: 5093 CLOSED line (5174): "the Hill fit has no limit (h to infinity
  is 0.147 lower), its 1.8e-3 short stop is filed"; item "LOW. The Hill
  body on step data stops 1.8e-3 short" filed. Findings row and 5093
  section corrected.
- F4: "HIGH. A fit clear of a band of poles is not looked at inside it,
  and stops 591 to 646 short in silence", with the review's numbers and
  script.

F2 to F4: done.

## 4. F6: findings record

`dev/nlstop2-findings.md`: no "IN PROGRESS", "PENDING" or "rerun the
final build" left; Status DONE; the fix round has its own section;
the 5093 row and section carry the F3 text. Three gaps, all record only:
- (a) Line 539: "58 against 59 on reference BLAS, a difference in the
  number of assertions run that is not explained here". Replace with:
  "58 against 59: the value check of the mm 55 test's at-limit ending
  (line 203, in a branch that is unchanged since 0.82.0), which OpenBLAS
  does not take; the explin assertions run on every BLAS
  (dev/reviews/2026-11-13-nlstop2-short.md)."
- (b) Line 423: "(69 against 71 on reference BLAS: the explin test's
  code-0 branch.)" accounts for one of the two. Should read "the explin
  test's code-0 branch (removed in the fix round) and the mm 55 test's
  at-limit branch".
- (c) Lines 561 to 564: the check is "still running at hand-back". It has
  now ended (section 5); add its Status line. The deletes it asks for
  (`dev/nlstop2-check/`, `dev/ciharden-out/`) are not done; this check did
  not do them.

F6: FIX (the three record edits above; no code change, no rerun).

## 5. The two runs open at hand-back

Neither was touched; this check waited on the process (PID 32380,
`dev/nlstop2-check.sh`), polling every 120 s.

- 0.3.32 explin draws (`dev/nlstop2-out/ob2-0.3.32-explintest.txt`) had
  already ended when this check began ("DONE"; `ob2-driver.log` "OB2
  DONE"). blas: line `dev/ciharden-out/Rob0.3.32-lapack`, lib: line the
  lane library. "explin test body: pass 12 of 12": k 7 at code 0 with
  `frmtmb_nl_limit`, the other 11 at code 1 with `frmtmb_conv_code`
  (6 with `frmtmb_conv_gradient` as well); `limit_message TRUE` on all 12.
- Core `R CMD check --as-cran` on the fix-round tree
  (`dev/nlstop2-out/check-driver2.txt`, `dev/nlstop2-check/check.log`),
  ended 16:16: "Status: 2 NOTEs". The NOTEs: example timing
  (frm_hazard_reads 10.20 s, ps 6.25 s, profile.frmtmb_fit 7.08 s,
  dharma_residuals 6.53 s elapsed; each 3.3 s CPU or less, so the
  machine's load) and "Skipping checking math rendering: package 'V8'
  unavailable". Both come from the environment, as in the first pass.
  Tests "[66m] OK"; vignettes rebuilt OK; PDF manual OK.

## 6. Reruns

Lane library, reference BLAS (`nlstop2-short-run.sh b`):
- `test-fit-restart.R`: pass 59, fail 0, err 0, skip 0, warn 0.
- `test-bracket-access.R`: pass 36, fail 0, err 0, skip 0, warn 0.
- `test-families.R`: pass 252, fail 0, err 0, skip 0, warn 0.

rellib-r21, `test-fit-restart.R`: pass 35, fail 17, err 1, the lane's
count. Failures are behavioral: explin no warning (245) and
`nl_limit_message()` NULL (246); binomial mm 11 and 44 without
`frmtmb_nl_pole_start` (264), 4.322 and 0.919 short (265 and 4 ulp draws
each at 272), K -0.2 and -0.1 (266); no `frmtmb_nl_pole` from K -2.5
(287). The error is the absent case's missing `nl_pole_walls()` (300),
the weak form the review accepted.
