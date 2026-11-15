# Short check of the refits fix round (round of 2026-11-13, base ec8178d3)

Scope: the fixes of the review `dev/reviews/2026-11-13-refits.md` (F1,
F2, F3, 1259 digit text, 3188 refresh, backlog filings). The lane was
not reviewed again. Lane library `wt-refits-lib` (frmtmb built
2026-11-13 19:39 UTC, after the last edit of `R/`, 19:36 UTC; each
probe prints that the installed functions hold the fix). brms 2.23.0,
`standata()` only. Scripts and logs: `dev/reviews/refits-short-*`.

## Verdict: FIX

Two small items. F1, F3, 1259, the backlog filings and the tests hold.

- **N1 (regression of the fix round, false refusal).** The 3188
  refresh writes the new response into the frame's `data_frame`, and
  `refit_data_matches()` compares every column with the call's data.
  So `vcov_cluster(refit(fm, Y2), ~ cl)` on a matrix-response fit is
  now refused ("those no longer hold the data the model was fitted
  to"), where rellib-r21 gives the right answer, equal to the explicit
  factor (`refits-short-extra-{before,after}.log`, M1, seed 5, n = 200,
  `multinomial()`):

  ```
  rellib-r21  vcov_cluster(rf2, ~ cl):   0.0523188 0.0514621 0.0493927 0.0617483
  lane        vcov_cluster(rf2, ~ cl):   ERROR ... no longer hold the data ...
  both        vcov_cluster(rf2, dd$cl):  0.0523188 0.0514621 0.0493927 0.0617483
  ```

  The same call on a univariate refit computes (its `data_frame` is not
  rewritten; `refits-short-3188.log` R14). Fix: leave the response
  columns out of the value comparison in `refit_data_matches()` (the
  rows and the other columns still prove the alignment), or compare a
  refit against the data before the refresh.

- **N2 (F2 incomplete, silent, pre-existing).** When the call's data
  are gone (fit made in a function on local data), `resolve_cluster()`
  still falls back to the formula environment, so a workspace `cl` of
  length n is used in place of the data's column. The F2 rule is "the
  formula's environment only when the data lack the variable"; gone
  data are not data that lack it, and `sandwich::vcovCL()`
  (`expand.model.frame()`) stops there. `refits-short-extra-after.log`
  X1 (probe2 data, seed 1): `vcov_cluster(fl, ~ cl)` = 0.0162685
  0.0645226 0.0215026 (the global `cl`), where the data's `cl` gives
  0.3053286 0.1815976 0.1108744. Identical on rellib-r21. Without the
  global the call is refused (probe2 V11), so the refusal path exists.
  Fix: in `R/sandwich.R` the last fallback
  `if (is.null(v2) && is.null(v)) v <- eval(ex, held, env)` runs only
  when `d0` was read or the call names no data
  (`!is.null(d0) || is.null(fit$call$data)`); else the V11 refusal.

## Item by item

### 1. F1 `refit_bernoulli_codes()`, `boot_refits()`: OK

`refits-short-f1.log` (seed 9, n = 40). `same()` is `identical()` of
`fixef()` against `refit(f01, ...)` on the 0/1 fit.

```
A  {0,1} fit: rep(1, n) refits (intercept sign +1); rep(0, n) -1;
   s, as.character(s), s == 1 all = fit; a 2 refused (neither)
B  {0,2} fit: s (codes) = refit(f01, s); 2s (values) = same;
   rep(1, n) = refit(f01, 1s); rep(2, n) = same; rep(0, n) = refit(f01, 0s);
   0/1/2 mix refused; residuals(newdata y = 1) refuses 1 (values only):
   no case where refit() and residuals() read one number differently
C  factor fit, levels "1","2": numeric s + 1 (values) and s (codes) = f01;
   rep(1, n) refused, names both readings and the two ways to say which;
   factor("1") = failure; rep(TRUE, n) = success; factor(s, 0:1) refused;
   factor with levels reversed = values by label; numeric fit {1,2}:
   factor("1") failure, TRUE success, 2s success, 0s failure (0 only a code)
   logical fit: rep(1, n) success, FALSE failure
D  NA in numeric, factor and the constant 1 vector: refused before the
   coding ("newresp has 1 NA value(s)")
E  {-1,0} fit: rep(0, n) refused (two ways); rep(-1, n) failure; s - 1 = f01
```

Bootstrap (seed 21, n = 12, nsim 40, seed 5), all-one-class draws on
each side: near-all-success data give 13 all-1 draws, near-all-failure
data 13 all-0 draws. Against the 0/1 coding, `frm_bootstrap()$t` is
`identical()` for 1/2, 0/2, -1/0, factor no/yes and logical, with the
same 15 of 40 non-converged in every coding, and no draw refused (0 NA
rows). `boot_refits()` maps draws with `as.logical()`, so they never
reach the guess.

`update(newdata = )` is a fresh read of the new data, not the refit
rule: on the 1/2 fit, `y` all 1 reads as a success (intercept +1), all
2 success, all 0 failure; a factor with levels 1, 2 all "1" failure.
brms 2.23.0 does the same: `update.brmsfit()` sets `data <- newdata`
and rebuilds the basis (`R/update.R:91-92, 271`), and a fresh
`standata()` of all 1s gives `Y = 1` (G6b). Non-constant vectors agree
with rules 1 to 3 (G1, G2: = f01). Not a defect.

`influence()` (seed 31, n = 30, one success at row 28): deleting it
leaves only the value 1 on the 1/2 coding; the 1/2 and 0/1 fixed tables
are `identical()` (row 28: -627.9464, -0.005399515 on both). The
reviewer's script prints `-Inf` here because it reads a column name
that does not exist (`fixed.effects[-case]`); the check above replaces
it.

Reruns, diffed against the lane's fix-round logs: `refits-review-bern.R`
(`refits-short-bern.log`), `refits-review-probe5.R`
(`refits-short-probe5.log`) and `dev/refits-f1-readers.R`
(`refits-short-readers.log`): all three identical.

### 2. F2 `resolve_cluster()`: OK on the listed rows; N2 open

probe2 rerun (`refits-short-probe2.log`): V10 = 0.3053286 0.1815976
0.1108744 = V0. V1, V2, V5, V6, V7, V8 equal the explicit factor; V3,
V4, V11, V12 refused; V9 silent as documented. Everything after the V
rows is identical to `dev/refits-f2-probe2-after.log`.

`refits-short-f2.log` (seed 1, 120 rows):

```
S1 cl2 only in the workspace, data lack it   = explicit cl2
S2 cl3 local to a calling function           = explicit
S3 workspace cl2 of length 120, 2 NA rows    refused (length 120 vs 118)
S4 workspace vector of length 118            = explicit
S5 cl a fixed term                           = explicit dd$cl
S6 cl a grouping factor (1 | cl)             = explicit dd$cl
S7, S8 the same with a global cl of length n = S5, S6 (data column wins)
S9 cl a fixed term, NA rows                  = explicit dd$cl[-c(3, 50)]
S11 global cl beside data cl, NA rows        = explicit (data column)
```

N1 and N2 above.

### 3. F3 NEWS: OK

New section, BREAKING on: character ordinal (digit text), `y ~ . - w`
/ `formula(fit)`, `refit()` bernoulli reading, `update()` arguments,
`getCall()` after a delta, `vcov_cluster(~ v)`, `influence()` with
`s()` / `t2()`: the six of F3 plus `update()`. No lane names, no
`dev/` paths. Version numbers appear only as facts: `brms 2.23.0`
(twice) and "saved by frmtmb 0.76.0 or earlier", as the review
accepted.

### 4. 1259: OK

`refits-short-1259.log` (seed 2501, n = 60). For `cumulative`,
`sratio`, `cratio`, `acat` and `hurdle_cumulative`, labels and digit
text: the frmtmb message starts with brms 2.23.0's sentence, equal
character for character ("Family 'X' requires either positive integers
or ordered factors as responses."; "non-negative" for
`hurdle_cumulative`), then names the response, its values and the fix.
Ordered factor and integer responses fit for all four. A character
response with NA rows, a single non-NA value, and the ordinal response
of an `mvbf()` (`refits-short-3188.log`) are refused the same way;
`categorical()` character still fits. Docs: `?influence.frmtmb_fit`
now says "integers or an ordered factor" (base: "character vector");
`test-thres-refit.R` drops the character coding; no other text in `R/`,
`man/`, `vignettes/` or the extensions claims a character ordinal
response fits.

### 5. 3188 refresh after `refit()`: OK for `Y`; N1; one LOW

`refits-short-3188.log` (seed 5, n = 200): after `refit(fm, Y2)` the
attribute's `Y.a`..`Y.c` and `data_frame$Y` equal `Y2`, `x` unchanged,
the original fit's attribute unchanged; `avg_predictions(rf)` = 
`colMeans(fitted(rf))` (3.434996 4.930001 3.635003).

LOW: with the `cbind(Ya, Yb, Yc)` spelling the attribute carries the
split matrix and the raw columns `Ya`, `Yb`, `Yc`; after `refit()` the
matrix columns are new and `Ya`..`Yc` old (R9 FALSE). No estimate
reads them (predictions ignore the response), so not blocking; file it.

Pre-existing, found on the way (identical on rellib-r21,
`refits-short-extra-*.log` U1/U2, seed 8, n = 30): `influence()` of a
univariate refit deletes from the ORIGINAL response, because its
`data_frame` is not rewritten: max |diff| against `influence()` of the
original fit 1.06e-05, against a fresh fit on the new response 1.85.
Silent wrong; outside this round's items. The refresh now covers
matrix responses only, so univariate and matrix refits differ here
too. Writing the response for every refit fixes it, and then N1's
comparison must leave the response out for univariate refits as well.
File it (MEDIUM or HIGH).

### 6. Backlog filings: OK

"Filed in the round of 2026-11-13" / "From lane refits", directly
above "## Filed at the 0.82.0 release": the `gp()` deletion refit item
(MEDIUM, dmax 2.928603 -> 2.919538) and the LOW `y ~ x + w - w` item
extended with the first fit of `y ~ . - w` / `.^2 - w` and `y ~ . - q`.

### 7. Tests on the lane library

One file per process, 3 at once, `NOT_CRAN=true` (brms fit gate not
set), every `lib:` the lane library (`refits-short-test-*.log`, driver
`refits-short-runtest.R`, list `refits-short-tests.jobs`). The extra
files are those a multiline grep finds with a character vector near an
ordinal family (core and extensions); none fits a character ordinal
response.

```
test-refit-plumbing.R      pass=57   fail=0 error=0 skip=0
test-thres-refit.R         pass=110  fail=0 error=0 skip=0
test-bernoulli-coding.R    pass=48   fail=0 error=0 skip=0
test-double-bar.R          pass=155  fail=0 error=0 skip=0
test-formula-dot.R         pass=27   fail=0 error=0 skip=0
latent test-hmm-starts.R   pass=95   fail=0 error=0 skip=0
test-thres.R               pass=67   fail=0 error=0 skip=0
test-pp-check-types.R      pass=216  fail=0 error=0 skip=0
test-ordinal-disc-thres.R  pass=179  fail=0 error=0 skip=0
test-ordinal-mixture.R     pass=99   fail=0 error=0 skip=0
test-nl-family.R           pass=333  fail=0 error=0 skip=0
test-diagnostics-ux.R      pass=130  fail=0 error=0 skip=0
test-cs-factor.R           pass=84   fail=0 error=0 skip=0
test-ce-bands.R            pass=170  fail=0 error=0 skip=0
test-ce-parity.R           pass=30   fail=0 error=0 skip=0
test-brms-families.R       pass=1257 fail=0 error=0 skip=0
test-brms-likelihood.R     pass=5    fail=0 error=0 skip=56 (brms gate)
sample test-simulators.R   pass=55   fail=0 error=0 skip=0
sample test-thres-draws.R  pass=4    fail=0 error=0 skip=0
sample test-sampling-ported.R pass=207 fail=0 error=0 skip=1
```

Warnings 0 in every file.

## To do for SHIP

1. N1: `refit_data_matches()` leaves the response out of the value
   comparison; test `vcov_cluster(refit(<multinomial fit>), ~ cl)` =
   the explicit factor.
2. N2: no fallback to the formula environment when the call names data
   that cannot be read; test the X1 case refused.
3. File: the univariate `influence(refit())` stale response (silent,
   pre-existing) and the `cbind()` attribute LOW.

## Re-check after fix round 2

Scope: "Fix round 2" of `dev/refits-findings.md`. Lane library frmtmb
built 2026-11-13 21:53 UTC, after the last edit of `R/` (21:52 UTC).
At most 2 R processes. Scripts `refits-short-r2.R` (`ARM=before` runs
only part 4b on rellib-r21), `refits-short-r2b.R`, the earlier
`refits-short-extra.R`; logs `refits-short-r2-*.log`,
`refits-short-r2b.log`.

### Verdict: SHIP

N1, N2 and the `influence()` of a refit hold on every case below. No
case is silent-wrong; no non-refit table changes.

### 1. M1, X1, U1/U2 (`refits-short-r2-extra.log`)

```
M1 vcov_cluster(refit(multinomial), ~ cl)  0.0523188 0.0514621 0.0493927 0.0617483
   = explicit dd$cl (M2), = rellib-r21     (was: refused)
X1 data gone, global cl of length n        refused "(d_local, which are gone)"
                                           (was: the global cl, 0.0162685 ...)
U1 influence(refit) vs fresh fit on ynew   3.33e-07  (was 1.850286)
U2 influence(refit) vs original fit        1.850289  (was 1.06e-05)
```

### 2. `ignore_response` accepts no other change

Only `resolve_cluster()` sets it (`R/sandwich.R:92`); `update()`
(`R/confint.R:3139`) keeps the full comparison. Seed 1, 120 rows:

```
N1a refit, data unchanged: ~ cl                  = explicit
N1b refit, dd$x changed                          refused
N1c fit, dd$x changed                            refused
N1d update(fit, . ~ . + z), dd$x changed         refused (New variables found)
N1e refit, two rows swapped, row names reset     refused
N1h fit, one row deleted                         refused
N1i log(y + 10) fit, raw dd$y changed            refused (only "log(y + 10)"
                                                 is skipped; raw y compared)
N1j matrix refit, data unchanged                 = explicit
N1k matrix refit, dm$x changed                   refused
N1l matrix refit, trials column dm$N changed     refused
N1f fit, only dd$y changed                       computes, = explicit (rows and
                                                 cl intact, so the map holds)
N1g update(fit, . ~ . + z), only dd$y changed    refused
```

### 3. N2 refuses only gone data

```
N2a top-level fit, data present: ~ cl                = explicit
N2b data present, variable only in workspace (cl_ws) = explicit
N2c' data = globalenv(), variables global: ~ clv     = explicit (0.3053286 ...)
N2e data passed by value through do.call()           = explicit
N2g fit with data2 (reads call[["data"]], not data2) = explicit
N2h ~ cl built inside a function, data global        = explicit
N2f fit in a function on its argument, data gone     refused, names `d`
```

`data = NULL` cannot be tested: `frm()` refuses it ("`data` is NULL"),
and without `data` it stops with "argument "data" is missing". The
no-data branch (`is.null(fit$call[["data"]])`) is reachable only by a
call edited by hand; `data = globalenv()` is the case a user has, and
it works.

### 4. `influence()` of refits; non-refit fits unchanged

4a, seed 8, n = 30, rows shuffled. Max |diff| of the fixed table of
`influence(refit(f, new))` against `influence()` of a fresh fit on the
new response, and against the original fit's table:

```
                                   vs fresh     vs original
R1  gaussian                       2.78e-08     2.27
R2  weights(w)                     1.56e-05     2.47
R3  binomial trials(N)             2.68e-07     1.33
R4  bernoulli 1/2 (values)         1.65e-08     0.945
R5  (1 | g), groups = "g"          5.84e-10     2.29
R6  ar(time = t, gr = g), unsorted 9.38e-08     2.31
R7  log(y + 10)                    4.70e-08     0.209
R8  multinomial matrix response    7.52e-06     1.80
R9  subset(keep) (24 of 30 rows)   3.14e-08     2.16
R10 na.action dropped 2 rows       3.25e-08     2.35
```

Every refit table is the fresh fit's to optimizer tolerance and none is
the original's. R5 theta differs by 0.57 on the log scale: both at the
boundary (sd 2.4e-5 vs 4.2e-5, max |diff| of sd 1.8e-05), and the full
fits already differ the same way (-10.65361 vs -10.08172,
`refits-short-r2b.log`); not the change. Multivariate: `refit()`
refuses an `mvbf()` fit ("refit() supports univariate models"), so no
multivariate refit exists; `influence()` of the `mvbf()` fit runs
(30 x 6) on the old path (`y_fit` is NULL for more than one response).
`influence(refit, data = du)` takes the given data, response included
(= the original fit's table, 1.3e-05).

4b, non-refit fits, rellib-r21 against the lane, seed 41, n = 60, rows
shuffled: `identical()` TRUE for all eight tables, fixed and theta:
`(1 | g)` by group, bernoulli 1/2, `weights()`, `subset()`, binomial
`trials()`, `ar()` by group, `log(y + 10)` and `mi()` with NA rows.

### 5. NEWS: accurate

The `vcov_cluster()` bullet adds gone data (refused) and the refit
response (not compared): matches X1, N2f, M1, N1j. The new bullet
"`influence()` of a `refit()` deletes rows of the refit's response"
matches R1 to R10 and U1/U2. LOW: `?influence.frmtmb_fit` does not say
that the default data of a refit carry the refit's response and that
`data =` overrides it (R12).

### 6. Tests on the lane library (2 at a time, warn 0)

```
test-refit-plumbing.R    pass=60  fail=0 error=0 skip=0
test-sandwich.R          pass=58  fail=0 error=0 skip=0
test-influence-plot.R    pass=21  fail=0 error=0 skip=0
test-refit-data.R        pass=122 fail=0 error=0 skip=0
```

The other files that call `influence()` (grep): test-thres-refit 110/0,
test-data2 29/0, test-mo-coding 166/0, test-se-reparam-inv 13/0,
test-v15 41/0, test-methods-audit 58/0, eam test-ndt-seam 83/0, eam
test-ndt-bound 79/0, learn test-rlddm-ndt 57/0 (skip 0, error 0 each).
