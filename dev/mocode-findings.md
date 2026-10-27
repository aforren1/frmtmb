# Lane mocode: `mo()` codes its integer as brms does

Base: 70476dc3 (frmtmb 0.75.0), branch `wt-mocode`, library
`C:/Users/adf44/source/r/wt-mocode-lib` in front of `rellib-r13`.
Started 2026-10-27. Nothing is committed.

## The defect, reproduced on rellib-r13

`dev/mocode-repro.R` rebuilds the brmsnames review's case (seed 5,
n = 300, `k` in 1..4, `k0 = k - 1`, `k2 = k + 1`, and `kn = k - 3`
added here). I read the review's `dev/review-brmsnames-mo.R` before its
worktree was removed and kept its data, its model and its sampler call
(2 chains, 1000 iterations, seed 11). Logs:
`dev/mocode-out/repro-base.txt`, `dev/mocode-out/repro-lane.txt`.

Base (0.75.0), every number equal to the review's table
(`dev/reviews/2026-10-24-brmsnames.md`, finding 7):

| term | range | logLik (df) | Intercept (SE) | slope (SE) |
|---|---|---|---|---|
| `mo(k0)` | 0..3 | -427.181869 (5) | 2.083383 (0.1137954) | 1.515906 (0.05474984) |
| `mo(k)` | 1..4 | -427.181869 (6) | 1.852026 (NaN) | 1.194768 (NaN) |
| `mo(k2)` | 2..5 | -427.181869 (7) | 1.961154 (NaN) | 0.9339891 (NaN) |
| `mo(kn)` | -2..1 | refused: "variable must be an ordered factor or non-negative integers" | | |

The `mo(k)` and `mo(k2)` fits warn "3 of 4 predictions move along a
direction the fit does not determine"; the fitted values at the four
categories agree with `mo(k0)` to 1e-6 and three of the four standard
errors are NaN. `fitted(mo(k))` at `k = 0`, below the fitted minimum,
returns 1.8520261 without a message. `frm_sample()` of `mo(k)`:
Intercept 0.912 (sd 1.330), `simo_mok1[1..4]` means 0.176, 0.466,
0.238, 0.119, two ESS warnings; of `mo(k0)`: Intercept 2.098 (0.117),
`simo_mok01[1..3]` 0.567, 0.285, 0.148. The review's numbers are the
same (0.912, 1.330, 0.176).

Lane: all four fits have logLik -427.181869 with df 5, Intercept
2.083383 (0.1137954) and slope 1.515906 (0.05474984), no warning; the
four columns of fitted values and of their standard errors are equal;
`k = 0` and `k = 5` are refused ("mo(k): new data outside the fitted
range 1..4"). `frm_sample()` of `mo(k)` gives what `mo(k0)` gives:
Intercept 2.098 (0.117), `simo_mok1[1..3]` 0.567, 0.285, 0.148, no
warning.

## brms 2.23.0's coding, measured

`dev/mocode-brms-coding.R` (`standata()`, `stancode()`, `brm(empty =
TRUE)`; nothing compiles), log `dev/mocode-out/brms-coding.txt`:

| input | Jmo | Xmo_1 |
|---|---|---|
| 0..3, 1..4, 2..5, -2..1 | 3 | 0 1 2 3 |
| gaps 1, 2, 5 | 4 | 0 1 4 (categories 2 and 3 empty) |
| gaps 0, 3 | 3 | 0 3 |
| two values (1, 2 or 0, 1) | 1 | 0 1 |
| one value 3 | 0 | 0 |
| 1.0..3.0 (double) | 2 | 0 1 2 |
| 1, 1.5, 2 | error "Monotonic predictors must be integers or ordered factors" |
| ordered a < b < c < d | 3 | 0 1 2 3 |
| ordered, an unused level (first, last or middle) | levels dropped first: unused last gives 2 for a, b, c; unused first gives 2 for b, c, d |
| unordered factor | error, as above |
| `mo(x) * z` on 1..4 | 3, 3 | 0 1 2 3 |
| `mo(x + 1)` on 1..4 | 3 | 0 1 2 3 |
| `(mo(x) \| g)` on 1..4 | 3 | 0 1 2 3 |
| 1..4 with `attr(x, "min") <- 0` | 4 | 1 2 3 4 |

New data (`standata(fit, newdata =)` on a 1..4 fit): `2..4` gives
`Xmo_1 = 1 2 3` with `Jmo = 3` (coded from the FITTED minimum, kept by
`validate_newdata()` in `attr(newdata$x, "min") <- min(mf$x)`, brms
`R/data-helpers.R` 541 to 557); `3` gives 2; `0` and `5` stop with
"Invalid values in variable 'x'"; `2.5` stops the same way. brms does
not store the minimum in the fit: `attr(fit$data$x, "min")` is NULL,
and `validate_newdata()` recomputes it from `fit$data` each time.
`fit$data$x` of an ordered factor with an unused level has the used
levels only, and a new-data value at the dropped level stops ("New
factor levels are not allowed").

frmtmb matches every row of the table except two (below, "Where
frmtmb still differs"). It already dropped unused ordered levels
(`dev/mocode-out/edge-base.txt`: an unused first, middle or last level
gives D = 2 or 3 as brms). It stores the minimum in the frame's `mo()`
entry (`min`), because the frame, not the data, is what every consumer
reads; the effect is brms's.

### brms is wrong in one corner (filed, not followed)

The undocumented `attr(x, "min")` that brms uses to carry the minimum
to new data is also read at the fit: 1..4 with `attr(x, "min") <- 0`
fits `Jmo = 4`, `Xmo_1 = 1 2 3 4`, but new data 1..4 are then coded
`0 1 2 3`, because `validate_newdata()` overwrites the attribute with
the minimum of the fitted values. The fit and its predictions use
different codes. frmtmb does not read the attribute. Filed as brms-40
in `dev/upstream-bugs.md` (not reported externally).

## The fix

- `R/frame.R`, the `mo()` loop of frame assembly: an integer variable
  is coded `v - min(v)`, `D = max(v) - min(v)`, and the frame's entry
  keeps `min`. Any finite whole number is taken, so negative integers
  fit; a non-numeric, non-finite or fractional value is refused ("must
  be an ordered factor or integers"). An ordered factor is coded as
  before, with `min = 0`. The refusal of fewer than three categories
  now names the term and the count ("mo() needs at least 3 ordered
  categories; mo(x) has 2").
- `R/predict.R`: `mo_min()` returns the minimum, and refuses a frame an
  older frmtmb built when its integer codes do not start at 0 (below).
  `mo_codes()`, which codes every `newdata` value for `predict()`,
  `fitted()`, `conditional_effects()`, `simulate()`, emmeans and
  marginaleffects routes and frmtmb.sample's `posterior_*()`,
  subtracts it and refuses a value outside the fitted range with that
  range in the message (`mo(k): new data outside the fitted range
  1..4`), as brms's `validate_newdata()` refuses it.
- `R/objective.R`: `mo_transform()` and `lp_eta_fixed()` call
  `mo_min()` for its refusal, so a 0.75.0 frame cannot be retaped
  (the objective of `frm_sample()`, `update()` of the objective and
  every numeric reader of the term).
- `R/fit.R` (and `man/frm.Rd`): the "Monotonic effects" section says how
  the categories are coded, with the gap and new-data rules.

### Every consumer of the codes and of D

Grep of `codes`, `$D`, `[["D"]]`, `mo_col_values`, `mo_codes`,
`mo_transform`, `lp[["mo"]]` over `R/` and `extensions/*/R`. None but
the three above needed a change: they read the codes and D from the
frame entry, which is now brms's.

| consumer | what it reads |
|---|---|
| `R/frame.R` `mo_indicator_basis()`, `check_cs_identified()` | codes, D |
| `R/objective.R` `lp_eta_fixed()`, `mo_transform()`, `mo_re_zterm()`, `mo_frame_terms()`, `simplex_frame_terms()`, `mo_chart_frame()` | codes, D, zeta |
| `R/predict.R` `mo_col_values()`, `lp_Z()`, `patch_mo_cols()`, `mo_codes()`, the newdata design paths (lines near 840 and 920) | codes, D, min |
| `R/autoscale.R` (unit of the column) | codes |
| `R/fit.R` `mo_search_terms()`, `mo_search()` (stereographic chart) | D |
| `R/priors.R` simo Dirichlet length and densities (near 2070, 3974, 4015 to 4075) | D |
| `R/se-check.R` `se_mo_gaps()`, `se_mo_null()`, labels | codes, D |
| `R/se-reparam.R`, `R/interop.R`, `R/brms-names.R`, `R/methods-fit.R`, `R/conditional-effects.R` | presence, column, label, D |
| `extensions/frmtmb.sample/R/sample.R` (`mo_sample_*`, `draws_mo_cols()`, `mo_simplex_nlp()`, check_laplace) | D, zeta, through core |

The `mo(x) * f` and `(mo(x) | g)` paths read the same entry (one per
column of the special-term model matrix), so the brms 2.22.0 behavior
kept on purpose (brms#1828) is unchanged. Measured: on the seed-5 data
with `z` and an 8-level `g`, `mo(%s) * z` and `mo(%s) + (mo(%s) | g)`
fitted on k0, k, k2 and kn give bitwise identical logLik
(-431.02567308214731 and -431.66083892540126) and coefficients, all
standard errors finite (`dev/mocode-out/ident.R`).

### Saved fits

A frame that 0.75.0 built has no `min`. It coded the integers
themselves, which is the same model when the data start at 0 and
another model otherwise. `mo_min()` takes 0 for an old frame of an
ordered factor or of integers whose codes start at 0, and refuses the
rest. `dev/mocode-saved.R` saves two fits with rellib-r13 and reads
them in the lane (`dev/mocode-out/saved-lane.txt`): the 1..4 fit is
refused by `fitted()`, `fitted(newdata)`, `predict(newdata)`,
`simulate()`, `residuals()`, `conditional_effects()` and
`frm_sample()`; `update()` refits it (Intercept 2.083383, slope
1.515906, brms's model). The 0..3 fit works in all of them, with the
numbers it had. `logLik()` and `fixef()` of the old 1..4 fit return
what the fit stored (-427.1819; 1.852026, 1.194768 with NaN standard
errors): they read stored values and compute nothing, so I left them.

## Agreement with brms on the sampler

`dev/mocode-agree.R`, the gated test in frmtmb.sample's
`test-mo-coding.R`: `y ~ mo(k)` on the seed-5 data, frm_sample() and
brms's own Stan program (compiled into `dev/stan-cache`), 2 chains of
2000 on seed 3 each. z is the difference of the means over the pooled
Monte Carlo standard error (`posterior::mcse_mean()`).

| par | frmtmb lane mean (sd) | brms mean (sd) | z | sd ratio | base z | base sd ratio |
|---|---|---|---|---|---|---|
| b_Intercept | 2.0914 (0.1140) | 2.0905 (0.1197) | 0.19 | 0.952 | -12.28 | 9.203 |
| bsp_mok | 1.5101 (0.0564) | 1.5116 (0.0561) | -0.65 | 1.005 | -6.57 | 4.935 |
| sigma | 1.0175 (0.0415) | 1.0146 (0.0423) | 1.97 | 0.981 | 1.96 | 1.054 |
| simo_mok1[1] | 0.5670 (0.0326) | 0.5673 (0.0331) | -0.29 | 0.985 | -43.16 | 3.856 |
| simo_mok1[2] | 0.2878 (0.0380) | 0.2863 (0.0376) | 1.13 | 1.010 | 35.66 | 2.052 |
| simo_mok1[3] | 0.1452 (0.0352) | 0.1464 (0.0331) | -0.95 | 1.063 | 31.97 | 1.533 |

brms's largest R-hat 1.0036. The base arm has a fourth weight
`simo_mok1[4]` (not in the table) and one bulk-ESS warning. The test
asserts |z| < 4 and a sd ratio in (0.8, 1.25), as
`test-default-priors-fit.R` does.

## Tests

New: `tests/testthat/test-mo-coding.R` (7 blocks, 71 expectations) and
`extensions/frmtmb.sample/tests/testthat/test-mo-coding.R` (2 blocks,
the second gated by FRMTMB_BRMS_FIT_TESTS). Each was seen to fail on
rellib-r13 first:

- core, base: pass 6, fail 16, error 6, warn 9
  (`dev/mocode-out/t-mocoding-base.txt`). Every block fails: the codes
  and D (1..4 coded 1..4, D = 4), the refusal message, the brms
  standata comparison (Jmo 4 against 3), the shifted fits (negative
  integers refused; NaN SEs), new data (`k = 0` predicted; "fitted range
  0..4"), the simo prior ("Expected input of length 4 ... got 3") and
  the saved-frame refusal (no refusal). Lane: pass 71, fail 0.
- frmtmb.sample, base, gated: pass 16, fail 12, error 1, warn 1
  (`dev/mocode-out/ts-mocoding-base-gated.txt`): 4 simo names, `k = 0`
  predicted, and both the z and the sd-ratio checks fail for the
  Intercept, the slope and the three weights shared with brms
  (10 of the 12 failures; sigma passes, z 1.96, ratio 1.054).
  Lane: pass 29, fail 0, warn 0.

The shift test asserts bitwise identity (`expect_identical()`). That is
an identity, not a tolerance: the four variables have the same codes,
so the fits run the same tape from the same start. It held in one
process on this machine and under OpenBLAS 0.3.26 with LAPACK (below).

## Suites

Every test file of core and of the eight extensions, one file per R
process, against the lane library with every extension installed from
this worktree (`dev/mocode-par.sh`, runner `dev/mocode-run-test.R`,
logs `dev/mocode-out/plain-files/`, `dev/mocode-out/gated-files/`).
The gated tier (FRMTMB_BRMS_FIT_TESTS=true) ran the 41 files that
mention `mo(`, a brms fixture or a brms fit; it skipped nothing. Every
log's `lib:` line names `wt-mocode-lib` for frmtmb and the package
(452 of 452 plain, 41 of 41 gated). The one-session runs are
`dev/release/run-sessions.sh mocode-plain 8` with FRMTMB_RELLIB set to
the lane library (log `dev/release/mocode-plain.log`); `files=270/280`
counts the files that recorded an expectation, and the forward pass
count equals the per-file total. OpenBLAS 0.3.26 with LAPACK
(`dev/ciharden-openblas.sh 0.3.26 lapack`; 74 Rblas and 616 Rlapack
symbols routed to OpenBLAS) ran the four core files that pin `mo()`
numbers and frmtmb.sample's new file (`dev/mocode-out/ob-*.txt`).
Generated by `dev/mocode-summary.R` and grep from those logs:

```
source: dev/mocode-out/plain.log 
plain ran 452 of 452 
             pkg files  pass fail error skip warn
          frmtmb   280 19831    0     0  204    0
 frmtmb.coupling    11   551    0     0    5    0
      frmtmb.eam    32  1838    0     0    5    0
   frmtmb.latent    12   681    0     0    2    0
    frmtmb.learn    15   430    0     0   13    0
      frmtmb.ode    11   549    0     0    1    0
   frmtmb.sample    68  3060    0     0   10    0
  frmtmb.spatial     8   453    0     0    5    0
 frmtmb.survival    15   552    0     0    0    0
files with fail, error or warn: 0 
files without a RESULT line: 0 

source: dev/mocode-out/gated.log 
gated ran 41 of 41 
             pkg files pass fail error skip warn
          frmtmb    33 4474    0     0    0    0
   frmtmb.sample     7  653    0     0    0    0
 frmtmb.survival     1   15    0     0    0    0
files with fail, error or warn: 0 
files without a RESULT line: 0 

SESSION frmtmb forward files=270/280 pass=19831 fail=0 err=0 skip=204 warn=0
SESSION frmtmb reverse files=270/280 pass=19981 fail=0 err=0 skip=204 warn=0
SESSION frmtmb preload files=270/280 pass=20074 fail=0 err=0 skip=206 warn=0
SESSION frmtmb.coupling forward files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling reverse files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.coupling preload files=11/11 pass=551 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam forward files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam reverse files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.eam preload files=32/32 pass=1838 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.latent forward files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent reverse files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.latent preload files=12/12 pass=681 fail=0 err=0 skip=2 warn=0
SESSION frmtmb.learn forward files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn reverse files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.learn preload files=15/15 pass=430 fail=0 err=0 skip=13 warn=0
SESSION frmtmb.ode forward files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode reverse files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.ode preload files=11/11 pass=549 fail=0 err=0 skip=1 warn=0
SESSION frmtmb.sample forward files=63/68 pass=3060 fail=0 err=0 skip=10 warn=0
SESSION frmtmb.sample reverse files=63/68 pass=3060 fail=0 err=0 skip=10 warn=0
SESSION frmtmb.sample preload files=63/68 pass=3060 fail=0 err=0 skip=10 warn=0
SESSION frmtmb.survival forward files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival reverse files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.survival preload files=15/15 pass=552 fail=0 err=0 skip=0 warn=0
SESSION frmtmb.spatial forward files=8/8 pass=453 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.spatial reverse files=8/8 pass=453 fail=0 err=0 skip=5 warn=0
SESSION frmtmb.spatial preload files=8/8 pass=453 fail=0 err=0 skip=5 warn=0
MOCODE-PLAIN ran 27 of 27
LEAK lines: 46, untagged: 0

RESULT test-mo-coding.R pass=4 fail=0 error=0 skip=1 warn=0
RESULT test-mo-coding.R pass=71 fail=0 error=0 skip=0 warn=0
RESULT test-mo-factor-re.R pass=34 fail=0 error=0 skip=2 warn=0
RESULT test-mo-terms.R pass=74 fail=0 error=0 skip=0 warn=0
RESULT test-optima.R pass=107 fail=0 error=0 skip=0 warn=0
```

## Ported brms-suite rows that read `mo()`

Correction from fix round 1: one row does flip, `brmsfit-methods:924`
(see "A ported row that flips" at the end). The first pass's gated
list missed the file that ports it.

`dev/brmsport-verdicts.tsv` and the bin-2 ledger have three kinds of
row that read `mo()`; none flips.

- `brm:96` (`brm(y ~ mo(g) * t2(x))` is refused as invalid): pass
  before and after; the refusal is a parse refusal before any coding.
- `brmsfit-methods:978` (`variables()`): a divergence by policy
  (`variables()` of an ML fit has no `simo_`), unchanged.
- Every `brmsfit-methods`, `emmeans`, `posterior_predict` and `log_lik`
  row on fixture 1, whose formula has `mo(Exp)`: `Exp` is an ordered
  factor (brms's `brmsfit_example1` data), whose coding the fix does
  not change. Measured: the frame entry is D = 4 with the same codes
  (md5 e754ba6e73690fc5fa8db9eb186eb835) on base and lane
  (`dev/mocode-fixture1.R`).

brms's own block "standata correctly prepares data for monotonic
effects" (`tests.standata.R` 379 to 431, `x1 = rep(1:4, 30)`) asserts
`Xmo_1 == x1 - 1` and `Jmo == max(x1) - 1`, the coding this lane
adopts. It is in bin 3S of `dev/brmssuite-classified.tsv`, not ported,
because its formulas are `mo(x1) * mo(x2) * y`, which frmtmb refuses
(`mo(x):mo(w)`). The core test `test-mo-coding.R` asserts the same
two facts on `mo(k)` alone.

## Where frmtmb still differs from brms (not fixed, filed here)

1. **Two categories.** brms fits `mo(x)` on two values (`Jmo = 1`, a
   simplex of one weight, which is a linear term in the code); frmtmb
   refuses fewer than 3 categories, as it did before. The fix moves one
   case across that line: values 1 and 2 fitted in 0.75.0 (D = 2, with
   an empty category 0, so not identified) and are refused now, while
   0 and 1 were refused before and still are. Fitting D = 1 needs a
   zero-length simplex through `mo_simplex()`, `mo_coords()`, the
   sphere basis (`contr.helmert(1)` stops), `mo_search()`, the Dirichlet
   prior and the sampler's softmax. I left it: it is a new case, not the
   coding defect. One value: brms gives `Jmo = 0`; frmtmb refuses.
2. **A fractional new-data value was rounded.** Fixed in fix round 1
   (item 1 there): refused, as brms refuses it; `conditional_effects()`
   holds the variable at `round(median())` (item 2 there).
3. `attr(x, "min")` is not read (brms-40 above).
4. `variables()` of an ML fit still has no `simo_` entries (filed by the
   0.75.0 brmsnames lane; not this lane's).

## Overlap with the other lanes

Rechecked after fix round 1 against both worktrees' `git status`.

- acprior: `R/fit.R` (its hunk is in the `prior` argument docs near line
  115, mine in "Monotonic effects" near 331), `R/compat.R` (its hunk
  near line 2021, mine at the `mo()` row near 1768), `man/frm.Rd`
  (regenerate), `NEWS.md` (both add a top development heading) and
  frmtmb.sample's `NEWS.md`. No shared function.
- loglik: `NEWS.md` and frmtmb.sample's `NEWS.md` only.
- The files only this lane touches: `R/predict.R`,
  `R/conditional-effects.R`, `R/interop.R`, `R/frame.R`,
  `R/objective.R`, `R/methods-fit.R`, `R/confint.R`,
  `R/extension-home.R`, `man/frmtmb-emmeans.Rd`, and the tests
  `test-mo-coding.R` (both), `test-emmeans.R`, `test-se-check.R`.
- `dev/upstream-bugs.md`: I added brms-40 and brms-41; neither other
  worktree has an entry numbered 40 or above.

## What needs the user

- Nothing blocks. The change is BREAKING for every fit whose integer
  `mo()` variable does not start at 0, since 0.15.0.

## R CMD check --as-cran (first pass)

- frmtmb.sample: `Status: OK` (`dev/mocode-check/frmtmb.sample/check.txt`).
- frmtmb: `Status: 1 NOTE`, the expected V8 math-rendering NOTE on the
  HTML manual (`dev/mocode-check/frmtmb/check.txt`); tests passed.

## Test files before and after

Core 279 files at 70476dc3, 280 now (`test-mo-coding.R` added);
frmtmb.sample 67, 68 now (`test-mo-coding.R` added). No test moved,
rewritten or dropped.

## Fix round 1 (review `dev/reviews/2026-10-27-mocode.md`, verdict FIX)

Numbered as the review's findings. Probe: `dev/mocode-fix1.R` (logs
`dev/mocode-out/fix1-lane2.txt`, `fix1-base2.txt`, seed 5, n = 300, and
seed 7, n = 200 for the gap case); brms 2.23.0 on the same data:
`dev/mocode-brms-fix1.R` (2 chains of 2000, seed 3, log
`dev/mocode-out/brms-fix1.txt`) and `dev/mocode-brms-subset.R` (copied
from the reviewer's script and rerun, log `dev/mocode-out/brms-subset.txt`).

### 1. Fractional new data are refused (major)

`mo_codes()` (R/predict.R) checks, in this order, that the new values
are numeric, have no missing value and are whole numbers, and then the
range, and codes `v - min` without rounding. On the 1..4 fit of
`y ~ mo(k)`:

| new `k` | 0.75.0 (rellib-r13) | lane | brms 2.23.0 (review's `brmsfit.txt`) |
|---|---|---|---|
| 0.6 | 2.130228 (the k = 1 value) | refused, "defined at whole numbers only, but new data give 0.6" | "Invalid values in variable 'k': '0.6'" |
| 2.4 | 4.781501 (k = 2) | refused, same | refused |
| 2.5 | 4.781501 (k = 2) | refused, same | refused |
| 4.4 | 6.46676 (k = 4) | refused, same | refused |
| 0 | 2.124897 (code 0, no row) | refused, "outside the fitted range 1..4 (0)" | refused |
| 5 | refused, "range 0..4" | refused, "outside the fitted range 1..4 (5)" | refused |

`marginaleffects` 1.0.0: `avg_slopes(fz, variables = "k")` on `y ~ mo(k)
+ z` was 0 on 0.75.0 (the finite difference 1.99985 rounded back); now
it stops, and marginaleffects' message carries frmtmb's reason ("In
addition, this error message was raised: mo(k) is defined at whole
numbers only, but new data give 1.99985, 2.99985, 0.99985, 3.99985").
`avg_comparisons(fz, variables = list(k = c(1, 4)))` works (4.31, SE
0.172); `variables = list(k = 1)` is refused, because k + 1 = 5 is past
the range. `avg_slopes(variables = "z")` is unchanged (0.5011905).
NEWS and `man/frm.Rd` now say exactly what is refused: a value that is
not a whole number, or that is outside the fitted range.

### 2. `conditional_effects()` and emmeans hold `mo()` as brms does (major)

`ce_ref_value(col, step)`, `ce_second_values(col, int_cond, step)` and
`ce_build_nd(step_vars =)` (R/conditional-effects.R), as the reviewer's
prototype, with one change: the step variables are the union of the
selected predictor's chain (`ce_step_vars()`) and every numeric `mo()`
variable of the model (`ce_mo_vars()`), as brms's `get_int_vars()`
lists the `mo()` variables of the whole model. Without the union a
`mo()` variable read only by another predictor (`sigma ~ mo(k)`) would
be held at its mean and then refused by item 1.

| `y ~ mo(k) + z`, seed 5 | brms 2.23.0 (`brms-fix1.txt`) | lane | 0.75.0 |
|---|---|---|---|
| `effects = "z"`: `k` held at | 3 | 3 | 2.486667 (coded as 2) |
| estimate at the first `z` (-3.498059) | 4.399008 | 4.403423 (= `fitted()` at k = 3) | 3.049038 |
| `effects = "z:k"`: values of `k` | 2 3 4 | 2 3 4, labels 2 3 4 | 1.366844 2.486667 3.606489, labels 1.37 2.49 3.61 |
| `effects = "k"`: grid; `z` held at | 1 2 3 4; 0.01882471 | 1 2 3 4; 0.01882471 | same |
| `emmeans(~ z)` | "Invalid values in variable 'k': '2.48666666666667'" | refused, "emmeans holds 'k' at 2.486667, but mo(k) is defined at the whole numbers 1 to 4 only (brms refuses this too). Give the value: at = list(k = 3), or cov.reduce = function(v) round(median(v))" | 4.811664 |
| `emmeans(~ z, at = list(k = 3))` | 6.157851 | 6.166052 | 6.16605 |

The brms and frmtmb estimates differ by the posterior mean against the
ML estimate, not by a category.

emmeans: I chose the refusal (`emm_check_mo_grid()`, R/interop.R),
which matches brms. Holding the variable at `round(median())` would be
more convenient, but emmeans documents that a numeric covariate is held
at its mean, and its output then labels the grid by the value it was
asked for; a silent `round(median())` would answer another question
than the one the label states. The refusal names both remedies, and
the second gives exactly `conditional_effects()`'s value. A grid with
`at =` or `cov.reduce =` that gives whole values works.

`tests/testthat/test-emmeans.R` (the block "a smooth or mo() term
reaches the marginal means", formerly line 168) asserted that
`emmeans(fm, "f")` on `y ~ f + mo(e)` evaluates at `e = mean(d$e)`,
the behavior this item removes. It now asserts the refusal without
`at =`, and the same design-basis comparison at `at = list(e = 2L)`.

### 3. The grid of a `mo()` effect in a `subset()` response (minor)

`ce_grid_values(step_range =)` takes the categories of the fitted
terms of the selected response that read the variable alone
(`ce_mo_range()`: `seq(min, min + D)`, the intersection over the
response's terms). On `bf(y | subset(k > 1) ~ mo(k)) + bf(y2 |
subset(s2) ~ mo(k))`: `resp = "y"` gives the grid 2 3 4 (estimates
4.781505 6.110253 6.466759), `resp = "y2"` gives 1 2 3 4. The lane's
first pass stopped with "outside the fitted range 2..4"; 0.75.0
stopped with "replacement has 3 rows, data has 4".

brms 2.23.0 is wrong here (rerun, `brms-subset.txt`): with `subset()`
in a multivariate model it codes new data from the new data's own
minimum. `y2` (fitted on 1..4, `Jmo_y2` 3) codes `k = 2 3 4` as
`0 1 2`, where its fit codes them `1 2 3`, and `y` (`Jmo_y` 2) codes
`1 2 3 4` as `0 1 2 3`, a code past its simplex. Without `subset()`
both give `1 2 3`. Filed as brms-41 in `dev/upstream-bugs.md`; neither
the acprior nor the loglik worktree has a brms-40 or brms-41 entry
(grep of their `dev/upstream-bugs.md`), so 40 and 41 are free.

### 4. One refusal for a saved 0.75.0 fit (minor)

`mo_check_fit()` (R/predict.R) calls `mo_min()` on every `mo()` entry.
It runs from `require_fit_packages()` (so `fixef()`, `anova()`,
`drop1()`, `conditional_effects()` and every other method that already
checks the fit's packages) and from `logLik()` (so `AIC()` and
`BIC()`), `vcov()`, `confint()`, `summary()` and `coef()`. The message
is `mo_min()`'s, the same as the prediction paths. `update()` calls
`require_structure_packages()` only, so it refits. On the real 0.75.0
fit (`dev/mocode-saved.R`, log `dev/mocode-out/saved-lane2.txt`):
`fitted()`, `predict()`, `simulate()`, `residuals()`,
`conditional_effects()`, `logLik()`, `fixef()`, `coef()`, `AIC()`,
`BIC()`, `vcov()`, `confint()`, `summary()`, `drop1()`,
`anova(fit, update(fit))` and `frm_sample()` are refused with:

    mo(k) of this fit was coded by frmtmb 0.75.0 or earlier, from 0
    rather than from the minimum of the data (1). Its simplex has a cell
    for each value below the minimum, so the fit is not brms's model and
    its coefficients are not identified. Fit the model again with this
    version of frmtmb.

`update()` gives Intercept 2.083383, slope 1.515906. The 0..3 fit
answers every one of them as before (logLik -427.1819, AIC 864.3637,
`anova()` df 5 and 5). The core test checks one method of each group
(the nine methods and `anova()`) on a constructed old frame, and that
`update()` of it gives the new fit's log-likelihood.

### 5. The build behind the record (minor)

The reviewer was right: the first pass's library was built before the
last edit of `mo_min()`'s message (the 452 files and the 27 sessions
ran the older text; both texts contain the phrase the tests match).
Every run below is on a build installed after the last source edit
(`dev/mocode-out/install5.txt`), and `dev/mocode-saved.R` was rerun on
it (the message above is the new one).

### 6. The standard error at an empty interior category (major, pre-existing)

The cause was larger than the gap. `lp_delta_A()`, which every
prediction standard error reads (`predict(se.fit =)`, `fitted()`,
`frm_lp_basis()`, `conditional_effects()` bands, emmeans' grid route,
`frm_curve()`), held a `mo()` term's VALUE column only: the delta
method left the simplex out of every prediction since 0.15.0. So the
standard errors were too small wherever the simplex is uncertain, and
finite at a category no row is in. `lp_delta_mo()` now adds the
columns `c * D * cumsum(dw / dz)` at each row's code, with `c` the
row's coefficient (the population-level one times its multiplier, plus
the group-level one of `(mo(x) | g)` when random effects are kept: from
`lp$Z` on the fitted rows, and from the unscaled columns that
`pred_design()` now keeps in `re_parts` for new data). With them
`jc_nonest()` finds the gap direction that `se_mo_null()` already put in
the null basis.

Gap case (seed 7, 200 rows, values 1, 2, 5; the fit warns that the
three simplex coordinates have no standard error):

| x | 0.75.0 Est.Error | lane Est.Error | brms 2.23.0 posterior sd |
|---|---|---|---|
| 1 | 0.1213724 | 0.122221 | 0.124137 |
| 2 | 0.09871537 | 0.1135415 | 0.115953 |
| 3 (no row) | 0.08586617 | NaN | 0.746853 |
| 4 (no row) | 0.09500675 | NaN | 0.782651 |
| 5 | 0.1212624 | 0.1212624 | 0.122501 |

`predict(se.fit = TRUE)` and the `conditional_effects()` band `se__`
give the same five numbers, with one warning "2 of 5 predictions move
along a direction the fit does not determine". On `y ~ mo(k) + z`
(seed 5) at k = 1..4, z = 0: lane 0.1193581 0.1259063 0.1178845
0.1243743; brms 0.118422 0.128666 0.11381 0.112717; 0.75.0 NaN NaN NaN
0.1243739 (its coding was not identified).

Correctness of the new columns, tested: `frm_lp_basis()`'s simplex
columns against central differences of `frm_linpred()` in the simplex
coordinates (h = 1e-5), for `mo(k)`, `mo(k) * z` and
`mo(k) + (mo(k) | g)` on new data with `re_formula = NULL`, and on the
fitted rows of the last; the largest difference is below 1e-6 of the
largest entry in each.

Two assertions of `tests/testthat/test-se-check.R` pinned the old
omission and were wrong: on `helper-mo-flat.R`'s data (`inc` in 0, 2,
so category 1 has no row) they asserted that every band of
`conditional_effects(fit, "inc:age")` is finite and that nothing warns.
The grid has `inc = 1`, whose prediction reads the lost weight. They
now assert NaN at `inc = 1`, finite bands at 0 and 2, and, for the
absent case, finite standard errors and no warning from `fitted()` at
`inc` 0 and 2 alone. The new expectations fail on rellib-r13 (4
failures, `dev/mocode-out/t-sec-base.txt`) and pass on the lane (pass
134, `t-sec-lane.txt`).

Doc sentence (`man/frm.Rd`, "Monotonic effects"): "The data do not
determine a prediction at such a category, so its standard error is
`NaN`."

### 7. Not done: two-category `mo()`

Backlog entry (not built, as the review asks). brms 2.23.0 fits
`mo(x)` on two values: `Jmo = 1`, codes 0 and 1, a simplex of one
weight fixed at 1 (Stan's `simplex[1]` has no free parameter), and a
`simo mox1` row with `dirichlet(1)` in `default_prior()` (review's
`brms.txt`); the term is then a linear term in the code. One value
gives `Jmo = 0`. frmtmb refuses both ("mo() needs at least 3 ordered
categories; mo(x) has 2", "has 1"). A build needs a zero-length
simplex: `extras[[zeta]] <- numeric(0)` in frame assembly; a `D = 1`
branch in `mo_simplex()` (return 1) and `mo_coords()` (return
`numeric(0)`), because `mo_sphere_basis(1)` calls `contr.helmert(1)`,
which stops; `mo_search()` and `mo_search_terms()` skipping the term;
the Dirichlet prior's length check and density; `se_mo_gaps()`;
`lp_delta_mo()` (no columns); frmtmb.sample's softmax chart, its
`dirichlet(1)` density and `draws_mo_cols()` (a constant
`simo_mox1[1]` of 1, as brms reports); and brms agreement tests on
`standata()` and on a fit. The 1, 2 case that 0.75.0 fitted had an
unidentified extra cell, so the refusal is no regression in what is
fitted correctly.

### 8. New-data refusals name their cause (nits)

| new data on a 1..4 integer fit | 0.75.0 | lane |
|---|---|---|
| `k = NA` | "missing value where TRUE/FALSE needed" | "mo(k): new data have missing values, and the model has no category for them" |
| an ordered factor | "'round' not meaningful for factors" (the first pass let the base warning "'-' is not meaningful for ordered factors" escape) | "mo(k) was fitted on integers, but new data give a factor. Give the integers 1 to 4" |
| `"2"` | "non-numeric argument to mathematical function" | "... but new data give a character vector. Give the integers 1 to 4" |
| `y ~ mo(mi(k))` at the fit | "could not find function \"mi\"" | "mo(mi(k)): mo() takes an ordered factor or integers, not the latent variable of mi() or me()" (brms: "Monotonic predictors must be integers or ordered factors") |

Each has a test in `test-mo-coding.R`. An ordered-factor fit still
takes `2L` or `"2"` as a level label, as brms does (review's table).

### 9. Style (nits)

- frmtmb.sample's gated block: `allow_warnings(..., "Effective Samples
  Size")` around both samplers instead of `suppressWarnings()`, and a
  comment that R-hat < 1.05 is the convention of
  `test-default-priors-fit.R` and only checks that brms's reference
  chains mixed.
- `man/frm.Rd`: the 37-word sentence is split ("An integer between the
  minimum and the maximum that no row has is also a category. For
  example, `1, 2, 5` gives `D = 4`.").
- `R/compat.R` (the `mo()` row): "coded from its minimum as brms codes
  it; new data must be whole numbers in the fitted range".

### Files touched in fix round 1

`R/predict.R` (`mo_codes()`, `mo_check_fit()`, `lp_delta_mo()`,
`pred_design()`'s `mo_rows` and `re_parts$mo_raw`, `lp_eta_design()`),
`R/conditional-effects.R`, `R/interop.R` (`emm_check_mo_grid()` and the
emmeans help), `R/frame.R` (`mo_call_heads()`, the `mo(mi())`
refusal), `R/extension-home.R`, `R/methods-fit.R`, `R/confint.R`,
`R/fit.R`, `R/compat.R`, `man/frm.Rd`, `man/frmtmb-emmeans.Rd`,
`NEWS.md`, `dev/upstream-bugs.md` (brms-41),
`tests/testthat/test-mo-coding.R`, `tests/testthat/test-emmeans.R`,
`tests/testthat/test-se-check.R`,
`extensions/frmtmb.sample/tests/testthat/test-mo-coding.R`, and the new
dev scripts `dev/mocode-fix1.R`, `dev/mocode-brms-fix1.R`,
`dev/mocode-brms-subset.R`; `dev/mocode-saved.R` gained the methods of
item 4.

### Tests seen to fail on rellib-r13

- `test-mo-coding.R` (core, now 112 expectations): base pass 12, fail
  36, error 8, warn 11 (`dev/mocode-out/t-mocoding-base3.txt`); lane
  pass 112, fail 0 (`t-mocoding-lane3.txt`).
- `test-emmeans.R`: base pass 91, fail 1 (the refusal); lane pass 92.
- `test-se-check.R`: base pass 130, fail 4; lane pass 134.

### Suites after fix round 1

One file per process on the build of `install5.txt` (every `lib:` line
names wt-mocode-lib), then the one-session runs of core and
frmtmb.sample, then OpenBLAS 0.3.26 with LAPACK. Generated from the
logs (`dev/mocode-out/plain2.log`, `gated2.log`,
`dev/release/mocode-plain2.log`, `dev/mocode-out/ob2-*.txt`):

```
source: dev/mocode-out/plain2.log 
plain2 ran 452 of 452 
             pkg files  pass fail error skip warn
          frmtmb   280 19870    3     0  204    1
 frmtmb.coupling    11   551    0     0    5    0
      frmtmb.eam    32  1838    0     0    5    0
   frmtmb.latent    12   681    0     0    2    0
    frmtmb.learn    15   430    0     0   13    0
      frmtmb.ode    11   549    0     0    1    0
   frmtmb.sample    68  3060    0     0   10    0
  frmtmb.spatial     8   453    0     0    5    0
 frmtmb.survival    15   552    0     0    0    0
files with fail, error or warn: 1 
    pkg            file pass fail error skip warn
 frmtmb test-se-check.R  126    3     0    0    1
files without a RESULT line: 0 

source: dev/mocode-out/gated2.log 
gated2 ran 100 of 100 
             pkg files pass fail error skip warn
          frmtmb    80 9380    1     0    0    0
      frmtmb.eam     1   23    0     0    0    0
   frmtmb.latent     2  429    0     0    0    0
    frmtmb.learn     1  146    0     0    0    0
   frmtmb.sample    13 1455    0     0    0    0
  frmtmb.spatial     2  160    0     0    0    0
 frmtmb.survival     1   15    0     0    0    0
files with fail, error or warn: 1 
    pkg                      file pass fail error skip warn
 frmtmb test-brms-suite-methods.R  161    1     0    0    0
files without a RESULT line: 0 

SESSION frmtmb forward files=270/280 pass=19878 fail=0 err=0 skip=204 warn=0
SESSION frmtmb reverse files=270/280 pass=20028 fail=0 err=0 skip=204 warn=0
SESSION frmtmb preload files=270/280 pass=20121 fail=0 err=0 skip=206 warn=0
SESSION frmtmb.sample forward files=63/68 pass=3060 fail=0 err=0 skip=10 warn=0
SESSION frmtmb.sample reverse files=63/68 pass=3060 fail=0 err=0 skip=10 warn=0
SESSION frmtmb.sample preload files=63/68 pass=3060 fail=0 err=0 skip=10 warn=0
MOCODE-PLAIN2 ran 6 of 6
LEAK lines: 34, untagged: 0

RESULT test-mo-coding.R pass=4 fail=0 error=0 skip=1 warn=0
RESULT test-emmeans.R pass=92 fail=0 error=0 skip=1 warn=0
RESULT test-mo-coding.R pass=112 fail=0 error=0 skip=0 warn=0
RESULT test-mo-factor-re.R pass=34 fail=0 error=0 skip=2 warn=0
RESULT test-mo-terms.R pass=74 fail=0 error=0 skip=0 warn=0
RESULT test-optima.R pass=107 fail=0 error=0 skip=0 warn=0
RESULT test-se-check.R pass=134 fail=0 error=0 skip=0 warn=0

test-se-check.R rerun after its two assertions were corrected:
RESULT test-se-check.R pass=134 fail=0 error=0 skip=0 warn=0
```

### A ported row that flips: brmsfit-methods:924

The gated tier's `test-brms-suite-methods.R` (pass 161, fail 1) reads
"brms brmsfit-methods:924 is recorded as pass and no longer holds:
STALE: reads 'up', whose last assignment failed: mo() needs at least 3
ordered categories; mo(Exp) has 1". brms's block "update has
reasonable outputs" refits fixture 1 with `newdata` whose `Exp` is 4
on every row, under `testmode = TRUE`, which never samples. 0.75.0
fitted that with D = 4 (codes 4, four empty categories); with the
coding from the minimum the term has one value and is refused. brms
cannot fit it either: `dev/mocode-brms-const.R` (log
`dev/mocode-out/brms-const.txt`) samples `y ~ z + mo(x)` with `x = 4`
on 40 rows, and Stan stops on the data, "Jmo[1] is 0, but must be
greater than or equal to 1", so `brm()` returns a fit with no draws.
The row flips from pass to a refusal, a policy divergence (brms
accepts it only without sampling). It flips with the first pass's
change (the D < 2 refusal of a coding from the minimum); the first
pass's gated list did not include this file, because its text has no
`mo(`. The verdict files are not edited here; consolidation
regenerates the ledger.


### R CMD check --as-cran after fix round 1

- frmtmb: `Status: 1 NOTE`, the V8 math-rendering NOTE on the HTML
  manual; tests in one session `[ FAIL 0 | WARN 0 | SKIP 232 | PASS
  19631 ]` (`dev/mocode-check/frmtmb/frmtmb.Rcheck/tests/testthat.Rout`).
  The first attempt of this pass stopped at "copying to build directory
  failed" because I deleted `dev/ciharden-out/` while R CMD build was
  copying the tree; the rerun is the one reported.
- frmtmb.sample: `Status: OK` (its tests changed).
- `dev/ciharden-out/` is deleted.

### What needs the user after fix round 1

- `brmsfit-methods:924` flips from pass to a refusal (above); the
  ledger is regenerated at consolidation.
- emmeans refuses a `mo()` variable held at its mean, as brms does
  (item 2). Holding it at `round(median())` instead is a one-line
  change in `emm_check_mo_grid()` if the user prefers convenience.
- Two-category `mo()` is not built (item 7).

## Fix round 2 (review's "Fix round 1 check", verdict FIX on R1 and R2)

### R1. marginaleffects carries the `mo()` simplex (major)

Done, not qualified. `interop_extra_comps(model, simplex = TRUE)`
(R/interop.R) is an ordinal fit's extras, then the coordinates of
every `mo()` simplex the outer vector holds. `interop_coef_names()`,
`interop_coef_vector()`, `set_coef.frmtmb_fit()`, `interop_vcov()` and
`coef_vcov_root()` read it; `hyp_par_cov(fit, extra =)` returns the
joint covariance of those extra components. The reviewer's pointer,
`simplex_extra_comps()`, lists the logistic-normal correlation charts
only, not the `mo()` coordinates, so a separate list was needed.

`dev/mocode-se3.R` (the reviewer's `dev/review-mocode-se3.R` plus the
round-2 checks; logs `dev/mocode-out/se3-lane.txt`, `se3-base.txt`).
Exact reference: `lm(y ~ factor(k) + z)` with the ML sigma. Seed 5,
n = 300, k = 1..4, z = 0:

| | k = 1 | k = 2 | k = 3 | k = 4 |
|---|---|---|---|---|
| exact | 0.1156304 | 0.1219742 | 0.1142026 | 0.1204897 |
| `predictions()`, lane | 0.1156304 | 0.1219741 | 0.1142026 | 0.1204896 |
| `predictions()`, lane before this round (review's `se3-lane.txt`) | 0.1156304 | 0.08460426 | 0.1035353 | 0.1204896 |
| `predictions()`, rellib-r13 (0.75.0's coding, not identified) | 40.62469 | 17.89227 | 6.252875 | 0.1204897 |

- `avg_comparisons(variables = list(k = c(1, 4)))`: 4.408083, SE
  0.1667724; exact 4.408082, 0.1667724.
- `avg_slopes(variables = "z")` (not a `mo()` variable): -0.008897168,
  SE 0.05860098; exact -0.008897571, 0.05861361 (base -0.008897429,
  0.05859927). These are marginaleffects' own finite differences.
- `hypotheses(hypothesis = "b2 = 0")`: -0.008897168, SE 0.05861362.
- `get_coef()` names: `(Intercept) z mok sigma_(Intercept) zeta1_1
  zeta1_2`.

At a gap (`se_mo_gaps()` not empty), `get_vcov()` sets the simplex's
rows and columns to NaN and warns once per fit
(`mo_gap_nan_warning()`, class `frmtmb_cov_products`, as the warning of a
factored covariance). marginaleffects then gives NaN for every standard
error of the fit, because a product with NaN is NaN even where the
Jacobian column is 0; the warning says so and names the routes that
give NaN only where the prediction moves along the direction. Gap fit
(seed 7), `predictions()` at x = 1..5: lane NaN NaN NaN NaN NaN with
the warning; base 0.1213724 0.09871537 0.08586617 0.09500675
0.1212624, no warning.

`dfbeta()`, `dfbetas()` and `cooks.distance()` (R/influence.R) and
insight's `get_parameters()` / `get_varcov()` (R/insight.R) read
`interop_coef_vector(simplex = FALSE)`, so their columns stay the
coefficients, as documented: `dfbeta()` on `y ~ mo(k) + z` has
`(Intercept) z mok sigma_(Intercept)` on lane and base alike.

Tests (`tests/testthat/test-mo-coding.R`): "marginaleffects' standard
errors carry the mo() simplex" (`predictions()` against the exact
reference to 1e-4 relative, and NaN with the warning at a gap) fails
on rellib-r13 at lines 313, 319 and 321; "dfbeta() keeps the
coefficients' columns on a mo() fit" is the absent case and passes on
both.

### R2. The emmeans remedy names the variable (minor)

The refusal, NEWS, `man/frmtmb-emmeans.Rd` and R/interop.R now give
`cov.reduce = list(k = function(v) round(median(v)))`, with the term's
own variable, and the help says that a bare function applies to every
numeric covariate. The test asserts the message and that the remedy's
emmean is identical to the one at `at = list(k = 3)`; the message
assertion fails on rellib-r13 (line 291).

### Nits

- `man/frm.Rd`: "so in a maximum-likelihood fit its standard error is
  `NaN`".
- `fixef()` and `variables()` of saved 0.75.0 `frm_sample()` draws:
  not done. frmtmb.sample may not call `frmtmb:::mo_check_fit()`
  (guardrail 3, no `:::`), so it needs a new exported seam, its docs
  and a sample-side test; more than a few lines. Every computed
  quantity on such draws is refused already (review, `draws-lane.txt`).
  Recorded here for a later round.

### Files touched in fix round 2

`R/interop.R`, `R/confint.R` (`hyp_par_cov(extra =)`),
`R/influence.R`, `R/insight.R`, `R/fit.R`, `man/frm.Rd`,
`man/frmtmb-emmeans.Rd`, `NEWS.md`, `tests/testthat/test-mo-coding.R`,
and `dev/mocode-se3.R`. frmtmb.sample is not touched in this round.

### Test runs after fix round 2

`test-mo-coding.R` (core): lane pass 120, fail 0
(`dev/mocode-out/t-mocoding-lane4.txt`); rellib-r13 pass 16, fail 40,
error 8, warn 12 (`t-mocoding-base4.txt`).

One file per process on the build of `dev/mocode-out/install7.txt`
(installed after the last R/ edit), the 67 files that reach interop,
marginaleffects, emmeans, dfbeta, insight or `mo(` (job list
`dev/mocode-out/jobs-r2.txt`), in both tiers; every `lib:` line names
wt-mocode-lib (134 of 134). Then the gated `test-brms-suite-methods.R`
alone, the core one-session runs, and OpenBLAS 0.3.26 with LAPACK
(74 Rblas and 616 Rlapack symbols routed; the copy is deleted again).

```
source: dev/mocode-out/r2plain.log 
r2plain ran 67 of 67 
             pkg files pass fail error skip warn
          frmtmb    61 5561    0     0  103    0
      frmtmb.eam     2  162    0     0    0    0
    frmtmb.learn     1   53    0     0    0    0
   frmtmb.sample     2   44    0     0    1    0
 frmtmb.survival     1   15    0     0    0    0
files with fail, error or warn: 0 
files without a RESULT line: 0 

source: dev/mocode-out/r2gated.log 
r2gated ran 67 of 67 
             pkg files pass fail error skip warn
          frmtmb    61 7101    0     0    0    0
      frmtmb.eam     2  162    0     0    0    0
    frmtmb.learn     1   53    0     0    0    0
   frmtmb.sample     2   69    0     0    0    0
 frmtmb.survival     1   15    0     0    0    0
files with fail, error or warn: 0 
files without a RESULT line: 0 

RESULT test-brms-suite-methods.R pass=161 fail=1 error=0 skip=0 warn=0

SESSION frmtmb forward files=270/280 pass=19886 fail=0 err=0 skip=204 warn=0
SESSION frmtmb reverse files=270/280 pass=20036 fail=0 err=0 skip=204 warn=0
SESSION frmtmb preload files=270/280 pass=20129 fail=0 err=0 skip=206 warn=0
MOCODE-PLAIN3 ran 3 of 3
LEAK lines: 20, untagged: 0

RESULT test-emmeans.R pass=92 fail=0 error=0 skip=1 warn=0
RESULT test-influence-plot.R pass=21 fail=0 error=0 skip=0 warn=0
RESULT test-interop.R pass=45 fail=0 error=0 skip=0 warn=0
RESULT test-mo-coding.R pass=120 fail=0 error=0 skip=0 warn=0
RESULT test-se-check.R pass=134 fail=0 error=0 skip=0 warn=0
```

The one failure of `test-brms-suite-methods.R` is the row
`brmsfit-methods:924` of fix round 1, unchanged.

R CMD check --as-cran of frmtmb for this round: Status: 1 NOTE (the V8 math-rendering NOTE), `dev/mocode-check/frmtmb/check.txt`; it finished after the hand-back. frmtmb.sample is not touched in this round, so its check of fix round 1 (Status: OK) stands.
