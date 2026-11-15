# Short check of lane silent2's fix round (round of 2026-11-13, base ec8178d3)

Scope: the fixes for F1 to F5, the nits and 5214 of
`dev/reviews/2026-11-13-silent2.md`, against the "Fix round" section of
`dev/silent2-findings.md`. The lane was not re-reviewed. Lane library
`wt-silent2-lib` (frmtmb built 2026-11-13 19:10 UTC; no file under `R/`
or `extensions/frmtmb.sample/R/` is newer than the install). Scripts
`dev/reviews/silent2-short-*`, logs `dev/reviews/silent2-short-out/`
(local). No lane file was edited. Free space on C: 144.34 GB at start,
144.41 GB at end. All runs ended before the coordinator's notice to use
at most 2 R processes. None failed on disk space or allocation, so none
was rerun.

## Verdict: FIX

One item: F2's key lets factor responses with equal codes and different
outcomes compare. Everything else holds.

## 1. F1, `re_term_contrasts()`: holds

Reviewer's scripts on the lane library (`rev-contr-lane.txt`,
`rev-probe-lane.txt`): every row 0 (`mm(g, g2)` in-sample vs newdata
2.22e-16, reset 0); "RE-only (1 + factor(cyl) | g) sum -> helmert, max
|diff|: 0".

New cases, `silent2-short-f1.R`, fit under `contr.sum`, fitted values of
9 rows of newdata against the in-sample values of the same rows, under
sum / helmert / treatment (`f1-lane.txt`, `f1-base.txt`):

```text
case                                    lane (sum/helm/treat)   base helmert   base treatment
(1 + factor(cyl):x | g)                 0 / 0 / 0                0              0
(1 + interaction(a, b) | g)             8.9e-16 (all three)      2.36           REFUSED
(1 + relevel(f, "b") | g)               0 / 0 / 0                4.29           REFUSED
(1 + cut(x, 3) | g), all 120 rows       0 / 0 / 0                2.57           REFUSED
(1 + cut(x, 3) | g), 9 rows             REFUSED (new levels)     REFUSED        REFUSED
(1 + factor(cyl) | gr(g, by = byv))     0 / 0 / 0                0.674          REFUSED
(1 + f | gr(factor(gnum)))              0 / 0 / 0                1              REFUSED
(1 + factor(cyl) | gr(factor(gnum)))    0 / 0 / 0                0.379          REFUSED
(1 + factor(cyl) | mm(gm1, gm2))        0 / 0 / 0                0.324          REFUSED
(1 + interaction(a, b) | mm(gm1, gm2))  4.4e-16 (all three)      3.54           REFUSED
(1 + factor(cyl) || g)                  0 / 0 / 0                1.22e-08       REFUSED
sigma ~ (1 + factor(cyl) | g)           0 / 0 / 0                0              0
```

The 8.9e-16 and 4.4e-16 are the same under `contr.sum`, the fit's own
option: rounding of newdata against in-sample, not contrasts.
`mm(factor(a), factor(b))` is refused at the fit on base and lane
("mm(): each membership variable must be a bare column name"), so it
cannot reach prediction; `mm()` was checked with a transformed factor on
the left side instead (above). `cut(x, 3)` on 9 rows refuses because
`cut()` rebuilds its breaks on newdata ("factor cut(x, 3) has new
levels"), on both builds. No silent difference on the lane.

## 2. F2, `loglik_rows_key()`: FIX

`silent2-short-f2.R`, lane (`f2-lane.txt`):

```text
## must run
y ~ 1 vs I(y) ~ x                                          ran
y ~ 1 vs y2 ~ x (copy)                                     ran
cumulative int vs ordered factor                           ran
bernoulli int vs logical                                   ran
## must refuse
rows 1..72 vs 9..80                                        REFUSED
NA rows of different predictors                            REFUSED
y vs rev(y)                                                REFUSED
y vs y + 1e-12                                             REFUSED
categorical fa (a..d) vs fb (w..z), equal codes            ran
cumulative oa (1,2,3,4) vs ob (1,2,3,5), equal codes       ran
cumulative oa (1..4) vs oc (2..5), equal codes             ran
bernoulli factor no/yes vs lo/hi, equal codes              ran
```

The four pairs compare and the old refusals refuse, but a factor
response is keyed by `as.integer(v)` alone (`R/frame.R:2238`), so two
factors of different outcomes whose codes agree pass as the same data
(4 cases above). Fix: key a factor by its labels (`as.character(v)`),
not its codes; where one fit's key is numeric and the other's is labels,
compare `as.character()` of the number with the labels. That keeps
`o` (integer 1..4) against `factor(o, ordered = TRUE)` and refuses all
four. Add `fa` against `fb` to the test.

Side notes, not FIX (false refusals, the safe side): the same labels
with the levels reversed refuse; an integer 2..5 against the ordered
factor of it refuses.

## 3. F3: holds

`silent2-short-f3.R` (`f3-lane.txt`), cs_sr_acat seed 1, logLik
-404.9915602:

```text
class frmtmb_ord_unplaced,frmtmb_warning,warning,condition
  mixture(sratio, acat): component 2 (a threshold the data do not place: it has no standard error) ended at a degenerate boundary ...
class frmtmb_warning,warning,condition
  Standard errors are not available for 4 of 15 parameters. tau_raw2_2, tau_raw2_3, bcs4_2, bcs4_3: ...
fixef NaN Est.Error: mu2_Intercept[2] mu2_Intercept[3] mu2_x[2] mu2_x[3]
```

No number and no "flat". The SE warning's four names are the four NaN
rows of `fixef()`.

`R/fit-end.R` (blob cf7be54d against base 20540993): the only change is
`if (inherits(w, "frmtmb_ord_unplaced")) return(invisible(NULL))` at the
top of `explained()`. Every other warning a `post$fit_check` hook raises
still sets `se_explained`, as on base; no other class changes. One
consequence: when the same hook also raises another warning (for
example `thres_fit_check()`), that warning still silences the SE list.
That is the base rule for those warnings, so not FIX.

Three-component test (seed 20261005 + 11, `starts = 0`, logLik
-530.6851732): degenerate warning (class `frmtmb_ord_unplaced`; component
1 "has no standard error", components 2 and 3 run off or collapsed) and
"Standard errors are not available for 6 of 14 parameters. mu1_x,
tau_raw1_1, tau_raw1_3, tau_raw2_2, tau_raw2_3, tau_raw3_3". True:
component 1's thresholds are -3974, 7029 and 9619 with NaN SEs and
`mu1_x` is -4468 with a NaN SE. `tau_raw2_2` and `tau_raw3_3` are the
log gaps of tied thresholds (mu2 [1] = [2] = -0.0722; mu3 [2] = [3] =
0.0629), flat toward minus infinity. The allowance is by its words
("Standard errors are not available") and hides no false warning.

## 4. F4 and F5: hold

F4: the CLOSED line under 5167 (backlog line 5259) says "in part", "the
1e-15 twin of cum_cum seed 5 only", "5 of the six flat moves (not
cs_sr_acat seed 4)"; the remainder is filed LOW (line 4834). The backlog
against base (CR stripped) has only added lines: 8 hunks, no item text
changed.

F5, `silent2-short-f5.R` (`f5-lane.txt`, `f5-base.txt`):

```text
case                                               lane                     base
vreal(rep(2,120)) fitted(newdata = data), 120      0 vs in-sample           0
vreal(rep(2,120)) fitted(newdata = data[1:10, ])   0 vs in-sample           REFUSED
vreal(rep(2,120)) predict(newdata 10 rows)         ran                      REFUSED
vreal(rep(2,120)) CE resolution = 120              0 vs fitted on grid      0
vreal(rep(2,120)) CE default                       ran                      REFUSED
vreal(factor(rep("a",120))) fitted(newdata 120)    0                        0
vreal(factor(rep("a",120))) fitted(newdata 10)     0                        REFUSED
vreal(factor(rep("a",120))) CE default             ran                      REFUSED
trials(rep(10L,60)) newdata 7 vs trials(10)        0                        REFUSED
vreal(rep_len(c(1,3),120)) newdata = data          REFUSED (new text)       ran
vreal(rep_len(c(1,3),120)) CE resolution = 120     REFUSED                  ran
vreal(c(rep(2,110), rep(NA,10))) fit               optimizer fails (NaN gradient), both builds
```

The value with one distinct element and NA cannot be fitted, so it never
reaches `aterm_scalar()`. By the code, `unique(c(2, NA))` has two
elements, so such a value would stay a vector and be refused, not made
a scalar. A factor constant is coded through its table before
`aterm_scalar()` sees it, so it is the number 1. A constant of length n
on 10 rows is a scalar and correct.

## 5. 5214: holds

(1) `f3-lane.txt`:

```text
RE, starts unset:  degenerate TRUE  names starts TRUE
RE, starts = 0:    degenerate TRUE  names starts FALSE
no RE, default:    degenerate TRUE  names starts FALSE
no RE, starts = 0: degenerate TRUE  names starts FALSE
```

The condition (`R/families.R:7153`) is no starts ran, random-effect
blocks present and `control$starts` NULL.

(2) Not shipped; filed LOW at backlog line 4821 with 1.31 and 1.36
times the fit (n = 1000, 50 groups, seeds 2 and 3; seed 1 builds it) and
`dev/silent2-recost.R`.

## 6. NEWS: holds, two nits

Top `# frmtmb (development version)` and `# frmtmb.sample (development
version)` sections; no lane names, no `dev/` paths, no version numbers;
the control-fit counts are gone; the `log_lik()` bullet is BREAKING.

- Nit: the contrasts bullet ("A prediction rebuilds a group-level term's
  design under the contrasts of the fit") changes values (base helmert
  diffs of 0.32 to 4.29 above), like the `log_lik()` bullet, and is not
  marked BREAKING.
- Nit: `R/families.R:7068` still says "flagged in three cases" over four
  bullets (the merge note's resolution sets "five").

## 7. Tests on the lane library: pass

One file per process, at most 3 at once (`t-*.txt`):

```text
RESULT test-newdata-aterms.R pass=67 fail=0 err=0 skip=0 warn=0
RESULT test-predict-newdata.R pass=26 fail=0 err=0 skip=0 warn=0
RESULT test-subset-rate.R pass=104 fail=0 err=0 skip=0 warn=0
RESULT test-ordinal-mixture.R pass=124 fail=0 err=0 skip=0 warn=0
RESULT test-laplace-draws.R pass=130 fail=0 err=0 skip=0 warn=0       (frmtmb.sample)
RESULT test-brms-draws-order.R pass=39 fail=0 err=0 skip=0 warn=0     (frmtmb.sample)
```

OpenBLAS 0.3.32 with LAPACK (`dev/ciharden-openblas.sh 0.3.32 lapack`:
Rblas 76 to OpenBLAS, Rlapack 617 to OpenBLAS), `OPENBLAS_CORETYPE=Haswell`,
`OPENBLAS_NUM_THREADS=4` (`ob332-ordmix.txt`):

```text
RESULT test-ordinal-mixture.R pass=124 fail=0 err=0 skip=0 warn=0
```

I built `dev/ciharden-out/` (463 MB) and deleted it afterwards.

## For the lane

Fix F2 as in item 2 (labels, not codes, for a factor response) and add
the `fa`/`fb` pair to the test. Optional: the two NEWS nits. Then SHIP.
