# Lane slopetrig: when autoscale's default engages

Round of 2026-11-16, base b255350e (frmtmb 0.83.0), branch
`wt-r84-slopetrig`, library `C:/Users/adf44/source/r/wt-slopetrig-lib`.
Items 5033 (HIGH) and 5058 (MEDIUM) of `dev/test-backlog.md`. The
"before" arm is rellib-r22 (0.83.0). Logs are local and gitignored, under
`dev/slopetrig-log/`; every number below is pasted from the log named
beside it. The prior lane's logs (slopescale, round of 2026-11-13) are
copied to `dev/slopetrig-prior/` because that tree is deleted after this
round starts.

Status: done (fix round 1 included).

## Step 1: item 5033 reproduced on rellib-r22

`dev/slopetrig-repro.R base`, log `dev/slopetrig-log/repro-base.log`
(the review's `dev/reviews/slopescale-review-uncentered3.R` on one
library; gaussian, 20 groups of 15, seeds 1 to 5; short is against the
same model on the centered, standardized column):

```
    c0    s                form n engaged short silent shortT
 5e-01 0.02 y ~ 1 + (1 + x | g) 5       5     0      0      0
 1e+00 0.05 y ~ 1 + (1 + x | g) 5       4     0      0      0
 5e+00 0.10 y ~ 1 + (1 + x | g) 5       0     4      4      0
 5e+00 1.00 y ~ 1 + (1 + x | g) 5       0     0      0      0
 2e+03 3.00 y ~ 1 + (1 + x | g) 5       0     5      5      0
 5e-01 0.02 y ~ x + (1 + x | g) 5       5     0      0      0
 1e+00 0.05 y ~ x + (1 + x | g) 5       4     0      0      0
 5e+00 0.10 y ~ x + (1 + x | g) 5       0     4      3      0
 5e+00 1.00 y ~ x + (1 + x | g) 5       0     0      0      0
 2e+03 3.00 y ~ x + (1 + x | g) 5       0     5      2      0
TOTAL fits 50 engaged 18 short 18 silent 14 shortT 0 maxshort 41.1 maxshortT 4.77e-09
```

The item's counts hold on 0.83.0: the year 10 of 10 short (7 silent),
5 + 0.1 z 8 of 10 (7 silent), up to 41.1 units; `autoscale = TRUE` is
within 4.77e-9 on all 50.

## Step 1b: item 5058 reproduced on rellib-r22

`dev/slopetrig-cor.R base`, log `dev/slopetrig-log/cor-base.log`
(`dev/slopescale-verdict.R`'s design, `y ~ 1 + (1 + x | g)`,
`autoscale = FALSE`, seeds 1 to 6; `d_*` is the change of the objective
when the correlation parameter alone moves to the clamp, to 1 and to 0;
negative is a gain):

```
  fam     s seed conv          t     short n_warn boundary_msg   d_clamp     d_one    d_zero
  gau 1e-04    4    0  2.656e+02  6.08e+00      0         TRUE -2.01e-04 -4.73e-01 -4.85e-01
  gau 1e-04    5    0  2.714e+02  3.90e+00      0         TRUE -9.19e-04 -1.57e+00 -5.05e-01
  gau 1e-05    1    1  1.349e+03  2.54e+00      0         TRUE -1.26e-03 -1.87e+00 -4.86e-01
  gau 1e-05    2    0  1.172e+03  4.51e+00      0         TRUE -2.78e-03 -3.35e+00 -8.82e-01
  gau 1e-05    3    1  1.605e+03  5.24e+00      1         TRUE -2.18e-03 -3.42e+00 -3.26e+00
  gau 1e-05    6    1 -2.281e+03  1.68e+00      1         TRUE -9.64e-04 -1.28e+00  5.75e-01
 pois 1e-05    5    0 -6.571e+02  3.51e+01      0         TRUE -9.95e-03 -1.61e+01 -2.36e+01
TOTAL fits 24 ; t past clamp 14 ; boundary message 7 ; short>1e-3 18 ; short and boundary message 7 ; short and silent 5 ; back move gains>1e-6 among past clamp 14
```

(The 7 rows with the boundary message; the full table is in the log.)
All 7 are short by 1.68 to 35.1 units. Five give no warning: the item's
four (gaussian, 1e-4 seeds 4 and 5, 1e-5 seeds 1 and 2) and a fifth the
item did not list, poisson 1e-5 seed 5, 35.1 units short. On all 14
fits whose parameter is past the clamp, the move to the clamp alone
gains (2.01e-4 to 0.0195).

## Step 2: calibration of the ratio threshold (rellib-r22)

Scripts `dev/slopetrig-calib.R` (coarse grid) and
`dev/slopetrig-calib2.R` (fine grid between ratios 2 and 7), summarized
by `dev/slopetrig-calib-sum.R` and `dev/slopetrig-calib-tab.R`. Each
case runs four arms in one process: `ref` (the same model on the
centered, standardized column), `F` (`autoscale = FALSE`), `T`
(`autoscale = TRUE`) and `D` (the default of the library). "Short" is
more than 1e-3 below the best of the four. The ratio in the tables is
over the nonzero rows (the designs below the dose rows have no zero
rows, so it is the ratio over every row there).

Coarse grid: x = c0 + r c0 z (c0 0.1 to 1000, r 1e-3 to 1), a calendar
year (2000 + s z, s 1 to 30), a dose on half the rows (0 or 1 + r z, r
0.01 to 0.3), the short check's two-value and ulp designs, and a 0/1
column on a level (c0 + trt + 0.01 z); gaussian and poisson; both
forms; seeds 1 to 5; 960 cases. Fine grid: target ratios 2, 3, 4, 5, 7
for c0 1, 10, 1000 and the dose, and the 0/1 column on levels 0.5 to 3;
seeds 1 to 10; 960 cases.

`dev/slopetrig-log/calib-base-sum.log`:

```
cases 960 ; NA ll: F 0 T 0 D 0 ref 0
short > 1e-3: F 366 T 0 D(base) 234 ref 0 ; max short T 2.77e-05
by spread rule alone (spread < 0.05): engaged 320 ; short F not engaged 234
among spread >= 0.05: short F 234 of 640 ; min ratio of a short F 3.555 ; max ratio of a not-short F 1007
```

`dev/slopetrig-log/calib2-base-sum.log`:

```
cases 960 ; NA ll: F 0 T 0 D 0 ref 0
short > 1e-3: F 23 T 1 D(base) 23 ref 0 ; max short T 0.00183
among spread >= 0.05: short F 23 of 960 ; min ratio of a short F 2.086 ; max ratio of a not-short F 7.627
```

The two grids together, the 1600 fits with spread >= 0.05 (those the
spread rule leaves alone), short under `autoscale = FALSE` by root mean
square (rows) and ratio (columns)
(`dev/slopetrig-log/calib-tab-base.log`):

```
         [1,1.8) [1.8,2.5) [2.5,3.5) [3.5,4.5) [4.5,6) [6,9) [9,15) [15,40) [40,Inf)
  rms<2  0/40    0/120     0/152     0/108     0/80    0/80  0/40   .       .
  2-5    .       .         .         .         0/40    2/40  .      .       .
  5-20   0/20    0/40      0/44      2/56      2/40    0/40  2/20   13/40   15/20
  20-200 0/20    .         0/4       0/16      .       .     5/20   8/20    62/80
  >200   0/20    1/40      1/44      7/56      7/40    3/40  1/20   7/20    119/140

short under autoscale = TRUE:
  every cell 0 but rms 2-5, ratio [6,9): 1/40 (1.83e-3)
by family, short F: 166 91 of 800 800   (gau, pois)
```

What this says:

- A short fit needs a column far from 0 against its spread: none at a
  ratio below 1.8 (0 of 100), 1 of 200 at 1.8 to 2.5 (ratio 2.086,
  c0 1000, 5.34 units, silent), rising to 196 of 240 (82 percent) at
  40 and above.
- It also needs a large column: none of the 620 fits with a root mean
  square below 2 was short, at any ratio up to 15 (the dose on half the
  rows, x = 1 + r z, the 0/1 column on levels 0.5 to 2). The ratio
  is a proxy for the mechanism, not the mechanism: the start puts the
  slope's sd at the response's unit per unit of x, so a column of root
  mean square 1000 starts with a slope effect 1000 times the
  intercept's, and a slope nearly collinear with the intercept is then
  stuck. autoscale divides the column by its root mean square, which
  removes both.
- On this one `us()` slope, `(1 + x | g)`, `autoscale = TRUE` is within
  1e-3 on 1599 of 1600 (and 960 of 960 in the coarse grid). The one:
  the 0/1 column on level 3, gaussian, `y ~ x + (1 + x | g)`, short
  1.83e-3. It does not hold beyond one `us()` slope: out of sample,
  `autoscale = TRUE` is short on two random slopes, nested terms,
  `cs()`, `gr(by =)` and a factor-by-numeric column (fix round 1,
  filed HIGH).

Choice: the threshold is the ratio below which nothing was short, 2,
as `autoscale_slope_sd` was set between where fits were short (0.03)
and where they were not (0.06). The ratio is over the rows the block
covers (every row, as first built; corrected in fix round 1), not the
nonzero rows as the spread is: a dose on half the rows (0, or 1 + r z)
is not collinear with the intercept (ratio about 1.4 over the rows its
`(1 + x | g)` block covers, which is every row) and none of its 160
fits was short; over its nonzero rows its ratio is 1 / r and would
engage every one of them. On every grid design the block covers every
row, so the correction changes no grid count (fix round 1). See "Needs
the user" for the alternatives the table supports.

## Step 2b: candidate rules scored on both grids (base arm)

`dev/slopetrig-calib-rules.R`, log `dev/slopetrig-log/calib-rules-base.log`
(every rule is "spread < 0.05 OR <rule>"; ratio_all is over every row;
short is under `autoscale = FALSE`, more than 1e-3):

```
cases 1920 ; spread rule alone engages 320 
none (spread rule only)    engaged  320 (by rule only    0) ; short F left 257 (max    213, silent 184) ; engaged not needed  188 ; short T where engaged 0
ratio_all > 2              engaged 1520 (by rule only 1200) ; short F left   0 (max      0, silent   0) ; engaged not needed 1131 ; short T where engaged 1
ratio_all > 2 & rms > 2    engaged 1200 (by rule only  880) ; short F left   0 (max      0, silent   0) ; engaged not needed  811 ; short T where engaged 1
ratio_all > 3              engaged 1384 (by rule only 1064) ; short F left   1 (max   5.34, silent   1) ; engaged not needed  996 ; short T where engaged 1
ratio_all > 5              engaged 1016 (by rule only  696) ; short F left  12 (max   25.9, silent  10) ; engaged not needed  639 ; short T where engaged 1
ratio_all > 10             engaged  704 (by rule only  384) ; short F left  25 (max   30.3, silent  21) ; engaged not needed  340 ; short T where engaged 0
ratio_all > 20             engaged  640 (by rule only  320) ; short F left  33 (max   30.3, silent  27) ; engaged not needed  284 ; short T where engaged 0
ratio_nz > 2               engaged 1748 (by rule only 1428) ; short F left   0 (max      0, silent   0) ; engaged not needed 1359 ; short T where engaged 1
dose-type designs (dose, twonz, ulpind): ratio_all range 1.25 1.41 
```

The spread rule alone leaves 257 short fits of 1920, 184 of them
silent, up to 213 units. ratio_all > 2 leaves none; it engages 1200
more fits than the spread rule, of which 943 (79 percent) did not need
it (`autoscale = FALSE` within 1e-3; review of lane slopetrig, finding
1). The 1131 in the table is of all 1520 fits engaged, the spread
rule's included. The one fit `autoscale = TRUE` leaves short (1.83e-3, the
0/1 column on level 3, seed 8) is engaged by every rule up to ratio 5.

## Step 3: the trigger (`R/autoscale.R`)

- `autoscale_block_cols()` also returns `ratio`, the root mean square
  over the sd of the column over the rows the block covers (over every
  row as first built; fix round 1), NA where the column does not
  qualify, and `autoscale_plan_z()` carries it in each plan entry.
- `autoscale_decide()` engages the default when a planned random-slope
  column has `spread < autoscale_slope_sd` (unchanged) or
  `ratio > autoscale_slope_ratio` (new, 2).
- No back-transform is new: the plan, the pre-fit and the reported
  fit are those `autoscale = TRUE` already used, and the reported fit
  is the model as written, so estimates, SEs, VarCorr(), ranef(),
  coef(), fitted values and draws are read off it as before. The new
  test checks the estimates (sd, ranef, coef, fixef, fitted, fitted at
  new data) against the fit on the centered, standardized column
  through the affine map (slope by 1 / sd, intercept effects by the
  slope times mean / sd) within 1e-4 relative, and the log-likelihood
  within 1e-8 relative. CORRECTION (fix round 1): as first written it
  checked standard errors only on the `y ~ x + ...` form (the x
  coefficient's Est.Error), not on `y ~ 1 + ...` and not on fitted
  values, so it could not see the skew of the review's finding 4. It
  now also checks the fitted values' standard errors on every case and
  the fixed intercept's standard error on the no-fixed-x form (the
  same parameter in both fits), within 1e-3 relative, except on the
  year with no fixed x, where the skew is real and filed MEDIUM. A
  maximum likelihood fit has no draws (`as_draws_df()` refuses);
  `predict()` intervals use the covariance the standard errors come
  from.
- Docs: `frmtmb_control(autoscale =)` names the rule and says that a
  random slope on a column far from 0 needs it (the old text said such
  a column does not, which holds for a fixed effect only).

### 5033's table after

`dev/slopetrig-repro.R lane`, log `dev/slopetrig-log/repro-lane.log`:

```
TOTAL fits 50 engaged 50 short 0 silent 0 shortT 0 maxshort 4.77e-09 maxshortT 4.77e-09
```

(before: `engaged 18 short 18 silent 14`, step 1). x = 5 + z (ratio
5.1) engages too and was not short before.

### The grids after

`dev/slopetrig-log/calib-lane-sum.log`, `calib2-lane-sum.log` (the "D"
arm is the lane's default; the F arm is unchanged autoscale code and
gives the same short counts as base):

```
coarse: short > 1e-3: F 366 T 0 D 0 ref 0 ; max short T 2.77e-05
fine:   short > 1e-3: F 23 T 1 D 1 ref 0 ; max short T 0.00183
```

The one: the 0/1 column on level 3 (ratio 7), gaussian, `y ~ x + (1 +
x | g)`, seed 8, 1.83e-3 short and silent, where `autoscale = FALSE`
reaches the maximum (filed LOW). Engaged per design over both grids
(column `as_D` of the lane tsv files): the year 80 of 80; x = c0 + r
c0 z 1140 of 1300 (not: the 100 at ratio 1.44 and 60 near ratio 2);
the 0/1 column on a level 220 of 220 (ratio 2.2 and up); the dose on
half the rows 40 of 280 (the 40 at r 0.01 and 0.03, by the spread
rule as before); the two-value and ulp designs 40 of 40 (spread rule).

Warnings of the default arm over the 1920 fits (awk over the four tsv
files of each arm; short = more than 1e-3 below the best arm):

```
base D: short 257 (warned 73, silent 184) ; not short but warned 25 (of them engaged 22)
lane D: short 1 (warned 0, silent 1) ; not short but warned 59 (of them engaged 58)
```

Of the 1200 fits the ratio rule engages (spread >= 0.05), 36 warn at
the maximum, all at a ratio of 200.5 or more (years with sd 1 or 3,
c0 >= 10 with r = 1e-3, the 0/1 column on level 100).
`autoscale = TRUE` warns on the same fits. These are the correlation
near -1 in the column's own units, filed MEDIUM; on base most of these
fits were short instead.

The F arm (unchanged autoscale code) shows item 5058's change: 18 fits
that were silent on base now warn, and every one of them is short, by
3.21 to 56.8 units; no F fit lost a warning (F fits with a warning:
base 107, lane 125).

### lme4 and glmmTMB

`dev/slopetrig-lme4.R`, the item's year and 5 + 0.1 z designs, both
forms, seeds 1 to 5 (20 data sets), ML, on the raw column and on the
centered, standardized column; shortfall against the best of the six
fits of each data set:

```
before (dev/slopetrig-log/lme4-base.log)
  frm_raw : 18 of 20 short > 1e-6 (max 41.1, NA 0)
  frm_c : 0 of 20 short > 1e-6 (max 9.98e-09, NA 0)
  lme4_raw : 16 of 20 short > 1e-6 (max 41.2, NA 0)
  lme4_c : 0 of 20 short > 1e-6 (max 9.64e-10, NA 0)
  tmb_raw : 2 of 20 short > 1e-6 (max 40.7, NA 10)
  tmb_c : 0 of 20 short > 1e-6 (max 6.01e-10, NA 0)
max |frm_raw - lme4_c|: 41.1 ; max |frm_raw - tmb_c|: 41.1
after (dev/slopetrig-log/lme4-lane.log)
  frm_raw : 0 of 20 short > 1e-6 (max 2.32e-09, NA 0)
  frm_c : 0 of 20 short > 1e-6 (max 1.14e-08, NA 0)
  lme4_raw : 16 of 20 short > 1e-6 (max 41.2, NA 0)
  lme4_c : 0 of 20 short > 1e-6 (max 2.44e-09, NA 0)
  tmb_raw : 2 of 20 short > 1e-6 (max 40.7, NA 10)
  tmb_c : 0 of 20 short > 1e-6 (max 2e-09, NA 0)
max |frm_raw - lme4_c|: 2.44e-09 ; max |frm_raw - tmb_c|: 2.32e-09
```

frmtmb on the raw column now agrees with lme4 and glmmTMB on the
centered column within 2.44e-9 on all 20. On the raw column lme4 is
short on 16 (up to 41.2), and glmmTMB errors on 10 and is short on 2:
neither rescales a random slope's column.

### Cost

`dev/slopetrig-cost.R lane`, log `dev/slopetrig-log/cost-lane.log`.
Each time is a block of fits grown past 1.2 s, the minimum of 5 rounds,
arms interleaved; "control" is `autoscale = FALSE` run again, divided
by the first:

```
sleepstudy           engaged FALSE n   32 default 0.0425 s FALSE 0.0434 s default/FALSE 0.98 control 1.01 ; short default -1.5e-11 FALSE -1.5e-11 ; identical par TRUE
unc c0=1 ratio3      engaged TRUE  n   32 default 0.0800 s FALSE 0.0506 s default/FALSE 1.58 control 0.94 ; short default 0 FALSE 1.5e-09 ; identical par FALSE
unc c0=10 ratio2.2   engaged TRUE  n   32 default 0.0816 s FALSE 0.0581 s default/FALSE 1.40 control 1.06 ; short default 0 FALSE 1.7e-11 ; identical par FALSE
binlev c0=1          engaged TRUE  n   32 default 0.0753 s FALSE 0.0544 s default/FALSE 1.39 control 0.93 ; short default 0 FALSE -1.1e-11 ; identical par FALSE
age 40+10z           engaged TRUE  n   16 default 0.0750 s FALSE 0.1013 s default/FALSE 0.74 control 0.98 ; short default 0 FALSE 6.4 ; identical par FALSE
year seed2           engaged TRUE  n   32 default 0.1534 s FALSE 0.0578 s default/FALSE 2.65 control 0.92 ; short default 0 FALSE 41 ; identical par FALSE
unc c0=1000 ratio4   engaged TRUE  n   16 default 0.0837 s FALSE 0.0831 s default/FALSE 1.01 control 0.98 ; short default 0 FALSE 10 ; identical par FALSE
pois 5000 x=10+3z    engaged TRUE  n    2 default 0.8650 s FALSE 0.8750 s default/FALSE 0.99 control 1.00 ; short default 0 FALSE -1.5e-11 ; identical par FALSE
gau 20000 x=10+3z    engaged TRUE  n    1 default 2.2000 s FALSE 1.4200 s default/FALSE 1.55 control 0.98 ; short default 0 FALSE 1.1e-09 ; identical par FALSE
```

Where the rule engages and is not needed, the default costs 0.99 to
1.58 times `autoscale = FALSE` (the control arm moves 0.92 to 1.06).
Where it is needed it costs 0.74 to 2.65 times, and the plain fit was
6.4 to 41 units short. A design from the field: a random slope on age
(40 + 10 z, ratio 4.1) was 6.4 units short under `autoscale = FALSE`
on seed 1, and the default is faster there. sleepstudy's Days (ratio
1.85 over every row) does not engage, and its fit is
`autoscale = FALSE`'s bit for bit. The 20000-row fit ran one fit per
block (n = 1), so each of its arms is the minimum of 5 single fits.

### Over roundings

`dev/slopetrig-ulp.R lane`, log `dev/slopetrig-log/ulp-lane.log`
(y * (1 + 2e-15 k), k = 0..20; x for the poisson counts):

```
engage c0 2000 seed 2 y ~ 1 + (1 + x | g) : engaged=21 loglik=21 sd=21 re=21 re_int=21 coef=21 fitted=21 warn_ok=15 of 21
engage c0 2000 seed 2 y ~ x + (1 + x | g) : engaged=21 loglik=21 sd=21 re=21 re_int=21 coef=21 fitted=21 warn_ok=21 of 21
engage c0 5 seed 1 y ~ 1 + (1 + x | g) : engaged=21 loglik=21 sd=21 re=21 re_int=21 coef=21 fitted=21 warn_ok=21 of 21
engage c0 5 seed 3 y ~ x + (1 + x | g) : engaged=21 loglik=21 sd=21 re=21 re_int=21 coef=21 fitted=21 warn_ok=21 of 21
shortcor 4 1e-04 gaussian : no boundary msg 21 short warning 21 lost shortcor 21 of 21
shortcor 5 1e-04 gaussian : no boundary msg 21 short warning 21 lost shortcor 21 of 21
shortcor 5 1e-05 poisson : no boundary msg 21 short warning 21 lost shortcor 21 of 21 (x, not the counts, moved by ulps)
warning texts seen in the engage cases:
  Optimizer did not report convergence: false convergence (8). A Newton step from the stop g   6
  Standard errors are not available for 1 of 5 parameters. theta_3: the likelihood is flat a  16
```

warn_ok was 15 of 21 on the first case: on 6 roundings the year fit
warns "false convergence (8)" at the maximum instead of the flat
correlation. The test now allows that text too (the filed MEDIUM
item); the rerun is below.
## Step 4: item 5058, a correlation run out at a short stop

Change (`R/se-check.R`, four hunks, all in the boundary verdict for a
correlation; nothing in the parts sparsese and gpselimit own):

- `se_edge_cor()`: a correlation parameter (dimension 2) or a Cholesky
  row (wider `us()` blocks) past the clamp is at its edge only if, also,
  no move back toward 0 with the rest held gains more than
  `se_edge_tol` (new `se_cor_gain_back()`: the row scaled so that its
  largest entry is at the clamp, at 1, and at 0; cached per row and
  point in `fit$cache$se_cor_back`). At a limit the log-likelihood with
  the rest held rises all the way out, so a gain proves the stop is
  short.
- `se_boundary_names()`: a flat correlation that the "further out"
  probe calls an edge, but whose row gains when moved back
  (`se_cor_short()`), gets kind `shortcor` instead of `limit`.
- `se_relabel()` maps kind `shortcor` to lost reason `shortcor`, and
  `se_lost_clauses()` words it: "the fit stopped short of the maximum:
  this correlation is at the limit of the range, where the likelihood
  is flat, but a value nearer 0 raises the log-likelihood, so it is not
  at a boundary. Refit from a value nearer 0 (start =; see
  par_template()), or, for a random slope on a badly scaled column,
  with frmtmb_control(autoscale = TRUE)". It is a warning (the lost-SE
  warning), not the boundary message.

After, on the item's design (`dev/slopetrig-cor.R lane`, log
`dev/slopetrig-log/cor-lane2.log`):

```
TOTAL fits 24 ; t past clamp 14 ; boundary message 0 ; short>1e-3 18 ; short and boundary message 0 ; short and silent 0 ; back move gains>1e-6 among past clamp 14
```

against base (`cor-base.log`): boundary message 7, short and silent 5.
Of the 7 that had the boundary message: 3 (gaussian 1e-4 seeds 4 and
5, poisson 1e-5 seed 5) now carry the new "stopped short" clause; the
4 at 1e-5 (gaussian seeds 1, 2, 3, 6) now carry the convergence
warning that was already there or newly reached ("the curvature at the
stop is unusable, so whether it is the maximum could not be
measured"), because with the correlation no longer held at the clamp
their Hessian row is noise.

Over roundings and with the guard's condition absent
(`dev/slopetrig-guard.R lane`, log `dev/slopetrig-log/guard-lane.log`;
y * (1 + 2e-15 k), k = 0..20; `bnd*` are the true boundary fits of
test-review-v29.R, no slope variance in the data, parameter past the
clamp on every draw):

```
bnd2   fits 21 past clamp 21 boundary msg 21 short warn  0 any warn  0 gain_back>1e-6  0 lost: theta_3=boundary 21
bnd7   fits 21 past clamp 21 boundary msg 21 short warn  0 any warn  0 gain_back>1e-6  0 lost: theta_3=boundary 21
bnd8   fits 21 past clamp 21 boundary msg 21 short warn  0 any warn  0 gain_back>1e-6  0 lost: theta_3=boundary 21
gau4   fits 21 past clamp 21 boundary msg  0 short warn 21 any warn 21 gain_back>1e-6 21 lost: theta_3=shortcor 21
gau5   fits 21 past clamp 21 boundary msg  0 short warn 21 any warn 21 gain_back>1e-6 21 lost: theta_3=shortcor 21
pois5  fits  1 past clamp  1 boundary msg  0 short warn  1 any warn  1 gain_back>1e-6  1 lost: theta_3=shortcor 1
```

On the 63 true boundary fits every back move loses (gain 0); on the 43
short stops the smallest gain is 0.485.

Not changed: the 7 other past-clamp fits of the item's design whose
correlation is flat but not at an edge (the move further out changes
the likelihood) keep base's "flat" warning, which is a warning but says
"the data do not determine it" where the fit is 1.68 to 30.5 units short.
Filed LOW.

## Tests

### Seen failing on rellib-r22

`dev/slopetrig-testfile.R base frmtmb tests/testthat/test-autoscale.R`.
The file as first written (log `dev/slopetrig-log/t-autoscale-base.log`):

```
RESULT test-autoscale.R pass=183 fail=43 err=2 skip=0 warn=1
  34 FAIL "a slope column far from 0 against its spread engages the default"
     (not engaged; 41.137548 units short on the year, fitted values
     0.717 relative off; 9.805252 short on the year's x form; ...)
   9 FAIL "a correlation run out at a short stop is not called a boundary"
     (the boundary message, no "stopped short" warning, lost reason not
     shortcor, on all 3 cases)
   1 ERROR "a slope column near its spread's scale leaves the fit alone"
   1 ERROR "a correlation at its limit keeps the boundary message"
```

The final file (log `t-autoscale-base4.log`):

```
RESULT test-autoscale.R pass=184 fail=38 err=3 skip=0 warn=0
  34 FAIL "a slope column far from 0 against its spread engages the default"
   4 FAIL "a correlation run out at a short stop is not called a boundary"
     (gaussian 1e-4 seed 4: the boundary message, no warning, a lost
     reason "boundary" or "limit", se_edge_cor() marks the parameter;
     then the error below ends the test)
   1 ERROR in the same test: se_cor_gain_back() absent (weak form)
   1 ERROR "a slope column near its spread's scale leaves the fit alone"
   1 ERROR "a correlation at its limit keeps the boundary message"
```

The last two are the guard-absent tests; on base they stop only at the
symbols this lane adds (`autoscale_slope_ratio`, `se_cor_gain_back()`).
What they assert of the fits (bit for bit `autoscale = FALSE`; the
boundary message kept) is base's behavior.

### OpenBLAS

`dev/ciharden-openblas.sh 0.3.26 lapack` and `0.3.32 lapack` (build
logs `dev/slopetrig-log/ob0326-build.log`, `ob0332-build.log`:
"Rlapack: 616 to OpenBLAS" and "617 to OpenBLAS"), run through
`dev/slopetrig-ob.sh` with 4 OpenBLAS threads.

The first OpenBLAS 0.3.26 run of test-autoscale.R failed 2
(`t-autoscale-ob0326.log`): on the poisson 1e-5 seed 5 fit the lost
reason was "nocurv", not "shortcor", and so there was no "stopped
short" text, although the fit warned and had no boundary message. The
verdict depends on the rounding, so the test now asserts what holds on
every rounding: the parameter is past the clamp, se_cor_gain_back() is
above se_edge_tol, se_edge_cor() does not mark it, there is no boundary
message, there is a warning, and no lost reason is "boundary" or
"limit". The same file on rellib-r22 still fails that test (the
boundary message and no warning).

Final runs:

```
reference BLAS  t-autoscale-lane5.log    RESULT test-autoscale.R pass=258 fail=0 err=0 skip=0 warn=0
OpenBLAS 0.3.26 t-autoscale-ob0326c.log  RESULT test-autoscale.R pass=258 fail=0 err=0 skip=0 warn=0
OpenBLAS 0.3.32 t-autoscale-ob0332.log   RESULT test-autoscale.R pass=258 fail=0 err=0 skip=0 warn=0
OpenBLAS 0.3.26 t-v29-ob0326.log         RESULT test-review-v29.R pass=187 fail=0 err=0 skip=0 warn=0
OpenBLAS 0.3.26 t-selimit-ob0326.log     RESULT test-se-limit.R pass=379 fail=0 err=0 skip=0 warn=0
```

Two edits came after some of these runs, and neither changes a fit or
a verdict. The 5058 test's four behavioral assertions moved ahead of
its call to se_cor_gain_back(), so that base fails them before it
stops at the missing symbol. `"shortcor"` moved to a line of its own
in se_lost_clauses()'s order vector (it now comes after "concave"),
off the lines that lanes gpselimit and alapse edit. The 0.3.32 run,
the suites and the sessions read the build before both edits. On the
final build: test-autoscale.R pass=258 fail=0 err=0 skip=0 warn=0
(`t-autoscale-lane6.log`; under OpenBLAS 0.3.26 before the second
edit, `t-autoscale-ob0326c.log`), and `dev/slopetrig-cor.R` gives the
same 24-row table as before (`cor-lane3.log`). R CMD check ran on the
final tree.

Over roundings (`dev/slopetrig-ulp.R`, 21 draws per case; the
allowed-warning list as in the test):

```
reference (ulp-lane2.log):  every engage case 21 of 21 on all 8 checks; shortcor 21, 21, 21 of 21
OpenBLAS 0.3.26 (ulp-ob0326.log): every engage case 21 of 21; shortcor 21, 21, 19 of 21 (poisson)
OpenBLAS 0.3.32 (ulp-ob0332.log): every engage case 21 of 21; shortcor 21, 21, 19 of 21 (poisson)
```

and `dev/slopetrig-guard.R` under OpenBLAS 0.3.26
(`guard-ob0326.log`): the three true boundary fits keep the boundary
message on 63 of 63 draws with gain 0; gau4 and gau5 are shortcor on
42 of 42; the poisson fit warns with lost reason "nocurv".

### Every suite, one file per process

`dev/release/run-par.sh` with 2 jobs, `FRMTMB_RELLIB` = the lane
library, `NOT_CRAN=true`, `FRMTMB_BRMS_FIT_TESTS=true`; logs
`dev/release/slopetrig-suite-core.log`, `slopetrig-suite-ext.log`. The
`lib:` line of all 549 files names the lane library for frmtmb.

```
frmtmb            349 files  pass=28823 fail=0 err=0 skip=15 warn=0
frmtmb.coupling    11 files  pass=551 fail=0 err=0 skip=5 warn=0
frmtmb.eam         33 files  pass=1857 fail=0 err=0 skip=3 warn=0
frmtmb.latent      13 files  pass=700 fail=0 err=0 skip=2 warn=0
frmtmb.learn       16 files  pass=545 fail=0 err=0 skip=2 warn=0
frmtmb.ode         11 files  pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample      92 files  pass=4651 fail=0 err=0 skip=1 warn=0
frmtmb.spatial      8 files  pass=634 fail=0 err=0 skip=0 warn=0
frmtmb.survival    16 files  pass=575 fail=0 err=0 skip=0 warn=0
```

Skips: core's are the drmTMB tier (run separately below), the scale
tier and test-fuzz.R's gate; every extension skip is the scale tier.
`test-drmtmb-agreement.R` with `FRMTMB_DRMTMB_FIT_TESTS=true`
(`t-drmtmb-lane.log`): `pass=131 fail=0 err=0 skip=0 warn=0`.

### Every suite, one session each (forward, reverse, preload)

`dev/release/run-sessions.sh slopetrig-sessions 2`, same environment,
log `dev/release/slopetrig-sessions.log`: 27 of 27 sessions, every one
`fail=0 err=0 warn=0` (frmtmb forward pass=28832, reverse 28997,
preload 29101; each extension's counts equal its one-file run). 218
LEAK lines, every one tagged `[with load of ...]` or `[extension
session state]`; no untagged line.

### R CMD check --as-cran

`dev/slopetrig-check.sh` on the final tree (core only; no extension
changed), `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, pandoc and TinyTeX
on PATH, `R_LIBS` = lane library, rellib-r22, user library; log copied
to `dev/slopetrig-log/check.log`: examples OK, tests OK (29 min),
vignettes OK, `Status: 1 NOTE`, the expected one ("Skipping checking
math rendering: package 'V8' unavailable"). A first check was stopped
during its build so that the check would read the final tree.

## Not done, and why

- The ported brms-suite ledger (dev/brmsport-*) is not regenerated
  (the consolidation does it once). A brms test whose data put a random
  slope on a column whose mean is more than 1.73 sd from 0 now runs the
  pre-fit, which can move its estimates within the optimizer's
  tolerance; no manual row was edited.
- The equivalence survey is not remeasured: its data (`dev/equiv-specs.R`,
  `eq_dat()`) draw every slope column as `round(rnorm(n), 3)`, a ratio
  near 1, so the trigger engages on none of its pairs; and the 5058
  change acts only on a fit whose correlation is past the clamp and
  gains when moved back.
- The start of a random slope's sd at unit / root mean square (Needs the
  user, 2) is not built or measured.
- The 7 past-clamp fits of the 5058 design that keep base's "flat"
  warning are filed LOW, not changed.

## Files touched

- `R/autoscale.R` (owned): `autoscale_block_cols()` (`ratio`),
  `autoscale_plan_z()` (carries `ratio`), new `autoscale_slope_ratio`,
  `autoscale_decide()` (the rule and its roxygen).
- `R/se-check.R`: `se_edge_cor()` (and its roxygen), new
  `se_cor_gain_back()` and `se_cor_short()` after `se_row_at_edge()`,
  `se_boundary_names()` (the `shortcor` kind), `se_relabel()` (one
  mapping), `se_lost_clauses()` (the `shortcor` text and its place in
  the order). Outside the brief's `se_edge_cor()`: the last three
  hunks, each a few lines, needed to give the short stop its warning.
- `R/fit.R`: the roxygen of `frmtmb_control(autoscale =)` only.
- `man/frmtmb_control.Rd` (roxygenised).
- `tests/testthat/test-autoscale.R`: five new tests (the one on
  `gr(by =)` added in fix round 1).
- `NEWS.md`: a new `# frmtmb (development version)` section, two
  bullets.
- `dev/test-backlog.md`: CLOSED lines under items 5033 and 5058; a new
  section "Filed in the round of 2026-11-16", subsection "From lane
  slopetrig", five items (one HIGH and two MEDIUM, rewritten or added
  in fix round 1; two LOW).
- New dev scripts: `dev/slopetrig-*.R`, `dev/slopetrig-*.sh`, this
  file. Logs under `dev/slopetrig-log/` (gitignored); the slopescale
  lane's logs copied to `dev/slopetrig-prior/` (gitignored).

## Merge notes (read from the other lanes' trees on 2026-11-17)

- `R/se-check.R`, `se_boundary_names()`: lane gpselimit adds a branch
  (`gp_basis_white()`, kind "gpwhite") after the `else if` of the gp
  ridge, near this lane's `shortcor` lines before it and its
  `kind <- c(...)` and return lines after it. Both can stand: they
  touch different kinds.
- `R/se-check.R`, `se_lost_clauses()`'s order vector: gpselimit edits
  the line with "boundary" (adds "gpwhite"), alapse the line with
  "nocurv" (adds `ald_se_reason`); this lane adds one line,
  `"shortcor",`, between the line with "concave" and the "nocurv" line.
  Expect a textual conflict there; the resolution is the union.
- No other lane's tree touches `R/autoscale.R`.

## Release

Bump: minor (two BREAKING changes of default behavior). Floors: none.
No new dependency, and no extension calls a function this lane added
or changed, so no extension needs a new core floor.


## Fix round 1 (review `dev/reviews/2026-11-17-slopetrig.md`: FIX)

### 1. The ratio over the rows the block covers (`R/autoscale.R`)

`autoscale_block_cols()` reads `ratio` over the rows where any of the
block's Z columns in that predictor, its intercept included, is
nonzero, no longer over every row of the predictor. A `gr(g, by = h)`
block covers only the rows of its own level's groups; over every row a
calendar year there read 1.41.

- The review's probe of 17 designs (`dev/reviews/slopetrig-review-probe.R`
  on the fixed build, log `dev/slopetrig-log/probe-review-lane2.log`):
  grby_year reads 737 and 652 (was 1.41 and 1.41) and engages; every
  other design reads what it read before (the review's `probe.log`).
- The calibration grids: on all 134 grid designs at seed 1 (both
  forms), the plan's ratio equals the every-row ratio exactly (max
  relative difference 0) and no decision changes
  (`dev/slopetrig-cover.R`, log `dev/slopetrig-log/cover.log`). Every
  grid design is one `(1 + x | g)` block, whose intercept covers every
  row, so no calibration count changes, and none was rerun.
- The dose on half the rows (0 or 1 + 0.2 z): under `(1 + x | g)` it
  reads 1.39, as before (every row 1.39, nonzero rows 5.51). Under
  `(0 + x | g)` the block covers only its nonzero rows and it now reads
  5.51 and engages; that design is outside the grids, and engaging
  costs the pre-fit only.
- grby_year, seeds 1 to 10 (`dev/slopetrig-grby.R`, log
  `dev/slopetrig-log/grby-lane.log`; before is the review's run of the
  previous build, `dev/reviews/slopetrig-review-oos-lane.tsv`):

  ```
  after: engaged 10 of 10 ; D short 7 (silent 2, max 15.6) ; T short 7 ; ref short 0
  before (review's run of the previous build): engaged 0 of 10 ; D short 10 (silent 9, max 41.4)
  ```

  Engaging is not a cure here: `autoscale = TRUE` is itself short on 7
  of 10 (filed HIGH). Short fits go from 10 to 7, silent ones from 9 to
  2, the largest shortfall from 41.4 to 15.6.
- New test "a gr(by =) block reads its ratio over the rows it covers"
  (the plan only, from `frm(dry_run = "frame")`): seen failing on the
  previous lane build (`t-autoscale-prefix.log`: the plan's ratio 1.41
  against about 690 expected, relative difference 0.998; the decision NULL),
  `pass=263 fail=3`; on the fixed build `pass=266 fail=0`.

### 2. Filed HIGH

"The default, and `autoscale = TRUE`, stay short where rescaling cannot
fix the geometry" (`dev/test-backlog.md`, section "Filed in the round
of 2026-11-16", lane slopetrig), with the review's counts: two random
slopes 10 of 10 (up to 43.7), nested 10 of 10 (up to 64.6), `cs()` 10
of 10 (up to 41.1), `gr(by =)` 7 of 10 under `autoscale = TRUE`,
`trt * year` 7 of 10 under `autoscale = TRUE`; the centered,
standardized fit reaches the maximum on all of them. The claim
"`autoscale = TRUE` within 1e-3 on 1599 of 1600" is scoped to one
`us()` slope in step 2 and in `autoscale_slope_ratio`'s roxygen. NEWS
and `frmtmb_control(autoscale =)` now say that with two such slopes, a
nested term, `cs()` or `gr(by =)` the pre-fit can stop short too, and
that centering the column fixes it.

### 3. The MEDIUM item rewritten, and a second MEDIUM

The item on warnings at the maximum now carries the review's counts
(28 "data do not determine it", 28 `VarCorr()`, 8 "Large maximum
absolute gradient ... may not have converged", 13 false convergence, 1
fit with both sd standard errors lost; all 36 within 3.2e-8 of the
centered fit, which gives 0 warnings; the out-of-sample counts). The
standard-error skew on the year with no fixed x (fixef SE 0.63 percent
low, fitted SE 2.87 percent off, `hypothesis()` Est.Error 0.17598
against 0.17710, the correlation's SE NaN; `autoscale = TRUE` the same)
is its own MEDIUM item.

### 4. Corrections to this record

- The test's standard-error coverage (step 3, marked CORRECTION). The
  new standard-error checks (fitted values' SEs on every case, the fixed
  intercept's SE on the no-fixed-x form, 1e-3 relative) leave out the
  year with no fixed x, where the skew is real (filed). On the other
  three cases they hold on 21 of 21 roundings with the reference BLAS
  and with OpenBLAS 0.3.26 (`ulp-lane3.log`, `ulp-ob0326b.log`,
  column `se`).
- Not needed: 943 of the 1200 newly engaged fits (79 percent), step 2b.
- Process: during the first OpenBLAS 0.3.26 build its own BLAS check
  ran as a third R process beside the two of the core suite, against
  the round's limit of two. In this round the build ran beside one
  other R process only.

### 5. Threshold

Kept at ratio > 2. "Needs the user" now carries the review's table of
the rules and the two alternatives it measured.

### Reruns on the fixed build

```
test-autoscale.R, reference BLAS (t-autoscale-lane7.log):   pass=266 fail=0 err=0 skip=0 warn=0
test-autoscale.R, OpenBLAS 0.3.26 (t-autoscale-ob0326d.log): pass=266 fail=0 err=0 skip=0 warn=0
test-gr-by.R, OpenBLAS 0.3.26 (t-grby-ob0326.log):           pass=71 fail=0 err=0 skip=0 warn=0
```

21 roundings (`dev/slopetrig-ulp.R`; reference `ulp-lane3.log`,
OpenBLAS 0.3.26 `ulp-ob0326b.log`): every engage case 21 of 21 on all
10 checks (engaged, log-likelihood, sd, ranef, ranef intercept, coef,
fitted, warnings allowed, the x coefficient and its SE, the new SE
checks); the 5058 cases 21, 21, 21 of 21 (reference) and 21, 21, 19 of
21 (OpenBLAS; the poisson's other 2 warn with reason "nocurv", which
the test accepts); the gr(by =) plan test 21 of 21 on both.

The core test files the change can reach, one per process (every file
that names `gr(... by =)`, `autoscale`, `subset(`, `mm(` or a random
term without an intercept, 86 files; `dev/release/run-par.sh
slopetrig-suite-reach`, 2 jobs, gated, log
`dev/release/slopetrig-suite-reach.log`, `lib:` the lane library on all
86): `ran 86 of 86`, pass=11211 fail=0 err=0 skip=0 warn=0.

The core suite in one session, three ways (`dev/release/run-sessions.sh
slopetrig-sessions2 2 frmtmb`, gated, log
`dev/release/slopetrig-sessions2.log`):

```
SESSION frmtmb forward files=349/349 pass=28840 fail=0 err=0 skip=15 warn=0
SESSION frmtmb reverse files=349/349 pass=29005 fail=0 err=0 skip=15 warn=0
SESSION frmtmb preload files=349/349 pass=29109 fail=0 err=0 skip=17 warn=0
```

135 LEAK lines, every one tagged; no untagged line.

R CMD check --as-cran of core on the fixed tree (`dev/slopetrig-check.sh`,
log `dev/slopetrig-log/check2.log`): tests OK (27 min), `Status: 1
NOTE`, the expected V8 math-rendering note. The check tree and
`dev/ciharden-out/` are deleted. The extension suites were not rerun
in this round (the fix round asked for core's). Their first-pass runs
read the every-row ratio, which equals the covered-row ratio on every
block with an intercept over all of a predictor's rows; a block without
an intercept, a `gr(by =)` block or a block over part of a predictor's
rows in an extension's tests is not covered by them.

## Needs the user

1. The threshold (kept at ratio > 2, fix round 1). The user's decision
   (2026-11-14) was to engage on a "large" ratio of root mean square to
   sd. The calibration puts the line at 2, which is not large: a column
   whose mean is more than 1.73 sd from 0 engages. Days 1 to 10, ages
   in years, visit numbers 2 to 6 and any calendar year all do; Days 0
   to 9 (ratio 1.85) does not. Where it engages and is not needed, the
   fit pays the pre-fit (a median 1.33 times `autoscale = FALSE`, 90th
   percentile 2.0; review finding 1). The review's table of the rules
   on the lane's grids, base arm (`dev/reviews/slopetrig-review-rules.R`,
   `.log`):

   ```
   rule                     newly engaged  not needed   short left (silent, max)
   ratio > 2                1200           943 (79%)    0
   ratio > 2 & rms > 2       880           623 (71%)    0
   ratio > 3                1064           808 (76%)    1 (1, 5.34)
   ratio > 5                 696           451 (65%)   12 (10, 25.9)
   ratio > 10                384           152 (40%)   25 (21, 30.3)
   ```

   - ratio > 2 is the only threshold that leaves nothing short on the
     grids. Its margin is thin: the nearest short fit reads 2.086.
   - ratio > 2 and root mean square > 2 covers the grids with 320 fewer
     engaged, but the root mean square is in x's units (a year and a
     decade answer differently), so it is not one geometric rule; it
     measures the mismatch of the column with the start's unit, which
     the start alternative below targets directly.
   - ratio > 10 leaves 25 short, 21 silent, up to 30.3 units.
   - Out of sample (17 designs, review finding 2) the threshold does not
     matter: every engaged design reads 47 or more, every missed one
     1.83 or less.

   Recommendation, and what is built: ratio > 2.
2. Two alternatives at the source, measured by the review and not built
   (`dev/reviews/slopetrig-review-start.R`, `.log`):
   - The start: each random slope's log sd moved by -log(root mean
     square), with `autoscale = FALSE` and no pre-fit. Item 5033's 10
     short cases (the year, seeds 1 to 5, and 5 + 0.1 z) reach the
     maximum, 10 of 10 within 9.1e-10, against 3.45 to 41.1 units short
     at the present start. It does not cure the designs where
     `autoscale = TRUE` fails: two slopes 0 of 3, nested 0 of 3, `cs()`
     0 of 3; `gr(by =)` improves on 3 of 3 but reaches none. It changes
     the start of every fit with a random slope on a column whose root
     mean square is not 1, and it was measured only on these 22 cases.
   - Centering: the centered, standardized fit reaches the maximum on
     all 12 of those cases and on 168 of 170 out-of-sample fits.
     Centering the column with the covariance mapped by
     A = [[1, -m/s], [0, 1/s]] (exact for `us()`) is what fixes every
     design; scaling by the root mean square leaves the column at the
     same angle to the intercept. That is the HIGH item filed in fix
     round 1.
