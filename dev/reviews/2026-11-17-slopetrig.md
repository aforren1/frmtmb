# Review of lane slopetrig (round of 2026-11-16)

Reviewer: adversarial review, 2026-11-17. Tree: `frmtmb-wt-slopetrig`,
library `wt-slopetrig-lib` (installed 2026-11-17 04:56, after the last
edit of `R/se-check.R` at 04:54), base `rellib-r22`. Code compared with
`git show b255350e:<path>` and `diff --strip-trailing-cr`; no `git diff`.
Every script and log named below is under `dev/reviews/`, prefix
`slopetrig-review-`. At most two R processes ran at once.

## Verdict: FIX

The lane does what it says on the designs it calibrated on, and nothing
measured got worse: on 170 out-of-sample fits the default goes from 114
short (75 silent) to 57 short (14 silent), with no fit where the lane's
default is more than 1e-3 below base's. Item 5058's rule held on 78 true
boundary correlations with a gain of exactly 0 every time. But the new
rule's population ("every row") misses `gr(by =)` blocks, where the
default stays short and silent; the cure it switches on (`autoscale =
TRUE`) is short itself on four designs outside the calibration; and the
filed item on the warnings it newly shows undercounts them and omits a
standard-error skew. Items:

1. **Ratio population (code, `R/autoscale.R`, owned).** Compute the ratio
   over the rows the block covers (rows where any Z column of the block,
   its intercept included, is nonzero), not over every row of the linear
   predictor. A `gr(g, by = h)` block reads 1.41 where the column is a
   calendar year (692 over its own rows); see finding 3. The dose on half
   the rows keeps its 1.4, because its block's intercept is nonzero on
   every row. Show `grby_year` before and after.
2. **File HIGH: the default is short where `autoscale = TRUE` is short.**
   Two random slopes, nested terms, `cs()`, `gr(by =)` and a
   factor-by-numeric column (finding 2). Scope the claim "`autoscale =
   TRUE` within 1e-3 on 1599 of 1600" in the findings and roxygen to one
   `us()` slope on `(1 + x | g)`.
3. **Rewrite the filed MEDIUM item** on the warnings at the maximum with
   all the kinds and counts and the standard-error skew (findings 4 and
   5).
4. **Correct the findings record**: the test checks standard errors only
   on the `y ~ x + ...` form (finding 4); newly engaged and not needed is
   943 of 1200 (finding 1); record the third R process (finding 8).
5. **Add to "Needs the user"** the measured alternatives (finding 7): the
   start at unit / root mean square fixed 10 of 10 of item 5033's cases
   with no pre-fit, and the centered fit reached the maximum on every
   design where `autoscale = TRUE` failed.

## Finding 1: the calibration regenerates, and its reference is sound

Script `sum.R` (my own summarizer; it regenerates x from the lane's
generator and seeds rather than using the lane's closed form), log
`sum.log`:

- The F, T and ref arms are identical on both libraries on 1920 of 1920
  fits; the regenerated nonzero-row ratio matches the tsv on 1920 of 1920.
- Left short by the spread rule alone: **257, 184 silent**, up to 213.4
  units. Matches.
- Short under `autoscale = FALSE` by ratio band (spread >= 0.05; the
  lane's bands are over nonzero rows): [1, 1.8) **0 of 100**; [1.8, 2) 0 of
  72 and [2, 2.5) 1 of 128, so **1 of 200** at 1.8 to 2.5. Matches. Over
  every row: 0 of 340, 0 of 60, 1 of 100, 1 of 196, 10 of 208, 13 of 312,
  36 of 144, 196 of 240.
- Newly engaged by the ratio rule: 1200; of them not needed
  (`autoscale = FALSE` within 1e-3): **943 of 1200 (79 percent)**. The
  lane's "1131" is of all 1520 engaged fits (its text says so); the brief
  for this review paraphrased it as "1131 of 1200", which is wrong.
- The engaged flag of the lane arm equals `spread < 0.05 | ratio_all > 2`
  on 1920 of 1920.
- Reference: the best of the four arms. The centered, standardized arm is
  never more than 5.88e-5 below that best, so "best" and "the centered
  fit" give the same counts at 1e-3 (both 389 short F, 257 short base D,
  1 short lane D). At 1e-6 the lane's default has 12 short (10 silent),
  all below 1e-3: the tolerance does not hide anything material.
- The margin is thin: the one short fit at 1.8 to 2.5 is at 2.086, 4
  percent above the threshold, and only 60 fits sit in [1.8, 2.0).
- Time of the lane's default over `autoscale = FALSE` where newly engaged
  and not needed: median 1.33, 90th percentile 2.0.

Weakness: every grid design is one `us()` block `(1 + x | g)`, 20 groups
of 15, one size of slope sd. The "dose on half the rows was never short"
evidence for "every row" uses a dose at level 1 only (root mean square
below 2, where nothing was short at any ratio), so it does not separate
"every row" from "large column". At a large level it is short (finding 2,
`fac_year`).

## Finding 2: out of sample

Scripts `designs.R`, `oos.R` (base and lane in parallel), summary
`oos-sum.R`, logs `oos-base.tsv`, `oos-lane.tsv`, `oos-sum.log`. Seeds 1
to 10; arms ref (affine twin of the column), F, T, D; short = more than
1e-3 below the best arm of both libraries.

```
design        engaged  short D base (silent)  short D lane (silent)  short T  max lane
binom_year      10       10 (10)                0                      0       4e-12
binom_dose      10        7 (5)                 0                      0       2e-12
nb_year         10       10 (6)                 0                      0       1.5e-08
nb_dose         10        0                     0                      0       1.4e-09
two_slopes      10       10 (6)                10 (2)                 10       43.7
two_slopes_b    10        0                     0                      0       1.7e-09
sigma_year      10        5 (4)                 0                      0       9.4e-10
nested_year     10       10 (5)                10 (1)                 10       64.6
mm_year         10        9 (9)                 0                      0       4.5e-09
outlier_dose     0        2 (2)                 2 (2)                  0       0.055
year_zero30      0        1 (0)                 1 (0)                  0       0.374
fac_year         0        8 (0)                 8 (0)                  7       10.2
fac_dose         0        6 (0)                 6 (0)                  5       8.30
grby_year        0       10 (9)                10 (9)                  7       41.4
ar_year         10        8 (6)                 0                      0       2.9e-05
cs_year         10       10 (10)               10 (0)                 10       41.1
subset_year     10        8 (3)                 0                      0       3.1e-05
TOTAL 170               114 (75)               57 (14)
```

- Fixed: binomial and negative binomial GLMMs, a slope in `sigma`, `mm()`,
  `ar()`, a `subset()` response. The `subset()` response's Z covers its
  own rows only, so its ratio is 717, not diluted.
- Engaged but not fixed, because `autoscale = TRUE` is short itself: two
  random slopes (10 of 10, up to 43.7), a nested term (10 of 10, up to
  64.6), `cs()` (10 of 10, up to 41.1). The centered fit reaches the
  maximum on 30 of 30. The change still helps here: silent short fits go
  from 21 to 3, because the reported fit now warns.
- Not engaged: `gr(by =)` (finding 3); a factor-by-numeric column
  `trt * year` (ratio 1.41 over every row, 8 of 10 short; `autoscale =
  TRUE` short on 7 of them, so engaging would not cure it); 3 outliers at
  50 on 5 + 0.1 z (ratio 1.57; 2 short, silent, 0.031 and 0.055); a year
  with 30 percent of rows coded 0 (ratio 1.82; 1 short, 0.37).
- No fit where the lane's default is more than 1e-3 below base's.

## Finding 3: "every row" misses `gr(by =)`

`probe.log`: `(1 + x | gr(g, by = h))` with x = 2000 + 3 z builds one block
per level of `h`; each block's Z column is 0 on the rows of the other
level's groups, and `autoscale_block_cols()` reads ratio 1.41 for both
blocks (692 over the block's own rows). The default does not engage, and
the fit is short on 10 of 10 seeds, **9 of them silent**, by 5.79 to 41.4
units (`oos-sum.log`), the failure of item 5033 itself. Those zeros are
structural (the group is not in the block), not data: the block's
intercept column is 0 there too. The fix in item 1 separates the two.
`autoscale = TRUE` is still short on 7 of 10 (finding 7), so item 1 only
makes it engage; item 2 files the rest.

No case was found where the rule engages wrongly apart from the cost of
the 943 fits that did not need it. Not tested: `mi()` (brms does not allow
it inside a random term, and a response `mi()` keeps every row in Z).

## Finding 4: back-transforms

Script `back.R`, log `back.log`; 8 newly engaged cases (year, 5 + 0.1 z,
age 40 + 10 z, ratios 2.26 and 2.44, a Poisson), both forms.

- Default against `autoscale = TRUE`: **identical** on 94 of 96 accessor
  comparisons (`opt$par`, `fixef()`, `vcov()`, `VarCorr()`, `ranef()`,
  `coef()`, `fitted()`, `predict(newdata)` with a seed, `fitted(newdata,
  re_formula = NA)`, `predict()` intervals, `hypothesis()`,
  `conditional_effects()`). The 2 others: `conditional_effects(fit, "x")`
  errors on both for `y ~ 1 + (1 + x | g)` (no fixed x). A maximum
  likelihood fit has no draws (`as_draws_df()` and `fitted(summary =
  FALSE)` refuse), so draws were checked through `predict()`'s intervals.
- Default against the centered fit mapped back: log-likelihood within
  2.4e-9; sd, correlation, ranef, coef, fixef and fitted values within
  5.6e-5 relative (optimizer tolerance, not 1e-8, since the two fits are
  separate optimizations); the correlation in x's units to 1e-12 (year)
  to 1.1e-5 (ratio 2.44).
- **Standard errors are skewed on the year with no fixed x** (`y ~ 1 + (1 +
  x | g)`, seed 2): fixef SE 0.63 percent below the centered fit's,
  fitted SE 2.87 percent off, `hypothesis("Intercept > 0")` Est.Error
  0.17598 against 0.17710, VarCorr correlation SE NaN. In exact arithmetic
  these are invariant to the reparameterization; the likely cause is the
  lost `theta_3` direction dropped from the joint covariance. The other 7
  cases agree within 2.7e-5 on fixef SEs and 2.4e-4 on fitted SEs. `autoscale = TRUE` gives the same numbers,
  so this is not a back-transform error of the lane, but the default now
  reports them. The lane's test checks `fixef` Est.Error only when x is a
  fixed effect and does not check `fitted()` SEs, so it cannot see this;
  the findings' "checks each but the draws" overstates it.

## Finding 5: the 36 warnings at the maximum are false alarms

Script `warn36.R`, log `warn36.log`, the 36 fits regenerated from the
lane's tsv:

- All 36 are within 3.2e-8 of the centered fit, which gives **0
  warnings** on all 36 (and 0 `VarCorr()` warnings).
- The correlation in x's units is determined: -0.99999953 and the like,
  with a delta-method SE from the centered fit of at most 3.3e-7.
- What the default prints on them: 28 "Standard errors are not available
  ... the data do not determine it" (false); 28 `VarCorr()` "move along a
  direction the fit does not determine" (false); **8 "Large maximum
  absolute gradient ... The fit may not have converged"** (false, and not
  in the filed item); 13 "false convergence (8)"; 1 fit (Poisson, c0 =
  1000) loses the SEs of both sds, theta_1 and theta_2 (not in the filed
  item).
- Out of sample the same happens: fits at the maximum that newly warn,
  lane against base: `binom_year` 4 against 0, `nb_year` 3 against 0,
  `sigma_year` 6 against 3, `ar_year` 9 against 1, `mm_year` 1 against
  0, `subset_year` 3 against 0.

On base all 36 were short (29 silent), so the default trades a wrong
value for a correct one with a wrong warning, which is the better side.
But the project rule says a warning on a correct fit is worse than none,
so the filed MEDIUM item must say how often, which kinds, and that a
gradient warning says "may not have converged" (item 3). Whether the
SE check should read curvature in the pre-fit's coordinates is the
user's call.

## Finding 6: item 5058's rule

Script `guard.R`, log `guard.log`: 7 true-boundary designs (v29's no
slope variance; a slope proportional to the intercept, true correlation
1, gaussian, Poisson and binomial; a 3-dimensional block with a Cholesky
row at its limit; 20000-row gaussian and 10000-row Poisson), seeds 1 to
5, each at y * (1 + 2e-15 k), k = 0 to 2 (x for counts).

- 105 fits, 78 with a correlation parameter past the clamp. **Every gain
  is exactly 0**; no `shortcor` kind; 57 boundary messages, all kept. The
  shortfall against the best of the default, `autoscale = TRUE` and a
  refit started at correlation 0 is at most 2.05e-4 (20000 rows, with a
  "nocurv" warning).
- Can it call a real boundary short? A gain with the rest held is a
  point of higher likelihood, so it proves the stop is short up to
  evaluation noise. The threshold is absolute, `se_edge_tol` = 1e-6, as
  for the sibling sd check (`se_sd_gain_up()`); measured noise was 0 on
  every fit here, including the large ones. A gain of a few ulps cannot
  pass 1e-6. With the gain 0 the code path equals base's by construction.
- `se_cor_short()` asks only `us()` blocks; on `cs()` the short fits warn
  through the convergence check instead (finding 2).

## Finding 7: the user's question, and the alternatives

Script `rules.R`, log `rules.log` (lane grids, base arm):

```
rule                     newly engaged  not needed   short left (silent, max)
ratio > 2                1200           943 (79%)    0
ratio > 2 & rms > 2       880           623 (71%)    0
ratio > 3                1064           808 (76%)    1 (1, 5.34)
ratio > 5                 696           451 (65%)   12 (10, 25.9)
ratio > 10                384           152 (40%)   25 (21, 30.3)
```

- 2: the only threshold that leaves nothing short on the grids; costs a
  median 1.33 times `autoscale = FALSE` where not needed. The margin to
  the nearest short fit is 4 percent (2.086).
- 2 and root mean square above 2: same coverage on the grids with 320
  fewer engaged; but the root mean square is in x's units (years and
  decades give different answers), so it is not one geometric rule. It
  measures the mismatch between the column and the start's unit, which
  is what the start alternative below targets directly.
- 10: leaves 25 short, 21 silent, up to 30.3 units.
- Out of sample the threshold does not matter: every engaged design reads
  47 or more, every missed one 1.83 or less.

The start the lane did not build (script `start.R`, log `start.log`):
`autoscale = FALSE`, start passed by the caller, each random slope's log
sd moved by -log(root mean square). The control (the fit's own initial
vector passed back) equals `autoscale = FALSE` on every case.

- Item 5033's 10 short cases (year seeds 1 to 5, 5 + 0.1 z): **10 of 10
  reach the maximum** (within 9.1e-10) with no pre-fit, against 3.45 to
  41.1 units short at the default start.
- Where `autoscale = TRUE` fails it does not cure either: two slopes 0 of
  3 (25.6, 43.7, 12.5), nested 0 of 3 (2.9, 3.0, 15.7), `cs()` 0 of 3;
  `gr(by =)` improves on 3 of 3 (16.7 to 4.11, 41.4 to 1.35, 18.3 to
  12.3) but reaches none.
- The centered, standardized fit reaches the maximum on all 12 of those,
  and on 168 of 170 out-of-sample fits. Centering, with the covariance
  mapped by A = [[1, -m/s], [0, 1/s]] (exact for `us()`), is what fixed
  every design; scaling by the root mean square leaves the column at the
  same angle to the intercept.

The start shift would change the start of every fit with a random slope
on a column whose root mean square is not 1; it was measured here only
on these 22 cases, not on the fits where the current start is fine.

## Finding 8: process

- The lane briefly ran a third R process (the OpenBLAS check beside two
  others), against the round's limit of two that counts the OpenBLAS
  build script's own R. Its findings do not record it.
- The lane stopped its own first `R CMD check` during the build so the
  check would read the final tree; the findings record this. It was the
  lane's own process. The final check is `Status: 1 NOTE` (V8), log
  `dev/slopetrig-log/check.log`.
- `dev/ciharden-out/` and `dev/slopetrig-check/` are gone, as required.
- Lane test results confirmed from its logs: `test-autoscale.R`
  pass=258 fail=0 on reference BLAS, OpenBLAS 0.3.26 and 0.3.32; on
  rellib-r22 fail=38 err=3, with behavioral failures before the missing
  symbols.

## Other checks

- `R/autoscale.R`: the rule adds `ratio` beside `spread` and one clause
  in `autoscale_decide()`; no back-transform changed. The roxygen's
  geometry is right: root mean square over sd is 1 / sin of the angle
  with the intercept (population sd), and ratio > 2 means a mean more
  than sqrt(3) sd from 0.
- `R/se-check.R`: four hunks as the findings say; `R/se-check.R` is CRLF
  in every worktree under `core.autocrlf = true`, so the line endings are
  not a lane change.
- NEWS numbers (3.3 to 41.1 units, 10 of 10, 7 silent; 5 of 7, 1.68 to
  35.1) match `repro-base.log` and `cor-base.log`.
