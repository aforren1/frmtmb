# Lane arstep: the exact gp() gap beside ar() on a step-shaped field

Round opened 2026-11-06, base 3cd2a0c8 (frmtmb 0.79.1). Worktree
`frmtmb-wt-arstep`, library `wt-arstep-lib`. "Base" is `rellib-r18`.
Logs are in `dev/arstep-out/` (local, gitignored). Every number below is
pasted from a log named beside it.

## Functions this lane changes (for lane fitedge)

`R/gp-start.R`: changed `gp_marginal_gate()`, `gp_marginal_setup()`,
`gp_marginal_data()` (new argument `ar`), `gp_marginal_nll()`,
`gp_marginal_cell()`, `gp_marginal_opt()`, `gp_start_working()`,
`gp_start_scan()`; new `gp_marginal_ar()`, `gp_ar_data()`, `gp_ar_left()`,
`gp_ar_both()`, `gp_marginal_alts()`, `gp_marginal_agrees()`,
`gp_second_descent()`. `R/fit.R`: one call to `gp_second_descent()` after
the first `optimize_obj()` in `fit_assembled()`, and the `start` text of
`?frm`. No change to `R/covstruct.R` or `R/gp-by.R`.

## 1. The review's row, located

`dev/arstep-seedscan.R` rebuilds `dev/review-ncgp-r3-ng.R`'s data (seed
`74000 + 31 * seed + n`) on the step curve, `gar` and `student`, and reads
the default fit against three short-length-scale starts. Seed 1 is the
review's: `gar` n 160 seed 1 gap -25.86 (review -25.9), `student` n 160
seed 1 -6.07 (review -6.07). Log `dev/arstep-out/seedscan-1-3.txt`.

The gap is not one row. On seeds 1 to 6 and 21 to 23 (base, three
starts only, so these are lower bounds on the gap):

```
ASS base gar 160 1 gap -2.586e+01 conv 0 w 0
ASS base gar 160 2 gap -2.681e+01 conv 0 w 0
ASS base gar 80 22 gap -8.181e+00 conv 0 w 0
ASS base gar 80 23 gap -1.375e+01 conv 0 w 0
ASS base gar 80 4 gap -8.537e-01 conv 0 w 0
ASS base gar 40 2 gap -4.344e-01 conv 0 w 0
ASS base student 80 1 gap -9.244e+00 conv 0 w 0
ASS base student 160 1 gap -6.072e+00 conv 0 w 0
ASS base student 160 22 gap -5.976e+00 conv 0 w 0
ASS base student 40 5 gap -8.771e+00 conv 0 w 1
ASS base student 80 6 gap -9.007e+00 conv 1 w 2
```

(Logs `dev/arstep-out/seedscan-{1-3,4-6,21-23}.txt`.)

## 2. The cause on the 25.9 row: not the marginal's choice of mode

Script `dev/arstep-diag.R base gar 160 1`, log
`dev/arstep-out/diag-gar-160-1.txt`; `dev/arstep-diag2.R`, log
`dev/arstep-out/diag2-gar-160-1.txt`. Outer parameters are `beta`,
`betad` (log sigma), `theta` (log sdgp, log lscale), `thetaac` (the
`ar` coefficient; under brms's default `cov = FALSE` it is the natural
coefficient, unconstrained).

The brief's account was: the working-response marginal has more than one
maximum and the start takes the wrong one. Measured, that is not what
happens.

- The working marginal (the response itself, an intercept and white
  noise; `ar()` is left out) has two local maxima over the length scale:
  0.150 (nll -124.28) and 0.208 (nll -128.62), on a profile over 121 log
  length scales. `gp_marginal_opt()` returns 0.2131 (nll -128.78): the
  better one. The start hands sdgp 0.742, lscale 0.213 to the fit, with
  `beta`, `betad` and `thetaac` at the template (0.118, -0.0046, 0).
- From that start the descent stops at ll 130.3068, conv 0, no warning.
- From the marginal's OTHER maximum (0.150), from 0.2076, and from all 15
  of the review's reference starts (sd 0.3, 1.5, 6 by lscale 0.1, 0.3, 1,
  3, 8; the rest at the template), the descent reaches 156.1689. Only the
  marginal's own best point, 0.2131, stops short. `start = list(theta =
  <working start>)` reproduces the default stop exactly, so the start is
  the whole cause; nothing else in the default route differs.

So the full model has two basins the working marginal cannot see, and
the working start lies just inside the wrong one. They differ in the
`ar` coefficient and the intercept, not mainly in the field:

| | short stop | full maximum |
|---|---|---|
| ll | 130.3068 | 156.1689 |
| sdgp | 0.7475 | 0.3166 |
| lscale | 0.2152 | 0.1260 |
| `ar` coefficient | 0.376 | 1.008 |
| sigma | 0.0457 | 0.0425 |
| intercept | 0.0104 | -1.0069 |

At the full maximum the `ar` term is a unit root: a residual that carries
its last value forward takes the step, the intercept sits on the left
level (-1), and the field shrinks to sdgp 0.32. At the short stop the
field takes the step and `ar` is near the truth (0.4). The full model's
profile over the length scale (each other coordinate optimized on the
fit's tape, chained from the short stop) shows the two basins side by
side: at lscale 0.215 the profile from the short stop reads 130.31 with
`ar` 0.376, and its neighbors at 0.168 and 0.276 read 147.68 and 129.12
with `ar` 1.007.

Note for the user: the generating `ar` is 0.4 and the maximum takes 1.008.
That is the likelihood of brms's `cov = FALSE` form, which frmtmb leaves
unconstrained "as in brms". This lane does not change it; see section 7.

## 3. The fix for `ar()`: the exact marginal reads a `cov = FALSE` AR term

brms's `cov = FALSE` AR, frmtmb's default for `ar()`, shifts each row's
mean by `sum_k ar_k (y - mu)` at lag `k` in its group, with zero before
the group's first row (`R/autocor.R`, the block comment above
`autocor_cond_positions()`). So `A (y - mu)` is white noise, with
`A = I - sum_k ar_k L_k` unit lower triangular in time order and
`det(A) = 1`. With the field integrated out, `A (y - X beta)` is normal
with covariance `A V A' + sigma^2 I`. That is a closed form, so the
gaussian marginal that the warm start already maximizes extends to it
exactly, with the `ar` coefficients as more parameters of the scan.

The scan adds AR coefficients 0, 0.5, 0.9 and 1 (first lag; other lags
0) to each cell, and starts the tape descent also from the best cell at
each coefficient. The start hands `thetaac` to the fit with `beta`,
`betad` and `theta`. Because this marginal restates the objective by
hand, it is used only where the objective reads the marginal's value at
the marginal's own maximum, to 1e-6 relative (`gp_marginal_agrees()`);
otherwise the working start runs as before.

Agreement, `dev/arstep-armarg.R lane 0`, log
`dev/arstep-out/armarg-lane-noref.txt`: `gp_marginal_nll()` against
`fit$obj$fn()` at 20 random outer points per row (seed 9100, sd 0.3
around the optimum), max relative difference:

```
AM gar step 160 1 marginal-vs-objective max rel diff 2.521e-14 over 20 points
AM gar step 160 2 marginal-vs-objective max rel diff 3.236e-14 over 20 points
AM gar step 80 22 marginal-vs-objective max rel diff 1.813e-13 over 20 points
AM gar step 80 23 marginal-vs-objective max rel diff 7.795e-14 over 20 points
AM gar step 80 4 marginal-vs-objective max rel diff 4.525e-14 over 20 points
AM gar step 40 2 marginal-vs-objective max rel diff 6.338e-14 over 20 points
AM gar mix 80 1 marginal-vs-objective max rel diff 2.488e-11 over 20 points
AM gar2 step 80 1 marginal-vs-objective max rel diff 6.276e-13 over 20 points
AM ggr step 120 1 marginal-vs-objective max rel diff 1.853e-13 over 20 points
```

(`gar2` is `ar(time = t, p = 2)`; `ggr` is `ar(time = t, gr = g)` with 3
groups of 40.)

The seed scan's short `ar()` rows on the lane build (same log):

| row | base ll | lane ll | base reference (3 starts) |
|---|---|---|---|
| step n 160 seed 1 | 130.3068 | 156.1689 | 156.1689 |
| step n 160 seed 2 | 100.4315 | 127.2413 | 127.2413 (15 starts) |
| step n 80 seed 22 | 7.2879 | 28.7882 | 15.4686 |
| step n 80 seed 23 | -2.7631 | 19.7057 | 10.9886 |
| step n 80 seed 4 | 4.9632 | 5.8169 | 5.8169 |
| step n 40 seed 2 | 8.1142 | 8.5486 | 8.5486 |

All at conv 0 with no warning. On two rows the lane is above every
reference start the scan tried.

## 4. The student rows: the objective itself depends on the inner start

`dev/arstep-diag3.R base student 160 1`, log
`dev/arstep-out/diag3-student-160-1.txt`: of the 15 reference starts, 11
end at conv 1 with warnings, at ll from 187.23 to 207.37, with sdgp 0.65
to 1.16 and lscale 0.37 to 0.58: no clean second basin. The default ends
at 201.2938, conv 0.

`dev/arstep-inner.R base student 160 1`, log
`dev/arstep-out/inner-student-160-1.txt`: at the reference's optimum
the Laplace objective reads -207.3655 from the inner mode the fit found
and -204.7744 from an inner start at zero. The Student-t likelihood is
not log-concave in the field, so the inner problem has more than one
mode and the outer objective is not a function of the outer parameters
alone. At the default's optimum all three inner starts read -201.2938.
So the -6.07 on this row is measured against a value that only a
particular inner path reaches. A second outer start does not address
this; section 5 measures it anyway.

## 5. The second start off the gaussian

Where the working route runs (another family, a non-identity link,
`sigma ~`, `ma()`, `arma()`, `ar(cov = TRUE)`, `mi()`, `me()`, `mo()`,
several responses), `gp_marginal_opt()` now keeps every local maximum
its descents reached, and `gp_start_working()` returns the best one at a
length scale more than a factor 1.1 away from the maximum
(`gp_marginal_alts()`). Where there is one, `frm()` descends a second
time from it after the first descent, and keeps the better optimum
(`gp_second_descent()`, called in `R/fit.R` right after the first
`optimize_obj()`). The first descent is the base build's own: the
objective's state is saved before the second start is evaluated and
restored after it (`obj_state_save()` / `obj_state_restore()`), because
the inner solve starts where the last evaluation left it. Without that
restore one row of the strategy comparison below came out 0.049 below
the base (`student 80 6 scan`, 74.6138 against 74.6624).

Strategies compared before choosing (`dev/arstep-alt.R`, log
`dev/arstep-out/alt-lane.txt`, on an earlier build that had a temporary
option `frmtmb.arstep_alt` and no AR marginal; 13 step rows from the
seed scan; wall is one fit, one process, loaded machine):

- `none`: base behavior.
- `local`: the working marginal's next local maximum (what shipped).
- `scan`: the objective scan of `gp_start_scan_units()` with no
  working start, as a second start on every working-route fit.
- `localscan`: `scan`, only where the marginal has a second maximum.
- `both`: `local` and `scan`.

```
student 160 1 none      ll 201.2938 wall 5.79
student 160 1 local     ll 201.2938 wall 12.59
student 160 1 scan      ll 201.2938 wall 52.56
student 160 22 none     ll 174.5739 wall 8.99
student 160 22 scan     ll 180.5503 wall 24.30
student 40 5 none       ll 25.6883 wall 0.47
student 40 5 scan       ll 36.3409 wall 1.47
student 80 6 none       ll 74.6624 wall 3.13
student 80 6 local      ll 84.0790 wall 7.09
gar 80 23 none          ll -2.7631 wall 0.67
gar 80 23 local         ll 6.7885 wall 0.89
```

`scan` closes three student rows that `local` does not, at 3 to 9 times
the fit's time on n 160 (52.6 s against 5.8 s on one row, for no gain).
On the student rows the objective depends on the inner start (section
4), so a gain there is a different inner mode, not a better start. I
kept `local`, which costs a second descent only where the working
marginal shows a second maximum. Both `ar()` rows in this list are now
the AR marginal's (section 3).

### The warning the brief proposed

The brief asked for a warning by class where the two starts disagree.
On the 150 rows the two descents disagree by more than 1e-3 on 29 fits
(section 9: the second higher on 6, lower on 23). After keeping the
better, 1 of the 29 is short of the reference (`pois_mid two 80`, -1.59,
where the second descent read 0.65 lower than the first) and 28 are at
it. On the 72 rows they disagree on 17, and 1 of those is short after
(`bern two 120`, -0.39, short by -0.416 on the base). A warning keyed on
disagreement would fire on 44 fits at the maximum to flag 2 short ones,
so none was
added: a check that fires on a correct model is worse than none
(dev/lane-rules.md). Generated from `dev/arstep-out/trace-ng-rows.txt`
and `summary-3.txt`.

## 6. The review's 150 off-gaussian rows, before and after

Script `dev/arstep-ng.R` (the review's `dev/review-ncgp-r3-ng.R` data and
15 reference starts, seed 1, i.e. `set.seed(74000 + 31 + n)`; seed 1
reproduces the review's rows: its 13 short and its -25.9, -11 and -6.07).
Logs: `ng-base-a.txt`, `ng-base-b.txt` (base with the 15 starts),
`ngt-base.txt` and `ngt-lane.txt` (default fit only, base and lane run
side by side in two processes for the wall time), and the two strategy
arms `ng-lane-none.txt`, `ng-lane-local.txt` (earlier build), all in
`dev/arstep-out/`. Summarizer `dev/arstep-summary.R`, output
`dev/arstep-out/summary-3.txt`. The reference per row is the best value
any arm or any reference start reached; short is more than 1e-3 under
it; silent is conv 0 and no warning. Pasted from the summary:

```
SUM ng all                     base       rows 150 short  13 silent  12 worst     -25.9 | wall on 150 common rows    450.5 s
SUM ng all                     lane       rows 150 short   5 silent   4 worst       -11 | wall on 150 common rows    506.8 s
SUM ng all                     lane       against base: worse 0, better 8
SUMS ng all                     base       short rows: bern mix 40 1 (-1.07); bern wig 80 1 (-0.258); bern step 40 1 (-0.402); pois_mid two 80 1 (-1.59); pois_mid smooth 40 1 (-0.54); binom5 two 40 1 (-0.2); negbin two 40 1 (-0.138); student step 40 1 (-0.203); student step 80 1 (-11); student step 160 1 (-6.07); gsigma step 40 1 (-1.45); gar two 40 1 (-0.262); gar step 160 1 (-25.9)
SUMS ng all                     lane       short rows: pois_mid two 80 1 (-1.59); student step 40 1 (-0.203); student step 80 1 (-11); student step 160 1 (-6.07); gsigma step 40 1 (-1.45)
SUM ng a: bern pois binom5     base       rows  60 short   6 silent   6 worst     -1.59 | wall on 60 common rows    131.2 s
SUM ng a: bern pois binom5     lane       rows  60 short   1 silent   1 worst     -1.59 | wall on 60 common rows    160.3 s
SUM ng a: bern pois binom5     lane       against base: worse 0, better 5
SUM ng b: the other six        base       rows  90 short   7 silent   6 worst     -25.9 | wall on 90 common rows    319.3 s
SUM ng b: the other six        lane       rows  90 short   4 silent   3 worst       -11 | wall on 90 common rows    346.4 s
SUM ng b: the other six        lane       against base: worse 0, better 3
SUM ng step curve              base       rows  30 short   6 silent   5 worst     -25.9 | wall on 30 common rows     81.9 s
SUM ng step curve              lane       rows  30 short   4 silent   3 worst       -11 | wall on 30 common rows     94.0 s
SUM ng step curve              lane       against base: worse 0, better 2
SUM ng student                 base       rows  15 short   3 silent   2 worst       -11 | wall on 15 common rows     96.6 s
SUM ng student                 lane       rows  15 short   3 silent   2 worst       -11 | wall on 15 common rows    106.9 s
SUM ng student                 lane       against base: worse 0, better 0
SUM ng gar                     base       rows  15 short   2 silent   2 worst     -25.9 | wall on 15 common rows     30.1 s
SUM ng gar                     lane       rows  15 short   0 silent   0 worst         0 | wall on 15 common rows     33.8 s
SUM ng gar                     lane       against base: worse 0, better 2
SUM ng reference starts ran on 150 of 150 base rows
```

So: 13 short to 5, worst 25.9 to 11; 0 rows worse; 8 rows better. Of
the 8, two are the AR marginal's (`gar step 160`, +25.9; `gar two 40`,
+0.262) and six the second descent's (`bern mix 40` +1.07, `bern wig
80` +0.258, `bern step 40` +0.402, `pois_mid smooth 40` +0.54, `binom5
two 40` +0.2, `negbin two 40` +0.138). The five still short: the three
student step rows (section 4), `pois_mid two 80` (-1.59) and `gsigma
step 40` (-1.45), all short on the base by the same amount.

Wall time of the default fit over the 150 rows, base and lane run side
by side on a loaded machine: 450.5 s against 506.8 s, 1.125 times. Over
the 60 low-information rows (set a) 131.2 s against 160.3 s, 1.22 times;
over the other 90, 319.3 s against 346.4 s, 1.085 times.

### The AR grid in the scan is needed

`dev/arstep-phigrid.R lane`, log `dev/arstep-out/phigrid.txt`: the scan's
AR coefficients (0, 0.5, 0.9, 1) against 0 alone (a namespace mutant of
`gp_marginal_opt()`), in one process, interleaved by row, on 35 `ar()`
rows (the seed scan's short rows and curves step and mix at n 40, 80,
160, seeds 1 to 5; 32 distinct rows). With 0 alone the fit is lower on 8
of 32 rows, by up to 13.32 nats (`step 80 22`: 15.4686 against 28.7882),
and higher on none. The grid costs: 80.8 s against 52.2 s over the 32
rows, 1.55 times (counted by an awk line over the log). The test file
now pins the grid on `step 80 22` (section 11).

The same log shows the AR marginal above the base's own reference on
rows the review's table did not flag: `step 160 3` reads 118.7502 where
the base read 75.3382 and its three reference starts found nothing
better (`dev/arstep-out/seedscan-1-3.txt`).

## 7. The review's 192 gaussian curves, before and after

Script `dev/arstep-gauss.R base 21:22` (the review's
`dev/review-ncgp-r3-gauss.R` data, seeds 21 and 22, `set.seed(73000 + 29
seed + n)`, reference the review's closed form) and `dev/arstep-gauss2.R
lane 21:22 0` (no closed form; the reference is the base arm's). Logs
`dev/arstep-out/gauss-base.txt`, `gauss-lane.txt`; summary
`dev/arstep-out/summary-4.txt`:

```
SUM gauss 192                  base       rows 192 short   0 silent   0 worst         0 | wall on 192 common rows    226.6 s
SUM gauss 192                  lane       rows 192 short   0 silent   0 worst         0 | wall on 192 common rows    240.5 s
SUM gauss 192                  lane       against base: worse 0, better 0
```

The log-likelihoods agree to the 6 printed decimals on all 192 rows
(join of the two logs by row, awk). By construction this route changed
only in bookkeeping: with no residual term `gp_marginal_ar()` returns
`NULL`, the scan's AR loop runs once with no coefficient, and the start
returns before any second start is read. The wall times, 226.6 s and
240.5 s, were not taken side by side (the base ran beside its closed
form), so they are not a cost measurement.

## 8. Cost of the AR marginal at size

`dev/arstep-cost.R <arm> 200,500,1000 1`, base and lane run side by side
in two processes; logs `dev/arstep-out/cost-base.txt`, `cost-lane.txt`.
`y ~ gp(x) + ar(time = t)` on n distinct positions, `sin(x) + 0.5
sin(3x)` plus AR(0.4) noise sd 0.1, seed `9300 + n`; the control is `y ~
gp(x)` on the same data, whose route this lane does not change, to read
the load on each process:

```
ACO base ar n 200 rep 1 wall 12.72 ll 145.506748 conv 0 w 0 marginal FALSE
ACO base control n 200 rep 1 wall 4.61 ll 134.616851 conv 0 w 0 marginal TRUE
ACO base ar n 500 rep 1 wall 319.65 ll 396.948574 conv 0 w 0 marginal FALSE
ACO base control n 500 rep 1 wall 158.94 ll 370.432388 conv 0 w 0 marginal TRUE
ACO lane ar n 200 rep 1 wall 14.31 ll 145.506748 conv 0 w 0 marginal TRUE
ACO lane control n 200 rep 1 wall 5.34 ll 134.616851 conv 0 w 0 marginal TRUE
ACO lane ar n 500 rep 1 wall 203.64 ll 396.948574 conv 0 w 0 marginal TRUE
ACO lane control n 500 rep 1 wall 150.64 ll 370.432388 conv 0 w 0 marginal TRUE
```

At n 1000 `gp()` refuses a dense covariance on both builds. Both builds
reach the same log-likelihood at both sizes. At n 200 the `ar()` fit
took 1.13 times the base, with the control at 1.16 times, so the
difference is load. At n 500 it took 0.64 times the base (203.6 s
against 319.7 s; control 0.95): the start lands at the maximum and the
descent is short. The machine was loaded by other lanes: the control at
n 500 took about 150 s, where NEWS for the last release gives 28 s
unloaded. One repeat only; I did not repeat it because the two arms
agree in direction with the 15 `ar()` rows of section 6 (30.1 s against
33.8 s, side by side).

## 9. What the second descent did, row by row, and what it cost

Traced runs on the lane build: `dev/arstep-ng3.R lane 1 <10 cases> 0`
and `dev/arstep-r2ng3.R lane 1 0`, traces in `dev/arstep-trace.R` (time
inside `gp_second_descent()`, the two descents' values, and the working
marginal's best value and second maximum). Logs
`dev/arstep-out/ng3-lane.txt`, `r2ng3-lane.txt`; summary
`dev/arstep-trace-sum.R`, output `dev/arstep-out/trace-sum-1.txt`:

```
TSUM ng3-lane.txt rows 150, second descent on 40, raised by >1e-3 on 6, lowered (kept the first) on 23; seconds 69.4 of 518.4 total
TSUM ng3-lane.txt marginal gap to the second maximum: rows that gained bern mix 40 1 8.66 (+1.07); bern wig 80 1 8.04 (+0.258); bern step 40 1 11.17 (+0.402); pois_mid smooth 40 1 0.05 (+0.54); binom5 two 40 1 0.02 (+0.2); negbin two 40 1 0.35 (+0.138); others: median 4.84, range 0.18 to 166.26
TSUM r2ng3-lane.txt rows 72, second descent on 25, raised by >1e-3 on 1, lowered (kept the first) on 16; seconds 34.7 of 226.2 total
TSUM r2ng3-lane.txt marginal gap to the second maximum: rows that gained bern two n 120 seed 1 0.41 (+0.0258); others: median 38.40, range 0.15 to 169.89
```

Every 150-set row with a second descent (`trace-ng-rows.txt`, an awk
line over the log; gap is the working marginal's second maximum minus its
best, in nats of the working marginal; gain is the second descent's
log-likelihood minus the first's):

```
bern mix 40 1 gap 8.66 gain 1.0746 secs 0.45
bern mix 80 1 gap 4.62 gain 0.0000 secs 1.66
bern mix 160 1 gap 13.79 gain -12.4389 secs 9.61
bern wig 80 1 gap 8.04 gain 0.2584 secs 0.50
bern two 160 1 gap 1.61 gain -1.9765 secs 6.46
bern step 40 1 gap 11.17 gain 0.4022 secs 0.17
pois_low mix 40 1 gap 0.18 gain 0.0000 secs 0.05
pois_low mix 80 1 gap 13.68 gain 0.0000 secs 0.49
pois_low wig 160 1 gap 5.88 gain -6.7156 secs 1.46
pois_low smooth 40 1 gap 0.20 gain -0.0151 secs 0.03
pois_low step 40 1 gap 0.19 gain -0.2581 secs 0.13
pois_low step 80 1 gap 1.79 gain 0.0000 secs 0.31
pois_low step 160 1 gap 1.30 gain -0.3067 secs 1.76
pois_mid wig 40 1 gap 2.60 gain -8.7734 secs 0.04
pois_mid wig 80 1 gap 0.19 gain 0.0000 secs 0.38
pois_mid wig 160 1 gap 16.04 gain -60.3064 secs 1.86
pois_mid two 80 1 gap 1.28 gain -0.6474 secs 0.38
pois_mid smooth 40 1 gap 0.05 gain 0.5395 secs 0.03
pois_mid smooth 80 1 gap 3.30 gain -5.6137 secs 0.25
pois_mid step 80 1 gap 0.63 gain 0.0000 secs 0.36
binom5 mix 40 1 gap 5.07 gain -12.4339 secs 0.23
binom5 mix 80 1 gap 15.91 gain -38.7366 secs 1.56
binom5 wig 80 1 gap 0.27 gain 0.0000 secs 1.89
binom5 wig 160 1 gap 11.86 gain -12.6757 secs 4.50
binom5 two 40 1 gap 0.02 gain 0.2004 secs 0.13
negbin two 40 1 gap 0.35 gain 0.1384 secs 0.55
negbin two 80 1 gap 4.43 gain 0.0000 secs 2.92
negbin two 160 1 gap 8.64 gain -8.0022 secs 14.28
gamma two 40 1 gap 18.81 gain -19.8969 secs 0.16
gamma two 80 1 gap 45.94 gain -45.6167 secs 0.55
gamma two 160 1 gap 166.26 gain -167.5078 secs 3.95
lognorm two 40 1 gap 22.92 gain -22.9187 secs 0.07
lognorm two 80 1 gap 63.98 gain -63.9834 secs 0.22
student two 40 1 gap 26.62 gain -26.6224 secs 0.24
student two 80 1 gap 53.18 gain -53.1830 secs 1.11
student step 160 1 gap 1.95 gain 0.0000 secs 10.22
gsigma mix 40 1 gap 0.72 gain 0.0000 secs 0.07
gsigma two 40 1 gap 25.42 gain -25.6993 secs 0.07
gsigma two 80 1 gap 72.10 gain -72.6331 secs 0.24
gsigma step 40 1 gap 1.06 gain 0.0000 secs 0.06
```

The cost of the second descent, timed inside the fit: 69.4 s of the
518.4 s the 150 fits took (13 percent), and 34.7 s of 226.2 s on the 72
(15 percent). The side-by-side wall ratios of sections 6 and 10 (1.125
and 1.31) carry the load difference of two processes; the in-process
share is the cleaner number.

A threshold on the marginal gap would cut most of that: every gain is
at a gap of 11.17 or less, and the second descents at gaps above 12
never gained. I did not add one. The largest gainer's gap (11.17) is
within 1 nat of where such a threshold would sit, chosen on these same
rows, and a row just past it would lose a gain silently, which is the
defect this lane is about.

Build note. `ngt-lane.txt` (section 6) ran on a build whose
`gp_marginal_data()` took the AR block as a tenth argument; the final
build reads it from the frame instead, so that its signature stays the
one lane fitedge also edits. `ng3-lane.txt` ran the same 150 rows on the
final build: the log-likelihoods agree to the 6 printed decimals on all
150 (join by row, awk). Every other lane run here is on the final build.

## 10. The 72 rows of the ncgp review's fix round 1 (the 48 among them)

Script `dev/arstep-r2ng.R` (copy of `dev/review-ncgp-r2-ng.R`: curves
mix, wig, two; n 60 and 120; seed 1, `set.seed(56000 + 13 + n)`; 18
reference starts). The brief's "48 smooth rows" are read as that
script's eight spellings other than the two Poissons, bernoulli and
1000 trials: the 48 rows behind the last release's "0 of 48 short" in
NEWS. Logs `dev/arstep-out/r2ng-base-ref.txt` (base, 18 starts),
`r2ng-base.txt` and `r2ng-lane.txt` (default fit only, side by side);
summary `dev/arstep-out/summary-6.txt`:

```
SUM r2ng 72                    base       rows  72 short   5 silent   5 worst     -1.77 | wall on 72 common rows    314.6 s
SUM r2ng 72                    lane       rows  72 short   4 silent   4 worst     -1.77 | wall on 72 common rows    412.6 s
SUM r2ng 72                    lane       against base: worse 0, better 2
SUMB r2ng 72                    lane       better: pois_low two 120 1 (+0.105); bern two 120 1 (+0.0258)
SUMS r2ng 72                    base       short rows: pois_low wig 120 1 (-1.77); pois_low two 60 1 (-0.342); pois_low two 120 1 (-0.105); bern mix 120 1 (-1.42); bern two 120 1 (-0.416)
SUMS r2ng 72                    lane       short rows: pois_low wig 120 1 (-1.77); pois_low two 60 1 (-0.342); bern mix 120 1 (-1.42); bern two 120 1 (-0.39)
SUM r2ng the 48                base       rows  48 short   0 silent   0 worst         0 | wall on 48 common rows    221.1 s
SUM r2ng the 48                lane       rows  48 short   0 silent   0 worst         0 | wall on 48 common rows    286.1 s
SUM r2ng the 48                lane       against base: worse 0, better 0
```

0 rows worse. On `pois_low two 120` the lane's first and second descents
both read the base's -170.2040 (trace, section 9), and the fit then
ends at -170.0985 by the steps that follow the descent in `frm()`; both
points are degenerate (the base's field sd is e^-10.2, the lane's length
scale e^-9.5). The side-by-side wall ratio on the 48 is 1.29; the
in-process share of the second descent on the 72 is 15 percent
(section 9).

## 11. The new test file and its evidence

`tests/testthat/test-gp-ar-start.R`, five tests:

1. The AR marginal reads the model's own density: the tape
   (`fit$obj$fn()`), `gp_marginal_nll()` and a reference written from the
   formula (dense `A`, `V = sdgp^2 K + sigma^2 A^-1 A^-T`, no frmtmb code)
   agree to 1e-8 relative at 3 points, on `ar()`, `ar(p = 2)` and
   `ar(gr = )`. Away from the optimum the test sets the AR coefficients
   below 1: past a unit root `A^-1` grows like `ar^n` and the formula's
   `solve()` lost up to 1.2e-3 at `ar` 1.74, n 40, while the tape and the
   marginal still agreed with each other to every printed digit
   (`dev/arstep-testdbg.R`, log `dev/arstep-out/testdbg.txt`).
2. The review's row reaches the maximum: the default fit's log-likelihood
   is at least that of a fit from `start = list(theta = c(log(1.5),
   log(0.1)))`, to 1e-6 relative.
3. The gate reads only pure AR with `cov = FALSE`: `ma()`, `arma()`,
   `ar(cov = TRUE)` and `weights()` beside `ar()` get no marginal, and a
   model with no residual term keeps the marginal it had.
4. A marginal the objective does not agree with is not used (mock: the
   marginal's value plus 1).
5. Off the gaussian, the second maximum is a start: the review's
   bernoulli row (n 40, base 1.07 short) reaches the reference fit, and
   the objective's state is left at the optimum the fit reports.

Seen to fail on the base (`dev/arstep-out/t-gp-ar-start-base.txt`,
`dev/arstep-run-tests.R` on rellib-r18):

```
DETAIL FAIL [the AR marginal reads the model's own density] Expected `is.null(s)` to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL ERROR [the AR marginal reads the model's own density] Error in `if (s$has_sigma) exp(2 * p$betad) else 0`: argument is of length zero
DETAIL FAIL [a step beside ar() reaches the maximum, not the field's basin] Expected `ll` >= `lr - 1e-06 * abs(lr)`. | Actual comparison: 130.3 < 156.2 | Difference: -25.9 < 0
DETAIL FAIL [a step beside ar() reaches the maximum, not the field's basin] Expected `is.null(fit$obj$env$.frmtmb_gp_marginal)` to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [the marginal reads only the AR forms it states] Expected `is.null(setup(y ~ gp(x) + ar(time = t)))` to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [an AR marginal the objective does not agree with is not used] Expected `is.null(honest$obj$env$.frmtmb_gp_marginal)` to be FALSE. | Differences: | `actual`:   TRUE  | `expected`: FALSE | 
DETAIL FAIL [off the gaussian, the working marginal's second maximum is a start] Expected `length(w$alts)` > 0L. | Actual comparison: 0 <= 0 | Difference: 0 <= 0
DETAIL FAIL [off the gaussian, the working marginal's second maximum is a start] Expected `ll` >= `lr - 1e-06 * abs(lr)`. | Actual comparison: -27.2 < -26.1 | Difference: -1.1 < 0
RESULT test-gp-ar-start.R pass=12 fail=7 err=1 skip=0 warn=0
```

The behavioral failures are `130.3 < 156.2` (the 25.9 row) and `-27.2 <
-26.1` (the 1.07 bernoulli row). On the lane
(`dev/arstep-out/t-gp-ar-start-lane.txt`):

```
RESULT test-gp-ar-start.R pass=39 fail=0 err=0 skip=0 warn=0
```

(That run predates dropping one assertion, `fit$opt$convergence == 0`
on the 25.9 row, which a platform could legitimately fail at the
maximum; the count there was 39 and is now 38 per run; the final count
is in section 12.)

Namespace mutants, `dev/arstep-mut.R lane <mutant>`, logs
`dev/arstep-out/mut-*.txt`:

```
MUT lane control pass 39 fail 0 error 0
MUT lane noagree pass 38 fail 1 error 0
MUTF lane noagree an AR marginal the objective does not agree with is not used
MUT lane noalt pass 37 fail 2 error 0
MUTF lane noalt off the gaussian, the working marginal's second maximum is a start
MUT lane noarphi pass 39 fail 0 error 0
MUT lane nogate pass 15 fail 4 error 1
MUTF lane nogate the AR marginal reads the model's own density
MUTF lane nogate a step beside ar() reaches the maximum, not the field's basin
MUTF lane nogate the marginal reads only the AR forms it states
MUTF lane nogate an AR marginal the objective does not agree with is not used
```

`noagree` (the agreement guard always true) fails test 4; `noalt` (no
second maximum) fails test 5; `nogate` (the 0.79.1 route for `ar()`)
fails tests 1 to 4. `noarphi` (the scan's AR grid at 0 only) passed the
first version of the file, so test 2 now also fits `step 80 22`, the
row where the grid matters most (section 6). Rerun on the final file
(`t-gp-ar-start-lane2.txt`, `t-gp-ar-start-base2.txt`,
`mut2-control.txt`, `mut2-noarphi.txt`):

```
RESULT test-gp-ar-start.R pass=39 fail=0 err=0 skip=0 warn=0   (lane)
RESULT test-gp-ar-start.R pass=11 fail=8 err=1 skip=0 warn=0   (base)
MUT lane control pass 39 fail 0 error 0
MUT lane noarphi pass 38 fail 1 error 0
MUTF lane noarphi a step beside ar() reaches the maximum, not the field's basin
```

On the base the new row fails behaviorally too: `7.3 < 28.8`.

## 13. What I did not do, and why

- The student step rows (11, 6.07, 0.20 short) are not fixed. Section 4:
  the objective depends on the inner start there, so no outer start
  settles it. A fix would be an inner solve that searches its own
  modes, which is a change to the Laplace step for every non-log-concave
  family, not to the start. Named in NEWS under Limits.
- `pois_mid two 80` (-1.59) and `gsigma step 40` (-1.45) stay short, as
  on the base. Neither descent reaches the reference there.
- No warning class for disagreeing starts (section 5: 44 false alarms
  for 2 short fits on the two tables).
- No threshold on the marginal gap to skip hopeless second starts
  (section 9: it would be tuned on the same rows it is judged on).
- The AR marginal covers pure `ar()` with `cov = FALSE` and no `se()` or
  `weights()`. `ma()`, `arma()` (the error recursion makes `A`
  nonlinear in the coefficients), `ar(cov = TRUE)` (a stationary
  covariance, another closed form) and `ar()` beside `se()`/`weights()`
  keep the working start and the second descent.
- No comparison with 0.78.0: this round has no 0.78.0 library
  (`rellib-r16` is gone); the brief's before is `rellib-r18`.
- Sessions three ways ran for core and frmtmb.sample only. By
  construction no other extension reaches the change: `gp_start_scan()`
  returns at its first line on a frame with no exact `gp()`, nothing
  else sets `.frmtmb_gp_alt`, so `gp_second_descent()` returns its input
  unchanged, and the gate is read only from `gp_start_scan()`. A grep of
  the extensions' `R/` and tests for `gp(` and `"gp"` finds exact
  `gp()` fits only in frmtmb.sample's tests; frmtmb.spatial's
  `helper-brms.R` names the covstruct but no spatial test uses it. Every
  extension's files still ran one per process (section 12).
- The ported brms-suite ledger was not regenerated (consolidation does
  that). Rows it can flip: those with an exact `gp()` off the gaussian
  or beside `ar()`.

## 14. For the user and the consolidator

- Merge overlap with lane fitedge in `R/gp-start.R`. fitedge adds a
  `max_rows` argument to `gp_marginal_setup()` and `gp_marginal_data()`
  and restructures the loop of `gp_start_working()` (a branch for
  predictors other than `mu`). This lane left both signatures alone (it
  reads the AR block from the frame inside `gp_marginal_data()`), but it
  edits `gp_start_working()`'s loop end and return line: keep fitedge's
  structure and add this lane's `alts` lines after a successful `m` on
  the `mu` branch only, `alts <- list()` before the loop and `alts =
  alts` in the returned list. fitedge's `gp_marginal_outer_gr()` refuses
  any outer vector with `thetaac`, so it does not meet the AR marginal;
  if it is ever widened, `gp_marginal_nll()` already reads `thetaac`.
- `R/fit.R`: one call after the first `optimize_obj()` in
  `fit_assembled()`, and the `start` paragraph of `?frm`. Lanes bigx,
  loglik, nestquad, nlstart and superassign also edit `R/fit.R`; I did
  not check their hunks against these two.
- At the maximum on the step rows beside `ar()` the coefficient is 1.008
  to 1.045, past a unit root, where the data were made with 0.4. That is
  the likelihood of brms's `cov = FALSE` form, which is unconstrained in
  brms 2.23.0 too (`vector[Kar] ar;`, comment in `R/autocor.R`), so it is
  parity. Whether frmtmb should say anything when it fits a
  non-stationary residual regression is a question for the user.
- `DESCRIPTION` changed only in `Config/roxygen2/version` (8.1.0 to
  8.1.1), because roxygenise ran with the user library's roxygen2
  8.1.1.

## 15. OpenBLAS 0.3.26 with LAPACK

`bash dev/ciharden-openblas.sh 0.3.26 lapack` (log
`dev/arstep-out/openblas-build.txt`: Rblas 74 symbols and Rlapack 616 to
OpenBLAS). That copy rounds differently from the reference BLAS
(`dev/arstep-out/ob-blascheck.txt`: `sum(chol2inv(chol(x)))` reads
499.9976969797616 against 499.99769697981083). The three gp test files
under it, on the lane library (`dev/arstep-out/ob-*.txt`):

```
RESULT test-gp-ar-start.R pass=39 fail=0 err=0 skip=0 warn=0
RESULT test-gp-start.R pass=220 fail=0 err=0 skip=0 warn=0
RESULT test-gp-noncentered.R pass=76 fail=0 err=0 skip=0 warn=0
```

0.3.32 was not run: no test of this lane sits on a rounding-decided
branch (each compares the fit with a reference fit and accepts any
point at or above it).

## 12. Test suites

Every file of core and of all eight extensions, one file per R process
(`dev/arstep-suite.sh`, which runs `dev/arstep-run-par.sh`, a copy of
`dev/release/run-par.sh` on rellib-r18), with `NOT_CRAN=true` and
`FRMTMB_BRMS_FIT_TESTS=true`, the lane library ahead of rellib-r18. Every
log's `lib:` line names `wt-arstep-lib` for frmtmb. Logs
`dev/arstep-out/core.log`, `ext.log` and the per-file `core-files/`,
`ext-files/`; tallies by awk over the logs:

```
frmtmb files 309 pass 25600 fail 0 err 0 skip 15 warn 0
frmtmb.coupling files 11 pass 551 fail 0 err 0 skip 5 warn 0
frmtmb.eam files 33 pass 1854 fail 0 err 0 skip 3 warn 0
frmtmb.latent files 12 pass 683 fail 0 err 0 skip 2 warn 0
frmtmb.learn files 16 pass 545 fail 0 err 0 skip 2 warn 0
frmtmb.ode files 11 pass 549 fail 0 err 0 skip 1 warn 0
frmtmb.sample files 85 pass 4326 fail 0 err 0 skip 1 warn 0
frmtmb.spatial files 8 pass 561 fail 0 err 0 skip 0 warn 0
frmtmb.survival files 15 pass 552 fail 0 err 0 skip 0 warn 0
```

Skips: core 15 (`test-drmtmb-agreement.R` 13, gated by
`FRMTMB_DRMTMB_FIT_TESTS`; `test-curve-scale.R` 1, `FRMTMB_SCALE_TESTS`;
`test-fuzz.R` 1, `FRMTMB_FUZZ`), extensions 14, all in their
`test-scale.R` (`FRMTMB_SCALE_TESTS`). None reaches an exact `gp()`.
`test-gp-ar-start.R` ran in its final form there: pass 39.

## 16. Check and sessions (started at the coordinator's request)

Started, not finished when this lane was made to hand back:
`bash dev/arstep-check.sh .` (output `dev/arstep-check/frmtmb/check.txt`,
Status line at its end) and `dev/arstep-run-sessions.sh sess 2 frmtmb
frmtmb.sample` with `NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS` and
`FRMTMB_RELLIB=wt-arstep-lib` (output `dev/arstep-out/sess.log`, one
SESSION line per session plus DETAIL and LEAK lines;
per-session logs in `dev/arstep-out/sess-sessions/`). Three R processes
at most. Results to be read from those files.

R CMD check --as-cran of core, `bash dev/arstep-check.sh .` (output
`dev/arstep-check/frmtmb/check.txt`): `Status: 1 NOTE`, the expected
V8 math-rendering NOTE on the HTML manual ("Skipping checking math
rendering: package 'V8' unavailable"). Tests: `Running 'testthat.R'
[41m] OK`; testthat inside it: `[ FAIL 0 | WARN 0 | SKIP 270 | PASS 22303 ]`
(the brms-fit gate is not set under check, so its files skip).

Sessions of the first round: the core forward and reverse sessions had
not finished (no SESSION line) when the fix round of the review
(`dev/reviews/2026-11-06-arstep.md`) began; I stopped them, since any
code change makes them stale. Their partial logs showed only LEAK lines
tagged `[with load of frmtmb.sample]` or `[with load of a Stan model
(rstan)]`. Sessions are rerun at the end of the fix round (section 17).

## 17. Fix round after the review (dev/reviews/2026-11-06-arstep.md)

### B1: the AR marginal now factors positions, not rows

The first build set every row as its own position under `ar()`
(`J <- seq_len(n)`) and factored `n x n` in each of 220 cells and 13
tape descents. Now `gp_marginal_data()` keeps the positions, and the AR
term enters through `A`:

- With `W = A P L`, `L L' = V` (the field's covariance on the `u`
  positions, `P` the incidence of rows on positions), `A (y - X beta)` is
  normal with covariance `W W' + r I`. With `R' R = G = P' A' A P`,
  `K = R'^-1 P' A' A D` (`D = (y, X)`), `Wd = D' A' A D - K' K` and `M = r I
  + R V R'`: the quadratic forms are `Wd / r + K' M^-1 K` and the log
  determinant `(n - u) log r + log det M`. Only `u x u` is factored.
- `G`, `P' A' A D` and `D' A' A D` are quadratic in the AR coefficients,
  so their lag pieces are summed once (`gp_ar_data()`) and combined at
  each coefficient (`gp_ar_proj()`). Where `G` does not factor (far past
  a unit root, where the columns of `A` shrink like `ar^-T`), a QR of `A
  P` on the rows takes over (`gp_ar_proj_qr()`).
- Where every row is its own position (`direct`), the marginal works on
  the rows with `A V A' + r I`, built by gathers, as the size is the same.
- The AR marginal is no longer taped: its descents run on numbers with
  `beta` and the common scale in closed form (`gp_ar_descend()`), from the
  same starting cells. `gp_marginal_nll()` reads it on numbers.

Agreement with the objective, `dev/arstep-armarg.R lane 0` (log
`dev/arstep-out/armarg-fix3.txt`, 20 random points per row, seed 9100):
the maximum relative difference is 1.4e-14 to 6.1e-12 on 8 rows. On
`step 160 2` one point read `Inf` before the QR fallback
(`dev/arstep-infdbg.R`); that row has 159 positions for 160 rows, so it
takes the position route. The fits reach the same log-likelihoods as the
first build on all 9 rows.

Cost, the reviewer's `dev/review-arstep-cost.R` (one process, lane against
the `noar` mutant that sends `ar()` down the base's route, interleaved),
logs `dev/arstep-out/rcost-fix3-panels.txt`, `rcost-fix3-a500.txt`:

```
COST p600 n 600 rep 1 lane     wall 1.18 ll 32.589404 conv 0 w 0 marg TRUE second - Rmax_MB 318
COST p600 n 600 rep 1 noar     wall 0.25 ll 32.589404 conv 0 w 0 marg FALSE second - Rmax_MB 249
COST p600 n 600 rep 2 noar     wall 0.33 ll 32.589404 conv 0 w 0 marg FALSE second - Rmax_MB 245
COST p600 n 600 rep 2 lane     wall 0.72 ll 32.589404 conv 0 w 0 marg TRUE second - Rmax_MB 232
COST p1000 n 1000 rep 1 lane     wall 0.52 ll 173.765935 conv 0 w 0 marg TRUE second - Rmax_MB 239
COST p1000 n 1000 rep 1 noar     wall 0.35 ll 173.765935 conv 0 w 0 marg FALSE second - Rmax_MB 244
COST p1000 n 1000 rep 2 noar     wall 0.34 ll 173.765935 conv 0 w 0 marg FALSE second - Rmax_MB 244
COST p1000 n 1000 rep 2 lane     wall 0.58 ll 173.765935 conv 0 w 0 marg TRUE second - Rmax_MB 239
COST a300 n 300 rep 1 lane     wall 20.64 ll 227.607924 conv 0 w 0 marg TRUE second - Rmax_MB 227
COST a300 n 300 rep 1 noar     wall 17.37 ll 227.607924 conv 0 w 0 marg FALSE second - Rmax_MB 231
COST a300 n 300 rep 2 noar     wall 16.59 ll 227.607924 conv 0 w 0 marg FALSE second - Rmax_MB 232
COST a300 n 300 rep 2 lane     wall 21.89 ll 227.607924 conv 0 w 0 marg TRUE second - Rmax_MB 239
COST pois240 n 240 rep 1 lane     wall 14.39 ll -459.715793 conv 0 w 0 marg FALSE second 1,0 Rmax_MB 238
COST pois240 n 240 rep 1 nosecond wall 8.44 ll -459.715793 conv 0 w 0 marg FALSE second - Rmax_MB 238
COST pois240 n 240 rep 2 nosecond wall 8.31 ll -459.715793 conv 0 w 0 marg FALSE second - Rmax_MB 238
COST pois240 n 240 rep 2 lane     wall 14.60 ll -459.715793 conv 0 w 0 marg FALSE second 1,0 Rmax_MB 238
COST bern240 n 240 rep 1 lane     wall 14.47 ll -134.893263 conv 0 w 0 marg FALSE second - Rmax_MB 238
COST bern240 n 240 rep 1 nosecond wall 15.02 ll -134.893263 conv 0 w 0 marg FALSE second - Rmax_MB 238
COST bern240 n 240 rep 2 nosecond wall 14.08 ll -134.893263 conv 0 w 0 marg FALSE second - Rmax_MB 238
COST bern240 n 240 rep 2 lane     wall 14.37 ll -134.893263 conv 0 w 0 marg FALSE second - Rmax_MB 238
COST a500 n 500 rep 1 lane     wall 91.36 ll 269.968071 conv 0 w 0 marg TRUE second - Rmax_MB 322
COST a500 n 500 rep 1 noar     wall 175.31 ll 262.192289 conv 0 w 0 marg FALSE second - Rmax_MB 279
```

| design | lane, s | base route, s | ratio | first build, s |
|---|---|---|---|---|
| p1000: 25 series by 40 times, 40 positions | 0.52, 0.58 | 0.35, 0.34 | 1.5 to 1.7 | 180.4, 175.7 |
| p600: 12 series by 50 times, 50 positions | 1.18, 0.72 | 0.25, 0.33 | 2.2 to 4.7 | 47.2, 45.4 |
| a300: 300 distinct positions | 20.64, 21.89 | 17.37, 16.59 | 1.19 to 1.32 | 22.7, 24.8 |
| a500: 500 distinct positions | 91.36 | 175.31 | 0.52 | 117.8 |

The lane's log-likelihood equals the base route's on p600, p1000 and
a300, and is 7.78 higher on a500 (269.968071 against 262.192289). On the
panels the remaining overhead is about 0.2 to 0.4 s, in the 220-cell scan
and the descents' R-level calls (`dev/arstep-prof.R p600`, log
`prof-p600.txt`: the marginal took 0.27 s of a 0.52 s fit). The first
p600 rep includes the process's first fit. p600 stays above the 1.5
times asked for; I did not find a cut that keeps every row (below).

Fewer starts were tried (`dev/arstep-trim.R`, log `trim.txt`, one process,
36 fits): starting from 2 length scales instead of 6 saved 104.8 s to
70.4 s over the rows, but read `mix 80 2` 0.107 lower (81.862259 against
81.969413) and did not speed the panels (0.64 against 0.66 s). Not taken.

### M1 and open question 3: a constant sigma

`dev/arstep-sigc.R base`, log `dev/arstep-out/sigc-base.txt`: on
`bf(y ~ gp(x), sigma = 0.3)` (seed 12, n 50) the base's marginal at the
fit's optimum read 1606500084705.9 against the objective's 18.778542,
because it read no noise. Its maximum (48.014998) was the start and the
stored `.frmtmb_gp_marginal`. The fit still reached the best of 6 starts
on both sigmas tried, and the proof did not quiet a stop there (the
objective was 18.78 and 48.09 against 48.01), so this was a wrong start
and a wrong proof value that happened to do no harm on these rows.

Fix: `gp_marginal_gate()` reads the constant (`sigc`) and
`gp_marginal_setup()` adds `sigma^2 / w` to the known noise, the route
`se()` already takes. Beside `ar()` the gate refuses a constant sigma,
so no scan runs and nothing warns. The marginal's own `nlminb()` calls
are also wrapped in `suppressWarnings()`: their NaN steps name nothing
the caller wrote.

### M2: the settle test

The reviewer's `dev/review-arstep-settle-test.R` is now the test "a second
descent that loses leaves the modes at the kept optimum".

### Minors

- m0: NEWS gives the per-design cost and the worst single second
  descent (about 2.2 times, the reviewer's rows 28 and 35 of part N, and
  1.73 times on `pois240` here).
- m1: `gp_second_descent()` sets `evals` and `nonfinite_trials` from the
  shared tally, so a lost second descent counts.
- m2: `gp_start_working()`'s `@return` names `alts`.
- m3: every line of the test file is 80 columns or less.
- m5, m6: NEWS bullet 1 is split; one bullet says the unit root is the
  usual answer on step-shaped data and that it is brms's likelihood.
- m4: the `on.exit()` restore in `gp_start_scan()` stays and is not
  pinned (the reviewer measured it harmless on 31 rows).

### Merge notes (for the consolidator)

- `R/fit.R`, lane nestquad: both lanes insert after the `}))` that closes
  the first `optimize_obj()` call in `fit_assembled()` (base line 1644).
  Keep this lane's two lines,
  `opt <- gp_second_descent(obj, opt, ctl_opt, bounds, par_units,` /
  `verbose = vb, tally = tally)`, inside the first branch, right after
  that `}))`, and before nestquad's `} else if (!is.null(qnest)) {`.
- `R/gp-start.R`, lane fitedge, `gp_start_working()`: take fitedge's
  loop. On its `mu` branch, after `m` is accepted and before
  `idx <- c(idx, from)`, add this lane's block:
  `at <- length(theta) + seq_along(from)` and
  `alts <- c(lapply(alts, function(th) c(th, m$theta)), lapply(m$alts,
  function(a) { th <- c(theta, m$theta); th[at] <- a; th }))`.
  On fitedge's non-`mu` branch, before its `idx <- c(idx, from)`, add
  `alts <- lapply(alts, function(th) c(th, m$theta))`, so every
  alternative keeps the length of `theta`. Keep `alts <- list()` before
  the loop and `alts = alts` in the returned list.
- `R/gp-start.R`, fitedge's `max_rows` argument to `gp_marginal_setup()`
  and `gp_marginal_data()`: this lane does not change either signature.
  Its hunks in `gp_marginal_data()` are the deleted `is.null(ar)` test
  (now the base's line again), the `ar`/`ard` lines after `u`, and the
  `ar = ard` element of the returned list. fitedge's
  `gp_marginal_outer_gr()` refuses `thetaac`, so it does not meet the AR
  marginal, which is not taped any more.

### No row worse: the reviewer's 168 fresh rows on the fixed build

`dev/review-arstep-fresh.R lane <part> 0` (the reviewer's script, its seeds), logs
`dev/arstep-out/fresh-{G,N,P}-fix.txt`; compared by `dev/arstep-fresh-cmp.R`
against the reviewer's base logs and references (output `fresh-cmp.txt`):

```
FC G rows 45 | base short 9 worst -32.6 | fixed short 0 silent 0 worst 0 | fixed worse than base 0, better 9 | wall base 266.2 fixed 415.3 (not side by side)
FC G fixed against the first lane build: lower by >1e-3 on 0, higher on 0
FC N rows 105 | base short 17 worst -6.18 | fixed short 11 silent 9 worst -5.1 | fixed worse than base 0, better 6 | wall base 325.6 fixed 375.3 (not side by side)
FC N fixed against the first lane build: lower by >1e-3 on 0, higher on 0
FC P rows 18 | base short 2 worst -0.0407 | fixed short 0 silent 0 worst 0 | fixed worse than base 0, better 2 | wall base 22.7 fixed 28.9 (not side by side)
FC P fixed against the first lane build: lower by >1e-3 on 0, higher on 0
FCSUM rows 168 base short 28 fixed short 11 fixed worse than base 0 better 17
```

Same counts as the review found on the first build, and no row moved by
more than 1e-3 from it. The wall sums were not taken side by side with
the base (three of my processes ran at once on a loaded machine), so
they are not a cost measurement; the cost table above is.

### The lane's sets on the fixed build

Default fit only, `ngt-lane-fix.txt` (`dev/arstep-ng.R`), `gauss-lane-fix.txt`
(`dev/arstep-gauss2.R`), `r2ng-lane-fix.txt` (`dev/arstep-r2ng.R`), against the
base logs and references of sections 6, 7 and 10; `dev/arstep-summary.R`,
output `summary-fix.txt` (walls not side by side, so not a cost):

```
SUM ng all                     base       rows 150 short  13 silent  12 worst     -25.9 | wall on 150 common rows    450.5 s
SUM ng all                     lane       rows 150 short   5 silent   4 worst       -11 | wall on 150 common rows    426.3 s
SUM ng all                     lane       against base: worse 0, better 8
SUM gauss 192                  base       rows 192 short   0 silent   0 worst         0 | wall on 192 common rows    226.6 s
SUM gauss 192                  lane       rows 192 short   0 silent   0 worst         0 | wall on 192 common rows    172.3 s
SUM gauss 192                  lane       against base: worse 0, better 0
SUM r2ng 72                    base       rows  72 short   5 silent   5 worst     -1.77 | wall on 72 common rows    314.6 s
SUM r2ng 72                    lane       rows  72 short   4 silent   4 worst     -1.77 | wall on 72 common rows    153.7 s
SUM r2ng 72                    lane       against base: worse 0, better 2
SUM r2ng the 48                base       rows  48 short   0 silent   0 worst         0 | wall on 48 common rows    221.1 s
SUM r2ng the 48                lane       rows  48 short   0 silent   0 worst         0 | wall on 48 common rows    107.0 s
SUM r2ng the 48                lane       against base: worse 0, better 0
```

The 150 rows read the same as the first build to 1e-3 on every row,
and the 192 gaussian rows the same as the base to the 6 printed
decimals on every row (joins by row, awk).

The first version of the constant-sigma fix broke `test-gp-start.R`'s
"the marginal reads the tape's objective at one point" on `y | se(s) ~
gp(x)` (6.99 apart; `dev/arstep-out/tfix-test-gp-start.R.txt`, an earlier
run of it): there the frame carries a constant sigma, but `se()` alone
replaces sigma (`resid_sd()` in `R/families.R`). The gate now keeps the
constant only without `se()` or with `se(sigma = TRUE)`, where the two
add. Final runs on the final build (`dev/arstep-out/tfix-*.txt`):

```
RESULT test-gp-ar-start.R pass=49 fail=0 err=0 skip=0 warn=0
RESULT test-gp-start.R pass=220 fail=0 err=0 skip=0 warn=0
RESULT test-gp-noncentered.R pass=76 fail=0 err=0 skip=0 warn=0
```

Mutants on `test-gp-ar-start.R` (`dev/arstep-mut.R`, logs
`dev/arstep-out/mut3-*.txt`, on the build before the `se()` change, which
does not reach these tests): control pass 49; `nosettle` (the reviewer's)
fails the settle test; `nosigc` (a constant sigma read as no noise)
fails the constant-sigma test; `rows` (every row its own position, the
first build's B1) fails the panel test. On rellib-r18 the file reads
pass 18, fail 11, error 1 (`t-gp-ar-start-base3.txt`), with the
behavioral failures `130.3 < 156.2`, `7.3 < 28.8`, `-27.2 < -26.1`, the
constant-sigma marginal 1.6e12 from the objective, and `s$u` NULL on
the panel.

### Still running at hand-back

- Core suite on the final build, one file per process, 3 at a time
  (`dev/arstep-suite.sh core2`): per-file logs
  `dev/arstep-out/core2-files/`, totals `dev/arstep-out/core2.log`
  ("CORE2 ran X of 309").
- Not started: `R CMD check --as-cran` (`bash dev/arstep-check.sh .`,
  output `dev/arstep-check/frmtmb/check.txt`) and the sessions
  (`dev/arstep-run-sessions.sh sess 2 frmtmb frmtmb.sample` with
  `NOT_CRAN`, `FRMTMB_BRMS_FIT_TESTS` and
  `FRMTMB_RELLIB=wt-arstep-lib`; output `dev/arstep-out/sess.log`). The
  first round's check read `Status: 1 NOTE` (V8) on the first build; it
  is stale.

Core suite on the final build finished after hand-back
(`dev/arstep-out/core2.log`, every `lib:` line names wt-arstep-lib):
`frmtmb files 309 pass 25610 fail 0 err 0 skip 15 warn 0` (the same 15
gated skips as section 12).
