# Lane nlband: a fit clear of a band of poles

Round of 2026-11-16, base b255350e (frmtmb 0.83.0). Worktree
`frmtmb-wt-nlband`, branch `wt-r84-nlband`, library `wt-nlband-lib`,
before arm `rellib-r22`. Brief: `round-briefs/r084-nlband.md`, items 5115
(HIGH) and 5099 (LOW).

Status: see the end (gated tier and check were running at hand-back).

Scripts are `dev/nlband-*.R` and `-*.sh`; logs are under
`dev/nlband-out/` (gitignored, local). The prior round's logs this lane
reads were copied from `frmtmb-wt-nlstop2` into `dev/nlband-prior/`
before that tree is deleted.

## What changed

- `R/fit-restart.R`: new `nl_pole_holes()` and its four sizes
  (`nl_hole_window`, `nl_hole_side`, `nl_hole_odds`, `nl_hole_max`);
  `nl_pole_walls()` returns `holes`; `nl_pole_escape()` adds the hole
  starts (two per hole), reads a hole run's end by distance only, and
  returns `where`; new `nl_pole_from()` and `nl_pole_where()` give the
  message and the warning their text. The 0.83.0 past-the-band starts,
  the sign probe, the saddle restart and `nl_starts = FALSE` are as they
  were.
- `R/fit.R`: `@param nl_starts` text only; `man/frmtmb_control.Rd`
  regenerated.
- `tests/testthat/test-fit-restart.R`: header line and four tests.
- `NEWS.md`: a `# frmtmb (development version)` section, one BREAKING
  bullet under `## Fitting`.
- `dev/test-backlog.md`: CLOSED lines under 5099 and 5115; section
  "Filed in the round of 2026-11-16", lane nlband, three LOW items.
- `dev/nlband-*` scripts and this file.

No extension file changed; no compat rule or brmsport row changed.

## 1. Reproduction on rellib-r22

`dev/nlband-repro.R` (`NLBAND_ARM=before`, log `repro-before.txt`):
`a / (x + b)`, y = 1 / (x - 2) + N(0, 0.2), x ~ U(0, 4) with 400 draws
and the rows in |x - 2| < 0.3 removed; the default fit, and the fit from
the review's start (1, -1.9).

```
ROW seed 1 n 334 default -543.816417 code 0 b 0.3865 a -0.5862 warn - msg - t 0.61 | gap-start 47.731985 code 0 b -2.0038 warn - msg - | short 591.548
ROW seed 2 n 348 default -525.219128 code 0 b 0.4811 a -0.4877 warn - msg - t 0.03 | gap-start 66.606612 code 0 b -2.0112 warn - msg - | short 591.826
ROW seed 3 n 355 default -589.073619 code 0 b 0.3709 a -0.4370 warn - msg - t 0.05 | gap-start 57.327072 code 0 b -2.0022 warn - msg - | short 646.401
```

The three values match the review's `pole-before.txt` to the printed
digits. 20 fresh seeds (101 to 120): every default stops at code 0 with
no warning and no message, b 0.27 to 0.61, 597.4 to 685.7 below the fit
from (1, -1.9) at b -1.9946 to -2.0071 (same log).

## 2. The rule (5115)

Cause. The walls of `b` in `a / (x + b)` are -x_i, one per row; they cut
the line of `b` into pieces, and the optimizer cannot cross a wall. The
default fit sits beyond the last wall (b 0.39, every denominator
positive) and is a maximum of its piece. The pole restart of 0.83.0
runs only when a band holds the fit (`nl_pole_pressed()`: denominators
of both signs, or the fit nearer a wall than the median gap) and only
starts past the band's ends, so no start reaches the gap at b -2.

Every piece between two adjacent walls has a maximum of its own, so a
start in every gap is out of reach on dense data (n - 1 starts).
(Corrected in the fix round: the first version said that a pole between two
adjacent rows fits no data. It can be the maximum: on 400 rows of U(0,
4) with no hole, y = 1 / (x - 2) + noise, both builds stop 2005.3,
2100.2 and 1574.7 short in silence; review of lane nlband, item 3,
filed HIGH.) The rule starts only in a hole: a gap between adjacent
walls that is wide against the gaps around it.

`nl_pole_holes(w)` (new, `R/fit-restart.R`), on the sorted walls `w`
(walls within 1e-9 of the range of the one before are merged):

- the local scale of gap i is the running median of the 41 gaps centered
  on it (`stats::runmed(d, 41, endrule = "constant")`, or the median of
  all gaps when there are 41 or fewer);
- gap i is a hole when it has 10 or more walls on each side and is wider
  than log2(1000 m) local medians, m gaps in all. Reason: on rows of a
  smooth density a gap is wider than t local means with probability
  exp(-t), so m gaps hold one wider than log(1000 m) means, which is
  log2(1000 m) medians (the median is ln 2 means), with odds of about 1
  in 1000;
- at most the 2 widest holes, by that ratio, ties by position.

The sizes are chosen by the data, never by a seed, and were picked on the
simulation below before any fit was run.

`nl_pole_escape()` keeps the 0.83.0 starts past the ends of a pressed
band, and adds one start at the middle of each hole the fit is not in,
whether or not the band holds the fit, the rest of the model at the fit.
A hole run, and its refit, is kept only if it ends no nearer a wall than
the walls' median gap (the "past" runs keep their 0.83.0 test, which
also refuses a mixed-sign end; a hole end is mixed-sign by design). The
gain rule (`nl_probe_gain()`), the refit through `run()`, the message
class `frmtmb_nl_pole_start` from the default start, the warning class
`frmtmb_nl_pole` from a given start, and `nl_starts = FALSE` are those of
0.83.0. The message and the warning name the coefficient and the hole's
walls when a hole start won.

### The rule's sizes, measured before use

`dev/nlband-holes.R` (log `holes.txt`): a fixed ratio against a running
median of 11 or 21 gaps fires on most draws at large n (the largest of n
exponential gaps grows as log n) and on most heavy-tailed draws (the
gaps of a tail grow). `dev/nlband-holes2.R` (log `holes2.txt`) then
compares four forms; seed 20261116, 400 draws per design for n 20 to
1000, 100 at 1e4, 20 at 1e5. The chosen form, window 41 (k 20), 10 walls
a side, odds 1000:
    NOHOLE k 20 side 10 a 1000 unif    n     20 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 unif    n     50 reps 400 | any 0.002 mean 0.00 max 1
    NOHOLE k 20 side 10 a 1000 unif    n    100 reps 400 | any 0.002 mean 0.00 max 1
    NOHOLE k 20 side 10 a 1000 unif    n    300 reps 400 | any 0.005 mean 0.00 max 1
    NOHOLE k 20 side 10 a 1000 unif    n   1000 reps 400 | any 0.015 mean 0.01 max 1
    NOHOLE k 20 side 10 a 1000 unif    n  10000 reps 100 | any 0.050 mean 0.05 max 1
    NOHOLE k 20 side 10 a 1000 unif    n 100000 reps  20 | any 0.050 mean 0.05 max 1
    NOHOLE k 20 side 10 a 1000 norm    n     20 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 norm    n     50 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 norm    n    100 reps 400 | any 0.007 mean 0.01 max 1
    NOHOLE k 20 side 10 a 1000 norm    n    300 reps 400 | any 0.010 mean 0.01 max 1
    NOHOLE k 20 side 10 a 1000 norm    n   1000 reps 400 | any 0.020 mean 0.02 max 2
    NOHOLE k 20 side 10 a 1000 norm    n  10000 reps 100 | any 0.050 mean 0.05 max 1
    NOHOLE k 20 side 10 a 1000 norm    n 100000 reps  20 | any 0.050 mean 0.05 max 1
    NOHOLE k 20 side 10 a 1000 lnorm   n     20 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 lnorm   n     50 reps 400 | any 0.020 mean 0.02 max 1
    NOHOLE k 20 side 10 a 1000 lnorm   n    100 reps 400 | any 0.007 mean 0.01 max 1
    NOHOLE k 20 side 10 a 1000 lnorm   n    300 reps 400 | any 0.020 mean 0.02 max 1
    NOHOLE k 20 side 10 a 1000 lnorm   n   1000 reps 400 | any 0.032 mean 0.04 max 2
    NOHOLE k 20 side 10 a 1000 lnorm   n  10000 reps 100 | any 0.070 mean 0.07 max 1
    NOHOLE k 20 side 10 a 1000 lnorm   n 100000 reps  20 | any 0.050 mean 0.05 max 1
    NOHOLE k 20 side 10 a 1000 exp     n     20 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 exp     n     50 reps 400 | any 0.010 mean 0.01 max 1
    NOHOLE k 20 side 10 a 1000 exp     n    100 reps 400 | any 0.020 mean 0.02 max 1
    NOHOLE k 20 side 10 a 1000 exp     n    300 reps 400 | any 0.010 mean 0.01 max 1
    NOHOLE k 20 side 10 a 1000 exp     n   1000 reps 400 | any 0.022 mean 0.02 max 1
    NOHOLE k 20 side 10 a 1000 exp     n  10000 reps 100 | any 0.040 mean 0.04 max 1
    NOHOLE k 20 side 10 a 1000 exp     n 100000 reps  20 | any 0.100 mean 0.10 max 1
    NOHOLE k 20 side 10 a 1000 t3      n     20 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 t3      n     50 reps 400 | any 0.002 mean 0.00 max 1
    NOHOLE k 20 side 10 a 1000 t3      n    100 reps 400 | any 0.010 mean 0.01 max 1
    NOHOLE k 20 side 10 a 1000 t3      n    300 reps 400 | any 0.020 mean 0.02 max 1
    NOHOLE k 20 side 10 a 1000 t3      n   1000 reps 400 | any 0.015 mean 0.01 max 1
    NOHOLE k 20 side 10 a 1000 t3      n  10000 reps 100 | any 0.080 mean 0.08 max 1
    NOHOLE k 20 side 10 a 1000 t3      n 100000 reps  20 | any 0.050 mean 0.10 max 2
    NOHOLE k 20 side 10 a 1000 unif_r3 n     20 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 unif_r3 n     50 reps 400 | any 0.002 mean 0.00 max 1
    NOHOLE k 20 side 10 a 1000 unif_r3 n    100 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 unif_r3 n    300 reps 400 | any 0.005 mean 0.00 max 1
    NOHOLE k 20 side 10 a 1000 unif_r3 n   1000 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 unif_r3 n  10000 reps 100 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 unif_r3 n 100000 reps  20 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 norm_r3 n     20 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 norm_r3 n     50 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 norm_r3 n    100 reps 400 | any 0.007 mean 0.01 max 1
    NOHOLE k 20 side 10 a 1000 norm_r3 n    300 reps 400 | any 0.005 mean 0.00 max 1
    NOHOLE k 20 side 10 a 1000 norm_r3 n   1000 reps 400 | any 0.005 mean 0.00 max 1
    NOHOLE k 20 side 10 a 1000 norm_r3 n  10000 reps 100 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 norm_r3 n 100000 reps  20 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 clust2  n     20 reps 400 | any 0.000 mean 0.00 max 0
    NOHOLE k 20 side 10 a 1000 clust2  n     50 reps 400 | any 1.000 mean 1.00 max 2
    NOHOLE k 20 side 10 a 1000 clust2  n    100 reps 400 | any 1.000 mean 1.00 max 2
    NOHOLE k 20 side 10 a 1000 clust2  n    300 reps 400 | any 1.000 mean 1.01 max 2
    NOHOLE k 20 side 10 a 1000 clust2  n   1000 reps 400 | any 1.000 mean 1.02 max 2
    NOHOLE k 20 side 10 a 1000 clust2  n  10000 reps 100 | any 1.000 mean 1.07 max 2
    NOHOLE k 20 side 10 a 1000 clust2  n 100000 reps  20 | any 1.000 mean 1.05 max 2
    HOLE k 20 side 10 a 1000 h 0.30 n    50 | detected 0.302
    HOLE k 20 side 10 a 1000 h 0.30 n   100 | detected 0.960
    HOLE k 20 side 10 a 1000 h 0.30 n   400 | detected 1.000
    HOLE k 20 side 10 a 1000 h 0.30 n  2000 | detected 0.998
    HOLE k 20 side 10 a 1000 h 0.10 n    50 | detected 0.007
    HOLE k 20 side 10 a 1000 h 0.10 n   100 | detected 0.030
    HOLE k 20 side 10 a 1000 h 0.10 n   400 | detected 0.995
    HOLE k 20 side 10 a 1000 h 0.10 n  2000 | detected 0.993
    HOLE k 20 side 10 a 1000 h 0.05 n    50 | detected 0.002
    HOLE k 20 side 10 a 1000 h 0.05 n   100 | detected 0.002
    HOLE k 20 side 10 a 1000 h 0.05 n   400 | detected 0.323
    HOLE k 20 side 10 a 1000 h 0.05 n  2000 | detected 0.998

"any" is the share of draws with at least one hole, "max" the most holes
on one draw. With no hole in the data, at most 3.2 percent of the
unimodal draws up to n 1000 get a start (lognormal at 1000), at most 2
of 20 at n 1e5 (exponential); clust2 (two clusters) has a real hole and
gets one start on every draw. The review's hole (half width 0.3, n 400)
is found on every draw, a hole of half width 0.1 on 398 of 400 at n 400.
Narrower forms (k 10, or 5 walls a side) start on up to 75 percent of
the draws without a hole (same log, `k 10 side 5`, unif n 1e5).

### The other coefficients of the body at a hole start

The first build started each hole run with the rest of the model at the
fit. `dev/nlband-guard.R` section 5 (`Vm * x / (K + x)` on y = 2 x / (x
- 2) + noise, the band data's x) showed that this fails where the fit is
on a ridge to a limit of the body: seeds 2 and 3 end, after the sign
probe, at Vm 3.0e8 and 3.4e5, logLik -953.61 and -1034.32, and the run
from the hole with Vm still at 3.0e8 ends at -5.5e9. From Vm 0 to 100,
or with Vm 2 from a `start`, every run reaches the hole's maximum, 64.94
and 81.18 (`dev/nlband-dbg2.R`, log `dbg2.txt`). Each hole now gets two
runs: the body's other coefficients at the fit, and at 0, where a
default fit starts them (one run where they are 0 at the fit). In the
fix round the second run became the whole model at the default start
(section 10). That build (`guard-after3.txt`): seeds 2 and 3 reach
64.943684 and 81.178905
with `frmtmb_nl_pole_start`; seed 1 was already there by the sign probe.
The starts past the ends of a band (0.83.0) keep the rest at the fit; I
did not change or measure them.

## 3. Guard-absent and neighbouring cases

`dev/nlband-guard.R`, one process per arm (`guard-before2.txt`,
`guard-after3.txt`), compared by `dev/nlband-cmp.sh`: SAME means
logLik and estimates to 17 digits, standard errors to 12 digits, and the
warning and message classes are identical. HOLES is the number of hole
runs the final build made.

```
SAME  ROW nohole clear 1 | holes 0 
SAME  ROW nohole clear 2 | holes 0 
SAME  ROW nohole clear 3 | holes 0 
SAME  ROW nohole clear 4 | holes 0 
SAME  ROW nohole clear 5 | holes 0 
SAME  ROW hole worse 1 | holes 2 
SAME  ROW hole worse 2 | holes 2 
SAME  ROW hole worse 3 | holes 2 
SAME  ROW hole worse 4 | holes 2 
SAME  ROW hole worse 5 | holes 2 
DIFF  ROW given clear 1 | holes 2 
SAME  ROW given inside 1 | holes 0 
DIFF  ROW given clear 2 | holes 2 
SAME  ROW given inside 2 | holes 0 
DIFF  ROW binomial band 1 | holes 2 
DIFF  ROW bernoulli band 1 | holes 2 
DIFF  ROW binomial band 2 | holes 2 
DIFF  ROW bernoulli band 2 | holes 2 
DIFF  ROW binomial band 3 | holes 2 
DIFF  ROW bernoulli band 3 | holes 2 
DIFF  ROW binomial band 4 | holes 2 
DIFF  ROW bernoulli band 4 | holes 2 
DIFF  ROW binomial band 5 | holes 2 
DIFF  ROW bernoulli band 5 | holes 2 
SAME  ROW Vm*x/(K+x) 1 | holes 0 
DIFF  ROW Vm*x/(K+x) 2 | holes 2 
DIFF  ROW Vm*x/(K+x) 3 | holes 2 
DIFF  ROW two holes 1 | holes 4 
DIFF  ROW two holes 2 | holes 4 
DIFF  ROW two holes 3 | holes 4 
SAME  ROW nl_starts FALSE | holes 0 
SAME  ROW mm clusters 1 | holes 2 
SAME  ROW mm clusters 2 | holes 2 
SAME  ROW mm clusters 3 | holes 2 
```

- No hole (x ~ U(0, 4), all 400 rows; truth b 0.5), seeds 1 to 5: no
  hole is found, no run, the same fit.
- A hole where the inside start is worse (the band data's x, truth b
  0.5), seeds 1 to 5: 2 hole runs each, both lose, the same fit.
- Two clusters of xp for `Vm * xp / (K + xp)` (0.1 to 1.5 and 3.5 to 5,
  a real hole), seeds 1 to 3: 2 runs each, both lose, the same fit.
- `nl_starts = FALSE`: no run, the same fit as 0.83.0.
- A given start clear of the poles, seeds 1 and 2: the same fit, now
  with the warning `frmtmb_nl_pole`, which names the hole:

      The fit from the given start reached log-likelihood -543.8164. The
      body has a pole (a denominator at 0) at one value of b_(Intercept)
      for each row of the data, and the optimizer cannot cross one. A
      start with b_(Intercept) between -2.313 and -1.681, a gap that the
      rows' poles leave, reached 47.73199, so the fit is not at the
      maximum. Give a `start` in that gap

- A given start inside the hole (1, -1.9): no hole run (the fit is in
  it), the same fit, no warning.
- binomial (10 trials, logit p = 1 / (x - 2)) and bernoulli (logit p = 2
  / (x - 2)) on the band data's x, seeds 1 to 5: 0.83.0 stops clear of
  the poles at b 0.03 to 0.49, silent; the final build reaches b -1.879
  to -2.202 inside the hole, code 0, no warning, each end 0.08 or more
  from the hole's walls, so none is a supremum on a wall (logLik moves
  from -1025.30 to -550.70 on binomial seed 1, from -204.49 to -120.81 on
  bernoulli seed 1; every row in the log).
- Two holes (|x - 1| < 0.25 and |x - 3| < 0.25 removed, pole at 3),
  seeds 1 to 3: 4 runs each; the fit moves from -644.16, -676.41,
  -682.12 to 59.53, 78.29, 68.07 at b -2.995 to -3.006.

## 4. Seeds 1 to 3 and 20 fresh seeds, after

`dev/nlband-repro.R` on the final build (`repro-after2.txt`): every one
of the 23 default fits now ends at code 0 with no warning and the
message `frmtmb_nl_pole_start`, at the fit from (1, -1.9): the two logLik
columns, printed to 6 decimals, agree on all 23 (`dev/nlband-ulp.R` below
puts the difference at 6.0e-13 or less):

```
ROW seed 1 n 334 default 47.731985 code 0 b -2.0038 a 0.9982 warn - msg frmtmb_nl_pole_start t 1.41 | gap-start 47.731985 code 0 b -2.0038 warn - msg 
ROW seed 2 n 348 default 66.606612 code 0 b -2.0112 a 0.9966 warn - msg frmtmb_nl_pole_start t 0.15 | gap-start 66.606612 code 0 b -2.0112 warn - msg 
ROW seed 3 n 355 default 57.327072 code 0 b -2.0022 a 0.9961 warn - msg frmtmb_nl_pole_start t 0.11 | gap-start 57.327072 code 0 b -2.0022 warn - msg 
```

The message (band data, seed 1, `guard-after3.txt` shows the same text
for the binomial rows):

    The fit of the nonlinear model reached log-likelihood -543.8164. The
    body has a pole (a denominator at 0) at one value of b_(Intercept)
    for each row of the data, and the optimizer cannot cross one. A
    start with b_(Intercept) between -2.313 and -1.681, a gap that the
    rows' poles leave, reached 47.73199, optimizer return code 0, and
    frm() returns that fit. Give `start` to choose your own

### Rounding

`dev/nlband-ulp.R` (`ulp-after2.txt`, one process, y * (1 + 2e-15 k), k
0 to 23):
- band data seeds 1 to 3: 72 draws; 72 end at code 0 with the message and no
  warning; the largest |fit from (1, -1.9) - default| is 6.04e-13.
- the hole loses (truth b 0.5, seeds 1 to 3): 72 draws; on 72 the fit is
  identical() (logLik, estimates, standard errors, condition classes) to
  the same fit in the same process with nl_pole_holes() returning NULL,
  which is the 0.83.0 escape; 72 have no message and no warning.

## 5. Cost

### The survey designs of lane nlstop2

The fresh rows (`dev/nlstop2-fresh.R`, `nlstop2-fresh2.R`), the scan
rows (`dev/nlstop2-scan.R`) and the limit battery (`dev/nlstop2-limit.R`),
the job list of `nlstop2-out/jobs-acc2.txt`, run through
`dev/nlband-survey.R` on rellib-r22 (`before`) and the final build
(`after2`); summary `dev/nlband-accsum.sh before after2`
(`accsum-after2.txt`). Every ROW line is compared with its wall-time
field dropped.

```
ALL rows 2555 differ 0
COSTSUM before jobs 25 escapes 1985 holes 0 past 11 evals 4098850 in escape 741 secs 881.2
COSTSUM after2 jobs 25 escapes 1985 holes 12 past 11 evals 4099581 in escape 1472 secs 1225.1
```

- 2555 rows, 25 jobs: no row differs (logLik, code, warning and message
  classes, estimates as each script prints them).
- `nl_pole_escape()` is called 1985 times. It makes 12 hole runs, on 6
  fits: the limit battery's mm and emax at n 200 seed 30 and n 2000
  seed 10, and Cprior's poismm at the same two, which share the x draws
  of runif(n, 0, 10) (emax: runif(n, 0, 20)). Each is a gap of uniform
  rows that the rule calls a hole, as the simulation predicts at about 1
  to 2 percent; every run loses and the fit is unchanged.
- Objective and gradient evaluations of every optimize_obj() run:
  4,098,850 before, 4,099,581 after, 731 more (0.018 percent), all
  inside the escape (741 to 1472). The 0.83.0 past-the-band runs are
  the same 11 on both arms.
- The `secs` sums are not comparable: the before arm ran two jobs at a
  time beside 10 other lanes, the after arm one at a time beside this
  lane's own suite.

### What a hole start costs one fit

`dev/nlband-cost.R` (`cost-after2.txt`), one process; "off" is the
same build with nl_pole_holes() returning NULL (0.83.0's escape), "off2"
the control. The wall times do not resolve it (controls 0.951 to
1.321), so the evaluation counts decide:

```
COST hole wins (band data, n 400): min s on 0.41 off 0.25 off2 0.28 | on/off 1.640 control 1.120
COST hole loses (truth b 0.5, n 400): min s on 0.29 off 0.25 off2 0.32 | on/off 1.160 control 1.280
COST no hole (n 400): min s on 0.29 off 0.28 off2 0.37 | on/off 1.036 control 1.321
COST hole loses (truth b 0.5, n 20000): min s on 2.15 off 2.03 off2 1.93 | on/off 1.059 control 0.951
COUNT hole wins (band data, n 400) row 1: evaluations on 398 off 261 ratio 1.525
COUNT hole wins (band data, n 400) row 2: evaluations on 454 off 296 ratio 1.534
COUNT hole wins (band data, n 400) row 3: evaluations on 433 off 286 ratio 1.514
COUNT hole loses (truth b 0.5, n 400) row 1: evaluations on 338 off 251 ratio 1.347
COUNT hole loses (truth b 0.5, n 400) row 2: evaluations on 385 off 289 ratio 1.332
COUNT hole loses (truth b 0.5, n 400) row 3: evaluations on 383 off 282 ratio 1.358
COUNT no hole (n 400) row 1: evaluations on 268 off 268 ratio 1.000
COUNT no hole (n 400) row 2: evaluations on 279 off 279 ratio 1.000
COUNT no hole (n 400) row 3: evaluations on 271 off 271 ratio 1.000
COUNT hole loses (truth b 0.5, n 20000) row 1: evaluations on 597 off 469 ratio 1.273
COUNT hole loses (truth b 0.5, n 20000) row 2: evaluations on 725 off 598 ratio 1.212
```

A fit with no hole pays three evaluations of its denominators and one
running median per coefficient and denominator, and no optimization. A
hole that loses costs 1.21 to 1.36 times the fit's evaluations (two
runs from the hole); one that wins, 1.51 to 1.53 (the runs and the
refit).

## 6. 5099: the denominators the restart does not read

Not extended; filed again in part. The two denominators the item names
have no zero to read:

- `1 + exp(u)` is above 1 for every real u.
- `EC50^h + x^h` with x >= 0: for EC50 > 0 both terms are at least 0 and
  EC50^h > 0, so the sum is positive; for EC50 < 0 and a non-integer h,
  EC50^h is NaN and the body is not defined. Read in u = EC50^h, the
  denominator is affine with walls at u = -x^h <= 0, which no EC50 > 0
  reaches. Only an integer h with EC50 < 0 has a zero (h = 1 is the
  affine Michaelis-Menten case the restart reads).

`dev/nlband-5099.R` (`5099-after3.txt`) checks this on four fits: each
coefficient over a grid of 4001 values from -100 to 100 times its
estimate's size, the others at the fit:

```
FIT hill logLik -143.262678 beta 9.8385 2.9379 2.0712 dens (EC50^h + x^h) 
  coef 1 (Emax_(Intercept)): walls NULL | grid -983.854 to 983.854: min finite 9.367, non-finite at 0, negative at 0
  coef 2 (EC50_(Intercept)): walls NULL | grid -293.788 to 293.788: min finite 0.04742, non-finite at 2000, negative at 0
  coef 3 (h_(Intercept)): walls NULL | grid -207.121 to 207.121: min finite 1.149e-97, non-finite at 0, negative at 0
FIT hillstep logLik -143.5856958 beta 10.015 4.7827 288.08 dens (EC50^h + x^h) 
  coef 1 (Emax_(Intercept)): walls NULL | grid -1001.48 to 1001.48: min finite 6.248e+195, non-finite at 0, negative at 0
  coef 2 (EC50_(Intercept)): walls NULL | grid -478.266 to 478.266: min finite 1.029e-74, non-finite at 3951, negative at 0
  coef 3 (h_(Intercept)): walls NULL | grid -28807.5 to 28807.5: min finite 0, non-finite at 3896, negative at 0
FIT logis logLik -141.2219891 beta 10.121 5.0457 1.0663 dens scal ; (1 + exp((xmid - x)/scal)) 
  coef 1 (Asym_(Intercept)): walls NULL | grid -1012.08 to 1012.08: min finite 1.01, non-finite at 0, negative at 0
  coef 2 (xmid_(Intercept)): walls NULL | grid -504.567 to 504.567: min finite 1, non-finite at 0, negative at 0
  coef 3 (scal_(Intercept)): walls NULL | grid -106.629 to 106.629: min finite -106.6, non-finite at 1, negative at 2000
FIT sig4 logLik -140.4547318 beta -0.13935 10.171 5.0192 1.105 dens scal ; (1 + exp((xmid - x)/scal)) 
  coef 1 (lo_(Intercept)): walls NULL | grid -100 to 100: min finite 1.012, non-finite at 0, negative at 0
  coef 2 (hi_(Intercept)): walls NULL | grid -1017.12 to 1017.12: min finite 1.012, non-finite at 0, negative at 0
  coef 3 (xmid_(Intercept)): walls NULL | grid -501.92 to 501.92: min finite 1, non-finite at 0, negative at 0
  coef 4 (scal_(Intercept)): walls NULL | grid -110.501 to 110.501: min finite -110.5, non-finite at 1, negative at 2000
```

`EC50^h + x^h` and `1 + exp(...)` are never negative; a non-finite
value comes only where EC50 < 0 (NaN), and the 0 of the hillstep h grid
is both terms underflowing at large negative h. The negative values are
of `scal`, the logistic's other denominator: it is affine in `scal` with
one wall (scal = 0 on every row), which nl_pole_walls() does not read
because a single wall leaves no band and no gap; the sign probe varies
the sign of `scal`. So reading these two forms would add no wall, and
the restart is right to skip them. A denominator affine in a monotone
function of a coefficient that does reach 0 (`a / (exp(b) - x)`, walls
at b = log(x_i)) is not read either; that is filed (LOW) with the
reason, not built: no row of the surveys has such a body.

The item's second half, binomial mm 22 of the fresh rows, is unchanged:
-643.501398 on both arms (`fresh-before/binomial.log`,
`fresh-after2/binomial.log`), Vm 7.2e6 and K 9.8e7 on the ridge to the
linear limit, warned. It is a fit far beyond one end of the band, which
the 0.83.0 rule does not restart (no pressure) and the hole rule does
not reach (no gap). Filed again (LOW).

## 7. OpenBLAS

`dev/ciharden-openblas.sh 0.3.26 lapack` and `0.3.32 lapack`
(`ob-build-0326.txt`, `ob-build-0332.txt`: Rlapack 616 routines to
OpenBLAS), then `dev/nlband-ob.sh` (OPENBLAS_CORETYPE=Haswell, 4
threads; this CPU has no AVX-512):

```
OpenBLAS 0.3.26 lapack:
RESULT test-fit-restart.R pass=106 fail=0 err=0 skip=0 warn=0
ulp band: 72 draws, 72 with the message and no warning, largest |difference| 6.18e-13
ulp worse: 72 draws, 72 identical to the escape without holes, 72 quiet
OpenBLAS 0.3.32 lapack:
RESULT test-fit-restart.R pass=106 fail=0 err=0 skip=0 warn=0
ulp band: 72 draws, 72 with the message and no warning, largest |difference| 6.18e-13
ulp worse: 72 draws, 72 identical to the escape without holes, 72 quiet
```

Reference BLAS reads pass 107 on the same file (`t-fit-restart-after2.txt`).
The one assertion fewer under OpenBLAS is the branch of the existing test
"a body on a ridge to its linear limit says so" that runs only when the
fit ends with the limit warning; its comment records that OpenBLAS 0.3.26
takes the other ending. Checked (`dev/nlband-obridge.R`): reference R ends
at -206.271197282 with `frmtmb_nl_limit`, OpenBLAS 0.3.26 at -206.27763956
with the standard-error warning, both code 0. No DETAIL line on either
copy.

## 8. Tests

`tests/testthat/test-fit-restart.R`, four new tests after the pole
tests:

- "a fit clear of its poles is started in a hole of the data": the band
  data, seeds 1 to 3, the default fit gives `frmtmb_nl_pole_start`,
  reaches the fit from (1, -1.9) within `grad_tol`, with b inside the
  hole; then 4 ulp draws of y on seed 1.
- "a given start clear of the poles warns of the hole": start (-0.5,
  0.4), the warning `frmtmb_nl_pole` (allowed by its text), the fit kept
  clear of the poles.
- "a hole whose start is worse, or no hole, keeps the fit": truth b 0.5
  with the hole (seeds 1, 2): no message, the hole is found, and the
  escape returns NULL without a refit (a recorder, since the escape
  catches a refit's error); every row kept: no hole; `nl_starts =
  FALSE` keeps b > 0; and the positive control: from that fit, with the
  switch on, the escape does refit, so the recorder can fire.
- "a hole is a gap much wider than the gaps around it": nl_pole_holes()
  on constructed walls: one hole between two even runs; none on an even
  run; none with fewer than 10 walls on one side; none with 20 or fewer
  gaps; of three holes the 2 widest, widest first.

Seen failing on rellib-r22 (`t-fit-restart-before1.txt`, the test file
before the recorder and the positive control were added):
`RESULT test-fit-restart.R pass=76 fail=20 err=1`. The failures are the
behavior: no `frmtmb_nl_pole_start` (3 seeds and 4 draws), logLik
591.548, 591.826, 646.401 and 591.548 (4 draws) below the reference, b
outside the hole, no `frmtmb_nl_pole` from the given start, `$holes`
NULL; the error is the missing nl_pole_holes() (the weak form, for the
unit test only). The final file on the final build
(`t-fit-restart-after2.txt`): pass 107, fail 0, err 0, skip 0, warn 0.

## 9. Suites

### Sessions

`dev/nlband-sessions.sh`: `test-fit-restart.R` and the ten `test-nl*.R`
files with the helpers and setup.R in one session each way
(`dev/nlband-session-tests/`, logs `session-<order>.txt`):

```
SESSION frmtmb forward files=11/11 pass=976 fail=0 err=0 skip=11 warn=0
LEAK lines: 0
SESSION frmtmb reverse files=11/11 pass=976 fail=0 err=0 skip=11 warn=0
LEAK lines: 0
SESSION frmtmb preload files=11/11 pass=976 fail=0 err=0 skip=11 warn=0
LEAK lines: 0
```
The 11 skips are test-nl-brms.R and test-nl-stan-funs.R behind
`FRMTMB_BRMS_FIT_TESTS`; both run in the gated tier.

### Plain tier, one file per process

`dev/nlband-tier.sh` over `dev/nlband-plain.jobs` (every test file of
core and the eight extensions, 549, split in two halves run at once),
NOT_CRAN=true, the lane's library first (`dev/release/nlband-plainA.log`,
`nlband-plainB.log`: "ran 275 of 275", "ran 274 of 274"). Every file's
lib: line names `wt-nlband-lib` for core. Per package
(`plain-sum.txt`):

```
frmtmb files=349 pass=25694 fail=0 err=0 skip=243 warn=0
frmtmb.coupling files=11 pass=551 fail=0 err=0 skip=5 warn=0
frmtmb.eam files=33 pass=1851 fail=0 err=0 skip=5 warn=0
frmtmb.latent files=13 pass=700 fail=0 err=0 skip=2 warn=0
frmtmb.learn files=16 pass=474 fail=0 err=0 skip=13 warn=0
frmtmb.ode files=11 pass=549 fail=0 err=0 skip=1 warn=0
frmtmb.sample files=92 pass=4286 fail=0 err=0 skip=18 warn=0
frmtmb.spatial files=8 pass=551 fail=0 err=0 skip=6 warn=0
frmtmb.survival files=16 pass=575 fail=0 err=0 skip=0 warn=0
```

### Still running at hand-back

- Gated tier (`dev/nlband-gated2.jobs`, 71 files: the union of the
  run-gated.ps1 set, the scale files and lane nlstop2's gated list;
  FRMTMB_BRMS_FIT_TESTS, DRMTMB and FUZZ on): 36 of 71 RESULT lines at
  hand-back, none with a failure or error. Log
  `dev/release/nlband-gated.log` when it ends.
- `R CMD check --as-cran` of core (`dev/nlband-check.sh`, log
  `dev/nlband-check/check.log`): running, no Status line yet.
- Not yet done: delete `dev/ciharden-out/` and `dev/nlband-check/`;
  rerun the final test file on rellib-r22 (the run recorded above was
  before the recorder and positive control were added).

## Bump and floor

Minor (one BREAKING behavior change in fits: a hole start can return a
different fit, with a message). No new dependency; floors unchanged.

Status: change done and measured; gated tier and check pending.

## 10. Fix round 1 (review dev/reviews/2026-11-17-nlband.md, verdict FIX)

State at the first round's end: gated tier (`dev/nlband-gated2.jobs`,
71 files) "NLBAND-GATED ran 71 of 71", pass 5393, fail 0, err 0, skip
15, warn 0; the 15 skips are the scale files (`test-scale.R`,
`test-curve-scale.R`), whose gate was not set. Core `R CMD check
--as-cran`: "Status: 1 NOTE" (V8 unavailable for the HTML manual's math).

Item 1 (separation). Neither existing check marks the separated hole
fit: `dev/nlband-dbg3.R` (`dbg3.txt`) on seeds 2 to 4 gives
nl_limit_message() NULL, NULL also with its one-standard-error gate
removed, and se_fit_analysis() tier 1 with nothing lost (SEs 8.6e133 and
1.4e131). New `nl_hole_located()`: the refit of a hole run is kept only
where the standard error of the hole's coefficient, from the SE check's
own analysis (se_fit_analysis(), shown covariance), is finite and below
half the hole's width. `dev/nlband-fix1.R` (`fix1-before.txt`,
`fix1-after1.txt`): SEP seeds 1 to 10 now identical to 0.83.0 (code 1,
`frmtmb_conv_code` + `frmtmb_conv_gradient`, same logLik and
estimates).

Item 2 (student). The second hole run now starts the whole model at
the cold start of make_start() (`nl_cold_par()`), the hole's coefficient
at the hole's middle. STUD seeds 1 to 3 now reach -42.621233,
-43.081863, -78.591654, equal to the fits from (1, -2) (0.83.0:
-549.88886, -524.53729, -599.66309). The pole path no longer says the
first fit's boundary message (`frmtmb_boundary_fit`), which on these
seeds said nu was at its limit beside a returned fit with log nu 0.95
(`R/par-template.R`, nl_fit_probed()).

Tests: two new tests (separated seeds 2 to 4; student seeds 1 to 3).
Seen failing on the lane's build with the two fixes taken out
(`dev/nlband-mockrun.R`, `t-fit-restart-mock.txt`): pass 113, fail 9,
warn 1 (the separated fits keep the message and give no warning; the
student fits 16.199, 20.326, 43.973 short). Final file, final build
(`t-fit-restart-after3.txt`): pass 122, fail 0, err 0. 21 roundings of
the new assertions (`dev/nlband-fix1ulp.R`): 126 of 126 PASS on
reference BLAS and on OpenBLAS 0.3.26 with LAPACK; the separated fits
end at code 1 with the two convergence warnings on all 63 draws of
each.

Item 6. The review's run of the final first-round test file on
rellib-r22: pass 79, fail 21, err 1 (`dev/reviews/nlband-review-out/
t-fit-restart-r22.txt`), the positive control among the failures. The
survey's "0 of 2555 rows differ" is at the printed precision of the ROW
lines (logLik to 6 decimals, estimates to 4 or 5 digits); the bitwise
evidence is the 72 identical() draws of section 4 and the review's 33
rows at 17 digits.

Items 3 to 5: roxygen of nl_pole_holes() and nl_pole_escape()
corrected; section 2 corrected; backlog filed (HIGH dense pole, MEDIUM
near-wall keep rule, LOW SE check at 1e133, LOW the review's misses);
NEWS cost corrected (1.47 to 1.58 with two holes, review finding 9).

Runs of the fix round so far: nl test files on reference BLAS, brms fit
tests on (`dev/release/nlband-nlref.log`): 11 of 11, pass 1122, fail 0,
err 0, skip 0, warn 0.

OpenBLAS 0.3.26 with LAPACK, the same 11 files (`dev/nlband-out/obnl/`):
11 of 11, fail 0, err 0, skip 0, warn 0, pass 1121 (test-fit-restart.R
121, one fewer than reference BLAS: the ridge test's branch, section 7).

Cost on the fixed build (`dev/nlband-cost.R`, `cost-after3.txt`,
evaluation counts on and off; the wall-time controls read 0.885 to
1.095, so the counts decide): one hole that loses, n 400: 1.371, 1.329,
1.309; n 20000: 1.273, 1.226; two holes that lose (the review's
jittered replicates): 1.589, 1.569, 1.523; a hole that wins: 1.586,
1.520, 1.531; no hole: 1.000 on 3 of 3. NEWS now gives 1.23 to 1.37
and 1.52 to 1.59.

Starting the one-session runs after the cost rerun was refused by the
permission system (auto mode classifier, "Destructive or Unwanted
Actions"); per the round's rules it was not worked around.

Survey on the fixed build (`accsum-after3.txt` against 0.83.0,
`accsum-after3-vs2.txt` against the first build): 2555 rows, 0 differ
in either comparison at printed precision; 12 hole runs on the same 6
fits, all losing; evaluations 4,099,561 (0.83.0 4,098,850, first build
4,099,581).

Owed runs, approved by the user on 2026-11-17 and run to the end:

- One-session runs of test-fit-restart.R and the ten test-nl*.R files
  (`dev/nlband-sessions.sh`, the copies in `dev/nlband-session-tests/`
  refreshed from `tests/testthat/` first; a run on the stale copies,
  976 pass, is not the record): forward, reverse and preload each
  "files=11/11 pass=991 fail=0 err=0 skip=11 warn=0", 0 LEAK lines.
- Plain core tier, one file per process (`dev/release/nlband-core2.log`):
  "ran 349 of 349", pass 25709, fail 0, err 0, skip 243, warn 0; every
  lib: line names `wt-nlband-lib`.
- `R CMD check --as-cran` of the fixed core: "Status: 1 NOTE" (V8
  unavailable for the HTML manual's math); tests OK.
- `dev/ciharden-out/` and `dev/nlband-check/` deleted.

Status: done. No job of this lane is running.

## 11. Fix round 2 (short check dev/reviews/nlband-short.md, verdict FIX)

Item 1 (the kept fit's own boundary message). nl_refit_held() now holds
the refit's `frmtmb_boundary_fit` messages beside its warnings, and
nl_fit_probed() says them with `pe$held`; every other message of the
refit stays muffled, as before. `dev/reviews/nlband-short-studg.R` on
the round-2 build (`studg-after4.txt`): student() on the gaussian band
data, seed 2, the default fit now says `frmtmb_nl_pole_start` and
`frmtmb_limit_fit` ("Boundary fit: nu_(Intercept) is at the end of its
range (20.4 ..."), as the fit from (1, -2) does.

Item 2 (the located bound). Design: the refit is kept where the
standard error of the hole's coefficient is below the span of the walls
(`nl_hole_se()`, `ck$span`), not half the hole's width. The bound is a
length the data set (the range of the rows' poles), and a supremum at
infinity has an unbounded standard error: the separated fits have 1e131
against a span of about 4, the short check's refused maxima at most 0.54
(1.4 half widths, 0.137 spans). Every run that beats the fit is now a
candidate, the highest first; a refit refused for any reason passes to
the next; the first hole refit refused by the bound is said in a message
of class `frmtmb_nl_pole_refused`, also when a later run is kept.

Sweep (`dev/nlband-sweep.R`, the short check's sweep with nl_hole_se()
traced, `sweep-after4.txt`), on the same 101 calls:

```
SWEEP fits 150 | calls 101 | refused by the round-1 rule 5 | by the round-2 rule 0 | max se/half 1.404 | max se/span 0.137
```

The five calls the round-1 bound refused, each now kept:

```
CALL amp 0.15 sd 2 seed 8 | se 0.5441 half 0.408 span 3.985 | refused round1 TRUE round2 FALSE | refit ll -122.565609 | returned ll -122.565609 code 0 b -2.1782 warn - msg frmtmb_nl_pole_start | GLOB 
CALL amp 0.15 sd 2 seed 13 | se 0.3827 half 0.2725 span 3.948 | refused round1 TRUE round2 FALSE | refit ll -127.942636 | returned ll -127.942636 code 0 b -1.9712 warn - msg frmtmb_nl_pole_start | GLO
CALL amp 0.20 sd 2 seed 8 | se 0.4494 half 0.408 span 3.985 | refused round1 TRUE round2 FALSE | refit ll -122.567193 | returned ll -122.567193 code 0 b -2.1399 warn - msg frmtmb_nl_pole_start | GLOB 
CALL amp 0.20 sd 2 seed 13 | se 0.2948 half 0.2725 span 3.948 | refused round1 TRUE round2 FALSE | refit ll -127.943073 | returned ll -127.943073 code 0 b -1.9799 warn - msg frmtmb_nl_pole_start | GLO
CALL amp 0.20 sd 2 seed 14 | se 0.3427 half 0.3339 span 3.898 | refused round1 TRUE round2 FALSE | refit ll -134.124191 | returned ll -134.124191 code 0 b -1.9979 warn - msg frmtmb_nl_pole_start | GLO
```

Each returns the refused refit's logLik with `frmtmb_nl_pole_start`,
code 0, no warning. None is the global maximum (a pole between rows,
the dense-pole item filed HIGH); each is 0.056 to 1.462 above the fit
the round-1 build returned.

Separated seeds 1 to 10 (`dev/nlband-fix1.R`, `fix1-after4.txt`): the
same logLik, code and warnings as 0.83.0 on all 10; seeds 2 to 10 add
`frmtmb_nl_pole_refused` (the refused refit is higher, about -1e-268);
seed 1 has no hole run that beats the fit. Student seeds 1 to 3 and band
seeds 1 to 3 unchanged from round 1.

Item 3. The student test's bound is now 1e-6 of the gap the hole start
closes (logLik of the fit from (1, -2) minus that of the fit with
`nl_starts = FALSE`, 481 to 521), not grad_tol.

Tests: two new tests, "a kept hole fit says its own boundary message"
(seed 2) and "a hole fit is kept where the data place its pole loosely"
(the sweep's amp 0.2, sd 2, seed 13); the separated test also asserts
the refusal message. Seen failing on the build with round 1's two rules
put back (`dev/nlband-mockrun2.R`, `t-fit-restart-mock2.txt`): pass 127,
fail 4 (no boundary message; the loose fit not kept: no message, b
outside the hole, logLik equal to the stuck fit's). Round-2 build
(`t-fit-restart-after4.txt`): pass 131, fail 0, err 0, skip 0, warn 0.
21 roundings of the new assertions (`dev/nlband-fix2ulp.R`,
`fix2ulp-ref.txt`): 105 of 105 PASS on reference BLAS (63 SEP, 21
STUDG, 21 LOOSE).

Consolidation note: `R/par-template.R`, nl_fit_probed()'s pole-path
filter, conflicts with lane gpselimit's change on the same line (base
line 1619). Resolution: keep both entries and both comments,

```r
      !inherits(cnd, c("warning", "frmtmb_nl_local", "frmtmb_nl_restart",
                       "frmtmb_boundary_fit", "frmtmb_gp_basis"))
```

(gpselimit's message carries both classes, so nlband's entry already
drops it on the pole path; short check section 6).

Line endings: an edit of this round had rewritten `R/fit-restart.R`
with LF endings; it and `tests/testthat/test-fit-restart.R`,
`dev/test-backlog.md` and this file were put back to CRLF, as the rest
of the tree is checked out (core.autocrlf true).

Reruns of fix round 2:

- `test-fit-restart.R` under OpenBLAS 0.3.26 with LAPACK
  (`ob2-0326-fit-restart.txt`): pass 130, fail 0, err 0 (one fewer than
  reference BLAS: the ridge test's branch); 21 roundings
  (`fix2ulp-ob0326.txt`): 105 of 105 PASS.
- One-session runs (copies refreshed first): forward, reverse and
  preload each "files=11/11 pass=1000 fail=0 err=0 skip=11 warn=0", 0
  LEAK lines.
- Survey on the round-2 build (`accsum-after4.txt` against 0.83.0,
  `accsum-after4-vs3.txt` against round 1): 2555 rows, 0 differ in
  either comparison at printed precision; 12 hole runs, all losing, the
  same as round 1; evaluations 4,099,561, equal to round 1.
- `R CMD check --as-cran` of the round-2 core: "Status: 1 NOTE" (V8
  unavailable for the HTML manual's math); tests OK.
- `dev/ciharden-out/` and `dev/nlband-check/` deleted.

Status at the end of fix round 2: done. No job of this lane is running.
