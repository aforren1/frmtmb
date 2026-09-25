# Phase 3a: real data in, items 3.1, 3.2 and 3.6

Lane `wt-phase3a`, worktree `frmtmb-wt-phase3a`, branch `wt-phase3a`, off
main at `51bfaa4e` (frmtmb 0.62.0; frmtmb.ode 0.6.0, frmtmb.coupling
0.5.0, frmtmb.spline 0.7.0). Nothing under core's `R/` changed. The base
build is `rellib-r3`; the lane build is `phase3a-lib`. Every script named
below takes `PHASE3A_ARM=base` or `lane` through `dev/phase3a-lib.R`, and
its output is in `dev/phase3a-log/`.

## The answer

* **3.6, the defect, is fixed.** `rp_floored()` now tests every censored
  row as well as every event row, and counts the censored ones apart as
  `n_nonmonotone_censored`. On the construction that filed the defect it
  fires on 4 of 6 seeds, exactly the 4 where five no-death centres land
  on a falling slope, and flags exactly those centres' 50 rows. On the
  designs sampled in round 0 it fired on **0 of 1600** simulated fits
  and **0 of 270** real-data fits at `df` 1 to 5. Those designs had no
  cure fraction: on a cure-fraction design with `gamma1 ~ arm` it fires
  on **27 of 80** fits (punch round 1, below). By the user's decision
  of 2026-09-24 it REFUSES only on event rows and WARNS on censored
  rows, with the size of the rise: those 27 now warn and answer, 53 are
  silent and none refuses (punch round 2, below). The
  inverse test that item 2.5 wrote to fail when this landed failed on
  the base build and is flipped.
* **3.1 ships, after punch round 1.** `frm_ode_records()` splits a
  NONMEM-shaped table into `data` and `events`. Round 0 validated it on
  schedules whose event times never coincided, and the review showed
  that was where it was wrong: a reset, steady state, replace or
  multiply at the same time as another event is read three ways by
  NONMEM, `frm_ode()` and rxode2, and round 0 passed it. It is refused
  now, `addl` repeats included, with the other holes the review found.
  Against `rxode2::et()`, its doses equal `etExpand()`'s on 126 of 126
  schedules; **those 126 contain no reset, because `etExpand()` errors
  on every table with `evid = 3`**, so reset schedules were compared on
  trajectories only. Against `rxode2::rxSolve()`, `frm_ode()` on its
  output agrees to 7.1e-13 of the trajectory's scale over 190
  continuous-time random schedules, to 1.0e-12 over the 132 of 300
  grid-time schedules it accepts, and to 4.2e-15 and 2.2e-13 on
  `nlmixr2data::warfarin` and `theo_md`. The Theoph vignette example is
  rewritten through it.
* **3.2 ships.** `frm_cross_spectrum()` takes lists of epochs of unequal
  length with a `group`, and a vector pair with a per-sample `group`.
  The epoch list equals the concatenated-and-split call bit for bit
  whenever they describe the same records. `frm_cross_pairs()` stacks
  channel pairs with a `pair` factor; a 4-pair frame fits with
  `(1 | pair)`, and with every parameter by pair it equals the four
  two-channel fits to 4.5e-04 of their own standard errors.

## 3.6: rp_floored() and a group with no events

### Reproduced on the released build first

`dev/phase3a-repro36.R` is `dev/frailty/frailty-floor3.R`'s construction
(400 subjects, 40 centres of 10, shape 0.6, centre `sd` 0.35, five
centres followed to a common administrative time with no deaths),
seeds 20260910 to 20260915, fitted with `gamma1 ~ (1 | centre)`.

Base build (`phase3a-repro36-base.txt`): on seeds 20260910, 11, 12 and
14 the minimum centre slope is -0.2852, -0.2583, -0.3056 and -0.2785,
five centres at or below zero each time, and every one of the three
gates is silent: `n_nonmonotone` 0, `rp_floored()` does not refuse,
`frm_curve()` does not refuse, no warning at fit end. These are the
frailty lane's numbers to four digits.

Lane build (`phase3a-repro36-lane.txt`), same seeds: `n_nonmonotone`
still 0, `n_nonmonotone_censored` **50** on each of the four, the
flagged rows `identical()` to the no-death centres' rows on all four,
`rp_floored()` and `frm_curve()` refuse, and the fit warns once. The 20
negative centre slopes run from -0.1773 to -0.3056 (the frailty lane's
"-0.21 to -0.31" was the per-seed minima; the help pages now say
-0.18 to -0.31). Seed 20260910 converges at code 0 with a maximum
gradient of 2.74e-05, as recorded. The six control fits with no such
centre report 0 on every count and refuse nothing.

### The charter question, and what was decided

The plan row asked whether the widened check answers a different
question and deserves its own name. It does answer a different
question, so it gets its own COUNT, not its own function:

* `n_nonmonotone`: EVENT rows with `d(eta)/d(log t) <= 0`. Their density
  is floored, so `logLik()` is a pseudo-likelihood. Unchanged.
* `n_nonmonotone_censored`: CENSORED rows with `d(eta)/d(log t) <= 0`,
  an interval row tested at both ends. Their `log S` is exact, so
  `logLik()` IS the model's; the model is what is wrong, because its
  survival function rises there.

`action = "error"` refused on either in rounds 0 and 1. Since punch
round 2 it refuses on the first and warns on the second; see "Punch
round 2" below. One function keeps `frm_curve()` and its two
companions gated through the one call they already make. A separate
function would have needed a second gate in each, and a user who calls
only `rp_floored()`, which every help page tells them to, would not see
it.

The row test is exact at `df = 1`, where the derivative is `gamma1 + u`
at every time. At `df >= 2` a dip between two observed times of one
group that recovers at both is not searched for. That is written on
both help pages.

### False-alarm rate

`dev/phase3a-falsealarm36.R` (simulated) and
`dev/phase3a-falsealarm36b.R` (real; the first file's real branch had a
data-loading bug and was never run to completion, see "Records"
below). One fit per row, and both counts come from ONE report, since
`n_nonmonotone` is the 0.7.0 quantity unchanged.

Simulated: four baseline hazards (Weibull 0.8 falling, Weibull 1.5
rising, log-logistic 2.5 rise-then-fall, lognormal sdlog 1), n = 200
and 1000, light and heavy censoring (administrative plus uniform
dropout, 17 to 23 percent and 63 to 71 percent), a treatment arm that
is proportional hazards only for the Weibull, 10 seeds each
(20260924 to 20260933), each fitted at `df` 1 to 5 both as
proportional hazards and with `gamma1 ~ trt`: 1600 fits.

Real: `flexsurv::bc`, and from survival `lung`, `veteran`, `colon`
(death), `rotterdam` (death), `gbsg` (recurrence), `pbc` (death),
`ovarian` and `mgus2` (death), each with one binary covariate, at `df`
1 to 5, on all three scales, proportional hazards and `gamma1 ~ x`:
270 fits. `bc` and `gbsg` are the same 686 patients with different
covariates, so this is nine designs from eight datasets.

The block below is generated by `dev/phase3a-fa36-report.R` and pasted
verbatim (full file `dev/phase3a-log/phase3a-fa36-report.txt`):

    <!-- generated by dev/phase3a-fa36-report.R -->
    simulated: 1600 rows, 1600 fitted, 0 fit errors
      nlminb code 0 on 1563 of 1600
      0.7.0 check fires (event rows): 0 of 1600
      widened check fires: 0 of 1600
      fires ONLY because of the widening: 0 of 1600, 95% upper bound 0.0023
      fit-end censored-row warning raised: 0 of 1600
      rp_floored(): refuse 0, warn 0, silent 1600 of 1600
      smallest d(eta)/d(log t) on any censored row: 0.04602
     df     form fits old new new_only
      1 gamma1~x  160   0   0        0
      2 gamma1~x  160   0   0        0
      3 gamma1~x  160   0   0        0
      4 gamma1~x  160   0   0        0
      5 gamma1~x  160   0   0        0
      1       ph  160   0   0        0
      2       ph  160   0   0        0
      3       ph  160   0   0        0
      4       ph  160   0   0        0
      5       ph  160   0   0        0

    censoring fraction by design (simulated):
           design    n cens_frac
     llogis_heavy  200     0.688
     llogis_light  200     0.216
      lnorm_heavy  200     0.665
      lnorm_light  200     0.190
     weib08_heavy  200     0.629
     weib08_light  200     0.173
     weib15_heavy  200     0.684
     weib15_light  200     0.217
     llogis_heavy 1000     0.706
     llogis_light 1000     0.225
      lnorm_heavy 1000     0.670
      lnorm_light 1000     0.193
     weib08_heavy 1000     0.633
     weib08_light 1000     0.184
     weib15_heavy 1000     0.687
     weib15_light 1000     0.226

    real: 270 rows, 270 fitted, 0 fit errors
      nlminb code 0 on 260 of 270
      0.7.0 check fires (event rows): 0 of 270
      widened check fires: 0 of 270
      fires ONLY because of the widening: 0 of 270, 95% upper bound 0.0136
      fit-end censored-row warning raised: 0 of 270
      rp_floored(): refuse 0, warn 0, silent 270 of 270
      smallest d(eta)/d(log t) on any censored row: 0.2631
     df     form fits old new new_only
      1 gamma1~x   27   0   0        0
      2 gamma1~x   27   0   0        0
      3 gamma1~x   27   0   0        0
      4 gamma1~x   27   0   0        0
      5 gamma1~x   27   0   0        0
      1       ph   27   0   0        0
      2       ph   27   0   0        0
      3       ph   27   0   0        0
      4       ph   27   0   0        0
      5       ph   27   0   0        0

    rows and censoring fraction (real):
              design    n cens_frac
                  bc  686     0.564
                lung  228     0.276
             veteran  137     0.066
         colon_death  929     0.513
     rotterdam_death 2982     0.573
            gbsg_rfs  686     0.564
           pbc_death  312     0.599
             ovarian   26     0.538
         mgus2_death 1384     0.304
    <!-- end generated -->

A zero is only evidence if the check can fire, and it does: on the
construction above it fires on 4 of 6 seeds and on exactly the rows it
should. The smallest censored-row derivative on any ordinary fit was
0.046 simulated and 0.263 real, so no fit came near the line. 37 of the
1600 simulated and 10 of the 270 real fits ended at a nonzero `nlminb`
code; they are counted, not dropped, and none fired.

### A cure fraction with a time-varying effect: it fires (punch round 1)

The review found the widened check firing on 15 of 80 `gamma1 ~ x`
fits of cure-fraction designs at `df` 3 to 6. Reproduced on this lane's
own construction (`dev/phase3a-cure36.R`, `phase3a-cure36.txt`): two
arms of 250, 40 and 55 percent cured, Weibull 1.3 for the rest with the
arm shifting the scale, administrative censoring at 6 plus uniform
dropout, `t | cens(censored) ~ x, gamma1 ~ x`, `df` 3 to 6, seeds
20260924 to 20260943. It fires on **27 of 80** (4, 7, 8 and 8 at `df`
3, 4, 5 and 6). On all 27 every flagged row lies past its own arm's
last event, and the fitted survival rises across the flagged rows by
5.2e-04 to 0.13. The count is right, since that survival function does
rise; the rise sits where the arm has no event to prevent it.

The round-0 sample missed this because none of its designs had a cure
fraction, so "0 of 1600" is a statement about those designs only. The
reason a fit WITHOUT a time-varying effect cannot fire past the last
event: every row shares one spline, and past the boundary knot its
derivative is the one at the last event row, which the event-row check
already covers.

This was held for the user in round 1 and decided in round 2, below.

### Punch round 2: refuse on event rows, warn on censored rows

The user's decision of 2026-09-24. `n_nonmonotone` (event rows, whose
density is floored) still refuses under `action = "error"`, and
`frm_curve()` and its two companions still refuse with it.
`n_nonmonotone_censored` (censored rows, interval rows at both ends)
now WARNS, under `action = "error"` too, and the curve functions
answer with that warning. The fit's own end-of-fit warning stays and
uses the same text. The warning gives the size: the new report field
`max_survival_rise` is, over the flagged rows, the largest amount the
fitted survival climbs above its running minimum on that row's own
coefficients, from the smallest observed time to the row's own (upper)
time on a 200-point grid.

The reasoning, as the Rd and NEWS give it: brms's `cox()` cannot
produce a rising survival (its baseline hazard is an M-spline basis
times non-negative weights); `flexsurv::dsurvspline()` clamps the
density to 0 where the derivative is not positive, so an EVENT row
there is `-Inf` and the fit avoids it while a CENSORED row is never
checked and a rise there passes silently; rstpm2 penalizes negative
hazards during the fit (`kappa.init = 1` up to `maxkappa = 1000`),
recorded as a possible future option in `dev/test-backlog.md` only.
frmtmb draws flexsurv's line and says so.

Reruns on the round-2 build:

| design | fits | refuse | warn | silent |
|---|---|---|---|---|
| cure fraction, `gamma1 ~ x`, `df` 3 to 6 (`phase3a-cure36.txt`) | 80 | 0 | 27 | 53 |
| round-0 simulated (`phase3a-fa36-report.txt`) | 1600 | 0 | 0 | 1600 |
| round-0 real data | 270 | 0 | 0 | 270 |
| no-death centres, seeds 20260910 to 15 (`phase3a-repro36-lane.txt`) | 6 | 0 | 4 | 2 |
| the same design's controls | 6 | 0 | 0 | 6 |

On the 27 warning cure fits the reported `max_survival_rise` runs from
5.23e-04 to 0.133, and agrees with the study's own independent
computation of the rise to the third digit (0.0390 against 0.0389 on
the first). The refusal is reached only by an event row on a
non-positive slope. None of this lane's designs reached it: 0 of 1870 on
the round-0 designs, 0 of 80 on the cure design with `gamma1 ~ x`, and 0
of 80 on the same cure design fitted with proportional hazards
(`dev/phase3a-cure36-ph.R`, `phase3a-cure36-ph.txt`, all 80 silent).
**An earlier draft said "no ordinary design refused in any round", and
that was false.** On the reviewer's cure design fitted with proportional
hazards at `df` 3 to 6, 43 of 120 fits refuse on EVENT rows. That is
unchanged 0.7.0 behavior and follows the user's rule; the 43 is the
reviewer's count, not re-measured here, and this lane's own cure design
does not reproduce it, so the difference is in the design.

**Tests, seen failing on the round-1 code.** The round-1 `rp-check.R`
and `royston-parmar.R` were copied aside before the edit and swapped
into the lane namespace (`dev/phase3a-round1-spline.R`,
`phase3a-p2-round1-spline.txt`): the new and changed tests fail 17
assertions in `test-rp-floored.R` and 2 in `test-frailty.R`, 0 errors.
The censored-only block (the cure design at seed 20260932, `df = 4`)
fails 7 there. The event-row block passes on round 1, as it must: that
half of the behavior did not change and the block pins it. The
silent-control block fails only its `max_survival_rise == 0`
assertion, because the field is new; its silence passed on round 1
too. On the lane build `test-rp-floored.R` passes 81 of 81 and
`test-frailty.R` 30 of 30.

Found on the way and filed, not fixed (`dev/test-backlog.md`):
`frm_curve(dpar = "mu")` on a `gamma1 ~ x` fit refuses with its own
covariance cross-check ("disagrees with frm_linpred(se.fit = TRUE) by 1
relative"), on the released 0.7.0 build too; `dpar = "gamma1"` on the
same fit answers, and `dpar = "mu"` on a proportional-hazards fit of
the same data answers. The new test reads `gamma1` for that reason.
That was one fit: on the reviewer's cure design `dpar = "gamma1"`
refuses too, on 3 of 4 fits, identically on the released build (the
reviewer's count). The backlog entry says so.

### Tests, seen failing on the released build

In `phase3a-seefail36-*` and `seefail/`, the new and flipped tests on
the base build:

* `test-rp-floored.R`, "a centre with no events and a falling slope is
  refused": the fitted construction at seed 20260910; the fit-end
  warning, the refusal, the `frm_curve()` refusal, the count of 50 and
  the row indices each fail on base.
* "the censored-row count reads an interval's upper end": a spline
  written so that only the upper ends of interval rows are
  non-monotone.
* "an ordinary fit reports zero censored rows where S rises":
  `flexsurv::bc` with `gamma1 ~ group` at `df = 3`, the inverse case.
  It fails on base only because the field does not exist there.
* `test-frailty.R`: the two assertions item 2.5 wrote to fail when this
  landed failed on base and are flipped to expect the count, with the
  likelihood count still asserted 0.
* The existing `expect_named()` in `test-rp-floored.R` gains the field.

Mutation (`phase3a-mutate36.R`): with the upper end not read, the
interval block fails 3 assertions; with no censored row tested (the
0.7.0 behavior), the no-events block fails 5 and the interval block 3.

## 3.1: frm_ode_records()

### What it does

One NONMEM-shaped table in, `list(data, events)` out, in `frm_ode()`'s
own `events` grammar. `evid` 0 is a data row; 1 an `"add"`; 3 a
`"reset"`; 4 a reset and then an add; rxode2's 5 and 6 `"replace"` and
`"multiply"`. A positive `rate` becomes `duration = amt / rate`; a
positive `dur` is the duration. `ii`, `addl` and `ss` are copied.
Column arguments match without regard to case. `data` keeps every
column except the dosing items. It does no arithmetic on a dose.

`dur` is an argument the plan's signature did not have. Without it,
every infusion rxode2's `et(dur = )` writes would be refused by name,
and rxode2 is the reference.

### Refused by name

Every refusal is a `frmtmb_ode_error` naming the column or the record:
`evid = 2` (a covariate-change record that would otherwise be lost or
become an observation), any other code including rxode2's classic 101,
`ss = 2`, a steady-state constant infusion, a negative (modeled) rate
or duration, a rate and a duration together, an infusion on a
non-dose record, a non-whole `addl`, a negative or missing dose, a
nonzero `amt` on an observation, a dose with no `cmt` or a `cmt` of
zero or below, an `mdv` that disagrees with `evid`, two `mdv` columns,
an NA `id`, the NONMEM items `date`, `dat1` to `dat3`, `cont`, `call`,
`pcmt`, `l1` and `l2`, an argument naming no column, two columns
differing only in case, an observation listed after an event at the
same time, and a reset inside a running infusion.

Two of those are there because the validation found them:

* **An observation at a dose time.** `frm_ode()` reads it before the
  dose. NONMEM reads records at one time in listed order. rxode2 reads
  it AFTER the dose whatever the order: `dev/phase3a-ode-rxode2.R`
  part D, depot at `t = 12` is 100.000185 in rxode2 with the sample
  listed before the dose and 100.000185 after, against 0.000185 in
  `frm_ode()`. So a sample listed after its dose is refused, a sample
  listed before it passes, and the help page says rxode2 differs.
  Both real datasets ship with their time-0 sample after the dose and
  are refused as shipped; reordered, they agree (below).
* **A reset inside a running infusion.** Among 200 random schedules, 10
  had an `evid = 3` while an infusion was running. There `frm_ode()`
  lets the infusion run on after the reset, `frm_lincmt()` refuses the
  schedule, and rxode2 5.1.7 returns NEGATIVE amounts on all 10, as low
  as -261.4. With three readings and one of them impossible, the
  records do not say one thing, and `frm_ode_records()` refuses. That
  `frm_ode()` runs the infusion on is a property of `frm_ode()`,
  pre-existing and not changed here; see "Filed, not fixed".

### Validation against rxode2

`dev/phase3a-ode-rxode2.R`, rxode2 5.1.7, a one-compartment oral model
written for rxode2 as an ODE, `ka = 1.1`, `ke = 0.23`, both integrators
at `atol = rtol = 1e-12`. Each difference is the largest absolute
difference over the observation rows divided by the largest absolute
state in rxode2's solution.

Nine hand-built schedules (oral `ii`/`addl`; IV infusion by `rate`; by
`dur`; oral plus IV in two ids; `evid` 3; `evid` 4; `evid` 5 and 6;
steady-state oral; steady-state infusion): expansion equal to
`etExpand()` on 8 of 8 that etExpand() expands; `frm_ode()` against
rxSolve at most 3.8e-13; `frm_lincmt()` at most 1.5e-12 (it refuses the
replace/multiply schedule by design).

200 random schedules (seeds 20260924 to 20261123; 1 to 3 ids, 1 to 5
doses each into either compartment, bolus, rate or dur, `ii`/`addl` on
30 percent, a reset on 20 percent of ids): 190 solved, 10 refused (all
10 the reset-in-infusion case). Over the 190, 4620 observations and
1217 dose events: `frm_ode()` central worst 7.12e-13, median 1.44e-14,
bitwise equal on 36; depot worst 2.56e-13; `frm_lincmt()` worst
3.05e-12, median 8.71e-13. Expansion equal to `etExpand()` on 126 of
126; etExpand() itself errors on the other 64, and
`dev/phase3a-etexpand-check.R` shows its failures are exactly the 74 of
200 schedules that carry an `evid = 3` record.

The bitwise agreement is not a construction error: against the
closed-form superposition on the first schedule both are 3.02e-13 of
the scale off, and they agree with each other to the last bit
(`identical()` TRUE), which reads as two LSODA implementations taking
the same steps on the same split intervals.
`frm_lincmt()`, an independent closed form, is the arm that matters.

Real records: `nlmixr2data::warfarin` (515 records, 483 observations,
32 doses; no `cmt` column, so doses are named into the depot) and
`nlmixr2data::theo_md` (348 records, 264 observations, 84 doses; its
classic `evid = 101` rewritten as 1), each with its time-0 samples
listed before its dose. `frm_ode()` central against rxSolve on the
records: 4.15e-15 and 2.20e-13; `frm_lincmt()`: 2.05e-13 and 2.32e-13.

### The Theoph example

`dev/phase3a-theoph.R`: the vignette's model fitted with the dose as
the initial condition, as the vignette had it, and through
`frm_ode_records()` as it has it now. Log likelihoods -191.19178194
and -191.19178193; the four fixed effects agree to 2.8e-06 absolute;
the two objectives at one parameter vector differ by 2.6e-08. But at the
default `atol = rtol = 1e-8` the records spelling stops at `nlminb`
code 1, false convergence, with a gradient of 1.67e-03, and restarts do
not change that (`phase3a-theoph-conv.R`: code 1 at restarts 1 and 3).
At `1e-10` it stops at code 0 with a gradient of 2.3e-05 and no
warning, at -191.19178184. The vignette uses `1e-10` and says why.

### Punch round 1: events at one time, and the holes the review found

The review's blocker: round 0 drew event times from a continuous range,
so no two records ever shared a time, and `frm_ode_records()` passed
every event-against-event conflict. `frm_ode()` applies a reset or a
steady state first at an instant whatever the listing, NONMEM reads
records at one time in listed order (from its documentation; no NONMEM
here), and rxode2 5.1.7 follows neither.

**The rule, decided by the user:** a reset (`evid` 3 or 4), a
steady-state dose (`ss = 1`), a replace or a multiply (`evid` 5 or 6)
that shares its time with another event of the same id is refused,
`frm_ode_records(): record k, a reset, shares time t with the event from
record j`, classed `frmtmb_ode_error`. Event times are written out over
`addl` first, so a repeat that lands on a reset is caught and the
message says it is a repeat. Two plain doses at one time still add and
pass. An observation at a steady-state record's time is refused in
either order, since `frm_ode()` reads the trough there (6.76) and
rxode2 the post-dose state (106.76); what NONMEM reads there was not
checked, and the Rd says so.

`dev/phase3a-ode-sametime.R`, central amount, atol = rtol = 1e-12,
`phase3a-ode-sametime-lane.txt`. A and B are the two listed orders:

| case | frm_ode | frm_lincmt | rxode2 A | rxode2 B | records |
|---|---|---|---|---|---|
| M2 bolus, then evid 3, t = 10, read at 11 | 79.45 | 79.45 | 79.45 | 79.45 | refused both orders |
| M3 bolus and ss = 1 dose at t = 10 | 112.03 | 112.03 | 32.58 | 32.58 | refused both orders |
| M4 4-hour infusion q4 x3, evid 3 at t = 4; t = 6, 13, 20 | 40.08, 72.65, 14.52 | same | 0.00, -22.33, -91.43 | | refused (addl repeat) |
| M5 bolus and replace 30 at t = 10 | refused by frm_ode | refused | 103.29 | 103.29 | refused both orders |
| M6 bolus and multiply 0.5 at t = 10 | refused by frm_ode | refused | 83.44 | 83.44 | refused both orders |
| M7 bolus listed before evid 4 at t = 10 | 91.13 | 91.13 | 91.13 | | refused |
| addl repeat at t = 12 on an evid 3 | 79.45 | 79.45 | 79.45 | | refused (addl repeat) |
| sample at an ss = 1 dose into central, t = 0 | 6.76 | 6.76 | 106.76 | 106.76 | refused both orders |
| control: two plain doses at t = 10 | 108.64 | 108.64 | 108.64 | 108.64 | accepted |

M2, M3 and M4 reproduce the review's numbers. M5 to M7 are this lane's
constructions, since the punch described M7 by its values only; M7 here
is 91.13 in all three, where the listed order would wipe the first
bolus. In M2, M7 and the repeat case the three programs agree with each
other and disagree with NONMEM's listed order, which is why agreement
with rxode2 is not the test here.

300 random schedules on a 6-hour grid (seeds 20260925 to 20261224; 1 to
3 ids, 1 to 5 events each, any of bolus, infusion, ss, evid 3, 4, 5, 6,
with `addl` on 30 percent of doses): 132 accepted, 156 refused by
`frm_ode_records()` (155 shared time, 1 observation after an event, 156
of 156 classed), 12 refused by `frm_ode()` itself for two `ss` rows in
one group. On the 132 accepted, `frm_ode()` against rxSolve is 1.04e-12
of the scale at worst, median 7.2e-15 (7 of them have every sampled
amount zero and are compared absolutely); none goes below zero in
rxode2 (lowest -8.2e-14), while 15 of the 156 refused do. The round-0
continuous-time run and both real datasets were rerun on the new build
with the same numbers as before.

**The other holes, all refused by name now** (tests in
`test-records.R`): `rate`, `ii`, `addl` or `ss` nonzero on an
observation; `amt`, `ii`, `addl` or `ss` nonzero on an `evid = 3`
record; a non-finite `rate`, `dur`, `ii`, `addl` or `ss` on an event
(`addl = Inf` used to reach a base coercion warning and be solved as
0); a non-integer `cmt` (1.5 was solved as state 1); time going
backwards within an id (the review's restart gave 116.74 against
58.37); an id split into two blocks (it was merged into one subject);
and columns named as a dose modifier `frm_ode()` does not apply:
PREDPP's `ALAGn`, `Fn`, `Rn`, `Dn`, `MTIMEn` and `XSCALE`/`TSCALE`, and
`lag`, `alag`, `tlag`, `tinf` with or without a number, which is where
the review's `lag`, `alag1`, `tinf` and `f1` fall. PREDPP's `Sn` and
`FO` are not refused: they scale or split the output, not a dose. The
PREDPP names come from its documented list of additional PK parameters
(`Sn, Fn, Rn, Dn, ALAGn, XSCALE, MTIME, MTDIFF, MNEXT, MPAST, FO`, NONMEM
Help, Appendix III); MTDIFF, MNEXT and MPAST are indicators, not data.
A covariate such as `WT` passes. Doses of an id with no observation now
raise a classed warning, which names the case where that is every dose
and `events` comes back NULL.

The boundary fix, `st <= tm[r]` in the reset-inside-infusion test, is
not reachable on its own through the function: an infusion that starts
at a reset's time is always an event sharing that time, which the new
rule refuses first. It is kept so that the loop is right by itself.

**Seen failing.** The new and changed `test-records.R` run against the
ROUND-0 `frm_ode_records()`, recovered from the round-0 check tarball
and put into the lane namespace (`dev/phase3a-round1-records.R`,
`phase3a-p1-records-round1.txt`): 66 assertions fail and 0 error. Every
new block fails there; the passes in them are the controls (two plain
doses, a covariate column, an `events` NULL that was already NULL). On
the released build, which has no `frm_ode_records()`, the file fails
156 of 160 with 0 errors; on the lane build it passes 160 of 160.

### Review of punch round 2: the last holes

The review returned MERGEABLE with one gap and records to close.

* **An item argument set to `NULL` while `d` has its column** is
  refused by name: `rate = NULL` with a `rate` column used to read an
  infusion as a bolus, `addl = NULL` dropped the repeats, `ss = NULL`
  the steady state. The default name is matched without regard to
  case, as the other name checks are. The Rd now says plainly that a
  dose item is read only through its argument and any other column is
  a covariate; a rate in a column called `RATE_MG_H` that is not named
  with `rate =` is data, and cannot be caught by name.
* **`F0` and `FO` pass alike.** Both are PREDPP's output fraction, which
  acts on an observation, not a dose, as `Sn` does. `Fn` for `n >= 1`,
  the bioavailability of compartment `n`, is still refused.
* **`CENS`, `LIMIT`, `BLQ` and `LLOQ` are refused**: rxode2 and nlmixr2
  read the first two as censoring of an observation, and the other two
  are the usual names for it; as data a censored concentration would
  be fitted as observed.
* **Times that differ only by rounding are one instant.** A repeat at
  `0 + 3 * 0.1` (0.30000000000000004) meets a reset written at 0.3 and
  is refused; a reset at 0.35 still passes. The first version of this
  used a tolerance of 1e-9 of the table's largest time, which the
  coordinator showed was a million times too wide: a table 1 ms apart at
  `t = 1e7` was refused, and at epoch seconds (1.7e9) a sample 1 s after
  its dose would be; and because it compared neighbors after sorting,
  closely spaced records chained into one instant. Both over-refused and
  neither answered wrong, since the clusters feed only refusals.

  **The bound now follows the rounding it guards against.** A repeat
  `t0 + k * ii` is one multiply and one add. Measured against the exact
  decimal instant over 183062 draws (`dev/phase3a-tol.R`, seed
  20260924; |t0| up to 2e9, ii up to 1e5, k up to 1000, 0 to 6
  decimals), the error is at most 1.96 `eps * max(|t0|, k * ii)` (99.99th
  percentile 1.83, median 0) and does NOT grow with `k`. The tolerance
  is `64 * eps` of that scale for each point (of `|t|` for a written
  time), a 32.7-fold margin over the worst draw: 1.4e-14 at `t = 1`,
  1.4e-7 at `t = 1e7`, 2.4e-5 at `t = 1.7e9`. Two points are one instant
  when they differ by at most the larger of their two tolerances, and a
  cluster is measured from its first point, so it cannot chain.

  Checked: `0 + 3 * 0.1` against 0.3 is refused; the 1e7 table with 1 ms
  spacing and a 1.7e9 table with 1 s spacing are accepted; a reset 1 ms
  after a dose at 1e7 is its own instant and passes. The 300 grid-time
  schedules give the same split as before (132 accepted, worst 1.04e-12
  against rxSolve; 156 refused by `frm_ode_records()`, all classed; 12
  refused by `frm_ode()`). The new test fails 3 assertions, 0 errors, on
  the build with the 1e-9 tolerance (`phase3a-p4-records-prev.txt`), and
  `test-records.R` passes 204 of 204 on the lane build.

**Seen failing.** The new tests run against the previous
`frm_ode_records()`, recovered from the check tarball built for punch
round 1 (`dev/phase3a-round2-records.R`,
`phase3a-p3-records-prev.txt`): 30 assertions fail, 0 errors (NULL
items 19, F0 1, censoring columns 8, floating-point repeat 2). On the
lane build `test-records.R` passes 199 of 199.

## 3.2: epochs, groups and pairs

### Epochs of unequal length

`frm_cross_spectrum(x, y, group = )` with two lists: the epochs sharing
a `group` label are one unit, read as the clean spans of one record, so
the segment length is the unit's usable count divided by `segments` and
no segment crosses an epoch boundary. A vector pair with a `group` of
one label per sample is `frmtmb::frm_periodogram()`'s spelling: each
label one record. The single-record path is unchanged: a vector pair
with no `group` goes through the same code with one epoch, and
`test-cross-spectrum.R` gives identical results on both builds.

Checked in `test-epochs.R`, all `identical()`:

* one epoch per unit, four unequal lengths (512, 700, 333, 901): the
  list equals the vector with a per-sample group, and each unit equals
  the single-record call on its epoch;
* two units of two 512-sample epochs at `segments = 4`, where the
  segment divides the epoch: the list equals the concatenated-and-split
  call. This is the row's "when lengths agree";
* one unit of 700, 1043 and 699 samples at `segments` 4, 8 and 16: the
  list equals the single-record call with an `NA` between epochs, the
  gap machinery already shipped; and it differs from the plain
  concatenation, whose segments cross the boundaries.

Mutation (`phase3a-mutate32.R`): with the epochs concatenated instead
of kept apart, the third block fails and the other four pass, which is
what the "lengths agree" identity says they should.

### Pairs

`frm_cross_pairs(X, pairs, ...)`: `X` a matrix or data frame whose
columns are channels, or a list of such epochs; `pairs` all pairs by
default, or a two-column matrix by name or number. Each block is
`frm_cross_spectrum()` on that pair, `identical()` in the test, with
`pair`, `ch1` and `ch2` in front. Refused by name: a channel paired with
itself, a pair repeated in either order, an unknown name or number, an
unnamed or unknown argument in `...` (so a misspelled `segmnets` is not
swallowed), one channel, a non-numeric column, epochs with different
channels, duplicate channel names.

`dev/phase3a-pairs-fit.R`, seed 45, 8192 samples, four channels with
loadings 1, 0.8, 0.5 and 0.25 on one white source, 16 segments, pairs
Fz-Cz, Fz-Pz, Cz-Pz and Pz-Oz, 1020 rows:

| pair | truth | separate fit | every parameter by pair | `coh ~ 1 + (1 \| pair)` |
|---|---|---|---|---|
| Fz-Cz | 0.1951 | 0.19320 | 0.19320 | 0.19280 |
| Fz-Pz | 0.1000 | 0.09299 | 0.09299 | 0.09285 |
| Cz-Pz | 0.0780 | 0.07409 | 0.07409 | 0.07404 |
| Pz-Oz | 0.0118 | 0.01073 | 0.01073 | 0.01140 |

On the logit scale, where `frm_coherence()`'s `.se` lives, the
by-pair fit is 4.5e-04 separate-fit standard errors from the separate
fits at worst and its standard errors agree to 7.2e-05 relative; that
is the identity the likelihood's sum structure predicts, to the
optimizer's tolerance. The `(1 | pair)` fit is 0.29 standard errors
away at worst, keeps the order of the pairs, and has `sd(coh | pair)`
1.13, so the pairs differ by far more than the data blur them and the
shrinkage is small. Both stacked fits end at code 0 with a gradient of
3.1e-03 and 2.1e-03, above core's 1e-03 warning line, at `nlminb`'s
relative tolerance on an objective of about 12,000; the four separate
fits end below 5e-04. The test lets that one warning through and fails
on any other.

## Nothing else moved

Every test file of the three packages, one file per R process
(`dev/phase3a-suite.sh`, `dev/phase3a-runtest.R`), in three arms: A,
the base build on the base tests (`git archive HEAD`); B, the lane
build on the same base tests; C, the lane build on the lane tests. The
comparison is one row per `test_that()` block on assertion count,
failures, errors, skips and warnings (`dev/phase3a-compare.R`,
`dev/phase3a-log/phase3a-compare.txt`), pasted verbatim:

    <!-- generated by dev/phase3a-compare.R -->

    ## frmtmb.ode
    A base build, base tests:  files 10, blocks 118, pass 341, fail 0, error 1, skip 1, warn 0
    B lane build, base tests:  files 10, blocks 118, pass 341, fail 0, error 1, skip 1, warn 0
    C lane build, lane tests:  files 11, blocks 138, pass 540, fail 0, error 1, skip 1, warn 0
    blocks that differ between A and B: 0

    ## frmtmb.coupling
    A base build, base tests:  files 9, blocks 66, pass 449, fail 0, error 0, skip 5, warn 3
    B lane build, base tests:  files 9, blocks 66, pass 449, fail 0, error 0, skip 5, warn 3
    C lane build, lane tests:  files 11, blocks 76, pass 542, fail 0, error 0, skip 5, warn 3
    blocks that differ between A and B: 0

    ## frmtmb.spline
    A base build, base tests:  files 14, blocks 93, pass 450, fail 0, error 0, skip 1, warn 0
    B lane build, base tests:  files 14, blocks 93, pass 448, fail 2, error 0, skip 1, warn 0
    C lane build, lane tests:  files 14, blocks 99, pass 491, fail 0, error 0, skip 1, warn 0
    blocks that differ between A and B: 2
                                                                                    key
     test-frailty.R :: a random effect on gamma1 is a per-centre shape, and it recovers
     test-rp-floored.R :: fitted() still refuses, and the count is named for what it is
     passed.a failed.a passed.b failed.b error.a error.b warning.a warning.b
           15        0       14        1   FALSE   FALSE         0         0
            6        0        5        1   FALSE   FALSE         0         0
    <!-- end generated -->

frmtmb.ode was rerun after the review of punch round 2 and frmtmb.spline
after punch round 2 (its later change is roxygen text only, which the R
CMD check below postdates); frmtmb.coupling's
code has not changed since round 0.

A and B are identical in every block of frmtmb.ode and
frmtmb.coupling. In frmtmb.spline exactly two blocks differ, and both
are the intended change: the `expect_named()` that lacks the new field,
and the frailty assertions item 2.5 wrote to fail when 3.6 landed
(`expect_silent(rp_floored(bad2))` now meets the warning; in round 1
it met a refusal and errored). Both are updated in the lane's tests.

The one error in frmtmb.ode, arms A, B and C alike, is
`test-ode.R`'s "a missing backend names the r-universe repository":
`local_mocked_bindings()` needs a package loaded by pkgload, which a
bare `test_file()` does not provide. It is the runner's, and it passes
under `R CMD check` below.

Arm D, the lane's new and changed test files on the BASE build
(`dev/phase3a-log/phase3a-armD.txt`, `seefail/`): `test-records.R`
fails 87 assertions, `test-epochs.R` 33, `test-cross-pairs.R` 50,
`test-rp-floored.R` 10 and `test-frailty.R` 2, with 0 errors in any of
them, so every assertion after the first ran and was seen to fail.

## R CMD check

`dev/phase3a-check.sh`: `R CMD build` with vignettes, then
`R CMD check --as-cran`, `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`,
against the lane library. Logs, `00check.log` and `testthat.Rout` are
copied to `dev/phase3a-log/check-<pkg>*`; each build log postdates
every file of its package.

| package | Status | in-check tests | NOTEs |
|---|---|---|---|
| frmtmb.ode (tolerance fix, 07:24) | 1 NOTE | FAIL 0, WARN 0, SKIP 2, PASS 543 | V8 math rendering |
| frmtmb.coupling | 1 NOTE | FAIL 0, WARN 3, SKIP 7, PASS 534 | examples timing (`cross_wishart` 2.29 s user against 9.91 s elapsed) |
| frmtmb.spline (review of punch round 2, 07:00) | 2 NOTEs | FAIL 0, WARN 0, SKIP 3, PASS 483 | examples timing (`frm_curve` 3.12 s user against 6.2 s elapsed, under the load of two checks and other lanes); V8 math rendering |

The timing NOTEs are on examples this lane did not touch, with user
time well under 5 s and elapsed several times that: the machine was
running about 30 R processes from the other lanes, which is the load
`dev/lane-rules.md` says this NOTE measures. The coupling WARN 3 are
`test-coherence.R`'s own convergence warnings, present identically in
arms A and B. The first frmtmb.ode build failed in the vignette; see
the prediction item under "Filed".

## What I did not do, and what is filed

* **Known over-refusal, left as is: a replace or multiply at the same
  time as an event on a DIFFERENT compartment.** The two commute, so the
  records do say one thing, but the same-time rule does not look at
  compartments and refuses them. It errs on the safe side; the Rd says
  so.
* **Unverified: an observation at an `addl` repeat time.** `frm_ode()`
  reads it before the implied dose, the trough. Whether NONMEM does
  depends on how PREDPP orders an implied dose against an observation
  record at the same time, which neither this lane nor the reviewer
  could check (no PREDPP guide was reachable). It is not refused.

* **Filed, not fixed: `frm_ode()` runs an infusion on through a
  reset.** A `"reset"` event sets every state to its value and leaves
  any infusion in progress running (`R/ode.R`, the rate loop reads
  every infusion with `time <= a < end` whatever reset came between).
  `frm_lincmt()` refuses the same schedule. Whether a reset should end
  an infusion is a semantic decision for `frm_ode()`, beyond this item.
  Round 0 said `frm_ode_records()` refuses the case "so it cannot be
  reached through records". **That was false**, and the review showed
  it: the round-0 test `st < tm[r]` was strict, so an infusion that
  starts at the reset's own time, an `addl` repeat included, passed.
  The reviewer's M4 gave 40.08, 72.65 and 14.52 in `frm_ode()` against
  0, -22.33 and -91.43 in rxode2, and 4 of the reviewer's 321 accepted
  schedules went negative in rxode2, down to -537. Punch round 1
  makes the test `st <= tm[r]`, and refuses any reset that shares a
  time with another event, which is the only way an infusion can start
  at a reset's time. On 300 grid-time schedules, 0 of the 132 accepted
  go below zero in rxode2 (lowest -8.2e-14), while 15 of the 156
  refused do. This is a count over one generator, not a proof that no
  accepted table reaches the path.
* **Filed, not fixed: a per-subject dosing table blocks prediction on a
  subset of subjects.** `frm_ode()` refuses an `events$group` that the
  data do not have, and on `newdata` that is every subject left out:
  the vignette's `frm_linpred(fit, newdata = <one subject>)` failed the
  first vignette build with "`events$group` names 11 groups that are
  not in `group`". The vignette now predicts for every subject and
  shows one. Dropping absent groups' events at predict time would be
  the fix, in `frm_ode()`, beyond this item.
* **Filed: the Theoph records fit at the default tolerance** ends at
  `nlminb` code 1 where the initial-condition spelling ends at code 0,
  6.3e-09 apart in log likelihood. The vignette tightens the
  integrator. The default itself is item 1.0d's territory.
* **Not measured: NONMEM itself.** No NONMEM is available here. Its
  record order rule and its reading of `ss`, `evid = 3` and
  `evid = 4` are taken from its documentation; rxode2 and nlmixr2data
  are the executable references.
* **Not done: an estimated `pcmt`, `date` handling, `evid = 2`
  covariate records, `ss = 2`, and rxode2's classic codes.** Each is
  refused by name, with what to do instead.
* **Not searched: a dip between observed times** of a group with no
  events at `df >= 2` (3.6), stated on both help pages.
* **Not modeled: the dependence between two pairs that share a
  channel** (3.2), stated in `?frm_cross_pairs` and the vignette.
* The `_pkgdown.yml` of frmtmb.ode and frmtmb.coupling had no
  `reference:` index, so pkgdown listed every export on its own. Each
  now has one naming every exported topic. `pkgdown::check_pkgdown()`
  could not run: it stops at "Pandoc not available" in this shell.

### Last item: times in refusal messages

A refusal named `format()`'s default seven digits, so a reset at
10000000.3 was reported as "shares time 1e+07". Every time in a
`frm_ode_records()` message now goes through `ode_rec_time()`, 15
significant digits: 10000000.3 prints as 10000000.3, 1700000001 as
1700000001, and `0 + 3 * 0.1` as 0.3 rather than 0.30000000000000004.
The large-time assertion in `test-records.R` now matches
"shares time 10000000.3"; it failed on the previous build
(`phase3a-p5-records-prev.txt`, 1 failure, 0 errors) and the file
passes 204 of 204 on the lane build (`phase3a-p5-records-lane.txt`).
No Rd changed, so the frmtmb.ode R CMD check at 07:24 was not rerun;
the change since it is `R/records.R` and `test-records.R` only.
