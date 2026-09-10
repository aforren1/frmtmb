# `frm_ode(n_ss =)`: the steady-state run-in

Item 1.0d of `dev/extension-gaps-plan.md`. Lane `nss`, worktree
`frmtmb-wt-nss` off `f9ee297`, frmtmb 0.55.2, frmtmb.ode 0.3.0,
RTMB 1.9, RTMBode 1.0 (r-universe), deSolve 1.42, R 4.6.1,
StanHeaders 2.32.10 from the pin.

The BEFORE arm is the round's shared reference build,
`C:/Users/adf44/source/r/rellib-0552`. The lane's own library is
`C:/Users/adf44/source/r/nsslib`; nothing was installed anywhere else,
except the eight throwaway variant builds of section 6, which went to
`C:/Users/adf44/source/r/nsslib-variant`.

Every script is in `extensions/frmtmb.ode/dev/nss/`. Each number below
names the script that produced it and, where one applies, the seed.

## What punch round 2 changed

Round 1 of `dev/reviews/2026-09-10-nss.md` returned **mergeable** with
three corrections, all of which were mine to make.

| finding | what was wrong | what it says now |
|---|---|---|
| NIT 9 | the C1 claim was wrong in both directions, and the low gate was **not C1 at all**: one-sided derivatives 0.05540 and 1.10804, a gap of `1 / (1 - rlow)` constant over four decades | the gate is the smootherstep `S` and not `S(x)/x`; no junction shows a corner, and the residual gap at the two curved junctions IS the factor's own second derivative (section 4a) |
| NIT 9 | the same shape made the gate OVERSHOOT the tail by 1.2487x below `rlow`, and that overshoot was the whole of the "agree to 9.3e-03" figure I attributed to the smoothing | the gate never exceeds the tail: largest ratio 1.0000, attained at `rlow` |
| NIT 10 | the tape figures were stale by half again | +258 and 255 nodes, re-counted (section 4b) |
| NIT 11 | the warning fell silent at `n_ss` of 650, 800 and 1000 on a 1653 h half-life, which is exactly where `?frm_ode` sends such a user | the report adds back the part of the tail the stand-down declined; those rows now warn at 1.01x to 1.02x of the true error (section 5) |

**And the lane's own see-it-fail harness failed open**, which belongs
in the record because it is the failure mode this project keeps
finding. After the shape changed under it, `nss-15-variant.py` matched
nothing, patched nothing, and all eight variants installed and tested
the SHIPPED code, reporting a clean `pass=142 fail=0` eight times. It
was caught only because eight identical clean rows are not a plausible
result. The script now cuts every anchor out of the current source and
raises when a patch is a no-op; the shell refuses to run a variant
whose `R/ode.R` still compares equal to the shipped one, and deletes
the previous log first, because a skipped variant had been leaving a
stale log that read as a pass.

## What punch round 1 changed

Reviewed in `dev/reviews/2026-09-10-nss.md`, which reproduced every
headline number here, several to four or five figures, confirmed the
record correction of section 1, and found one BLOCKER and eight NITs.
The reviewer's scripts are in `dev/rev-nss/` and were not edited.

| finding | what was wrong | what it says now |
|---|---|---|
| BLOCKER | the rewritten warning could not fire on `ss_extrapolate = FALSE`, because every guard drives the correction to exactly zero and the report was built from the correction; `n_ss` = 1 was silent on both arms | every branch falls back to the cycle-to-cycle movement, `stalled` is out of the `rel` gate, and all eight cells of the reviewer's table warn (section 5) |
| NIT 6 | the correction was built before the `isTRUE()` branch, so BOTH arms carried +147 tape nodes | built inside it; the truncated arm is the base commit node for node, the default arm costs +258 (section 4b) |
| NIT 4 | the objective was C0 and not C1 at the stand-down, which is the property used to reject the alternative design | a degree-7 shape, `ode_ss_rcap` moved to 0.9875 to hold the bound; round 2 then corrected the smoothness claim itself (section 4a) |
| NIT 1 | "4.2e-01 against 3.9e-03", two different quantities, called 108x | 4.96x like for like, direction unchanged (section 2) |
| NIT 2 | "auto" rejected for being worse than the status quo, which a floor repairs | rejected on dominance: the same error for 6.4x the solves (section 2) |
| NIT 3 | the only losing case named was "an oscillator, not pharmacokinetics" | flip-flop absorption, 8 of 338, worst 2.27x, with the bound (section 8) |
| NIT 5 | "manufactures a precision" read as general | true of `k21`, false of `ke`, which goes the other way in all twelve fits (section 7) |
| NIT 7 | "the integrator's tolerance" | amplified by about `1 / (1 - r)`, up to 62.8x (section 4) |
| NIT 8 | two shipped tables gave two numbers for one row | the Rd says which is which (section 4) |

Two shipped sentences were wrong as written and are corrected: NEWS
said `ss_extrapolate = FALSE` restores the run-in "bit for bit", which
was true of values and false of the diagnostic and of the tape, and it
said "no miss in either arm", which the blocker made false.

## What was shipped, in one paragraph

The reviewer's `n_ss = "auto"` was **measured and rejected**: chosen at
the starting values it is WORSE than today's fixed 20 on the start a
user gets when they pass no `start` at all. What is shipped instead is
`ss_extrapolate = TRUE`: the run-in's own successive differences give
the per-cycle contraction, so the geometric tail that truncation drops
is summed rather than dropped. It is arithmetic on tape variables, so
it runs DURING a fit at whatever parameters the fit has reached, which
is the property "auto" cannot have, and it costs **no extra solve**.
The numeric-path warning was rewritten to report a distance to the
limit rather than the cycle-to-cycle movement.

## 1. Reproducing the record first

`dev/nss/nss-09-verify.R` on the reference build, section A, against
`frm_lincmt()` at its default `n_ss = Inf`, which is the exact limit:

| terminal half-life | `ii` | measured here | `dev/reviews/2026-09-09-lincmt.md` |
|---|---|---|---|
| 23.2 h | 8 | 3.939e-03 | 3.94e-03 |
| 107.1 h | 24 | 1.182e-02 | 1.18e-02 |
| 264.6 h | 24 | 9.024e-02 | 9.02e-02 |
| 669.9 h | 24 | 2.092e-01 | 2.09e-01 |

`dev/nss/nss-13-warntext.R` on the reference build reproduces the
warning's understatement to two figures: 0.79x, 3.08x, 6.56x and
10.81x against the recorded 0.8x, 3.1x, 6.6x and 10.8x.

**One correction to the record.** The review's section C prose says the
`k21` factor of 3.16 was measured "terminal half-life 107 h, `ii` = 24",
and the plan's 1.0d row and this lane's brief repeat it. The script it
names, `dev/rev-lincmt-nss.R`, calls
`run("slow peripheral, ii = 24", 0.1, 0.2, 0.008, 1.0, 24, 10, 101)`,
which is `ke = 0.1, k12 = 0.2, k21 = 0.008`; that is `lambda_z` =
0.00262 and a terminal half-life of **264.6 hours**, the table's next
row down. The 107 hour row is `ke = 0.15, k12 = 0.3, k21 = 0.02`. The
`k21` factor of 3.16 therefore belongs to the 264.6 hour design.
Section 7 fits both: 107 hours, because that is what the acceptance
criterion names, and 264.6 hours, because that is where the recorded
number came from and where it reproduces to five figures.

## 2. The reviewer's proposal, and why it is not shipped

`n_ss = "auto"` would run the run-in numerically at tape-build time,
read `r` off the last differences, and take
`n_ss = ceiling(log(1e-9) / log(max r))`.

**The estimator itself is fine.** `dev/nss/nss-03-ratio.R` section A,
on exact matrix exponentials so that nothing is confounded with solver
error: read from three iterates, five, eight or twenty, the ratio
agrees with the dominant eigenvalue of the cycle map to six figures on
every row of the compartment grid, and the `n_ss` it implies at 1e-9 is
identical from three iterates and from twenty on seven of eight rows
and differs by one cycle on the eighth.

**The stated unmeasured risk is what kills it.**
`dev/nss/nss-08-rmove.R`, three designs x three starting points x six
seeds = 54 fits, seeds 101 to 106, `frm_lincmt()` as the vehicle
because `r` is a function of the parameters alone and the closed form
reaches the same optimum as `frm_ode()` to seven significant digits.
`n_start` is what "auto" would choose where the tape is built; `n_grp`
is what the worst SUBJECT needed at the optimum, which is the run-in
that is actually performed; `short@n_s` is the shortfall "auto" would
have shipped.

| design | start | `n_start` | `n_grp` at the optimum | shortfall shipped |
|---|---|---|---|---|
| 107 h, `ii` 24 | at the truth | 134 | 100 to 337 | 8.7e-13 to 2.6e-04 |
| 107 h, `ii` 24 | a plausible guess | 23 | 100 to 337 | 8.5e-03 to 2.4e-01 |
| 107 h, `ii` 24 | frmtmb's cold start | 3 | 10 to 11 | 1.9e-03 to 3.1e-03 |
| 265 h, `ii` 24 | at the truth | 330 | 150 to 6554, one none | 1.4e-20 to 1.0 |
| 265 h, `ii` 24 | a plausible guess | 23 | 150 to 6554, one none | 4.1e-02 to 1.0 |
| 23 h, `ii` 8 | frmtmb's cold start | 7 | 72 to 169 | 1.3e-01 to 4.2e-01 |

Read the last row. On the schedule the package's own examples use,
started where `frm()` starts a nonlinear model when the user passes no
`start` (every beta at 0, `R/fit.R` `make_start()`), "auto" picks 7
cycles where the worst subject needs 169.

**Punch round 1 corrected the size of that, against me.** Round 0 wrote
"4.2e-01 where today's blind 20 ships 3.9e-03" and called it a
hundredfold regression. Those are two different quantities: 4.2e-01 is
the theoretical shortfall at the worst SUBJECT at the optimum, and
3.9e-03 is a measured trajectory error at the POPULATION truth. Like
for like on the same subject (`dev/rev-nss/rev-nss-01-auto.R` B), the
blind 20 ships **8.50e-02** on that subject and "auto" ships 4.22e-01,
so the regression on the seed quoted is **4.96x, not 108x**. Across the
six seeds the ratio runs 4.96 to 42.1 and "auto" is worse on every one,
so the direction survives and only the magnitude was overstated. A
number that overstates in its author's own favour is the thing this
project does not do.

Even started AT THE TRUTH, which no real fit is, "auto" targets 1e-9
and lands between 8.7e-13 and 3.5e-01 on the five 265 h datasets whose
optimum has a steady state, because the population value at the start
is not the worst subject at the optimum: every random effect is 0 when
the tape is built. The sixth (seed 103) reached `r = 1.00000` at its
optimum, where no finite `n_ss` is enough at all and the shortfall is
1.

### The ground the rejection actually stands on

**Round 0's stated ground does not survive.** `n_ss = max(20, auto)` is
by construction never worse than today, and round 0 did not try a
floor: it is identical to today on the 23 h cold rows and 8.67e-13 to
2.62e-04 at 107 h started at the truth. So "worse than the status quo"
is a statement about an unfloored rule and does not decide the item.

**The rejection stands on DOMINANCE instead**, which the floor cannot
repair. `dev/nss/nss-25-dominance.R`, 107.1 h terminal half-life,
`ii` 24, `atol = rtol = 1e-8`, error worst over one dosing interval
against the exact limit, `deSolve::lsoda` entries counted per group:

| arm | error | `lsoda` |
|---|---|---|
| today, `n_ss` = 20, truncated | 1.182e-02 | 21 |
| floored auto, `n_ss` = 134, truncated | 2.854e-09 | 135 |
| the same at `n_ss` = 266 | 2.618e-09 | 267 |
| shipped, `n_ss` = 20, tail summed | **2.623e-09** | **21** |

They reach the same place. 2.6e-09 is the integrator's own floor at
this tolerance, which the `n_ss` = 266 row shows, and the correction
reaches it for **6.4x fewer solves**. A floored "auto" is dominated on
both axes at once, and no choice of floor changes that, because the
floor only ever raises `n_ss`.

Two things a floor cannot fix either. At the 265 h design, seed 103,
the optimum has `r` = 1.00000, where no finite `n_ss` is enough at all;
and "auto" is fixed where the tape is built while the correction is
recomputed at whatever parameters the fit has reached.

**Recomputing `n_ss` at every objective evaluation** was considered and
is not possible in this design: `n_ss` sets the NUMBER OF SOLVES, which
is the tape's structure, and the tape is built once. A tape whose shape
depended on a value would have to be rebuilt inside the optimizer, and
the objective would then be discontinuous in the parameters wherever
the chosen `n_ss` changed, which the Laplace approximation's inner
Newton solve cannot carry. That argument put a burden on the shipped
design too, and section 4a is where it is discharged.

## 3. What is shipped

`ss_extrapolate = TRUE`. Three successive cycle-start states give
`d1 = y[n-1] - y[n-2]` and `d2 = y[n] - y[n-1]`; the per-cycle
contraction is `r`, so what truncation drops is `d2 * r / (1 - r)`, and
that is added. This is Aitken's extrapolation and for linear kinetics
it is exact up to the second-slowest mode.

Three properties are not optional, and each has a construction that
breaks the rule without it (`dev/nss/nss-04-hazards.R`,
`nss-05-guard.R`, `nss-06-rule.R`):

1. **the ratio is read per state.** A state with no steady state at
   all, an AUC compartment being the ordinary example, has a ratio of
   1. Share one ratio across states and that state sets it: on a
   one-compartment oral model at `ke * ii` = 0.12 with an AUC state
   riding along, one shared ratio leaves the CENTRAL compartment at
   9.072e-02, exactly where truncation leaves it, against 4.735e-15
   read per state (`nss-14-absent.R`). With a hard cap and no
   stand-down the same shape is **2.4e+01**, a 2400 percent error,
   against truncation's 3.97e-14 (`nss-04-hazards.R` H2).
2. **the denominator carries the integrator's noise floor.** An oral
   depot empties to the same trough every cycle once the run-in has
   settled, so its successive differences are zero to the last bit and
   an undamped ratio is 0 / 0. Measured directly: `d1 = (0, 9.16)`,
   `d2 = (0, 8.13)`, undamped `r = (NaN, 0.887)`, undamped answer
   `(NaN, 792.3)`, damped answer `(6.1e-04, 792.3)`. This is not a
   corner: EVERY oral model reaches it.
3. **the correction stands down as the ratio reaches 1.** A ratio above
   1 says the component is moving further each cycle, and capping it at
   the top of its range reads that as "nearly stopped", which is the
   worst reading available because it multiplies the last difference by
   `1 / (1 - cap)`. On an already-converged system, `ke * ii = 48`, a
   cap turns a 0.000e+00 error into **7.879e-01**
   (`nss-05-guard.R`); on a driven oscillator it turns 4.089e-01 into
   3.596e+01 (`nss-06-rule.R`).

The rule sweep in `nss-06-rule.R`, nineteen systems, scores the shipped
rule against truncation on its WORST row, because a rule is only as
good as that: worst row 5.677e-01 against truncation's 6.087e-01,
worst ratio to truncation 3.19, median ratio 0.065. `ode_ss_noise = 4`
comes from that sweep.

**Punch round 1 added a fourth property and moved a constant.** The
stand-down was a linear ramp against a hard cap, which is C0 and not
C1; it is now a degree-7 shape whose first three derivatives vanish at
both ends, and the low end is gated by the same shape rather than
clamped. Section 4a has the measurement and what it costs.
`ode_ss_rcap` moved from 0.99 to **0.9875** as part of that: the smooth
shape is wider than the ramp it replaced, so at 0.99 the largest factor
the correction could apply would have risen from 99.0 to 124.0. At
0.9875 that bound is **99.02**, which is the bound the hard cap gave,
99.00 (`dev/nss/nss-22-kink.R`).

## 4. What it buys, through `frm_ode()` itself

`dev/nss/nss-09-verify.R` section A, worst over one dosing interval
against `frm_lincmt()` at `n_ss = Inf`, at the shipped `n_ss = 20` and
the shipped `atol = rtol = 1e-8`:

| terminal half-life | `ii` | truncated | extrapolated | gain |
|---|---|---|---|---|
| 23.2 h | 8 | 3.939e-03 | 9.879e-09 | 4.0e+05 |
| 107.1 h | 24 | 1.182e-02 | 2.623e-09 | 4.5e+06 |
| 107.1 h | 12 | 8.763e-02 | 8.935e-09 | 9.8e+06 |
| 264.6 h | 24 | 9.024e-02 | 5.733e-09 | 1.6e+07 |
| 669.9 h | 24 | 2.092e-01 | 1.351e-08 | 1.6e+07 |

The extrapolated column is the INTEGRATOR's tolerance, not the run-in's:
`nss-03-ratio.R` section C runs the same rows on exact matrix
exponentials and gets 1.2e-16 to 1.1e-14.

**Round 0 said that too broadly.** The correction multiplies the last
difference by up to 99, and that difference carries the solver's own
error, so what is left is the integrator's tolerance **amplified by
about `1 / (1 - r)`**. In units of `atol = rtol`
(`dev/rev-nss/rev-nss-11-amplify.R`): 0.3 to 1.4 on the five rows
above, 21.6 and 67.9 at a 670 hour half-life at 1e-6 and 1e-10, and
**62.8** at a 1655 hour half-life at 1e-8. The sentence is right on the
published rows and wrong as a general statement, and `?frm_ode` now
carries the amplification.

**Two shipped tables gave two numbers for the same row.** `?frm_ode`
gives 107 h, `ii` 24 as 4.5e-02 truncated, worst over the run-in
STATES; NEWS.md gives 1.2e-02, worst over one dosing interval on the
observable. Both were labelled and neither pointed at the other. The Rd
now says which is which and names the other number.

**The gradient is what makes this a wrong answer during a fit rather
than a slow one.** `dev/nss/nss-11-grad.R`, two-compartment oral,
107 hour half-life, `ii` = 24, `n_ss = 20`, `atol = rtol = 1e-10`,
against the exact steady state's derivative:

| | value vs exact | AD vs its own numeric path | AD vs the EXACT limit | `d/dlog(k21)` |
|---|---|---|---|---|
| truncated | 2.45e-02 | 1.38e-07 | 9.47e-02 | **+10.372** |
| extrapolated | 2.10e-11 | 1.11e-05 | 2.70e-10 | **-1.302** |

The truncated tape differentiates its own answer correctly, to 1.4e-07,
and that answer's derivative with respect to `log(k21)` **has the wrong
sign**. An optimizer handed +10.37 where the truth is -1.30 is pushed
the wrong way in `k21`, which is how a fit at `n_ss = 20` moved `k21`
by a factor of three. The extrapolated arm's larger "AD vs its own
numeric path" figure is the finite difference, not the tape: at
`h = 1e-5` with solver noise at 1e-10 amplified by `1 / (1 - r)` = 6.9,
the difference quotient carries about 7e-05 of its own, and the same
gradient matches the EXACT limit to 2.7e-10.

At `n_ss = 5` the truncated gradient is 5.34e-01 wrong and the
extrapolated one 2.64e-10.

## 4a. The objective's smoothness, which round 0 shipped unchecked

The reviewer found that the correction was continuous in the parameters
and its derivative was not, at the stand-down point, and ranked it a
NIT. It is not a NIT here, for a reason the reviewer names and I
accept: **it is the property round 0 used to reject the alternative
design.** Section 2 rejects recomputing `n_ss` inside the optimizer
because that would make the objective non-smooth; shipping a
non-smooth objective and then holding smoothness against the
alternative is not an argument, it is a double standard.

Measured, `dev/nss/nss-22-kink.R`, two-compartment oral, `ii` 24,
`n_ss` 20, `atol = rtol = 1e-12`, sweeping `k21` across the value that
puts `exp(-lambda_z * ii)` exactly at the stand-down point. The test is
not whether the objective jumps, which it never did, but whether the
two ONE-SIDED derivative limits converge as the bracket tightens:

| bracket | round 0, hard cap | round 1, smooth shape |
|---|---|---|
| eps 1e-3 | -0.0991 and +126.78 | -0.112456 and -0.112689 |
| eps 1e-4 | -0.0909 and +147.74 | -0.112566 and -0.112588 |
| eps 1e-5 | -0.0901 and +149.89 | -0.112576 and -0.112578 |
| eps 1e-6 | not measured | -0.112577 and -0.112577 |

Round 0's two limits were about -0.09 and +150 and did not converge:
a genuine derivative discontinuity, a sign flip and a factor of 1660.
Round 1's gap falls 2.33e-04, 2.23e-05, 2.23e-06, 2.23e-07, linearly
with the bracket, which is what a C1 function does. The same holds at
`r` = 1 (gap 3.8e-09 falling to 3.8e-12) and at the low gate
(7.9e-04 falling to 7.9e-07). `ss_extrapolate = FALSE` was smooth
before and still is.

**Is it reachable?** Yes, and the lane's own run proves it: at the
265 h design, seed 103, the exact-limit fit ran `lk21` from its start
at -4.83 to **-31.05**, and the path from one to the other crosses the
stand-down point. `r` is per group and depends on the random effects,
so the Laplace inner Newton solve sees it in `b` as well as the outer
optimizer seeing it in `beta`.

**What the smoothing costs.** A cubic smoothstep, which is what the
reviewer suggested, is C1 at the stand-down point and leaves the kink
at `r` = 1: near there the correction behaves like `3(1-r)/h^2` whose
derivative does not vanish. A degree-7 shape whose first three
derivatives vanish at both ends fixes both. It is wider than the ramp
it replaces, so the largest factor the correction can apply would have
gone from 99.0 to 124.0 at the old `ode_ss_rcap`; moving that constant
to 0.9875 puts the bound back at **99.02** against the hard cap's
99.00. The price paid for that is 0.0025 of the ratio range over which
the correction is exact, which at `ii` = 24 moves the boundary from a
1655 hour terminal half-life to 1322 hours.

### Round 2: the claim was wrong in both directions, and the probe was inert

Round 1's sentence, "C1 and not C2, the second derivative still jumps
at the two boundaries", was wrong three ways and the probe behind it
could not have seen the worst of it. The measurement has to be made on
the SHAPE, not on an objective built over it
(`dev/nss/nss-26-corners.R`, the reviewer's construction: with `ya` = 0,
`yb` = 1/r and `yc` = 1/r + 1 the damping is negligible and the ratio
the shipped function reads is `r` to 4.4e-16, so what comes back is the
factor itself). There are **four** junctions, not two:

| junction | round 1 gap at eps 1e-2, 1e-3, 1e-4, 1e-5 | round 2 |
|---|---|---|
| `r` = 0 | 1.68e-01, 2.67e-04, 2.73e-07, 0 | 3.37e-02, 5.34e-06, 0, 0 |
| `r` = rlow = 0.05 | 8.91e-01, **1.054, 1.053, 1.053** | 1.16e-01, 2.06e-03, 2.33e-04, 2.33e-05 |
| `r` = rcap = 0.9875 | 1.01e+04, 9.30e+02, 1.02e+02, 1.02e+01 | unchanged |
| `r` = 1 | 9.57e+03, 1.18e+03, 1.41e+01, 1.43e-01 | unchanged |

Round 1's low gate was **not C1 at all**: its two one-sided derivatives
were 0.05540 and 1.10804, a gap of 1.0526 that is exactly
`1 / (1 - rlow)` and did not fall over four decades. Round 1's own
probe swept `k21` through a solved model at `n_ss` = 20, where a ratio
of 0.05 leaves the run-in `0.05^20` = 1e-26 from its limit, so the
damping shut the gate on both sides and the corner was never in the
code path. **A probe that cannot reach the thing it certifies is not a
measurement**, and this is the second time in this lane that a check
was written where the case is absent rather than present but quiet.

**The cause, and the one-character repair.** `step(x)` is the degree-7
smootherstep `S` divided by `x`, so `step'(1)` is -1. That -1 is
exactly right for the STAND-DOWN, because the branch below the cap is
`r / (1 - r)` and -1 is the slope that joins it, which is why that
junction comes out smooth rather than merely continuous. It is exactly
wrong for a GATE, which has to saturate flat. Using `S` itself, `x^4`
in place of `x^3`, fixes it.

The same -1 made the gate **overshoot the geometric tail by 1.2487x**
below `rlow`, peaking at `r` = 0.035, and that overshoot was the whole
of the "the two shapes agree to 9.3e-03" figure round 1 attributed to
the smoothing: below the cap the two stand-down shapes are identical.
Measured now, largest factor over tail below `rlow` is **1.0000**,
attained at `rlow`, where a gate should reach it.

**How to read what is left.** A one-sided first difference carries an
error of `f''(r) * eps` of its own, so at a junction where the factor
has curvature the probe reads the curvature and cannot separate C1 from
C3. The control is the factor's own second derivative, `2 / (1 - r)^3`:

| junction | gap / eps at 1e-5 | `2 / (1 - r)^3` |
|---|---|---|
| `r` = rlow | 2.333 | 2.333 |
| `r` = rcap | 1.024e+06 | 1.024e+06 |

To four figures. So no junction shows a corner, `r` = 0 falls to
exactly 0 and `r` = 1 falls like the square of the bracket, and at the
two curved junctions there is no discontinuity left in the reading to
find. `?frm_ode` now says that rather than naming a class.

## 4b. The tape, which round 0 did not count

The claim "no extra solve" is verified and unchanged: `n_ss + 1`
`deSolve::lsoda` entries per group, 6 at `n_ss` = 5 and 21 at
`n_ss` = 20, on both arms and both builds. The tape is the other half
and round 0 said nothing about it. The reviewer counted it and found
**+147 nodes in BOTH arms**, because the correction was built before
the `isTRUE(ss_extrapolate)` branch, so the arm documented as
unchanged put the whole correction on the tape and discarded the value.

The call is now inside the branch. `dev/nss/nss-21-tape.R`, node counts
from the printed operation stack, which is load-independent:

| `n_ss` | base commit | lane, `FALSE` | lane, `TRUE` |
|---|---|---|---|
| 4 | 64 | **64** | 322 |
| 5 | 73 | **73** | 331 |
| 20 | 208 | **208** | 466 |
| 40 | 388 | **388** | 646 |

The truncated arm is the base commit node for node. The default arm
costs **+258 nodes**, a fixed cost independent of `n_ss`, which is
**+124 percent** of the base tape at `n_ss` = 20 on this design and
+66 percent at `n_ss` = 40. One call to the correction on three states
is **255** nodes measured on its own, against the **144** the reviewer
measured for the round-0 shape with the same probe, so **111** of them
are the price of section 4a's smooth shape. The solve nodes dominate
wall time, so this is small, but it is a real cost on a tape a fit
evaluates tens of thousands of times and it belongs in the record.

**Round 1 recorded +174 and 171, and those were stale**: they were
taken before the low gate was added and before `ode_ss_rcap` moved, and
were not taken again after. Corrected here, in the round-1 summary
table and in section 12's plan row. A count is not a count if it is
stale.

The absolute node counts here differ from the reviewer's (64 against
68 at `n_ss` = 4) because the two probes use different observation
grids. The DELTAS are what is comparable, and those agree exactly.

## 5. The diagnostic

The old warning reported the movement between the last two cycles. What
it reports now is the distance from THE VALUE RETURNED to the limit:
with the tail summed, the disagreement between successive extrapolants;
with it truncated, the size of the tail that was dropped.

`dev/nss/nss-13-warntext.R`, one compartment, `ka` = 1, `ii` = 12,
`n_ss = 20`, truncated, against the true error over one interval:

| half-life | `ke * ii` | base prints | lane prints | true | base understates by | lane by |
|---|---|---|---|---|---|---|
| 13.9 h | 0.60 | 5.05e-06 | 6.14e-06 | 3.99e-06 | 0.79x | 0.65x |
| 34.7 h | 0.24 | 2.25e-03 | 8.30e-03 | 6.93e-03 | 3.08x | 0.84x |
| 69.3 h | 0.12 | 1.27e-02 | 9.98e-02 | 8.33e-02 | 6.56x | 0.83x |
| 138.6 h | 0.06 | 2.67e-02 | 4.31e-01 | 2.89e-01 | 10.81x | 0.67x |

The lane's figure is a state-space distance and the "true" column is a
concentration error over one interval, which is why the ratio is 0.65
to 0.84 rather than 1.00. It errs toward saying too much, which is the
safe direction; the base figure errs toward saying too little, and by
more as the error grows.

### The BLOCKER: round 0's version of this could not fire

The reviewer's `dev/rev-nss/rev-nss-08-failopen.R`. On
`ss_extrapolate = FALSE` round 0 reported the size of the correction it
had DECLINED to make, and all three guards in section 3 work by driving
that correction to exactly zero. So on every state a guard exists for,
the reported distance was 0 and the check could not fire, however far
the run-in was from its limit. `stalled` sat inside the same
`rel > ss_tol` gate, so it was suppressed with it. Separately, at
`n_ss` = 1 the run-in returned before any check at all, on both arms.

That is the same defect round 0's own code comment claimed to have
fixed for the other arm, which is the part I mind: the absent case was
constructed and the present-but-quiet case was not.

Fixed, and re-measured on the reviewer's three constructions
(`dev/nss/nss-20-failopen.R`), `ss_tol` = 1e-6, `atol = rtol = 1e-10`:

| construction | base | round 0 `FALSE` | round 1 `FALSE` | round 1 `TRUE` |
|---|---|---|---|---|
| AUC state read as the output, q12, `n_ss` 20 | 0.0534 | **silent** | 0.0534 | 0.0534 |
| a run-in that oscillates, `ii` 1 | 0.101 | **silent** | 0.101 | 0.106 |
| `n_ss` = 1, true error 0.2701 | 1 | **silent** | 1 | 1 |
| `n_ss` = 2, true error 0.2313 | 0.461 | 0.461 | 2.75 | 0.733 |

All eight cells warn. What closed the hole is that every branch of the
report now falls back to `move`, the movement between the last two
cycle-start states, which is the one quantity the run-in has always
been able to compute and the only one with no hole in it. On the
extrapolated arm the residual between successive extrapolants already
degrades to `move` when the guards stand the correction down; on the
truncated arm the dropped tail does not, so a state whose movement is
not already negligible and whose ratio is outside (0, 1) puts its own
movement back as a floor. The movement test has to come first: an oral
depot's differences are zero to the last bit, which reads as a ratio of
0, and without it every oral model would carry that floor. The last row
over-reports, 2.75 against a true 0.2313, because at `n_ss` = 2 the
geometric model is being read off the first two cycles; over-reporting
is the safe direction and it is documented as conservative.

**The false-alarm measurement, rerun and widened.**
`dev/nss/nss-12-warn.R`, now **31** schedules: the original 22 ordinary
therapeutics, three flip-flop rows, which the reviewer pointed out the
sweep could not previously reach because every oral row ran `ka` at 1.0
or 1.1 per hour, three short run-ins at `n_ss` of 1, 2 and 4, and three
LONG ones at `n_ss` of 650, 1000 and 2000 on a 1653 hour half-life,
which round 1's set did not reach either.
Each row is classified against the TRUTH, `frm_lincmt()` at
`n_ss = Inf`, not against itself:

| arm | correct silence | correct warning | FALSE ALARM | MISS |
|---|---|---|---|---|
| base commit (22 rows) | 5 | 17 | 0 | 0 |
| lane, `ss_extrapolate = FALSE` | 6 | 25 | **0** | **0** |
| lane, default | 24 | 7 | **0** | **0** |

The seven rows the default still warns on are the three flip-flop rows,
`n_ss` of 1 and 2, and the 1653 hour half-life at `n_ss` of 650 and
1000, and it is right to warn on all seven.

The warning also names the case where no `n_ss` is enough: when a state
is not contracting between cycles at all, it says so, and that is now
reported on its own evidence rather than through the distance.

### Round 2: it went quiet exactly where the documentation sends you

The reviewer found three misses: one compartment, terminal half-life
1653 h, `ii` 24, default arm, at `n_ss` of 650, 800 and 1000, with true
errors 4.68e-05, 1.04e-05 and 1.39e-06 and no warning at all. That
`r` is 0.98999, ABOVE `ode_ss_rcap`, so the stand-down deliberately
declines part of the tail; successive extrapolants then agree with each
other while both sit short of the limit, and a report built only from
their disagreement cannot see the difference. It ranks above a nit
because `?frm_ode` tells a long-half-life user to raise `n_ss`, so the
documentation was walking them into the region where the check goes
quiet.

The report now adds back the part of the tail the stand-down declined,
which is computable from the same factor and is zero wherever the
correction was applied in full, so the ordinary case is untouched.
Measured, `dev/nss/nss-27-largenss.R`, `atol = rtol = 1e-12`:

| `n_ss` | true error | warned | said | said / true |
|---|---|---|---|---|
| 200 | 4.374e-03 | yes | 4.44e-03 | 1.02 |
| 400 | 5.844e-04 | yes | 5.91e-04 | 1.01 |
| 650 | 4.721e-05 | yes | 4.77e-05 | 1.01 |
| 800 | 1.042e-05 | yes | 1.06e-05 | 1.02 |
| 1000 | 1.393e-06 | yes | 1.41e-06 | 1.01 |
| 2000 | 1.861e-09 | no | | |
| 4000 | 4.905e-11 | no | | |

0 misses and 0 false alarms on both arms across the nine `n_ss` values
of that sweep. **That is not a claim about the whole range**, and the
review measured where it fails: in the gate band, at `n_ss` of 3 to 6,
13 false alarms in 82 runs on true errors of 1.8e-10 to 4.95e-07.
There `gate` falls under `ode_ss_trust`, so `modelled` is FALSE and
both terms fall back to `move`, which OVERSTATES by about `1 / r`,
which is the mirror of the defect punch round 1 fixed. Eleven of the
thirteen predate punch round 2, measured by running the `noundone`
variant over the same 60 gate-band runs, so this is pre-existing and
slightly worsened rather than a regression.

**A second limit that raising `n_ss` cannot pass**, found while adding
those rows to the sweep and worth saying because the documentation
sends people there. The cycles are chained solves, so the integrator's
error accumulates over all `n_ss + 1` of them: at the shipped
`atol = rtol = 1e-8` and `n_ss` = 1000 it contributes about **2e-05**,
which is larger than the run-in shortfall the extra cycles were bought
to remove. No run-in check can see that, because it is not the run-in.
The three long rows in the sweep therefore carry 1e-12, and `?frm_ode`
now says to tighten the tolerances alongside `n_ss` past a few hundred
cycles.

**And the number is a detector, not a measurement.** Over the
reviewer's 64 fresh runs the ratio of what the warning prints to the
true error ran from 0.0104 to 9.17e+09. It is built out of the same
geometric model the correction is, so where that model is poor the
number is poor with it, and it is bounded below by the cycle-to-cycle
movement rather than being an estimate in its own right. Round 1's
section 5 read as though the 0.65x to 0.84x it measured on four rows
were a property. `?frm_ode` and the warning's own text now say plainly
that it detects rather than measures.

**The refusal that was considered and not built.** The brief lists
"refusing, by name, a model whose `r` at the optimum implies a shortfall
above some threshold". At the optimum there is nothing left to refuse
on an ordinary schedule: the shortfall is at the integrator's
tolerance, so such a refusal would fire only on the cases of section 8,
and those are already named by the warning. A refusal at the STARTING
values would inherit exactly the defect that kills "auto" in section 2:
it would be measured where the fit is not.

## 6. Every test was seen to fail

`dev/nss/nss-15-seefail.sh` builds EIGHT variants of the package, each
keeping the whole API and removing one property, and runs
`test-ode-events.R` against each from a throwaway library.

**Read the harness's own guard first.** In round 2 this script matched
nothing after the shape changed under it, patched nothing, and reported
`pass=142 fail=0` on all eight variants, every one of which had
installed the shipped code. It now cuts each anchor out of the current
source, raises when a patch is a no-op, refuses to run a variant whose
`R/ode.R` still compares equal to the shipped one, and deletes the
previous log before starting so that a skipped variant cannot leave a
stale one behind. Eight identical clean rows is what gave it away, and
that is a thin thread to have caught it by.

The shipped build is `pass=142 fail=0 err=0` over 52 tests.

| variant | what is absent | result | tests broken |
|---|---|---|---|
| `none` | the correction is not applied at all | 132 / 10 fail | the tail is summed; a state with no steady state; a zero difference; the low gate; a long run-in; the warning's distance; the gradient |
| `joint` | one ratio shared across states | 141 / 1 fail | a state with no steady state |
| `nodamp` | the noise floor in the denominator | 139 / 3 fail | a zero difference; a growing state |
| `hardcap` | the smooth stand-down, a hard cap instead | 139 / 3 fail | the low gate; the objective stays differentiable |
| `oldgate` | the flat-saturating gate, `S(x)/x` instead | 139 / 3 fail | the low gate takes tail away and never adds to it |
| `outside` | the correction built before the branch | 140 / 2 fail | `ss_extrapolate = FALSE` does not build the correction |
| `oldreport` | round 0's diagnostic, built from the correction | 135 / 3 fail, 1 err | the guards cannot silence the warning; a long run-in |
| `noundone` | round 1's report, without the declined tail | 137 / 1 fail, 1 err | a long run-in in the stand-down band still warns |

Every one of the round's code fixes, across all three rounds, has a
test that fails without it, and the failures are behavioural rather
than a missing symbol.

`none` is the behavioural failure the brief asks for:
`ss_extrapolate = FALSE` reproduces the base commit **bit for bit** on
**26** schedules (`dev/nss/nss-24-backcompat.R`), `identical()` TRUE on
26 of 26 and the largest absolute difference exactly 0, not a printed
zero. The set is deliberately awkward: an oral depot whose differences
are exactly zero, a state with no steady state read both ways,
flip-flop absorption, an oscillator, infusions, three compartments,
`addl`, `n_ss` of 1, 2, 3, 4, 5, 20, 40 and 137, and tolerances at
1e-6, 1e-8 and 1e-12.

## 7. The acceptance criterion, and what it does and does not settle

> a two-compartment oral fit at a 107 h terminal half-life recovers
> `k21` inside its own interval at the shipped default, or the default
> is raised and the cost of raising it is stated.

**Met.** `dev/nss/nss-17-accept-lin.R`, six seeds, 101 to 106. Design:
two-compartment oral, `ke` 0.15, `k12` 0.3, `k21` 0.02, `ka` 1.0,
`V` 10, `ii` 24, terminal half-life 107.1 h, 30 subjects x 7 samples,
random effects on `lke` and `lka`, data simulated from the EXACT steady
state, `se = TRUE`, intervals from `confint()`. The vehicle is
`frm_lincmt()`, where `n_ss = Inf` is the limit the shipped default now
reaches and `n_ss = 20` writes out exactly the cycles the old default
truncated at.

**The substitution is licensed by measurement, not by the record.**
`dev/nss/nss-19-objgrad.R` builds both engines' tapes on the SAME
30-subject dataset, the same 66-element parameter vector (five
structural parameters, a residual scale and 60 random effects), and
compares the objective and its whole gradient rather than fitting
twice:

| arm | objective, relative | gradient, max relative over 66 |
|---|---|---|
| `frm_ode(ss_extrapolate = TRUE)` vs `frm_lincmt()` | 6.88e-09 | 3.40e-07 |
| `frm_ode(ss_extrapolate = FALSE)` vs `frm_lincmt(n_ss = 20)` | 3.37e-09 | 1.20e-08 |

A fit is determined by its objective and its derivatives, so two
engines that agree on both to eight figures reach the same optimum.
The same script gives the contrast at population scale: on this
dataset, truncating the run-in moves the objective by **60.66 units**
and turns `d/dlk21` from **+2.915 into -294.63**. That is the sign
flip of section 4 again, at the scale of a real fit and on the
acceptance design itself.

A population `frm_ode()` fit of the same pair was attempted three times
(`dev/nss/nss-10-accept.R`) and each run was killed by the session
before it finished; at 30 subjects x 21 solves per group per evaluation
under four other lanes' load it did not fit inside a background slot.
The objective-and-gradient comparison above is stronger evidence than
that fit would have been, and costs one tape build instead of an hour.

| seed | limit `lk21` [95%] | covers | truncated `lk21` [95%] | covers |
|---|---|---|---|---|
| 101 | -4.09 [-4.94, -3.25] | yes | -3.93 [-4.51, -3.34] | yes |
| 102 | -3.33 [-3.67, -2.98] | no | -3.33 [-3.67, -2.99] | no |
| 103 | -4.30 [-5.37, -3.23] | yes | -4.39 [-5.40, -3.38] | yes |
| 104 | -4.00 [-4.73, -3.27] | yes | -3.91 [-4.49, -3.32] | yes |
| 105 | -4.07 [-4.88, -3.26] | yes | -4.01 [-4.69, -3.34] | yes |
| 106 | -4.56 [-5.93, -3.20] | yes | -4.30 [-5.06, -3.54] | yes |

Truth is `log(0.02)` = -3.9120. The shipped default covers 5 of 6; the
one miss, seed 102, misses in BOTH arms by the same 0.28, so it is the
dataset and not the run-in.

**It does not discriminate at 107 h, and saying so is the point.**
Truncating also covers 5 of 6 there, because the interval on `lk21` is
**0.69 to 2.74 wide** on the log scale across the six seeds and the
run-in moves the estimate by at most 0.263. Round 0 wrote "0.6 to 1.4
wide", which is the half-width, so it undersold its own point by a
factor of two; the widths are in the table above. What discriminates at
107 h is the trajectory, 1.2e-02 to 2.6e-09, and the gradient's sign in
section 4. A criterion written on coverage of a wide interval is
satisfied by a biased estimator.

**What the criterion should have said**, and any of the three costs
nothing extra on the same six fits:

1. the trajectory rather than the estimate: the worst error over one
   dosing interval against `frm_lincmt()` at `n_ss = Inf` must be below
   `ss_tol`. Measured, 1.182e-02 against 2.623e-09;
2. the gradient's SIGN: `d/dlog(k21)` from the taped objective must
   agree in sign with the exact limit's. Measured, +10.37 against
   -1.30;
3. the reported precision: the standard error of `lk21` must be within
   a factor of 2 of the limit's on a design where `k21` is weakly
   identified. Measured at 264.6 h, 41.6x on seed 106.

The third is the one a user cannot check from their own output, so it
is the one worth writing into a plan row.

**The design where it bites is 264.6 h, not 107 h**, per section 1.
`dev/nss/nss-18-accept-265.R`, same shape with `ke` 0.1, `k12` 0.2,
`k21` 0.008. Seed 101 reproduces the review to five figures: the limit
gives `lk21` = -5.5455 against its recorded -5.54547 and truncation
-4.3923 against -4.39225, a `k21` factor of **3.168**.

| seed | limit `lk21` (se) | truncated `lk21` (se) | se ratio |
|---|---|---|---|
| 101 | -5.545 (1.724) | -4.392 (0.318) | 5.4x |
| 102 | -3.746 (0.272) | -3.748 (0.265) | 1.0x |
| 103 | -31.05 (0.266) | -4.886 (0.265) | 1.0x |
| 104 | -4.624 (0.667) | -4.305 (0.365) | 1.8x |
| 105 | -5.239 (1.267) | -4.673 (0.472) | 2.7x |
| 106 | -7.538 (12.81) | -4.599 (0.308) | 41.6x |

Coverage there is 4 of 6 for the limit and 5 of 6 for truncation, and
that is not an argument for truncating. At a 264.6 hour terminal
half-life sampled over one 24 hour interval, `k21` is barely identified:
solved correctly, the model SAYS SO, with standard errors from 0.27 to
12.8 and one fit (seed 103) running to the boundary at `k21` = 3e-14
with `nlminb` reporting false convergence. Truncating the run-in
removes the slow mode from the likelihood, and the model then reports
`k21` as well determined in all six fits, with a standard error 1.0x to
**41.6x** smaller.

**That claim holds for `k21` and does NOT generalize, which round 0 did
not say.** Standard error ratio, limit over truncated, above 1 meaning
truncation reports more precision, from the reviewer's twelve fits
(`dev/rev-nss/rev-nss-05.log` and `rev-nss-05b.log`):

| design | `lke` | `lka` | `lk12` | `lk21` | `lV` |
|---|---|---|---|---|---|
| 264.6 h | 0.12 to 0.81 | 0.99 to 1.07 | 0.47 to 1.46 | **1.00 to 41.64** | 0.94 to 1.31 |
| 107.1 h | 0.29 to 0.98 | 1.00 to 1.04 | 1.01 to 1.46 | 1.01 to 1.80 | 1.00 to 1.10 |

`lke`'s ratio is below 1 in all twelve, so truncation INFLATES the
standard error of `ke` by 1.02x to 8.3x while deflating that of `k21`
by up to 41.6x. The honest statement is that **truncation removes the
slow mode from the likelihood and redistributes the information**: the
parameter that controls that mode looks well determined and the fast
parameters absorb the misfit. It manufactures a precision for `k21`. It
is not a general property of the fit, and it is still not visible in a
coverage count.

## 8. What a user still does not get

- **Flip-flop absorption, which IS pharmacokinetics.** Round 0 wrote
  that the only losing case was a driven oscillator and "it is not
  pharmacokinetics". The reviewer showed that is wrong and gave the
  mechanism: the ratio the run-in reads is a weighted mean of the cycle
  map's modes, and when two weights have opposite signs it lies outside
  their range and the correction overshoots. Opposite signs are the
  NORMAL arrangement in an oral model, which is why a concentration
  rises before it falls, so the hazard is `ka` near `lambda_z`, which
  is what extended-release and depot formulations are built to produce.

  Re-measured against the shipped shape, `dev/nss/nss-23-flipflop.R`:
  338 finite two-compartment oral cycle maps, `lambda_z * ii` from 0.02
  to 3 and `ka / lambda_z` from 0.1 to 50. The correction loses on
  **8**, worst by a factor of **2.27**, and every losing case has
  `lambda_z * ii` = 0.05 with `ka / lambda_z` of 1.5 or 2. The median
  ratio over all rows is 0.222.

  **The bound is what makes the default defensible.** The smallest
  error truncation leaves on a losing case is **4.02e-01** and the
  largest it leaves on a winning one is **9.61e-01**, so the losing
  region is entirely inside "truncation is at least 40 percent wrong",
  and the warning fires on every losing case in both arms. The
  correction never loses where truncation was usable, and where it
  loses the user is told. Through `frm_ode()` on the reviewer's worst
  row, a 332.7 h terminal half-life with `ka / lambda_z` = 1.5 dosed
  daily: **5.2286e-01 truncated against 7.9513e-01 extrapolated**, both
  warned. The reviewer measured 5.229e-01 and 6.414e-01 against the
  round-0 shape, so the C1 shape of section 4a costs 1.24x on this row.
  `?frm_ode` and NEWS name flip-flop and carry the bound.
- **A run-in that oscillates rather than converging.** A damped
  oscillator dosed every `ii` has complex cycle-map eigenvalues, and a
  single-mode extrapolation does not describe it. At `ii = 1` the
  truncated answer is 7.167e-01 wrong and the extrapolated one
  9.950e-01, so the correction is **1.39x worse** there
  (`nss-09-verify.R` section B). No `n_ss` in this range is enough for
  that system either. It is not pharmacokinetics, but `frm_ode()` takes
  an arbitrary `dyn`.
- **A strongly nonlinear system.** Michaelis-Menten elimination deep in
  its saturated regime, `Km = 50` against a concentration of order 1:
  9.635e-02 truncated, 6.226e-03 extrapolated, a 15x gain and no more,
  because the cycle map is not affine.
- **Two modes of nearly the same rate.** The correction removes their
  combination and leaves what is left of the other: 5.713e-04 to
  9.958e-06 when the two rates differ by 0.5 percent.
- **A state with no steady state still has no steady state.** Reading
  an AUC compartment on an `ss` row gives 9.479e-01 relative error
  either way. The warning now says so; the value is no better.
- **A check during a fit.** RTMB still refuses comparison on the tape,
  so the warning is still numeric-path only. What changed is that the
  answer the fit uses is right, not that the fit can check it.

## 9. Defects found and NOT fixed

- **A derivative written `0 * y[k]` cannot be solved.**
  `dev/nss/nss-16-nonconsec.R`, identical on the lane build and on the
  reference build, so it is pre-existing and unrelated to this item:

  | second derivative | result |
  |---|---|
  | `y[1]`, depends on a state | solves |
  | `p[1]`, a plain constant | solves |
  | `0 * y[2]` | "Tape has Non-consecutive outputs" |

  It surfaces as a failed solve, so on the default `on_error` the
  group's rows carry the penalty and the user gets the generic
  failed-solve warning rather than anything naming the cause. Writing a
  held state as `0 * y[k]` is the natural spelling of "this compartment
  does not change", so this is worth a refusal that names it. It is a
  `dynamics` validation defect on `frm_ode()`'s general path with no
  connection to `n_ss`, so it belongs in a plan row of its own rather
  than folded into an item that is already a behavioural break. The row
  text is in section 13.
- **`frm_lincmt(n_ss = )` still truncates.** It is exact at its default
  `n_ss = Inf`, so this is not a wrong answer anyone gets by accident,
  but a user who sets a finite `n_ss` there to match `frm_ode()` no
  longer matches it. `test-lincmt.R` now compares the two DEFAULTS,
  which agree because both reach the limit, and keeps a matched-
  truncation assertion beside it.
- **`frm_ode()` absorbs unknown arguments into `...`.** `...` is passed
  to the integrator, so a misspelled or unsupported argument is
  silently ignored: the reference build accepted `ss_extrapolate` and
  did nothing with it, which is how section 4's BEFORE arm was measured
  without a separate script. That is convenient here and is a trap in
  general.

## 10. What was run

`frmtmb.ode`'s whole suite, one test file per R process
(`dev/nss/nss-suite.sh`), `NOT_CRAN=true`, against the lane library:

| file | assertions | baseline | fail | err | skip |
|---|---|---|---|---|---|
| `test-bracket-access.R` | 1 | 1 | 0 | 0 | 0 |
| `test-compat.R` | 15 | 15 | 0 | 0 | 0 |
| `test-lincmt.R` | 65 | 63 | 0 | 0 | 0 |
| `test-ode.R` | 72 | 72 | 0 | 0 | 0 |
| `test-ode-events.R` | 142 | 108 | 0 | 0 | 0 |
| `test-ode-nlf.R` | 4 | 4 | 0 | 0 | 0 |
| `test-ode-tv.R` | 27 | 27 | 0 | 0 | 0 |
| `test-scale.R` | 0 | 0 | 0 | 0 | 1 |
| **total** | **326** | **290** | **0** | **0** | **1** |

The baseline column is `dev/suite-baseline.tsv` at the base commit. The
one skip is the gated scale tier, which skips without
`FRMTMB_SCALE_TESTS`. `+34` in `test-ode-events.R` and `+2` in
`test-lincmt.R` are this lane's: +15 in round 0, +9 in round 1 and +10
in round 2.

`R CMD check --as-cran` on a tarball built WITH vignettes:
**Status: 1 NOTE**, the expected "Skipping checking math rendering:
package 'V8' unavailable" on the HTML manual. Tests 71 s under
`test_check()`, vignette rebuild included. The tarball was unpacked and
diffed against the working tree before the result was believed: `R/`,
`man/`, `NEWS.md`, the vignette and both changed test files match.

The eight guard-absent variant builds of section 6 each got their own
throwaway source directory and library, `nss-var-<v>` and
`nsslib-var-<v>` under `C:/Users/adf44/source/r/`, deleted after each.
They share nothing, because two copies of the whole-set driver sharing
one work directory deadlocked twice. Nothing was installed into the
user library or into `rellib-0552`.

## 11. The version bump this needs

**Minor, 0.3.0 to 0.4.0.** It changes the numbers a shipped model
returns for `ss` rows, which is a behavioural break even though it is a
correction, and it adds an argument. Choosing it is the consolidating
session's, not the lane's; the NEWS bullet says plainly what stops
working and how to get the old numbers back.

## 12. What `dev/extension-gaps-plan.md`'s 1.0d row should say

> **1.0d DONE.** `frm_ode()`'s steady-state run-in truncated a
> geometric tail whose size is `exp(-n_ss * lambda_z * ii)`, and the
> tail is now summed instead: `ss_extrapolate = TRUE` reads the
> per-cycle contraction off the run-in's own last differences and adds
> `d * r / (1 - r)`. It runs on the tape, so it tracks the parameters
> during a fit, and it costs no extra solve (21 per group at
> `n_ss = 20`, either way) and +258 tape nodes on the default arm. The
> default `n_ss` is unchanged at 20. On a two-compartment oral model at
> a 107 h terminal half-life dosed daily the error over one dosing
> interval falls from 1.2e-02 to 2.6e-09, and the gradient with respect
> to `log(k21)` from +10.37, the WRONG SIGN, to -1.30 against the exact
> limit's -1.302; at population scale on a 30-subject dataset the
> objective moves 60.66 units and `d/dlk21` goes from +2.915 to
> -294.63. The correction is applied in full up to a measured ratio of
> 0.9875, gated to nothing below 0.05 and stood down to nothing at 1,
> so the objective has no corner at any of the four junctions, which a
> hard cap did (-0.09 against +150, not converging). It LOSES to truncating on flip-flop absorption, 8 of
> 338 cycle maps and worst by 2.27x, bounded: the smallest error
> truncation leaves on a losing case is 0.40, so it never loses where
> truncation was usable, and the warning fires in both arms there.
> The numeric-path warning now reports a distance to the limit rather
> than the cycle-to-cycle movement it used to, which understated the
> error by up to 10.8x and did so most where the error was largest,
> though it is a DETECTOR and not a measurement and `?frm_ode` says so;
> 0
> false alarms and 0 misses on 31 schedules in both arms, including
> three where the correction's own guards fire, two where `n_ss` is too
> small for a tail to be read, and three where it is in the hundreds
> and inside the stand-down band.
>
> `n_ss = "auto"` was measured and REJECTED on **dominance**: floored
> at 20 it is never worse than today, so "worse than the status quo" is
> not the ground; but on the 107 h design it needs 135 `lsoda` entries
> per group to reach 2.85e-09 where the correction reaches 2.62e-09
> with 21, the same place for 6.4x the solves, and it is fixed where
> the tape is built while the correction tracks the fit. Not fixed: a
> run-in that oscillates, a strongly nonlinear cycle map, and a state
> that has no steady state, all named in `dev/nss-findings.md`
> section 8.
>
> **The acceptance criterion as written was met and worthless**: the
> interval on `lk21` is 0.69 to 2.74 wide on the log scale and the
> run-in moves the estimate by at most 0.263, so both arms cover 5 of 6
> with the same miss. It should have asked for the trajectory error
> against `frm_lincmt(n_ss = Inf)` to be below `ss_tol`, or the
> gradient's SIGN to agree with the exact limit's, or the standard
> error of `lk21` to be within a factor of 2 of the limit's on a design
> where `k21` is weakly identified (measured 41.6x). The third is the
> one a user cannot check from their own output.

## 13. A plan row that is not this lane's

> **`frm_ode()` cannot solve a derivative written `0 * y[k]`, and says
> nothing useful about it.** `list(c(-ka * y[1], 0 * y[2]))` fails with
> "Tape has Non-consecutive outputs" from RTMBode, while `y[1]` and a
> plain constant both solve; identical on 0.3.0 and on the lane build,
> so it is pre-existing. It surfaces as a failed solve, so the group's
> rows carry the penalty and the user gets the generic failed-solve
> warning with nothing pointing at the cause. Writing a held
> compartment as `0 * y[k]` is the natural spelling of "this state does
> not change", so a user cannot get from the symptom to the fix. Needs
> a refusal that names the state and the spelling, with its own
> see-it-fail evidence. It is a `dynamics` validation defect on
> `frm_ode()`'s general path with no connection to `n_ss`, which is why
> 1.0d filed it rather than widening itself. Construction:
> `dev/nss/nss-16-nonconsec.R`. Half a day.
