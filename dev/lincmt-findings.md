# frm_lincmt(): the analytic linear compartment model

Item 1.0c of `dev/extension-gaps-plan.md`, moved to Phase 1 by Phase 0
because one `frm(se = TRUE)` on the plan's ode design took 3948 s.
Lane `lincmt`, worktree `frmtmb-wt-lincmt` off 92e9330, frmtmb 0.55.1,
frmtmb.ode 0.2.0, RTMB 1.9, RTMBode 1.0, deSolve 1.42, R 4.6.1.

Every measurement script is in
`extensions/frmtmb.ode/dev/lincmt/`, copied there rather than left in a
session scratchpad. Each number below names the script that produced
it and, where a number depends on one, the seed.

## What punch round 1 changed

Reviewed in `dev/reviews/2026-09-09-lincmt.md`, which found no blocker
in the code and four wrong or unsupported claims in the prose. Every
correction below is arithmetic on numbers already measured, or a short
targeted run; the Phase 0 design was not re-run.

| finding | what was wrong | what it says now |
|---|---|---|
| F1 | the vignette's "57.3 s" existed nowhere in the repository | 55 s at the default, 119 s at `n_ss = 20`, both from the kept log |
| F2 | the vignette's "0.2 percent short" was the cycle-to-cycle movement | 0.4 percent short, with the movement named as the other number |
| F7 | the speedup was derived "within a pass", and the arms are sequential | at least 67x and at least 31x, with 95x and 47x cited as corroboration from a design where the pairing IS sound |
| F10 | the vignette's "the same estimates to every digit" | seven significant digits, with the digits shown |
| F3 | `?frm_lincmt` quoted only the coalescence sweeps | the bound, which holds everywhere, plus the wide-box figure |
| F4 | "right to eps of its own scale" | 3.6e-10 of its own scale on a seven-decade box, and the bound rather than an adjective |
| F5 | "why the textbook form is not used" read as covering both cancellations | narrowed to the absorption term, with the other one named and bounded |
| F6 | the gradient sweep's degenerate arm was misnamed | renamed, and the real eigenvalue collision added, which found a `NaN` |
| F8 | pass 1 has no log | said, with the reason it is load-bearing only against the lane's own interest |
| F9 | two labels on the solve-count table | corrected |
| N1 | a `NaN` at all rates underflowed to zero | documented, and the obvious fix measured and rejected |
| S1 | I declined the `n_ss` item on two grounds | both grounds fail on measurement; `?frm_ode`'s wrong sentence is corrected and the item text is written |

## What punch round 2 changed

Round 2 of the same review, section 0, findings R1 to R5. Verdict
mergeable, no blocker in the code. Nothing here changed behavior and
nothing here needed a run: every correction is a sentence or a division.

| finding | what was wrong | what it says now |
|---|---|---|
| R1 | the Boundary section said the derivative at a double root is "genuinely unbounded", and that clamping would cost the value's 5e-13 | both false, and both were the stated reasons for not fixing: the derivative is -9.2779361 to eight significant digits across three decades of step size, and a 1e-14 clamp costs 5.177e-13 against the shipped 5.181e-13 with no change at all on 400 ordinary draws |
| R2 | "a bit-level knife edge rather than a neighborhood" | 0 NaN in 4000 draws on a PK box and 4000 on 1e-4 to 20, 0.07 percent at an eight-decade spread, 1.33 percent at ten |
| R3 | the silent VALUE error in that region was not reported at all | up to 1.3e-04 of the trajectory's peak where the gradient is finite and nothing fires, now stated as a domain with the clean side measured |
| R4 | "at least 68x" in four places; 3716.6 / 55.0 is 67.575 | at least 67x in all four, with the arithmetic shown |
| R5 | N1's silent factor was 1.5, from a script that never applies the offset | 1.125: coefficients 4/9 not 1/3, height eight-ninths. The rejection stands |

**A note on this document's own provenance.** While making the R1 edit
I replaced a span that reached further than intended and truncated the
sections between "The gradient identity" and "The suite". They are
restored above. Every number in the restored text is one the review
independently reproduced in its own tables, or one that a named script
in `dev/lincmt/` re-runs; nothing was recalled without a source. A
reader wanting to check the record rather than trust it should re-run
the script named beside each number.

Following F6 turned up something round 1 did not have: adding the
eigenvalue collision it asked for showed that `lincmt_disp()`'s
gradient is `NaN` at a double root. **Round 2 then showed that my
description of THAT was wrong in both halves** (R1, R2, R3): the
derivative exists and is finite, the region is a neighborhood rather
than a knife edge, and inside it the VALUE is silently wrong by up to
1.3e-04 of the trajectory's peak where nothing fires at all. The
section "The gradient identity" carries the corrected account.

## What was built

`extensions/frmtmb.ode/R/lincmt.R`. `frm_lincmt()` evaluates a one- to
three-compartment linear model in closed form, with an optional depot:
the amount in a compartment is the superposition of the impulse
responses of the doses that precede each observation, and a
steady-state record is the geometric limit of that superposition rather
than a simulated run-in.

The schedule grammar is `frm_ode()`'s and is NOT re-implemented. The
`events` table goes through the same `ode_split_events()`, so `time`,
`value`, `state`, `method`, `duration`, `ii`, `addl`, `ss` and `group`
mean what `?frm_ode` says they mean, including that an observation at a
dose time reads the trough and an observation at a reset time reads the
state before the reset. What the closed form cannot express it refuses
by name.

The whole model is evaluated in ONE vectorized pass over the whole
data. There is no per-group loop on the tape, because there is no
solver whose per-group ceiling the loop exists to respect.

## The identity, per schedule kind

`dev/lincmt/lincmt-ident.R`. 118 cases: 1, 2 and 3 compartments crossed
with depot and no depot, crossed with ten schedule kinds (init bolus,
single dose, `addl`/`ii`, `ss`, `ss` plus `addl`, single infusion,
repeated infusion, infusion at steady state, `reset`, `t0` not zero),
plus a depot-amount output, plus 54 degenerate-rate cases, plus a
three-group model with a per-group schedule and an estimated
`event_scale`. The reference is `frm_ode()` at `atol = rtol = 1e-12`
rather than at its 1e-8 default, so that what is measured is the closed
form and not the integrator.

| measure | worst over 118 cases |
|---|---|
| `max abs difference / max abs value` | **2.62e-12** |
| errors | 0 |

Per schedule kind, worst over the shapes that carry it:

| schedule kind | cases | worst difference / scale |
|---|---|---|
| infusion at steady state | 6 | 2.62e-12 |
| `t0` not zero | 6 | 2.39e-12 |
| `ss` | 6 | 2.29e-12 |
| `ss` plus `addl` | 6 | 2.29e-12 |
| single dose | 6 | 2.21e-12 |
| `addl`/`ii` | 6 | 2.01e-12 |
| init bolus | 6 | 1.82e-12 |
| `reset` | 6 | 1.82e-12 |
| repeated infusion | 6 | 1.77e-12 |
| single infusion | 6 | 9.25e-13 |
| depot amount as the output | 3 | 3.14e-13 |
| the 54 degenerate-rate cases | 54 | 1.95e-12 |
| three groups, per-group schedule, `event_scale` | 1 | 1.13e-13 |

The order tracks how many segments the SOLVER has to chain, not
anything about the closed form: the arms the solver finds hardest are
the ones the two disagree on most.

**Pointwise relative error is the wrong measure here and the run says
why.** Its worst value over the same 118 cases is 1.0, in ONE case:
1 compartment with a depot, `ka == ke` exactly, a steady-state row, at
an observation where `frm_ode()` at `atol = 1e-12` returns 1.89e-16 of
the trajectory's own maximum. At a value that small the reference is
the solver's absolute tolerance and nothing else. The next largest
pointwise value in the table is 7.1e-07, at a point that is likewise
below 1e-10 of the trajectory's scale. Where the reference carries
signal, pointwise and scale-relative agree.

### Schedule semantics, where a disagreement would be a defect

`dev/lincmt/lincmt-edge.R`. Nine constructions where the two could
disagree about MEANING rather than about arithmetic, at
`atol = rtol = 1e-12`:

| case | difference / scale |
|---|---|
| `ss` row at t = 12, observations either side | 1.26e-13 |
| `reset` at t = 12 with an observation there | 3.21e-13 |
| `reset` at t0, `init` discarded | 0 (identical) |
| observations exactly on every `addl` dose | 1.26e-13 |
| duplicated and unsorted observation times | 3.21e-13 |
| two `add` rows at one instant | 3.29e-13 |
| overlapping infusions | 1.51e-12 |
| infusion at steady state with `duration == ii` | 1.19e-12 |
| a group with an empty schedule beside groups with one | 3.21e-13 |

**No defect in `frm_ode()` was found.** That is the answer to the
task's "the most valuable thing you could find": the closed form is an
independent implementation of the same contract and it agrees with the
solver on all 127 constructions above, including every corner of the
event grammar. What it did find is an APPROXIMATION, below, and a
defect in a MEASUREMENT helper, further below.

### The gradient identity, across every shape and schedule

`dev/lincmt/lincmt-ad2.R`. The value agreeing is half the claim. For
each of the six model shapes crossed with six schedules (single dose,
`addl`/`ii`, `ss`, `ss` plus `addl`, repeated infusion, infusion at
steady state), the tape's Jacobian with respect to the log rate
constants is compared with `frm_ode()`'s adjoint Jacobian at
`atol = rtol = 1e-12`, at three parameter points.

**F6, corrected, and following it found something.** The second point
is every rate constant exactly equal, and an earlier version of this
document called that a point "where every pair of them coalesces at
once". It is not. The review's algebra is right: at
`ke = k12 = k21 = k13 = k31 = r` the characteristic polynomial factors
as `(L - r)(L^2 - 4 r L + r^2)`, so the eigenvalues are `r`,
`(2 - sqrt 3) r` and `(2 + sqrt 3) r` and are distinct. Confirmed
against the code at `r = 0.3`: 1.1196, 0.3000, 0.0804. What that point
DOES coalesce is `ka` with one eigenvalue, which is `lincmt_diff()`'s
singularity, so the arm is real but was misnamed. It is now called
"ka on an eigenvalue".

A third point was added for the collision the sweep lacked: `k13` at
1e-300, which decouples the third compartment, with `k31` on the slow
root of the reduced quadratic, which makes that root double.
`lincmt_disp()` reaches it through `acos()`, whose derivative is
infinite where a double root puts its argument.

| arm | cases | worst gradient difference | worst value difference |
|---|---|---|---|
| ordinary | 36 | 5.83e-13 | 7.18e-13 |
| `ka` on an eigenvalue | 36 | 1.89e-13 | 7.67e-13 |
| a double eigenvalue | 12 | **2.74e-09** | 7.75e-13 |
| all | 84 | 2.74e-09 | 7.75e-13 |

The double eigenvalue is three decades looser than the rest, which is
what an infinite eigenvalue sensitivity buys, and it is still finite
and still agrees.

**And there is a boundary. My round 1 description of it was wrong in
both halves, and the review measured both.**

What I wrote was that the gradient is `NaN` at a bit-level knife edge,
and that the `NaN` is correct because "at a double root the eigenvalues
have a square-root branch point, so the derivative with respect to any
parameter that SPLITS the root is genuinely unbounded". The branch
point is real and it is in the EIGENVALUES. It is not in what
`frm_lincmt()` returns.

**R1: the derivative exists and is finite.** The central compartment's
response is a SYMMETRIC function of the three roots, hence a function
of the characteristic polynomial's coefficients, which are polynomials
in the rate constants; the residue formula's poles at a double root are
removable, so the trajectory is analytic in the rate constants at the
tangency. Measured at my own construction (`dev/rev-lincmt-f6.R`):

| | d/d(log ke) | d/d(log k31) | d/d(log ka) |
|---|---|---|---|
| the tape | **NaN** | **NaN** | 3.5969724 |
| central difference, h = 1e-4 | -9.2779361 | 1.95e-10 | 3.5969724 |
| central difference, h = 1e-5 | -9.2779361 | -1.95e-09 | 3.5969724 |
| central difference, h = 1e-6 | -9.2779361 | 0 | 3.5969724 |
| the tape one part in 1e8 away | -9.2779361 | -2.32e-08 | 3.5969724 |

Two components are `NaN`, not one, as I reported. And `d/d(log ke)` is
-9.2779361 to eight significant digits across three decades of step
size. A derivative a central difference resolves to eight digits is not
unbounded. **My own `dev/lincmt/lincmt-f6.R` printed the evidence
against my conclusion**: its "FD spread" column, the difference's own
step sensitivity, reads 9.573e-10 at the tangency, which is the column
saying the difference knows its answer there. I read the NaN in the
tape column and stopped.

So the correct statement is: the derivative exists, and the
IMPLEMENTATION returns `NaN` because it routes through the
eigenvalues, where `acos()` takes an argument of exactly one.

**R1, second half: the clamp I rejected without measuring costs
nothing at that point.** I wrote that padding the `acos()` argument
"would cost the 5e-13 the value currently holds". Measured
(`dev/rev-lincmt-f6c.R`): at a clamp of `1 - 1e-14` the gradient is
finite and right and the value error is **5.177e-13 against the
shipped 5.181e-13**, with 0.000e+00 change over 400 ordinary draws.
The padding is free there because it perturbs the two colliding roots
SYMMETRICALLY, and the trajectory is a symmetric function of them. The
assertion was wrong.

**R2: it is not a knife edge, and the region is wider than the
construction.** Measured incidence of a `NaN` gradient, 4000 draws per
box and 3000 per spread cell:

| box | draws | `NaN` |
|---|---|---|
| 1e-2 to 10 per hour, a PK box | 4000 | 0 |
| 1e-4 to 20 | 4000 | 0 |
| 1e-8 to 50 | 4000 | 38 |
| 1e-20 to 50 | 4000 | 749 |

| decades of spread | `NaN` rate |
|---|---|
| 2, 3, 4, 5, 6 | 0.00% |
| 8 | 0.07% |
| 10 | **1.33%**, one draw in 75 |
| 12 | 3.53% |

A `NaN` at `k13 = 1.8e-05` with no root within 1e-10 of another is in
that region and is not in my construction. "A bit-level knife edge" is
withdrawn.

**R3, and this is the part that matters more than the `NaN`: in the
same region the VALUE is silently wrong.** On 25 ten-decade draws where
the gradient comes back FINITE, so nothing fires at all, the worst
error against a 300-bit reference is **1.33e-04 of the trajectory's own
maximum**, four decades outside the item's bar. The `NaN` is the loud
end of a soft region, not a boundary guarding a right answer. A silent
value error outranks a loud `NaN`, and I did not report it because I
never looked outside the construction.

**Which is also why the clamp is not the fix.** Applied to 40 draws
from that region it makes every gradient finite and changes the VALUE
by a factor of 2.8e+02 to 2.5e+06. It is free where the collision is
genuine and destructive where the roots were already garbage, and
telling those apart needs a test on whether `pp` and `qq` still carry
relative precision, which is a comparison on an AD value and is the
constraint this whole file is written around. So the rejection stands
and both of the reasons I gave for it were wrong.

**The domain, which is the thing to ship.** Up to a six-decade spread
of rate constants the closed form's worst disagreement with `frm_ode()`
at `atol = rtol = 1e-12` is BELOW what the solver's own tolerance costs
on the same draws, in every cell, for two compartments and for three
(9.9e-11 against 1.3e-07 at two decades; 1.4e-07 against 1.9e-05 at
six). The direction an optimizer actually travels does not degrade at
all: driving the four peripheral rates from one to twelve decades below
the absorption rate leaves the worst error at **6.9e-11**. Past eight
decades `frm_ode()` stops being a reference at all, hitting DLSODA's
step limit and disagreeing with itself by 1e+05 and more. The region
needs the FAST rates extreme too, at a terminal half-life of seconds.

**Reachability in a fit.** Fitting three compartments to
two-compartment data from six starting points, including `log(k13)` at
-20, -30, -35 and -700, all six fits returned normally and `nlminb`
left `log(k13)` exactly where it began, because the gradient in that
direction is zero once `k13` stops mattering. An optimizer cannot walk
INTO the region; a user can start inside it, and then the fit stops
with `NA/NaN gradient evaluation`, which is loud and unactionable.

All of this is now in `?frm_lincmt`, as a stated domain plus a Boundary
section that says the derivative exists and the implementation does not
find it.

## The numerics, against a 240-bit reference

`dev/lincmt/lincmt-mpfr.R`. `frm_ode()` cannot settle whether the
closed form is right in a region where the closed form's textbook
spelling has no digits left, because the solver has its own error
there. The reference is therefore `Rmpfr` at 240 bits, exponentiating
the rate matrix by scaling and squaring and integrating it with the
augmented matrix `[[M, b], [0, 0]]` for an infusion. It never forms an
eigenvalue and never divides by a difference of rate constants. The
review checked that reference against an INDEPENDENT 300-bit one
sharing no line with it, built two ways including the textbook partial
fraction with eigenvalues found by bisection and Newton; the two agree
to 2.46e-86.

| sweep | points | worst relative error |
|---|---|---|
| `ka` down to `ke`, 1e0 to 0, two lags | 36 | 1.52e-15 |
| the same at a steady state, exact limit | 14 | 3.15e-16 |
| 2 cmt double root, `k12` to 0 with `k21 == ke` | 18 | 5.05e-16 |
| 3 cmt double root, `k31` to `k21` | 18 | 7.24e-15 |
| 3 cmt triple root, every rate equal | 24 | 4.47e-15 |
| infusion, `duration` 0.001 to 4, four lags each | 16 | 3.28e-16 |
| random, rates in [1e-3, 5], lag in [0.05, 72], 1 cmt | 60 | 6.87e-15 |
| the same, 2 cmt + depot | 60 | **7.07e-11** |
| the same, 3 cmt + depot | 60 | 3.44e-12 |

Every sweep includes the point where the rate constants are EXACTLY
equal, not merely close. Seed 11 for the random sweeps. The review
re-ran all nine against its own build and reproduced every one to the
printed digits.

The worst point in this lane's own box is 7.07e-11, at `ke = 2.393`,
`k12 = 0.001689`, `k21 = 0.004481`, `k13 = 0.2153`, `k31 = 0.001134`,
`ka = 0.4100`, lag 34.88. That is the one cancellation left in the
formulation and it is understood: the two-compartment disposition
coefficients are `(D +/- g) / 2D` with `D = sqrt(g^2 + k12 k21)`, and
when `k12 k21` is negligible against `g^2` the SMALLER coefficient is a
difference of two numbers that agree to nine digits.

It was left rather than removed, because removing it needs a selection
on the sign of `g`, which is a comparison on an AD value, and the
branch-free spellings of that selection are all wrong AT `g == 0`,
which `k10 + k12 == k21` reaches exactly in floating point (0.1 + 0.1
is 0.2 to the bit). A wrong derivative on a reachable set is worse than
a bounded relative error on an unreachable one.

**F4, corrected.** An earlier version of this document said the
ABSOLUTE error "stays at machine epsilon and the trajectory is right to
eps of its own scale". That is true on this lane's box and false on a
wider one. The review's searches (`dev/rev-lincmt-bound.R` seed 909,
`dev/rev-lincmt-bound2.R` seed 5150) put the worst error relative to
the trajectory's own maximum at 1.88e-13 for two compartments and
**3.57e-10** for three, over a seven-decade rate box, three to five
decades above eps. The review says plainly that its hill climb on that
quantity did not finish, so 3.57e-10 is a random search's worst and not
a converged maximum; it is not presented as one here either.

The right statement is a BOUND rather than an adjective, and it is the
review's:

    pointwise relative error <~ (error / peak) x (peak / value)

The absolute error is held at a few parts in 1e10 of the trajectory's
own maximum. The POINTWISE relative error at a point sitting `10^-d`
below `Cmax` is therefore at most about that times `10^d`: within four
decades of `Cmax` about 1e-11, within six about 1e-9, and the item's
1e-8 bar is crossed only beyond seven decades below `Cmax`, which is
far below any assay's limit of quantification. The pointwise error
alone IS unbounded, and the review drove it to 1.00, no correct digits,
by taking the two-compartment decoupling to where the small coefficient
is exactly zero; the error relative to the peak over that same sweep
never left 4.6e-16 to 2.3e-15. The two facts are one fact.

So: 3.57e-10 of the trajectory's scale on a seven-decade box, 1.4
decades inside the item's 1e-8, and not eps.

### Why the textbook form is not used, for the ABSORPTION term

`dev/lincmt/lincmt-mpfr.R` and the test file. At `ka - ke = 1e-12` the
textbook `D ka (exp(-ke t) - exp(-ka t)) / (V (ka - ke))` is wrong by
3.4e-05 relative where `frm_lincmt()` is wrong by 6.6e-17, and at
`ka == ke` the textbook form is `NaN`. `ka == ke` is not a corner: it
is the flip-flop case of an oral model, and it is the first objective
evaluation of any fit that starts two log rates at the same value.

**F5, narrowed.** That comparison is about the ABSORPTION cancellation
and only about it, and an earlier version of this section did not say
so. The review put the two spellings side by side in the OTHER
cancellation this file contains, the near-decoupled two-compartment
disposition coefficient (`dev/rev-lincmt-cancel.R` section E), and
there they lose the same digits within a factor of a few in either
direction: at `k12 = k21 = 1e-5` the shipped form is wrong by 8.28e-08
pointwise and the textbook by 1.23e-06; at 1e-7, 2.30e-02 against
9.73e-03; at 1e-8, 1.00 against 2.00. `coef2` is spelled `1 - c1`, and
that subtraction is itself a difference of nearly equal numbers.

So the claim is: the absorption singularity is removed COMPLETELY, and
the disposition cancellation is not removed but is bounded in the
measure that matters, as above. `?frm_lincmt` now says both.

### Which of the two is the more accurate

The review settled a claim this document could not. Saying the 2.6e-12
is "the solver's tolerance rather than the closed form's error" is
testable by tightening the solver and seeing which of the two moves.
`dev/rev-lincmt-solverbar.R`, three multi-dose schedules, `frm_ode()`
at 1e-9, 1e-12 and 1e-14: the disagreement falls from 1.65e-10 to
8.48e-13 to 1.70e-15 on the first, and from 4.02e-10 to 2.22e-13 to
4.48e-15 on the second. It does not stall. So the closed form is not
merely within the solver's tolerance of the solver, it is what the
solver converges TO.

## The tape, and the defect this lane found in its own code

`dev/lincmt/lincmt-ad.R`, `dev/lincmt/lincmt-seenfail.R`.

The first version of `lincmt.R` handled the removable singularity by
OFFSET alone: `phi(x) = -expm1(-(x + 1e-300)) / (x + 1e-300)`. Its
VALUE is correct to the last bit at `x = 0`. Its GRADIENT is not.
RTMB differentiates `n / z` by the chain rule into `n' / z - n / z^2`,
two numbers of size 1e300 whose difference is `-1/2`, so every digit is
lost:

| at `ka == ke` exactly | d/d(log ka) | d/d(log ke) | d/d(log V) |
|---|---|---|---|
| offset only | -1.139e+286 | 1.139e+286 | -31.401 |
| finite difference | 4.596 | -26.805 | -31.402 |
| shipped | 4.596 | -26.805 | -31.401 |

The value was 31.4012667478 in BOTH, to twelve digits. This is exactly
the failure the item warned about: right in R arithmetic, wrong in its
gradient. An optimizer that starts at `lka == lke` takes a first step
of size 1e286.

The fix is `lincmt_phi()`: the Taylor series where it converges, the
`expm1()` form where the series would need too many terms, blended by a
smoothstep whose first and second derivatives vanish at both ends.
Nothing branches: both arms are evaluated for every element, the series
argument is capped so that it cannot overflow in the arm whose weight
is zero, and the weight is built from `abs()`.

`test-lincmt.R`'s "the gradient is finite and right where two rate
constants meet" pins it. **I have seen it fail**: run against the
offset-only primitive it reports a ratio of 3.63e284 against a bar of
2. `dev/lincmt/lincmt-seenfail.R` reproduces both arms in one process,
and the review re-ran that construction and got 3.626e+284.

The review also checked that the fix is complete rather than local:
`lincmt_phi()` is reached through eight distinct branches of
`lincmt_resp()`, and each was driven to exact coalescence with its
gradient compared against a central difference, between 1.78e-11 and
3.67e-11 on all eight, with Hessians between 1.59e-09 and 2.56e-09.

With the fix, on the Phase 0 schedule at 12 observations:

| quantity | ordinary parameters | `ka == ke` exactly |
|---|---|---|
| value, closed form vs `frm_ode()` (1e-12) | 1.17e-14 | 4.21e-13 |
| gradient, tape vs `frm_ode()`'s adjoint | 1.68e-13 | 3.43e-13 |
| gradient, tape vs central difference | 1.49e-10 | 1.08e-10 |

With an `ss` row instead: value 5.00e-14 and 3.71e-14, gradient
1.51e-13 and 4.36e-13.

The Laplace approximation needs a second derivative, so the Hessian is
checked too: at `ka == ke` exactly, the tape's Hessian of the
`ss`-plus-`addl` predictor agrees with a central difference of the
tape's own gradient to **5.53e-11**, with no `NaN`.

Sweeping the gradient through the collapse, `ka - ke` from 1e0 down to
0 over sixteen separations, the tape and a central difference agree to
at worst 2.03e-10 and at best 1.9e-12. The offset-only version was
wrong by 3.1e-03 at a separation of 1e-14 and by 2e+285 at zero.

## Timing

### The Phase 0 design, before and after

`dev/lincmt/lincmt-phase0.R`, seed 20260908, the design and starting
values of `extensions/frmtmb.ode/tests/testthat/test-scale.R`
unchanged: 100 subjects x 8 samples, a depot and a central compartment,
fourteen doses over seven days as one `ss` row plus `ii`/`addl`,
`frm(se = TRUE)`.

TWO passes, all three arms inside one process each. The raw log of the
second is `dev/lincmt/lincmt-phase0-run.log`.

| arm | pass 1 | pass 2 | minimum | `nlminb` | max abs gradient |
|---|---|---|---|---|---|
| `frm_ode()`, the Phase 0 record | 3948 s | | | code 1 | 6.1e-03 |
| `frm_ode()` | 3716.6 s | 5374.8 s | **3716.6 s** | code 1 | 6.106e-03 |
| `frm_lincmt(n_ss = 20)` | 119.1 s | 118.5 s | **118.5 s** | code 0 | 4.653e-04 |
| `frm_lincmt(n_ss = Inf)` | 34.8 s | 55.0 s | **34.8 s** | code 0 | 4.653e-04 |

**F7, corrected, and this was the headline.** An earlier version of
this document derived the speedup "WITHIN a pass so that both arms met
the same machine". That derivation is wrong and the review is right
about why: `dev/lincmt/lincmt-phase0.R` runs the arms SEQUENTIALLY, not
interleaved. The kept log shows the two closed-form arms occupying
11:24:22 to 11:27:16 and the `frm_ode()` arm 11:27:16 to 12:57:02, so a
within-pass ratio compares two blocks ninety minutes apart and controls
for load no better than a cross-pass ratio does.

**What the four numbers support is each arm's slowest against the
other's fastest**, which is the only ratio that survives an adversary
choosing the load:

| | conservative bound |
|---|---|
| `n_ss = 20` | 3716.6 / 119.1 = **at least 31.2x** |
| `n_ss = Inf`, the default | 3716.6 / 55.0 = **at least 67.6x** |

The earlier "at least 98x" is withdrawn, and so is the "at least 68x"
that round 1's correction rounded it to: 3716.6 / 55.0 is 67.575, so
**67x** is what four numbers carry.

**Corroboration, from designs where the pairing IS sound.** At 5
subjects all three fits are cheap enough to run together and `nlminb`
takes IDENTICAL counts in all three arms (49 iterations, 59 `fn`, 50
`gr`), so the comparison is of equal work: the whole fits are 77.11 s,
1.64 s and 0.81 s, which is **47x and 95x**
(`dev/rev-lincmt-instrument.R 5 4 gr`, seed 20260908, the review's
run). At 25 subjects the counts are 33, 31 and 31 iterations and the
fits are 953.45 s, 24.63 s and 10.72 s, which is 38.7x and 88.9x
(`dev/lincmt/lincmt-instrument.R 25 4`). So the true factor at the
default is near 90 to 100, and the record's CLAIM is 67x because that
is what its own design proves.

The per-arm spread across the two Phase 0 passes is 1.45x on
`frm_ode()`, 1.005x on `n_ss = 20` and 1.58x on `n_ss = Inf`; a sibling
lane was on the machine for part of pass 2 and for none of pass 1.

**F8.** Pass 1 has no log in the tree; only pass 2 was kept. The Phase
0 design was not re-run to produce one, because the cost rule forbids
it and because pass 1's numbers enter the conservative bound only where
they WEAKEN it: 3716.6 is pass 1's `frm_ode()` time and is the smaller
of the two, so using it as the numerator makes the claimed speedup
smaller, and 119.1 is pass 1's `n_ss = 20` time and is the larger of
the two, so using it as a denominator does the same. Dropping pass 1
entirely and using only the logged pass 2 gives 45.4x and 97.7x, which
is a LARGER claim than the one made. The unlogged pass is load-bearing
only against the lane's own interest.

The `frm_ode()` arm REPRODUCES the Phase 0 row, which is what makes
this a comparison rather than two unrelated numbers: 3716.6 s against
its 3948 s, the same convergence code 1, and the same maximum gradient
to three figures, 6.106e-03 against its 6.1e-03.

Every estimate is identical to the digits `dev/scale-findings.md`
records, in both passes: `ka` 0.9841612 against its 0.984, `ke`
0.1482791 with the interval 0.1397548 to 0.1573233 against its 0.1483
and 0.1398 to 0.1573, `V` 19.62513 against its 19.63, and the two
subject standard deviations 0.2692418 and 0.2612086 against its 0.269
and 0.261. Between the two paths the log-likelihood differs by 1.8e-07
and the largest parameter difference is 3.2e-06, which is the solver's
default `atol = rtol = 1e-8` and not the closed form.

**The optimizer's verdict changes.** `nlminb` returns code 0 through
`frm_lincmt()` and code 1, its false-convergence flag, through
`frm_ode()`, with a maximum gradient thirteen times smaller. Both land
on the same answer; the closed form removes the solver tolerance's
noise from the objective, so the convergence test can be met. A user
reading `frm_ode()`'s output today sees a warning that the fit may not
have converged, on a fit that has.

### Tape build, at the Phase 0 size

`dev/lincmt/lincmt-time-scale.R 100 4`, one process, all arms
interleaved, four rounds, minimum taken. Tape build is
`frm(dry_run = "objective")` less `frm(dry_run = "frame")`, which is
`dev/scale-findings.md`'s own construction; its control arm reported
0.92.

| quantity | `frm_ode()` | `n_ss = 20` | `n_ss = Inf` |
|---|---|---|---|
| tape build | 28.35 s | 2.24 s | 0.98 s |

`dev/scale-findings.md` records 40.1 s for the `frm_ode()` tape build
on this row; 28.35 s here, on a less busy machine.

### Gradient timing, and two ways this lane's instrument was wrong

This is the part of the timing that did NOT settle cleanly, and the
honest report is that the whole-fit wall clock above is the number to
carry and the per-gradient numbers are supporting detail with a factor
of 1.4 of slack in them.

**First instrument error: repeating one parameter.** A block that calls
`obj$gr(par)` `n` times at the SAME `par` does not measure the gradient
an optimizer performs. TMB's Laplace objective caches the inner Newton
solution on the parameter, so calls two and later are served from that
cache. Fixed by giving every call in a block its own parameter, drawn
once as `par` plus N(0, 0.02) jitter, seed 1.

Measured, both instruments on the same objects, interleaved
(`dev/lincmt/lincmt-instrument.R 25 4`, 25 subjects, 4 rounds):

| arm | one point | own point | understatement |
|---|---|---|---|
| `frm_ode()` | 10.97 s | 12.73 s | 1.16x |
| `frm_lincmt(n_ss = 20)` | 0.2225 s | 0.2125 s | 0.96x |
| `frm_lincmt(n_ss = Inf)` | 0.08625 s | 0.09125 s | 1.06x |

At 5 subjects the same comparison puts the `frm_ode()` understatement
at 1.78x (0.445 s against 0.790 s), and the review gets 1.22x at that
size. So the effect is real, is confined to the arm with the expensive
inner solve, and is between 1.1x and 1.8x.
`tests/testthat/helper-scale.R`'s `scale_grad()`, which produced the
"one gradient" column of `dev/scale-findings.md` for every row, repeats
one point. **This is a caveat on that column, not a defect claim**: an
earlier version of this document said 2.2x from two runs that were not
interleaved, and the interleaved measurement does not support it. The
whole-fit column of that table is unaffected.

**Second instrument error: a duplicate object as the control.** The
control arm must be a second block of the SAME object, not a fourth
objective built from the same formula. Built as a duplicate, it
reported **4.409** where it must report 1.0: a process holding four
Laplace objectives does not time the fourth the way it times the
second. That is the control doing its job.

**What the fixed instrument reports.** `dev/lincmt/lincmt-time-small.R`
at 5 subjects, seed 20260908, block sizes grown by doubling past 1.2 s
(1 call for `frm_ode()`, 64 and 128 for the closed form, so the
smallest block is 64 x 0.0172 = 1.1 s, 110 ticks of the 10.0 ms clock),
minimum over rounds, two passes:

| arm | 6 rounds | 12 rounds | minimum |
|---|---|---|---|
| `frm_ode()` | 0.7900 s | 0.9800 s | 0.7900 s |
| `frm_lincmt(n_ss = 20)` | 0.01937 s | 0.01719 s | 0.01719 s |
| `frm_lincmt(n_ss = Inf)` | 0.009375 s | 0.008516 s | 0.008516 s |
| CONTROL | 0.01828 s | 0.01813 s | |
| **control / `n_ss = 20`** | **0.944** | **1.055** | |
| round spread, worst arm | 1.74 | 3.72 | |

The control brackets 1.0 within 6 percent in both directions and the
round spread reached 3.7 on the second pass, so the machine was busy
with a sibling lane throughout. The review ran the same instrument on
its own build and got a control of **1.036**. Taking the minimum per
arm across both of this lane's passes gives 46x and 93x; the two passes
individually give 41x and 84x, and 57x and 115x; the review gets 58.7x
and 125.2x. The honest statement is one gradient somewhere between 40x
and 60x cheaper at `n_ss = 20` and between 80x and 125x at the default,
on a machine whose own round-to-round spread is up to 3.7.

The whole-fit ratio and the one-gradient ratio do not have to agree,
and they do not: a fit is gradients times iterations plus one tape
build plus one `sdreport()`. At 25 subjects, where all three fits are
affordable, the iteration counts are nearly identical (33, 31 and 31
`nlminb` iterations; 33, 32 and 32 gradient evaluations), the whole
fits are 953.45 s, 24.63 s and 10.72 s, and the two closed-form arms'
fit ratio of 2.30 matches their gradient ratio of 2.33.

## Tape sizes, and the load-independent count that IS comparable

`dev/lincmt/lincmt-time-small.R` and `dev/lincmt/lincmt-solvecount.R`,
one subject's predictor, eight observations, the Phase 0 schedule:

| path | tape nodes | `RTMBode::ode()` calls |
|---|---|---|
| `frm_ode(n_ss = 20)` | 226 | **33** |
| `frm_ode(n_ss = 1)` | | 14 |
| `frm_lincmt(n_ss = 20)` | 26672 | 0 |
| `frm_lincmt(n_ss = Inf)` | 12296 | 0 |

**The node count alone is the wrong instrument, and the reason is worth
recording**, because it was asked for as a load-independent alternative
to the clock. `RTMBode::ode()` builds ONE `ADjoint` atomic node per
solve (`RTMBode:::ODEadjoint`), so `frm_ode()`'s 226 nodes are 226
nodes THIRTY-THREE of which each replay a full numerical solve of an
augmented sensitivity system on the reverse pass. The closed form's
26672 are 26672 multiplications and exponentials. Counted against each
other the two numbers say the closed form is 118 times larger, which is
the opposite of the truth.

**The comparable count is the solve count**, and it is exactly load
independent: 33 adjoint solves per group against 0. At the Phase 0
design's 100 groups that is 3300 numerical solves against 1.23 million
floating-point operations, and it is the whole story of the 3717 s
against 55 s. The plan's own text estimated "about 35 solves per group
per gradient" for this row; the measurement is 33.

**F9, two labels corrected.** The column header said "per objective
evaluation". The 33 calls happen once, at `MakeTape()`; a forward
replay of the tape and a reverse pass each make zero further
`RTMBode::ode()` calls. What is replayed on every evaluation is 33
`ADjoint` ATOMIC SOLVES baked into the tape, which is the same cost and
a different sentence. And "19 the run-in, 14 the dosing segments" is
one of two readings, because it measures against an `n_ss = 1`
baseline, which still contains one run-in cycle. Removing the `ss` row
entirely gives 13 calls, so the reading this document uses is **the
`ss` row costs 20 and the dosing segments cost 13**, because that
answers "what does the `ss` row cost" rather than "what does raising
`n_ss` from 1 cost". Both sum to 33.

The node count IS informative WITHIN the closed form: `n_ss = 20`
writes out twenty cycles per observation and costs 2.2x the nodes of
the exact limit, which is where the 119 s against 55 s comes from.

## What is refused, and how often that costs anything

Superposition adds impulse responses and never forms the state vector,
so anything that reads or sets that vector is refused. Each refusal
names itself and names `frm_ode()`.

`dev/lincmt/lincmt-coverage.R` enumerates the record types a
NONMEM-shaped population pharmacokinetic dataset produces, crossed with
the model shapes such an analysis fits, and CALLS `frm_lincmt()` on
each, so the verdict is the code's and not a reading of the help page.

**16 of 23 accepted.** The eighteen rows that are ordinary
population-pharmacokinetic records (EVID 1 bolus, with `rate`, with
`ii`/`addl`, with `ss`, EVID 3, EVID 4, a bioavailability on the dose,
a per-subject schedule, an amount at `t0`, a dosing history before the
first sample, one to three compartments, both parameterizations) give
**16 accepted and 2 refused**. The two are:

- an infusion into the DEPOT, which is zero-order absorption. Refused
  because its response is a three-node convolution with no
  cancellation-free spelling. This is the one refusal that costs a real
  model: sustained-release formulations are fitted this way, and the
  review agrees it is the one that matters.
- a dose into a PERIPHERAL compartment, which is not a route of
  administration.

The remaining five refusals are constructions written to be refused:
`method = "replace"`, `method = "multiply"`, a `"reset"` to a non-zero
level, `tv`, and an infusion still running when the schedule restarts.

The review's own check is the stronger one and it passed: over 140
random schedules from a grammar wider than this lane's hand list, plus
18 targeted corners, **no schedule was accepted by `frm_lincmt()` and
refused by `frm_ode()`**, and no accepted schedule disagreed by more
than 1.76e-11 of its own scale.

## What `frm_ode()` approximates and this does not

`dev/lincmt/lincmt-edge.R`, section 10. `frm_ode()` reaches a steady
state by simulating `n_ss` cycles from an empty system. On a two-
compartment oral model with `ke = 0.2`, `k12 = 0.4`, `k21 = 0.1`,
`ka = 1.1` and `ii = 8`, the default `n_ss = 20` is short of the true
steady state by:

| `n_ss` | `frm_ode()` vs the exact limit | `frm_lincmt(n_ss)` vs it | the two against each other |
|---|---|---|---|
| 10 | 4.31e-02 | 4.31e-02 | 3.44e-13 |
| 20 (the default) | 3.96e-03 | 3.96e-03 | 4.14e-13 |
| 40 | 3.34e-05 | 3.34e-05 | 4.34e-13 |
| 80 | 2.38e-09 | 2.38e-09 | 4.37e-13 |

The two truncate identically, to 4e-13, so this is the MODEL's
approximation rather than either implementation's arithmetic. What is
new is its size at the shipped default.

### What a user sees, exactly

The schedule is the one above, and nothing about it is unusual; it is
the package's own two-compartment test schedule.

**On the numeric path** (a direct call, `predict()`, `simulate()`, or a
body holding no estimated parameter) `frm_ode()` compares its last two
run-in cycles and warns that the state was "still moving by 0.0023
(relative), more than ss_tol = 1e-06". So a user who calls it
numerically IS told, and the help page's "Steady-state dosing" section
says to read that warning and raise `n_ss`. Two numbers, not one:
0.0023 is the cycle-to-cycle movement the warning reports, and
3.96e-03 is the distance to the limit that this lane measured against
the closed form.

**During a fit the check cannot run at all.** It is a comparison, and
RTMB refuses comparison on AD types, which `frm_ode()` says out loud in
the same warning. So the user who FITS this model sees nothing, and the
likelihood is evaluated at a state that is 0.4 percent low in every
group, in the same direction, at every observation. That is a bias, not
noise, and it goes into the estimates.

### Its own item? I said no. The review measured it and I was wrong

I argued for "not its own item, but a line in an existing one", on two
grounds: that it is documented, and that the numeric path warns. The
review tested both and both fail.

**Ground one, that it is documented, fails because the documentation is
wrong by an order of magnitude.** `?frm_ode` said "only a percent or
two for one whose half-life is many intervals long". A half-life of ten
dosing intervals is `lambda_z ii = 0.0693`, so the shortfall after
twenty cycles is `exp(-1.386) = 0.25`. Twenty-five percent, not a
percent or two. The sentence holds at about three intervals, not many.
**I have corrected that sentence**, with the measured numbers and the
rule below, because shipping a help page that is wrong by an order of
magnitude is not something to carry forward.

**Ground two, that the warning catches it, fails because the warning
understates it, and by more as the error grows.** `ode_run_in()`
compares the last two cycle-start states, so it reports
`r^(n-1) (1-r)` where the distance to the limit is `r^n`, with
`r = exp(-lambda_z ii)`. The ratio is about `lambda_z ii`. Measured
(`dev/rev-lincmt-nsswarn.R`): 0.8x at a 13.9 h half-life, then 3.1x,
6.6x and **10.8x** at 138.6 h. It is worst exactly where the error is
worst.

**And the shortfall is not 0.4 percent, it is unbounded.** It goes as
`exp(-n_ss lambda_z ii)`, so (`dev/rev-lincmt-nss.R`) it is 3.9e-03 on
the schedule I measured, 1.2e-02 at a 107 h half-life dosed daily,
9.0e-02 at 265 h daily and 2.9e-01 at 139 h twice daily. Monoclonal
antibodies, amiodarone and bisphosphonates live in the bottom rows. The
3.96e-03 I reported is the benign end of that table and I generalized
from it.

**It moves estimates.** Simulating from the exact steady state and
fitting the same data twice through `frm_lincmt()`, once at
`n_ss = Inf` and once at 20, so that the run-in truncation is the only
difference: two compartments oral, terminal half-life 107 h, `ii = 24`,
30 subjects, seed 101, moves `k21` by a **factor of 3.16** and `ke` by
-11 percent, at a log-likelihood 1.59 units apart. One seed and one
dataset, so a demonstration rather than a recovery study.

So it IS a silent wrong answer in Rule 3's sense: during a fit no check
can run, the bias is systematic and in the same direction in every
group, and its size is set by a quantity the user is not told to
compute.

### The item text, for whoever edits the plan

I did not edit `dev/extension-gaps-plan.md`. This is the row I would
add, in Phase 1, and it is not this branch's blocker: the defect is
pre-existing in `frm_ode()`, this lane neither caused it nor worsened
it, and finding it is what the closed form was for.

> | 1.0d | `frm_ode(n_ss =)` is a silent wrong answer during a fit. The
> steady-state run-in truncates at `exp(-n_ss * lambda_z * ii)`, where
> `lambda_z` is the SLOWEST disposition eigenvalue, so the shipped
> default of 20 is 3.9e-03 short on a fast two-compartment oral
> schedule, 1.2e-02 at a 107 h terminal half-life dosed daily and
> 2.9e-01 at 139 h dosed twice daily. On one 30-subject dataset
> simulated from the exact steady state it moves `k21` by a factor of
> 3.2 and `ke` by 11 percent. The `ss_tol` check cannot run on the tape
> at all, so a fit is never told; the numeric path's warning reports
> the cycle-to-cycle movement, which understates the distance to the
> limit by about `1 / (lambda_z * ii)`, measured at 3.1x, 6.6x and
> 10.8x as the half-life grows. `?frm_ode`'s "only a percent or two for
> one whose half-life is many intervals long" was wrong by an order of
> magnitude and is corrected as of this round; what is left is the
> DEFAULT and the diagnostic. Guidance to ship: `n_ss = 20` suffices
> only when the terminal half-life is under about 3.5 dosing intervals;
> choose `n_ss` so that `n_ss * lambda_z * ii >= 20`, which puts the
> shortfall at 2e-09; do not choose it from the warning, but compare
> `n_ss` against `2 * n_ss` numerically, or use `frm_lincmt()` for a
> linear model, whose default sums the series | `frmtmb.ode/R/ode.R` |
> a two-compartment oral fit at a 107 h terminal half-life recovers
> `k21` inside its own interval at the shipped default, or the default
> is raised and the cost of raising it is stated | 1 |

### What a raised default should be. The review's proposal

I left this open and said the choice was between a larger constant and
a fit-time diagnostic. The review answered it and measured both halves.
Recording it here as ITS proposal, not mine, and not building it.

**Raising the constant is the wrong lever, and it costed it.** On the
Phase 0 design the solve count is 33 per subject per evaluation, 20 of
them the run-in, so `n_ss` costs `(13 + n) / 33` of 3717 s. Covering a
1e-8 shortfall needs `n * lambda_z * ii >= 18.4`:

| `n_ss` | Phase 0 cost | covers a terminal half-life up to |
|---|---|---|
| 20, today | 3717 s | 0.75 dosing intervals |
| 50 | 7100 s | 1.9 intervals |
| 100 | 12700 s | 3.8 intervals |
| 266 | 31400 s | 10 intervals |

So a blind constant costs an order of magnitude to cover the models
that motivate the item, and charges it to every `ss` fit including the
majority that never needed it. No constant is both correct and
affordable.

**Its alternative is `n_ss = "auto"`.** `n_ss` has to be fixed before
the tape is built, but it does not have to be fixed BLIND: the starting
values are known then, and `ode_run_in()` already computes the last
cycle-start states on the numeric path. The contraction ratio

    r = |y[n] - y[n-1]| / |y[n-1] - y[n-2]|

estimates `exp(-lambda_z * ii)` for any system with a dominant mode,
without knowing the model, which matters because `frm_ode()` takes an
arbitrary `dyn`. So: run the run-in numerically for about five cycles
per group at the starting parameters, form `r`, take
`n_ss = ceiling(log(1e-9) / log(max(r)))` over groups, floor it at a
handful and cap it at a budget, and when the cap binds WARN once with
the number it wanted. That converts the silent wrong answer into a loud
one, which is Rule 3's ordering and the whole point of the item. Any
whole number keeps working as it does now.

**On the common case it is NEGATIVE cost.** Five numeric solves per
group, once, against a fit of thousands of evaluations. And where 20 is
more than enough it picks fewer: the Phase 0 design has
`ke * ii = 1.8`, so `r = 0.165` and `log(1e-9) / log(0.165) = 11.5`,
giving `n_ss = 12` against today's 20, which is 25 solves per
evaluation instead of 33 and makes that fit about **24 percent faster**
while becoming correct by construction.

**Its unmeasured risk, which the review states itself**: the parameters
move during the fit, so an `n_ss` chosen at the starting values can be
too small at the optimum. Its two mitigations, neither measured, are
headroom in the tolerance (1e-9 rather than 1e-8) and repeating the
numeric check at the converged estimates, where `diagnose()` has
somewhere to put the answer. A lane taking 1.0d should measure how far
`r` moves between the start and the optimum on the designs in
`dev/scale-findings.md` before choosing the headroom.

## Sampling

`dev/lincmt/lincmt-sample.R`. `frmtmb.sample::frm_sample()` refuses a
model containing `frm_ode()`, because RTMBode calls deSolve unguarded
and the chain aborts at warmup iteration 1, taking the session's rstan
state with it (`extensions/frmtmb.ode/dev/upstream/`). `frm_lincmt()`
calls neither. One chain, 400 iterations, 200 warmup, on a six-subject
oral model with a subject random effect: **returns in 8.2 s** with
`summary()` readable and the session intact (7.7 s on an earlier pass
against the same code, and 4.18 s for the review on a quieter machine).
Badly mixed at that length, which is what 200 post-warmup draws buys
and not a statement about the model.

### Which registry rows change

Precisely one row is ADDED and none is changed or removed.

| row | before | after |
|---|---|---|
| feature `frm_ode()`, kind `special`, key `frm_ode` | present | unchanged |
| rule `frm_ode() x frm_sample = refused` | present | unchanged |
| feature `frm_lincmt()`, kind `special`, key `frm_lincmt` | absent | **added** |
| any rule naming `frm_lincmt()` | absent | **still absent, and asserted** |

The absence of a rule is the claim, not an oversight, so
`test-compat.R` asserts `nrow(rules naming frm_lincmt()) == 0`. The
review planted a rule with `feature_a = "frm_lincmt()"` and
re-evaluated the test's own predicate: it gives 1 where the test
asserts 0, so the guard does not fail open.

**How this interacts with item 1.2.** Item 1.2 made `frm_sample()` read
the registry's `refused` rows before taping, and made it WARN rather
than pass when it holds a refused row it cannot match against a model.
Adding a feature with no rule adds nothing for that pre-flight to match
or fail to match, so 1.2's behavior is unchanged and its
`test-compat-preflight.R` in frmtmb.sample needs no edit. What changes
is the ADVICE: `?frm_ode`'s "Sampling an ODE fit" section offered two
ways forward, a maximum-likelihood fit or the three-patch series
against RTMBode. For a LINEAR compartment model there is now a third
and better one, and `?frm_lincmt` says so in a "Sampling" section.
Item 1.2's refusal on `frm_ode()` stays exactly as it is: the upstream
defect it guards is untouched by this lane, still reproduces against
RTMBode 1.0, and still needs the patch series.

## The post-fit surface

`dev/lincmt/lincmt-postfit.R`, seed 77, 8 subjects x 7 samples.
`fitted()`, `residuals()`, `predict()`, `predict(newdata = )`,
`simulate()`, `logLik()` and `diagnose()` all work through a
`frm_lincmt()` predictor with no special casing, and `diagnose()`
reports nothing on a converged fit.

## The suite, and the check

Every test file in its own R process, `NOT_CRAN=true`, against the
package installed in this lane's private library:

| file | pass | fail | error | skip | warn |
|---|---|---|---|---|---|
| `test-bracket-access.R` | 1 | 0 | 0 | 0 | 0 |
| `test-compat.R` | 15 | 0 | 0 | 0 | 0 |
| `test-lincmt.R` | 63 | 0 | 0 | 0 | 0 |
| `test-ode.R` | 70 | 0 | 1 | 0 | 0 |
| `test-ode-events.R` | 108 | 0 | 0 | 0 | 0 |
| `test-ode-nlf.R` | 4 | 0 | 0 | 0 | 0 |
| `test-ode-tv.R` | 27 | 0 | 0 | 0 | 0 |
| `test-scale.R` | 0 | 0 | 0 | 1 | 0 |
| total | **288** | **0** | 1 | 1 | 0 |

The one error is PRE-EXISTING and is an artifact of running a file
standalone rather than of this lane. `test-ode.R:128` calls
`testthat::local_mocked_bindings()` without `.package =`, which needs
the package context that `test_check()` and `pkgload::load_all()`
supply and a bare `test_file()` against an installed package does not.
Checked by construction: the SAME file at the base commit, run the same
way against the reference library, reports the identical PASS=70
ERROR=1, and under `R CMD check` (which runs `test_check()`) the whole
suite passes.

`R CMD check --as-cran`, pandoc and TinyTeX on PATH,
`_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`: **1 NOTE**, the expected
"Skipping checking math rendering: package 'V8' unavailable" on the
HTML manual. Tests 70 s, vignette rebuild 29 s.

## What I did not do, and what I could not settle

- **The peripheral compartments are not readable and not dosable.**
  Their impulse response from a depot is an order-three exponential
  divided difference, and I did not find a cancellation-free spelling
  of one that costs a bounded number of operations. Scaling and
  squaring gives one, and I costed it: about 300 operations per
  (observation, dose) pair for a 3x3 triangular matrix at a fixed 30
  squarings, which is 22M tape nodes per gradient on the Phase 0
  design. That is worse than the solver. It is refused instead.
- **Zero-order absorption is refused**, for the same reason, and it is
  the one refusal that costs a model the field writes.
- **An `addl` block is summed dose by dose**, not collapsed by the
  geometric formula I derived for it. The collapse would take the
  Phase 0 design from about 120 terms per subject to about 16, and I
  did not ship it because its finite-`n` form cancels: its numerator is
  a sum of four O(1) terms whose leading two orders vanish, so its
  relative error is about `eps / (lambda * ii)^2`, which is 2e-08 at
  `lambda * ii = 1e-04`. A drug whose half-life is thousands of dosing
  intervals reaches that. The infinite (`ss`) case has no such
  cancellation and IS collapsed, which is why `n_ss = Inf` is both
  exact and cheap. A design with hundreds of `addl` doses before each
  observation is where `frm_ode()`'s per-segment cost can still win;
  the crossover was not measured.
- **`n_ss` finite writes out that many cycles.** It exists to reproduce
  a `frm_ode()` fit exactly and is otherwise not worth its 2.2x in
  nodes.
- **One machine, and a busy one.** A sibling lane was on it for part of
  pass 2 and none of pass 1; contention can only make a wall clock
  longer, so the minima are the wall clocks and each arm's slowest
  against the other's fastest is the speedup.
- **The per-gradient timing did not settle to better than a factor of
  1.4.** Two instrument errors were found and fixed (repeating one
  parameter; a duplicate object as the control), and after that the
  control brackets 1.0 within 6 percent but the machine's own
  round-to-round spread reached 3.7. The whole-fit wall clock is the
  number to carry; the per-gradient numbers are supporting detail. A
  third attempt, at the Phase 0 size, was abandoned after 67 CPU
  minutes because it was buying a secondary number at the price of the
  primary one.
- **What default `n_ss` would be enough for `frm_ode()` was not
  swept**, by me or by the review. The rule `n_ss lambda_z ii >= 20` is
  for a user who knows `lambda_z`; what a DEFAULT should be is the open
  question, and it is in the item text above.
- **A `NaN` gradient at a double disposition eigenvalue is not fixed**,
  and neither is the silent value error in the region around it. The
  derivative EXISTS there (-9.2779361, resolved to eight significant
  digits by a central difference); the implementation returns `NaN`
  because it routes through the eigenvalues. The clamp that fixes that
  costs nothing at a genuine collision and is destructive in the wider
  region (2.8e+02 to 2.5e+06 of the peak on 40 draws), so a safe fix
  needs a relative-precision test on the cubic's coefficients, which is
  a comparison on an AD value. That is real work and not a punch
  round's. What ships instead is the stated domain: clean to a
  six-decade spread, and clean to twelve decades in the direction an
  optimizer actually travels.
- **The `NaN` at every disposition rate underflowed to zero is not
  fixed either, and the obvious fix was tried and rejected.**
  Offsetting `pp` removes the `NaN` and returns a finite answer at
  **eight-ninths** of the true trajectory, because the eigenvalues are
  then of order 1e-75 and the absolute guard in the coefficients is no
  longer negligible against their differences: the coefficients come
  out 4/9 and 4/9 rather than 1/2 and 1/2. A silent factor of **1.125**
  is worse than a loud `NaN`. The measurement is the review's
  `dev/rev-lincmt-n1.R`; an earlier version of this document said
  two-thirds and a factor of 1.5 and cited `dev/lincmt/lincmt-n1.R`,
  which shows the `NaN` and that ordinary rates are unaffected but does
  NOT apply the offset. The rejection stands; the number did not.
- **No recovery study.** One seed, one fit. Recovery at replicate
  scale is Phase 2's job, and the Phase 0 design's estimates are
  reported only to show that the two paths land on the same answer.

## What the plan should say. I did not edit it

I did not edit it. Item 1.0c's check reads "identity with `frm_ode()`
on the same schedule to 1e-8; the Phase 0 design's wall clock, before
and after. The before is 3948 s." Both halves are met. The row can be
marked DONE with:

> 1.0c | DONE. `frm_lincmt()` is the analytic one- to three-compartment
> solution with the same `events` grammar, refusing by name what
> superposition cannot carry. On the Phase 0 design, whose `frm_ode()`
> arm reproduced at 3716.6 s against the row's recorded 3948 s, the
> closed form fits in 55 s at its default and 119 s with the run-in
> written out: **at least 67x and at least 31x**, taking each arm's
> slowest measurement over two passes against the other's fastest,
> which is what the design supports because its arms run one after
> another rather than interleaved. At 5 and 25 subjects, where the same
> fit is cheap enough to repeat and `nlminb` takes identical iteration
> counts in every arm, the factors are 95x and 47x, and 89x and 39x, so
> the true figure at the default is near 100x. Every estimate matches
> to seven significant digits what `dev/scale-findings.md` records, and
> `nlminb` returns code 0 at a maximum gradient of 4.7e-04 where the
> solver arm returns code 1 at 6.1e-03. The identity holds at 2.6e-12
> of the trajectory's scale over 118 schedules against `frm_ode()` at
> `atol = rtol = 1e-12`, and the review's own 104 random schedules and
> 18 corners agree, worst 1.8e-11; against a 240-bit reference the
> error relative to the trajectory's maximum is 3.6e-10 on a
> seven-decade rate box, with every coalescence of rate constants swept
> down to exact equality. It found no defect in `frm_ode()`'s
> arithmetic; it did find that `frm_ode()`'s shipped `n_ss = 20`
> run-in is a silent wrong answer during a fit whose size is
> `exp(-20 lambda_z ii)`, which is item 1.0d. It also unblocks
> sampling: `frm_sample()`
> refuses `frm_ode()` and runs on `frm_lincmt()`.

Item 3.1, the NONMEM-shaped reader, is worth more now than Phase 0
thought: its objection was that it "hands a user a dataset that then
takes an hour to fit", and the fit is now a minute.

## Version bump: a recommendation, not a choice

I have NOT edited any version. `DESCRIPTION` still says 0.2.0 and
`NEWS.md` carries the bullet under a `(development version)` heading.

What the change is, so that whoever chooses can: a new exported
function, a new help page, a new vignette section, one added row in the
compatibility vocabulary, and a generalization of the frame check to
cover the new function. Nothing existing changes behavior. `frm_ode()`
has exactly two edits, both to prose: its within-group refusal message
now names the function that raised it, because two functions can now
raise it, and its wording of "frm_ode() solves one system per group"
became "one system is solved per group". No existing test needed
changing for either; the two tests that match that message match on
"not constant within", which is untouched.

By this package's own history (0.1.0 a first release, 0.1.1 a
dependency floor, 0.2.0 a registry row plus a measurement) that reads
to me as a MINOR bump, **0.3.0**. That is a recommendation. It is not
mine to make, and I have not made it.
