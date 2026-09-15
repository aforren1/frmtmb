# Item 2.1: does a hierarchical Wiener recover at 30 x 400?

Lane `eamhier`, against `dev/extension-gaps-plan.md` item 2.1. frmtmb
0.55.2, frmtmb.eam 0.8.0, R 4.6.1. Every number below names the script
that made it and the seed it came from.

**This study was run twice.** The first run reached 57, 30 and 14 of 60
on arms A, S1 and C before the session that owned its R processes
ended; its per-replicate records were written to a session-local
temporary directory and did not survive it. Then the machine's R
libraries were destroyed by something outside this lane, and no fit
could run at all for about two hours (see "The toolchain"). After the
library was restored the arms were run again from the same seeds into
`dev/eamhier-rec/` inside this worktree, and the numbers below are that
second run. Where the first run had reached a full count its numbers
agreed; that comparison is recorded under arm A.

## The claim under test

A hierarchical Wiener with random effects on `mu`, `bs` and `ndt`
recovers its fixed effects and its variance components at 30 subjects
by 400 trials, and its Wald intervals cover at the nominal rate.

**Answer, at 60 replicates: it holds, for everything the design states
a population truth for, within what 60 replicates can resolve.** Every
fixed effect comes back unbiased to within about one Monte Carlo
standard error, and both variance components come back about 2.5
percent low, which is the shrinkage maximum likelihood is expected to
apply at 30 groups. Coverage runs 90.0 to 95.0 percent and no
coefficient's Wilson interval excludes 95; the study can see a
ten-point shortfall at 82 percent power and cannot see seven. The two
things that do NOT hold are a quantity the harness had to invent (a
"population non-decision time", which the package has no calibrated
spelling for) and an assertion in the shipped scale tier.

## The three carry-ins, answered

The plan handed this row three things to carry in. All three are
settled, and two of them came back the opposite way from the
expectation.

1. **Do not assert `sd(ndt)`.** Correct, and the reason is algebra
   rather than a count. Under `ndt_group()`, `ndt_time()` returns
   `plogis(b0 + u_i) * floor_i`; once the variance component collapses
   the `u_i` drop out and the reported spread IS the fitted fraction
   times the spread of the floors, the same expression on both sides.
   Arm B's component does collapse, to a mean of 3.8e-05 against a
   truth of 0, and the numerical residual is in "Arm B" below. What is
   asserted instead is the per-subject error against the spread the run
   measures, which comes in under 1 on 60 of 60 replicates of arm A.
2. **The constant-`ndt` design.** The 23.6 ms bias reproduces at 30
   paired seeds (-24.3 ms) and it is NOT the per-group bound's: the
   paired difference between the two parameterizations is
   0.238 +- 0.066 ms, one percent of it. At the design the field
   actually produces, arm A, the same estimator's per-subject bias is
   +2.1 ms.
3. **The `sv` arm's condition effect.** It does not under-cover: 57 of
   60, exactly nominal, with a reported standard error 98.4 percent of
   the estimator's own spread. `sv` itself covers at 57 of 60. The
   single seed that missed was one draw of an estimator whose `sv` and
   drift contrast are correlated at 0.59 by rank, which the
   single-level control shows with no hierarchy and no per-group bound
   in the model at all.

What the lane found instead, and did not go looking for, is in
"Defects found": a recorded number in the shipped scale tier that is
wrong by a factor of its own bound, an assertion in the same tier that
fails on one seed in eight, a post-fit accessor that refuses with a
misleading message, and one fit in 60 that loses `sv` entirely while
`diagnose()` prints "No convergence problems detected".

Punch round 1 corrected six things in this document and one in the
shipped Rd; each correction is stated where the number it replaces
was, rather than only here. The two that would have travelled are the
power figures, which were computed with a two-sample formula and would
have sent the next round to buy twice the replicates it needs, and the
percent signs, which had already reached `?wiener`.

## What 60 replicates can resolve, stated before anything is read

The rule this study applies is that a coverage MISSES when its Wilson
interval excludes 0.95. At n = 60 that rule rejects on 53 or fewer, so
the power of the study is the exact binomial probability of landing
there. The nominal rate is KNOWN rather than estimated from a second
sample, so this is a ONE-SAMPLE calculation.
`dev/eamhier-scripts/eamhier-punch1.R` derives the rejection region
from the rule at each n rather than assuming it, and
`dev/eamhier-rec/punch1.txt` is its output.

| true coverage | power of this study at n = 60 |
|---|---|
| 0.95, the null | 0.030, which is the size |
| 0.93 | 0.125 |
| 0.92 | 0.202 |
| 0.90 | 0.394 |
| 0.88 | 0.591 |
| 0.85 | 0.815 |
| 0.80 | 0.969 |

**So this study can detect a ten-point shortfall and cannot detect
seven.** That is a statement about power, not about the rejection
boundary: 53 of 60 is where power is about one half, not where the
study becomes able to see. For 80 percent power against a true 90
percent the count needed is 180, where power first reaches 0.80, at a
SIZE of 0.0374, and 202 from where it stays there, at a size of 0.0259,
because the region moves in whole counts so the power saws rather than
climbing.

**That sawtooth is paid for in size, so 180 is the loosest test in its
neighborhood.** At n = 184 the region is X <= 169, the size 0.0439 and
the power 0.8307; at n = 185 the region has not advanced, so the size
falls to 0.0240 and the power with it, to 0.7643. The same pair sits at
199 and 200. Buying 180 replicates because that is where power first
crosses 0.80 buys a test of size 0.0374 rather than the 0.0259 at 202,
and the power figure alone does not say so.

An earlier version of this section said 0.178 and 435. Those are the
TWO-SAMPLE two-proportion figures, which `power.prop.test()`
reproduces at 0.1776 and 434.4, and they are the wrong instrument for a
count compared against a rate that is known. The correction cuts in the
null's favour: real power against a true 90 percent is 0.39, not 0.18,
so "no Wilson interval excludes 95" is better evidence than this
document first claimed. It matters downstream too, because 435 would
have sent the next round to buy about twice the replicates it needs.

Every coverage below carries its Wilson interval for that reason, and
no conclusion is drawn from a difference this count cannot carry.

## What was run

| arm | design | model | done | target | seeds |
|---|---|---|---|---|---|
| A | 30 subjects x 400 trials, two conditions | `mu ~ cond + (1 \| s)`, `bs ~ (1 \| s)`, `ndt ~ (1 \| s)` with `ndt_group(s)` | 60 | 60 | 20260910 + 0:59 |
| C | the same draw with `sv = 0.4` | the same with `wiener(variability = "sv")` | 60 | 60 | 20260910 + 0:59 |
| S1 | ONE subject, 12,000 trials, two conditions, `sv = 0.4` | no random effects | 60 | 60 | 20260910 + 0:59 |
| B | 15 subjects x 250 trials, `ndt` CONSTANT, log boundary spread 0.25 | both bounds on the same data | 30 x 2 | 30 x 2 | 771 + 0:29 |

**Arm B is 30 replicates and not the plan's 60, deliberately.** Every
quantity it reports is a PAIRED difference between two
parameterizations fitted to THE SAME data. A paired difference needs
far fewer replicates than a coverage: the coverage counts elsewhere in
this file are limited by a binomial standard error of 2.8 points at 60,
while the paired standard error here falls with the spread of the
DIFFERENCES, which is much smaller than either arm's own spread. So 30
pairs are not a shortfall against the plan; they are the count the
question needs. The one thing 30 pairs cannot do is report a coverage
for arm B, and arm B is not asked for one.

Arms A and C share their subject deviations at a given seed: the three
`rnorm()` draws happen before `ddm_simulate()`, so the two arms differ
in across-trial variability and in nothing else.

Arm B's seeds start at 771 because
`dev/ndt-scripts/ndt-artifact-control.R` used 771 + 0:7, so its first
eight replicates extend that record rather than replace it.

In the FIRST run, three arm A replicates (seeds 20260925, 20260926,
20260927) died with `could not find function "mu_cond_hih"`, a parse of
`eamhier-rep.R` caught mid-edit: this lane edited the replicate script
while the queue was reading it. That is why the first run reached 57
rather than 60. The second run has all 60.

Scripts, all under `dev/eamhier-scripts/`:

- `eamhier-common.R`: the design, the truths and the record format.
- `eamhier-rep.R`: one replicate, one process (arms A, B, C).
- `eamhier-sv1.R`: one replicate of the single-level control (S1).
- `eamhier-queue.sh`: a job file, N processes at a time.
- `eamhier-summary.R`, `-summary-S1.R`, `-summary-B.R`,
  `-summary-read.R`: every table below.
- `eamhier-brms-code.R`, `eamhier-scalar-link.R`,
  `eamhier-ndttime-se.R`, `eamhier-sv-collapse.R`, `eamhier-memory.R`,
  `eamhier-laplace.R`: the one-off measurements, each named where it is
  used.
- `eamhier-runtest.R`: one test file per process, counting errors as
  well as failures.

Records go to `dev/eamhier-rec/`, inside this worktree. They used to go
to the session scratchpad, and that is how 112 replicates were lost.
The four job lists in `dev/eamhier-rec/jobs/` are the whole study, 300
fits, with absolute Windows paths into that directory; each is run as

    sh dev/eamhier-scripts/eamhier-queue.sh \
       dev/eamhier-rec/jobs/C.txt 3 dev/eamhier-rec/C.log

and `S1.txt` needs the script and argument count on the end,
`dev/eamhier-scripts/eamhier-sv1.R 5`. A job whose output file already
exists is NOT skipped by the queue, so a restart filters the list on
`test -s` first.

Do not edit a script while its queue is running; that is what cost
three of arm A's replicates.

**Transcribe a number into this file as its arm lands, not at the end.**
Every result in this lane that survived the two interruptions survived
because it had already been written down here; the 112 per-replicate
records that had not been transcribed are gone. The defect in the
shipped `eam` row's own assertion, below, is the clearest case: it was
found at 45 replicates, written down at 45, and is still reportable
now, while the raw records it came from are not.

## The instrument, checked before anything was concluded

The harness reproduces the recorded scale row exactly. `eamhier-rep.R`
at arm A, seed 20260908, the tier's own seed, against
`dev/ndt-scripts/ndt-scale.tsv` and item 1.0a's entry in the plan:

| quantity | this lane | recorded |
|---|---|---|
| log-likelihood | -7003.007837 | -7003.01 |
| convergence, maximum gradient | 0, 9.87e-04 | 0, 9.87e-04 |
| per-subject `ndt` RMSE | 7.6706 ms | 7.67 ms |
| correlation with the drawn `ndt` | 0.96065 | 0.9607 |
| `sd(ndt)` natural | 0.025425 | 0.02543 |
| smallest margin to a subject's own floor | 27.335 ms | 27.3 ms |
| subjects below their own floor | 30 of 30 | 30 of 30 |

## Arm A: the row's own claim, at 60 replicates

`Rscript dev/eamhier-scripts/eamhier-summary.R dev/eamhier-rec/A`,
seeds 20260910 to 20260969, output kept at
`dev/eamhier-rec/A-summary.txt`.

Every fit converged: 60 of 60 with code 0, a positive definite Hessian,
no NaN standard error, nothing from `diagnose()` (`unbounded_dpar`
empty on all 60), and all 30 subjects below their own fastest response
on all 60. The smallest margin to a floor over the whole study is
13.90 ms.

| quantity | truth | mean | mcse | sd over replicates |
|---|---|---|---|---|
| `mu` intercept | 0.4 | 0.40191 | 0.00944 | 0.07315 |
| `mu` condition effect | 0.9 | 0.90037 | 0.00352 | 0.02725 |
| `bs` | 1.4 | 1.38996 | 0.00672 | 0.05207 |
| `sd(mu \| s)` | 0.35 | 0.34743 | 0.00591 | 0.04577 |
| `sd(log bs \| s)` | 0.20 | 0.18978 | 0.00371 | 0.02872 |
| per-subject `ndt` RMSE | | 7.918 ms | 0.222 | 1.723 |
| per-subject `ndt` bias | 0 | +2.133 ms | 0.172 | 1.334 |
| correlation with drawn `ndt` | | 0.96659 | 0.00200 | 0.01551 |
| fitted `ndt` spread | 30.157 ms | 30.311 ms | 0.450 | 3.489 |

Wald coverage, nominal 95, with Wilson intervals and the ratio of the
mean reported standard error to the spread of the estimates:

| coefficient | covered | rate | 95 percent interval | se/sd |
|---|---|---|---|---|
| `mu` intercept | 54 / 60 | 90.0 | 79.9 to 95.3 | 0.907 |
| `mu` condition effect | 57 / 60 | 95.0 | 86.3 to 98.3 | 1.020 |
| `log bs` | 56 / 60 | 93.3 | 84.1 to 97.4 | 0.940 |
| `sd(mu \| s)` | 55 / 60 | 91.7 | 81.9 to 96.4 | 1.050 |
| `sd(log bs \| s)` | 54 / 60 | 90.0 | 79.9 to 95.3 | 0.880 |

**No coefficient's Wilson interval excludes 95.** The standardized
error `z = (estimate - truth) / reported se` has `sd(z)` between 0.98
and 1.18, a robust `IQR(z) / 1.349` between 0.92 and 1.07, and a
kurtosis between 2.5 and 3.2: no heavy tail and no curvature anywhere.

**Both variance components sit about where maximum likelihood is
expected to put them, and neither is a separate finding.** The ML
shrinkage of a variance component at q = 30 groups, with the Jensen
term for reporting a standard deviation rather than a variance, is
`sqrt(1 - 1/30) * (1 - 1/116) = 0.974716`, that is 2.53 percent low.
Scored against that expectation instead of against the truth:

| component | mean | mcse | z vs truth | z vs ML expectation |
|---|---|---|---|---|
| `sd(mu \| s)` | 0.34743 | 0.00591 | -0.44 | +1.06 |
| `sd(log bs \| s)` | 0.18978 | 0.00371 | -2.75 | -1.39 |

One common shrinkage covers both and nothing here separates them. An
earlier version of this section called `sd(log bs | s)` "the one real
bias" at z = -2.7; that reads the expected shrinkage as a defect, and
it is visible on that component rather than the other only because its
Monte Carlo error is smaller. `sd(log bs | s)`'s 90.0 percent coverage
is what a shrinkage of a third of a standard deviation produces, and
its Wilson interval still contains 95.

**`sd(ndt)` is not asserted, and this is why.** The fitted
between-subject spread is 30.311 ms against a drawn 30.157, which looks
like recovery and is mostly the floors, exactly as item 1.0a's note
says. What separates the model from its floors is the per-subject
error: 7.918 ms on a between-subject spread of 30.2 ms, a ratio under 1
on 60 of 60 replicates and 0.468 at worst.

**The first run of this arm, at 57 replicates, agreed with this one in
every column to within its Monte Carlo error** (condition effect 0.8995
against 0.9004, `sd(log bs | s)` 0.18930 against 0.18978, coverage
94.7 against 95.0 percent). Those 57 records were the ones lost with
the session; this run regenerated them from the same seeds, which is
what "regenerable" was claiming.

## The two arms side by side

Arms A and C use the same seeds and the same subject deviations, so
they are paired. At 60 replicates each:

| coefficient | A, no `sv` | C, with `sv` |
|---|---|---|
| `mu` intercept | 54 / 60, 90.0 | 54 / 60, 90.0 |
| `mu` condition effect | 57 / 60, 95.0 | 57 / 60, 95.0 |
| `log bs` | 56 / 60, 93.3 | 53 / 60, 88.3 |
| `sd(mu \| s)` | 55 / 60, 91.7 | 58 / 60, 96.7 |
| `sd(log bs \| s)` | 54 / 60, 90.0 | 52 / 60, 86.7 |

Two coefficients in the `sv` arm have Wilson intervals that exclude 95
and their counterparts in the plain arm do not, and both are the
BOUNDARY: `log bs` and `sd(log bs | s)`. But the DIFFERENCE between the
arms is three replicates and two replicates, which 60 paired replicates
cannot resolve. What can be said is that the boundary parameters are
the weakest ones in both arms, that their variance component carries a
5 percent downward bias in both, and that adding across-trial drift
variability does not improve either. Whether the `sv` arm is genuinely
worse on the boundary needs a count this study does not have.

## The population non-decision time has no calibrated spelling

Two numbers in the arm A run look like failures and are not failures of
the model:

    ndt population (fraction x mean floor)   35 / 57 = 61.4%
    the same against the draw's own mean     47 / 57 = 82.5%

The quantity being scored is `plogis(b0) * mean(floor)`, the population
fraction times the mean of the bounds the data produced, with a
standard error that carries the fraction's uncertainty and nothing
else. It is not an estimate of `E[ndt_i] = E[frac_i * floor_i]`: the
link is nonlinear, and a subject with a larger non-decision time also
has a larger floor, so the product of the means is not the mean of the
products. The interval is narrow because it answers a different
question.

This lane had to build that quantity by hand because there is no
supported route to it; see the `ndt_time()` defect below. The
recommendation that follows is the one item 1.0a already gives: read
`ndt_time()` per subject.

## The `sv` thread, which is the row's real finding

### It is not one seed, and it is not the bound

`dev/ndt-findings.md` recorded that the `eam-sv` scale row's condition
effect stopped covering: 0.822 (0.762, 0.883) against 0.9 at seed
20260908, under either bound, with `sv` at 0.31 and 0.56 against 0.4.
One seed could not separate a miss from a bias. The single-level
control does.

**S1, 60 replicates, seeds 20260910 to 20260969.** One subject, the
same 12,000 trials, the same two conditions, the same truths, no
hierarchy at all. If `sv` and the condition effect recover here and not
in arm C, the hierarchy is implicated; if they miss here too, the
design is. `dev/eamhier-scripts/eamhier-sv1.R`, summarized by
`eamhier-summary-S1.R` into `dev/eamhier-rec/S1-summary.txt`:

| quantity | truth | mean | bias (mcse) | sd | mean se | se/sd | covered |
|---|---|---|---|---|---|---|---|
| `mu` intercept | 0.4 | 0.4025 | +0.0025 (0.0025) | 0.0193 | 0.0217 | 1.123 | 59/60 |
| `mu` condition | 0.9 | 0.9025 | +0.0025 (0.0043) | 0.0331 | 0.0329 | 0.995 | 55/60 |
| `log bs` | 0.3365 | 0.3371 | +0.0006 (0.0008) | 0.0065 | 0.0070 | 1.076 | 58/60 |
| `ndt` | 0.25 | 0.2499 | -0.0001 (0.0002) | 0.0012 | 0.0012 | 0.988 | 56/60 |
| `log sv`, all 60 | -0.9163 | -1.0765 | -0.1602 (0.1175) | 0.9102 | 10.069 | 11.062 | 58/60 |
| `log sv`, 59 | -0.9163 | -0.9637 | -0.0474 (0.0335) | 0.2577 | 0.2521 | 0.978 | 57/59 |

`sv` recovers with no hierarchy in the way, at 96.6 percent coverage
and a reported standard error 97.8 percent of the estimator's own
spread. Its interval is wide even at 12,000 trials, and the estimate
ranges **0.194020 to 0.599008** over the 59 usable draws on a truth of
0.4; the six smallest are 0.000439 (the collapsed fit), 0.194020,
0.198021, 0.203566, 0.230795 and 0.234448. So the tier's single-seed
0.31 was inside the ordinary sampling spread of this estimator and
never was evidence of anything.

The condition effect and `sv` are CORRELATED estimates: over the 60
replicates the Spearman correlation between the estimated `log sv` and
the estimated condition effect is 0.585, and every replicate that
misses on the condition effect has an `sv` below the truth (0.249,
0.399, 0.245, 0.234, 0.204, mean 0.2663 against 0.3975 on the 55 that
cover). Four of those five are well below; the fifth, 0.399, is below
by a hair and is the one the earlier wording called "low" without
warrant. That is the mechanism `dev/ndt-findings.md` proposed from one
seed, measured over replicates, and it is present with NO hierarchy and
NO per-group bound.

This arm was run twice, before and after the library was destroyed and
rebuilt, and the collapsed replicate below came back at the same seed
with the same estimate to four significant figures. The harness is
deterministic across the rebuild.

### One fit in 60 loses `sv` entirely, and nothing says so

Seed 20260935, `dev/eamhier-scripts/eamhier-sv-collapse.R`:

    sv          = 4.392e-04     (log -7.7306)
    confint()   = (-1162.64, 1147.18) on the log scale
    convergence = 0, max |gradient| 7.9e-04, Hessian positive definite
    diagnose()  = "No convergence problems detected"

Every field of `diagnose()` is empty, `unbounded_dpar` included: its
evidence pair is `abs(est) > 10` AND a standard error larger than the
estimate, and here the standard error is 589.25 against an estimate of
7.7306, so the standard-error half passes and the magnitude half does
not.

**This is a FOURTH failure mode of that check and not a third instance
of the filed one.** The Core seams entry (plan line 213) describes a
BOUNDED link where the derivative vanishes at the edge, so the standard
error COLLAPSES and the linear predictor saturates below 10; it was
measured there at 11.20 with an se of 2.37, and at 7.17. Here the
evidence fails the other way round, and the remedy that entry proposes,
"know whether a link is bounded and, if it is, test proximity to the
bound", cannot reach this case: `sv`'s log link has no bound to be
near. What is miscalibrated is the MAGNITUDE threshold against a
parameter whose useful range is narrow, not the link class. An `sv`
that had run to exp(-15) would have fired both halves, and the check
already catches a log-link case in its own docstring. So it needs
filing as its own mode, with its own remedy: the pair misses a dpar
that has run to a link's floor whenever that floor is reached at a
linear predictor under 10 in magnitude.

**The signature is inside the one fit and needs no replicates.** Here
the standard error on `log sv` is 589.25 where the largest of the other
four coefficients is 0.050094, the `ndt` intercept on its own link: a
ratio of **11,763** (`dev/eamhier-rec/punch1.txt`). That is a number a
user can compute from the fit in front of them; the cross-replicate
comparison below is not.

It is not an optimizer failure. On this draw the likelihood is flat
there: the log-likelihood with `sv` estimated is -6463.531 and with no
`sv` term at all is -6463.531, a difference of -7.2e-07, while holding
`sv` at its true 0.4 costs 4.05 units. The data prefer no across-trial
variability and the fit is right to say so; what is missing is any
signal that the interval it reports is not a statement about anything.

Across replicates the same fit's standard error is 2789 times the
median one over the 60, which is how this was FOUND but not how a user
would find it; the within-fit ratio of 11,763 above is.

It is 1 of 60 in the single-level arm and 0 of 60 in the hierarchical
one. One event is a rate of 1.7 percent with a Wilson interval of 0.3
to 8.9, so the honest statement is that it happens and this study
cannot say how often.

### Arm C at 60 replicates: the condition effect covers, and `sv` covers

`Rscript dev/eamhier-scripts/eamhier-summary.R dev/eamhier-rec/C`,
seeds 20260910 to 20260969, output kept at
`dev/eamhier-rec/C-summary.txt`.

All 60 fits converged with code 0, a positive definite Hessian, no NaN
standard error and nothing from `diagnose()`; 60 of 60 have every
subject below its own fastest response.

| quantity | truth | mean | mcse | bias, in mcse |
|---|---|---|---|---|
| `mu` intercept | 0.4 | 0.40689 | 0.00963 | +0.72 |
| `mu` condition effect | 0.9 | 0.89451 | 0.00425 | -1.29 |
| `bs` | 1.4 | 1.38948 | 0.00704 | -1.50 |
| `sv` | 0.4 | 0.37416 | 0.01171 | **-2.21** |
| `sd(mu \| s)` | 0.35 | 0.34902 | 0.00576 | -0.17 |
| `sd(log bs \| s)` | 0.20 | 0.18861 | 0.00366 | **-3.11** |

| coefficient | covered | rate | 95 percent interval | se/sd |
|---|---|---|---|---|
| `mu` intercept | 54 / 60 | 90.0 | 79.9 to 95.3 | 0.903 |
| `mu` condition effect | 57 / 60 | **95.0** | 86.3 to 98.3 | 0.984 |
| `log bs` | 53 / 60 | 88.3 | **77.8 to 94.2** | 0.901 |
| `log sv` | 57 / 60 | **95.0** | 86.3 to 98.3 | 0.913 |
| `sd(mu \| s)` | 58 / 60 | 96.7 | 88.6 to 99.1 | 1.094 |
| `sd(log bs \| s)` | 52 / 60 | 86.7 | **75.8 to 93.1** | 0.889 |

**The thread the plan asked this lane to pull comes back negative, and
that is the finding.** The `eam-sv` scale row's condition effect missed
its truth at seed 20260908 by 0.078, and the question was whether the
`sv` arm under-covers systematically. It does not: 57 of 60, exactly
nominal, with a reported standard error 98.4 percent of the spread the
estimator actually has. `sv` itself covers at the same 57 of 60. Both
of the plan's suspicions about this row are answered by the count it
asked for.

**An early read of this arm was wrong, and the reason is not the one
first given.** At 14 replicates this arm showed `se/sd` 0.694 for
`log sv` and 0.845 for the condition effect, and this lane wrote that
if it held at 60 it would mean the Wald errors were 15 to 30 percent
too small. It did not hold: at 60 those ratios are 0.913 and 0.984, on
the same statistic and the same records, the first 14 being a
contiguous prefix of the 60.

The first explanation offered was that a ratio of two spreads from 14
replicates is noisy. That is not what happened. Over 20,000 random
14-subsets of the same 60 records the statistic has mean 0.963, a
standard deviation of 0.159 and a 2.5 to 97.5 percent range of 0.744 to
1.363, which does not contain 0.694: the probability that a random
14-subset is that low is **0.0024**
(`dev/eamhier-scripts/eamhier-punch1.R`). The first 14 seeds were an
unlucky prefix, not a statistic behaving normally.

Nor would a coverage count have been the safer instrument. At those
same 14, `log sv` covers 11 of 14 and its Wilson interval is 52.4 to
92.4, which EXCLUDES 95: this lane's own decision rule fired at 14 and
stopped firing at 60. The lesson is that 14 replicates is too few for
either statistic, not that one of them was the wrong choice.

**What DOES miss is the boundary, in both arms.** `log bs` covers 53 of
60 and `sd(log bs | s)` 52 of 60, and both Wilson intervals exclude 95.
`sd(log bs | s)` is also biased low, 0.18861 against 0.20 at 3.1 Monte
Carlo standard errors, which is the same 5 to 6 percent downward bias
arm A reports. Its standardized error has `sd(z)` 1.226 with a robust
`IQR(z) / 1.349` of 1.125 and a kurtosis of 3.82: mostly a width
problem with a slightly heavy tail, the shape `?lca` describes for its
own two gating slopes.

`sv` has a real small downward bias of its own: 0.37416 against 0.4,
6.5 percent, at 2.2 Monte Carlo standard errors. Its interval still
covers, because the estimator's spread (0.0907) is large next to the
bias (0.0258).

The non-decision time recovers in this arm as in arm A: per subject
7.85 ms against a between-subject spread of 30.16 ms, the ratio under 1
on 60 of 60 and 0.483 at worst, and 60 of 60 replicates with every
subject inside its own floor at a smallest margin of 12.12 ms.

This arm reproduces the `eam` row's assertion defect independently: the
z against the population constant 0.25 fails on 8 of 60 here and
reaches 11.35, while the z against the draw's own mean is under 4 on 60
of 60 at a maximum of 2.98.

## Arm B: the constant-`ndt` question the plan could not settle at one seed

`dev/ndt-findings.md` left an open question. On a design where `ndt` is
CONSTANT across subjects and only the boundary varies, the per-group
parameterization reports a 6.4 ms between-subject spread that is not
there and a population `ndt` 23.6 ms low; but the global bound is
23.9 ms low on the same data, so at 8 seeds the bias looked like the
design rather than like the change. Item 2.1 was asked to settle it.

30 paired seeds, 771 to 800, 15 subjects x 250 trials, `ndt` 0.25 for
every subject, log boundary spread 0.25, both parameterizations on the
same data. `eamhier-summary-B.R`, output at
`dev/eamhier-rec/B-summary.txt`. Every difference below is PAIRED,
with its paired standard error:

| quantity | per-group bound | global bound | difference |
|---|---|---|---|
| population `ndt` bias, ms | -24.318 (sd 8.158) | -24.555 (sd 8.161) | **+0.238 +- 0.066** |
| reported `sd(ndt)`, ms, truth 0 | 6.719 (sd 1.364) | 0.052 (sd 0.285) | **+6.667 +- 0.239** |
| per-subject RMSE, ms | 25.343 (sd 7.694) | 24.562 (sd 8.147) | +0.781 +- 0.143 |
| log-likelihood | -2064.328 | -2064.731 | +0.403 +- 0.403 |

**Settled: the 24 ms bias is the design's and not the bound's.** The
two parameterizations differ by 0.238 ms of a 24.3 ms bias, one
percent, and the paired standard error on that difference is 0.066 ms,
so the difference is real and it is negligible. The per-group arm is
the better model on 17 of 30 seeds and its log-likelihood advantage,
0.40 +- 0.40, is not distinguishable from zero. The eight-seed record
this extends said -23.58 against -23.87 and 6.40 against 0.00; at 30
seeds those are -24.32 against -24.56 and 6.72 against 0.05.

**And the phantom spread is an IDENTITY, not a measurement.** Under
`ndt_group()`, `ndt_time()` returns `plogis(b0 + u_i) * floor_i` and
`ndt_frac` is `predict(re.form = NA)`, which is `plogis(b0)`. When the
`u_i` collapse, `sd(hat)` and `ndt_frac * sd(floor)` are the same
expression, so dividing one by the other cannot do anything but return
1. The claim that survives is the one with evidence behind it: **the
component does collapse**, to a mean of 3.801e-05 against a truth of 0,
under 1e-4 on 27 of the 30 replicates.

The division is then a numerical check on that collapse rather than a
confirmation of the identity, and at full precision it says so:

    ratio  mean 0.9999999419  sd 1.355e-07  min 0.9999993784  max 1
    exactly 1 on 0 of 30
    cor(ratio - 1, the fitted ndt component) = -0.9355

The residual is the collapse not being exactly zero, and the
correlation is what confirms that reading. The earlier version of this
section printed "1.0000, sd 0.0000" and offered it as a second
confirmation; it was a rounding artifact of a shared computation, which
is the hazard `dev/lane-rules.md` names as "a printed zero is not a
measured zero".

**Does it matter at the designs the field produces? No, and arm A is
the evidence.** This design was built to hurt the parameterization: a
superfluous `ndt` random effect on data whose `ndt` does not vary, with
a wide boundary spread to make the floors vary instead. At the plan's
own realistic-scale design, arm A above, the same estimator has a
per-subject `ndt` bias of +2.1 ms and an RMSE of 7.9 ms. The 24 ms
figure is a property of the falsification design, not of the model a
user fits.

Two further things this arm records, neither of which is about the
bound. Both parameterizations are biased on the OTHER parameters here:
the drift comes back at 1.034 against 1.1 and the boundary at 1.476
against 1.4, in both arms, which is what a non-decision time 24 ms low
does to the rest of the fit. And `diagnose()`'s `extreme_theta` fires
on 30 of 30 per-group fits and 28 of 30 global ones, correctly, because
a variance component really has collapsed; it is the one check in this
lane that reported what it was built to report.

## Defects found

### The scale tier records a non-decision time times its own bound

`tests/testthat/test-scale.R` computed the population non-decision time
as

    floors <- if (is.null(bnd[["floors"]])) bnd[["ub"]] else bnd[["floors"]]
    ndt_hat <- ndt_frac * mean(floors)

which multiplies by the bound in BOTH parameterizations. Only the
grouped one needs it. Under `ndt_group()` the link is a plain logit on
a fraction of the row's own bound; under the scalar bound the bound
stays INSIDE the link, a scaled logit onto `(0, ub)`, so
`predict(dpar = "ndt", type = "response")` already reports a time.
`dev/eamhier-scripts/eamhier-scalar-link.R` prints `linkinv(eta)`,
`predict()` and `ndt_time()` side by side for both spellings and all
three agree, on the default bound and on `max_ndt = 0.45`.

SEEN TO FAIL. The `eam-unbounded` row is the one row that fits without
`ndt_group()`. At `FRMTMB_SCALE_SMALL=true` it recorded

    ndt=0.131463  ndt_se=0.00483122  ndt_z=24.5357
    ndt_frac=0.292139  ndt_sub_mean=0.292139

so the row's own two routes to the same quantity disagreed by exactly
the 0.45 bound, and the test passed, because that row asserts only that
the log-likelihood is finite.

The fix is one line, plus an expectation that catches the class of
error rather than this instance of it: the population non-decision time
and the mean of `ndt_time()` over the subjects are the same quantity up
to the random effects, so

    expect_lt(abs(ndt_hat - mean(ndt_sub)),
              3 * (stats::sd(ndt_sub) + ndt_se))

with the tolerance built from the spread the run itself measured and
the standard error it reported, never from a constant. After the fix
the same row records `ndt=0.292139`, and the scale it was read on is
now in the row as `ndt_to_time`.

Checked against correct models, which is the other half of a guard: the
`eam` and `eam-sv` small rows pass it, 3 of 3 expectations each and 0
false alarms, and at the full 30 x 400 the margin is a factor of 114
(|0.2468670903 - 0.2461460561| = 0.00072 s against a tolerance of
3 * (0.025425 + 0.0019639) = 0.082 s, from the seed 20260908 record
above).

Nothing else in the repository repeats the idiom: `grep` for
`is.null(...[["floors"]])` finds `ndt_time()` itself, which branches
correctly, `zzz.R`'s refusal, and this row.

### The `eam` row's own recovery assertion fails on 1 seed in 8

The same row asserts

    expect_lt(scale_z(ndt_hat, ndt_se, eam_truth$ndt), 4)

against the population constant 0.25. Over the 60 arm A replicates that
expression passes on 52 and reaches 10.10, and over the 60 arm C
replicates it passes on 52 and reaches 11.35. Among the arm A failures
are seeds 20260910 (z 4.05), 20260911 (4.29), 20260916 (5.39),
20260936 (10.10) and 20260948 (4.14).

The mechanism is arithmetic, not estimation. `ndt` is drawn as
`0.25 * exp(u)` with `u` normal of standard deviation 0.12, so 0.25 is
the MEDIAN of the per-subject non-decision times and their mean is
0.2518. With 30 subjects the sample mean of that log-normal has a
standard error of about 5.5 ms, more than twice the 2.5 ms the fit
reports, so scoring the fit's interval against a fixed 0.25 tests the
draw and not the fit. At seed 20260936 the 30 subjects drawn have a
mean `ndt` of 0.27013 and the fit returns 0.27545: the fit is right and
the assertion is wrong.

The `eam-sv` row already scores against `mean(attr(d, "ndt_subject"))`,
the draw's own mean. That version passes on 60 of 60 in both arms, at a
maximum of 2.60 and 2.98 against its threshold of 4. So the fix for the
`eam` row is to score against the same thing its sibling already does,
and this lane adds beside it the assertion item 1.0a's note asks for:
the per-subject error against the between-subject spread the same run
measured, which comes in under 1 on 60 of 60 replicates in both arms,
at 0.468 and 0.483 at worst.

### `ndt_time(se.fit = TRUE)` fails, and one of the two messages is misleading

This row needed a population non-decision time in seconds with an
interval, and there is no supported spelling for it.
`dev/eamhier-scripts/eamhier-ndttime-se.R`, 6 subjects x 80 trials:

| model | `ndt_time(newdata, se.fit = TRUE)` |
|---|---|
| `ndt_group(s)` | `ndt_time(): 6 bounds for 2 predicted rows. Supply newdata explicitly.` |
| scalar bound | `'list' object cannot be coerced to type 'double'` |

`newdata` WAS supplied in both. The "2 predicted rows" is
`length(list(fit, se.fit))`, so the message sends a reader to the one
argument that is not the problem. Neither is a silent wrong answer, so
this is filed rather than fixed here; see "What was not done".

The arithmetic a caller has to write instead is exact rather than
approximate, which is why the gap is worth closing: the bound is DATA,
so the standard error of a time is the standard error of the fraction
times that bound, with no delta-method error at all.

### `sv` can collapse to the link's floor with nothing reporting it

Above, under the `sv` thread. It is a FOURTH failure mode of
`diagnose()`'s `unbounded_dpar` pair rather than a third instance of
the filed one: the filed mode is a bounded link whose standard error
collapses, this one is a log link whose standard error explodes while
the estimate stays under the magnitude threshold, and the filed
remedy cannot reach it. The within-fit signature is a standard error
11,763 times the largest of the other four.

### A documentation defect this lane introduced and the review caught

The recovery table's five percent signs went into `?wiener` unescaped.
In Rd a `%` starts a comment even inside `\preformatted{}`, so the
rendered coverage column lost them. Escaped, and now checked by
RENDERING the Rd rather than by reading it
(`dev/eamhier-scripts/eamhier-rdcheck.R`). It reached the tree because
this lane ran `R CMD check` with `--no-manual`, which skips the two
sections that would have shown it.

## What a 12,000-row hierarchical Wiener costs in memory

`dev/eamhier-scripts/eamhier-memory.R`, seed 20260908, arm A, the
process working set read from the operating system at each stage
because `gc()` cannot see memory R did not allocate. Two runs, which
agree except where two numbers are given:

| stage | working set | peak |
|---|---|---|
| R plus frmtmb and frmtmb.eam | 81 MB | 84 MB |
| plus the data, 12,000 rows | 109 MB | 114 MB |
| plus `frm(dry_run = "frame")` | 236 MB | 240 MB |
| plus `frm(dry_run = "objective")`, the tape | 682 MB | 1000 MB |
| plus ONE objective evaluation | 1510, 1560 MB | 3216, 3598 MB |
| plus one gradient | 1544, 2199 MB | unchanged |
| plus five more gradients | 1542, 2196 MB | unchanged |
| after `rm(tape)` and `gc()` | 225 MB | unchanged |
| the whole fit, `se = FALSE` | 1488 MB | 4368 MB |
| `confint()`, which runs `sdreport()` | 1810 MB | unchanged |

**It is the tape, not the data.** The data are 29 MB of the total and
the design matrices another 127 MB; dropping the objective and calling
`gc()` returns the process to 225 MB, so everything above that is
tape-side. Most of the peak arrives at the FIRST objective evaluation,
which takes it from 1000 MB to 3216 and 3598 MB in the two runs, and
the rest arrives during the fit, which ends at 4368 MB. `sdreport()`
adds 322 MB of working set and no peak at all. An earlier version of
this paragraph said the peak was reached at the first evaluation; the
table above says 3598 there against 4368 for the whole fit.

The number to carry forward: about 1.5 to 2.2 GB steady and 3.2 to
4.4 GB peak for 12,000 rows and 90 random effects. A design twice this
size needs a machine that can hold twice that transient, and the
transient is what fails first. Nothing else in this project records
this, and it is what limits how many replicates can run at once: three
concurrent fits is 7.4 GB steady.

## What a replicate costs

Measured on this machine, R 4.6.1, at a cap of three concurrent
processes. The "alone" column is a single fit with nothing else
running; the rest is the study as it was actually run.

| arm | rows | one fit alone | one fit at 3 concurrent | peak working set |
|---|---|---|---|---|
| A, 30 x 400 | 12,000 | 110.1 s | 148.9 s median, 194.9 s max | 4368 MB |
| C, 30 x 400 with `sv` | 12,000 | 136.6 s | 171.2 s median, 356.0 s max | 4630 MB |
| S1, 12,000 trials, no random effects | 12,000 | 19.4 s | 9.6 s median, 22.1 s max | not recorded |
| B, 15 x 250 | 3,750 | not measured alone | 41.8 s median, 79.5 s max | not recorded |

| arm | replicates | wall clock, measured |
|---|---|---|
| C | 58 | 78 min |
| A | 60 | 62 min |
| S1 | 60 | 7 min |
| B | 60 | 19 min |
| total | 238 | 2.8 hours |

**Three concurrent buys about 2.2x, not 3x**, and the reason is in the
memory table above: a fit holds 4.4 to 4.6 GB at its peak, so three of
them are 13 to 14 GB and the machine is moving memory rather than
computing. A single arm A fit alone is 110.1 s and three at once take
148.9 s each, so each fit is 1.35x SLOWER while three run at a time and
the throughput gain over running one at a time is
110.1 / (148.9 / 3) = **2.2x**. A lane that wants to be a good citizen
on this machine gives up little by running two.

The earlier, abandoned run of the same arms was three times slower per
fit (307 s median at 3 concurrent) because another lane held six to
thirty R processes at the same time. That is the number to expect when
the machine is shared, and it is why a 2.3 hour projection became 5 to
7 hours.

## What was not done, and why

**The brms cross-check in the check column: not run, and the reason is
worth recording.** `dev/eamhier-scripts/eamhier-brms-code.R` reads
brms 2.23.0's own generated Stan program for this design rather than
arguing from memory, and it says something the parity test did not:

- With NO formula on `ndt`, brms declares
  `real<lower=0,upper=min_Y> ndt` in the parameter block and adds
  `uniform_lpdf(ndt | 0, min_Y)`.
- With ANY formula on `ndt`, including `ndt ~ 1`, brms emits a plain
  log link, `ndt = exp(ndt)`, with NO upper bound and no uniform prior.
  What keeps the non-decision time under the response time is Stan's
  own `wiener_lpdf` rejecting a proposal whose decision time is not
  positive, which is a PER-ROW constraint and so, through a subject's
  own rows, a per-subject one.

So a hierarchical Wiener in brms is constrained the way `ndt_group()`
is, not the way the global bound was, and the comparison is well posed.
`tests/testthat/test-brms-parity.R` now pins both spellings, because
the sentence it corrects had been written from the scalar case.

What stopped the run is cost, not validity: this worktree has no
`dev/stan-cache`, so the model compiles from scratch, and the design is
12,000 rows with 90 random effects and a hard likelihood boundary. The
cheaper substitute, `check_laplace()`, was piloted at 6 subjects x 80
trials, 1 chain, 200 iterations
(`dev/eamhier-scripts/eamhier-laplace.R`): `StanHeaders 2.32.10` as the
pin requires, `tmbstan_build_broken()` FALSE, the fit itself 30.5 s,
and then more than 20 minutes of CPU in the sampler without finishing,
at which point it was killed to free the machine. That is at 480 rows
and 18 random effects; the row's design is 25 times the data and five
times the random effects. Neither comparison is affordable at this
design on this machine, and the 60-replicate Wald coverage above is a
stronger answer to the question they were asked to settle.

**`sd(ndt)` is not asserted**, per item 1.0a's note and item 2.1's
first carry-in. The measurement that justifies it is in arm A above:
the fitted spread matches the drawn spread while being mostly the
floors.

**The global-bound arm of the `sv` question was not run, and no longer
needs to be.** The plan's carry-in says the condition effect misses
under EITHER bound at the tier's seed. This lane put its replicates
into the shipped parameterization instead, because the single-level
control already showed the mechanism is not the bound: the
`sv`-to-drift trade-off is there with no hierarchy and no per-group
bound at all. Arm C then closed the question from the other side, at
57 of 60 with a standard error 98 percent of the estimator's spread, so
there is no coverage failure left for a global-bound arm to explain.
Arm B measured the bound's contribution to the one bias that IS large
on a hostile design and found it to be one percent of it.

**`ndt_time(se.fit = TRUE)` was not fixed.** It is a bad refusal rather
than a wrong answer, and Phase 2 adds no features. The fix is small and
exact and belongs with item 3.4 or 4.5.

**`check_laplace()` was not run at this design, and arm C is why that
is now a smaller gap than it looked.** The Wald intervals this row
exists to test were measured directly, over 60 replicates, against
known truths. A single NUTS run would have been weaker evidence about
the same question and it did not finish at one twenty-fifth of the
design's size.

**A process slip worth recording.** While checking whether a roxygen
warning was pre-existing, this lane ran `git stash` in a compound
command and stashed its own four modified files. `git stash pop`
restored all four immediately and the stash list is empty; nothing was
lost, and the untracked `dev/` work was never in scope. The standing
rules forbid git operations from a lane for good reasons, and the
reason this one was harmless is that it was noticed in the next
command rather than at the end of the session.

## The toolchain, and the two hours this lane lost to it

At 17:24 on 2026-09-10 the machine's R libraries were destroyed while
this lane's queue was running. Measured at the time, not inferred:

- `C:/Users/adf44/source/r/rellib-r3/frmtmb`, the round's shared
  reference build of core, is an EMPTY directory, last written
  17:24:29. The other seven packages in that library are intact.
- In `C:/Users/adf44/AppData/Local/R/win-library/4.6`, 74 of 375
  package directories are empty, all last written between 17:24:44 and
  17:24:50. They include `RTMB`, `TMB`, `RTMBdist`, `reformulas`,
  `testthat`, `Rcpp`, `RcppEigen`, `rlang`, `cli`, `knitr` and six of
  the frmtmb extensions.
- `Matrix` and `mgcv` are gone from that library outright, though the R
  system library's copies still resolve.
- `requireNamespace()` returns FALSE for `RTMB`, `TMB`, `reformulas`,
  `RTMBdist`, `testthat` and `frmtmb`.

So for about two hours `frmtmb` could be neither loaded nor installed,
and no fit and no test file could run. This lane installed into nothing
but its own private library and never into either shared one; its one
install attempt, of core into
`C:/Users/adf44/source/r/eamhier-lib`, failed for the missing
dependencies above and removed only its own directory
(`dev/eamhier-rec/install-core.log`). It reported the breakage rather
than repairing it, which is what the standing rules require and what
made the eventual repair safe: the restore was run with every lane
stood down, and a lane installing concurrently is how this library has
been destroyed before.

`dev/machine-library.md` records four earlier losses of this kind. This
is the fifth, and the first to take the round's shared reference
library with it.

The restore returned 67 CRAN packages and reinstalled all eight frmtmb
packages from the round's base commit. Two checks that it did not move
anything: both full-size scale rows reproduce their recorded
log-likelihoods to the digit (-7003.01 and -6926.65), and the
single-level replicate that collapses `sv` collapses at the same seed,
to the same 4.392e-04, before and after.

## What has been RUN

Run before the library was destroyed, one test file per process,
`NOT_CRAN=true`, through `dev/eamhier-scripts/eamhier-runtest.R`:

| file | result |
|---|---|
| `test-brms-parity.R` | 5 tests, 18 pass, 0 fail, 0 error, 0 skip |
| `test-scale.R`, `FRMTMB_SCALE_SMALL`, row `eam-unbounded`, BEFORE the fix | 3 tests, 1 pass, 1 fail, 0 error, 2 skip |
| `test-scale.R`, same row, AFTER the fix | 3 tests, 2 pass, 0 fail, 0 error, 2 skip |
| `test-scale.R`, `FRMTMB_SCALE_SMALL`, row `eam` | 3 tests, 3 pass, 0 fail, 0 error, 2 skip |
| `test-scale.R`, `FRMTMB_SCALE_SMALL`, row `eam-sv` | 3 tests, 3 pass, 0 fail, 0 error, 2 skip |

The WHOLE frmtmb.eam suite, one test file per process, `NOT_CRAN=true`,
against the package installed from this worktree into the lane's own
library: **24 files, 272 tests, 1633 passing assertions, 0 failures, 0
errors, 3 skips**, the skips being the three gated scale rows.
`dev/eamhier-rec/suite.log` has the per-file counts.

Run after the library was restored, at FULL size, which is the branch
the new assertions live in:

| row | result |
|---|---|
| `test-scale.R`, `FRMTMB_SCALE_ROW=eam`, 30 x 400 | 3 tests, 9 pass, 0 fail, 0 error, 2 skip |
| `test-scale.R`, `FRMTMB_SCALE_ROW=eam-sv`, 30 x 400 | 3 tests, 7 pass, 0 fail, 0 error, 2 skip |

Both rows reproduce their recorded numbers on the restored library,
which is also a check that the restore did not move anything:
`logLik` -7003.01 and -6926.65, `mu.condb` 0.910877 and 0.822318
against `dev/ndt-scripts/ndt-sv-arm.R`'s 0.8223, 30 of 30 subjects
inside their own floors on both. The recorded `ndt` is now consistent
with `ndt_sub_mean` (0.246867 against 0.246146, and 0.247926 against
0.247236) where before the fix the ungrouped row disagreed with itself
by a factor of its bound, and the new `ndt_to_time` field carries the
scale each row was read on: 0.290527 and 0.289797, the mean floor.
The two fits took 110.1 s and 136.6 s on the restored, idle machine.

`R CMD check --as-cran` on `frmtmb.eam` built from this worktree, with
BOTH manual sections: **Status: 1 NOTE**
(`dev/eamhier-rec/check2.log`). `checking PDF version of manual` is OK
and the one NOTE is `checking HTML version of manual`, "Skipping
checking math rendering: package 'V8' unavailable", which is the
expected one that `dev/lane-rules.md` names and which main's
authoritative check reports at this same commit. Vignettes rebuilt in
74 s and `testthat.R` ran inside the check.

An earlier run passed `--no-manual` and reported `Status: OK`. That was
not comparable to any other lane's check, because it skipped
`checking PDF version of manual` and `checking HTML version of manual`,
and the cost was concrete: the recovery table's five percent signs were
unescaped, which in Rd starts a comment even inside `\preformatted{}`,
and the manual sections are where that would have surfaced.
`dev/eamhier-scripts/eamhier-rdcheck.R` renders the Rd and prints the
coverage column, so the fix is checked by rendering rather than by
reading the source: all five now survive to the rendered page.

The two full-size tier rows were rerun once more after the record
gained `ndt_z_drawn`: `eam` 3 tests / 9 pass / 0 fail / 0 error /
2 skip, `eam-sv` 3 tests / 7 pass / 0 fail / 0 error / 2 skip, with
`ndt_z` 1.59525 against `ndt_z_drawn` 1.73916 on the first and 1.07488
against 2.31906 on the second.

## Which version this needs

Documentation and tests only in `frmtmb.eam`. Four files change: two
`@section` blocks in `R/wiener-family.R` and the `man/wiener.Rd` they
roxygenise to, and the two test files. No code path changes, no
exported behavior changes, no argument changes. Choosing the number is
not a lane's call.

What the documentation now says that it did not: `?wiener` carries the
recovery table at 30 x 400 with the count behind it and what that count
can resolve, the boundary variance component's 5 percent downward bias,
and two paragraphs on how well `sv` is identified, including the fit in
60 that loses it with nothing reported.
