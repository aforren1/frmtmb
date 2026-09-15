# Item 2.6, coupling: what the condition contrast recovers

Date: 2026-09-10. Worktree `frmtmb-wt-coh`, branch `wt-coh`, based on
the round 3 release (`0f94605`). frmtmb 0.55.2, frmtmb.coupling 0.3.0,
R 4.6.1, Windows 11. Reference library
`C:/Users/adf44/source/r/rellib-r3` read only; everything built here
went into the private library `C:/Users/adf44/source/r/coh-lib`, which
is first on `.libPaths()`.

Item 2.6 of `dev/extension-gaps-plan.md`:

> coupling: `coh ~ cond + s(freq, by = cond) + (1 | id) + (1 | id:cond)`,
> promoted from the scratchpad benchmark to a recovery assertion on the
> condition contrast, with its standard error captured this time.

Phase 0 had already fit that model and the four rungs below it at one
seed and found the contrast covered, so the survey's 0.77 was noise.
What one seed cannot say is whether the correct model covers at the
NOMINAL RATE and what the cheaper rungs cost. This lane measured both
over replicates.

## What punch round 1 changed, for a re-check

The review is `dev/reviews/20260910-coh.md`. The headline table, the
rule, the guard and the calibration all reproduced and none of them
moved. What moved:

1. **The bias attribution, which was wrong twice.** It called the two
   offset rungs "the two rungs that also drop the smooth" when only one
   of them drops the smooth, and it named the frequency term as what
   moves the estimate when paired rung by rung the smooth is the
   SMALLEST of the three. Worse, it read a different estimand as an
   error. Rewritten in this document, in `?cross_wishart`, in the
   vignette, in `NEWS.md` and in the plan row; the section "The offset
   in the low rungs is an estimand" carries the replacement, and every
   figure in it was re-derived here with `dev/coh-marginal.R` rather
   than copied from the review.
2. **The false-alarm rate quoted in the plan row**, which belonged to
   the rejected unpaired statistic. The row now gives the shipped
   paired quotient's own rate and says which statistic each number
   belongs to. Round 2 then found that rate quoted too precisely, and
   item 4 below is what replaced it.
3. Six nits: the guard's discrimination floor is now a comment beside
   it, the pairing property it rests on is now asserted and was seen to
   fail against the one spelling that breaks it, the `id` rung's
   residual is named, the intercept half of the Rd's claim is now
   measured, the examples NOTE has its second arm, the selection check
   on the six dropped replicates is recorded, and the CI workflow's
   expectation counts are refreshed.

Punch round 2, which was the last, moved one number and four sentences:

4. **"About 1 in 2,200" was one shape's answer given to two
   significant figures.** It is demoted everywhere it appears, behind
   the observation it came from: 0 of 40, bounding the four-seed rate
   below 0.26. See "What the guard's own numbers are" below, which
   re-derives the review's attack in `dev/coh-falsealarm2.R`.
5. Four prose corrections: the null arm is a second `V` the formula
   SURVIVES rather than a confirmation, and `0.346` is an
   approximation good to about 0.001 here; the residual arithmetic is
   an identity and now says so; and two Rd sentences, one that flattened
   conservative against anti-conservative and one that confessed to a
   draft users never saw.

## The claim, and what happened to it

> The correct model's interval on the condition contrast covers at the
> nominal rate, and the misspecified rungs' intervals do not, because
> they are too narrow rather than because they are in the wrong place.

**Held, with one correction and one addition.**

- The correct model covers 141 of 148, **0.953 (0.906, 0.977)**, which
  contains 0.95. Its standardized error `(estimate - truth) / se` has a
  spread of **1.062** over replicates, where 1.0 is an honest interval.
- The survey's model, `... + (1 | id)` with no `(1 | id:cond)`, covers
  **78 of 148, 0.527 (0.447, 0.606)**, with a standardized error of
  **2.774**: its interval is two and three quarter times too narrow.
  Paired against the correct model on the same 148 data sets it sits
  **+0.00360** away, on a truth of 0.5, and coverage predicted from the
  location and the width alone is 0.519 against 0.527 observed. It is
  the width. This is the rung that keeps a grouping and keeps the wrong
  one, and it is the only rung where a wider interval is the remedy.
- **The correction, and it is not a smaller version of the same
  failure. The two rungs with NO random effects estimate a DIFFERENT
  QUANTITY, and they estimate it correctly.** A logit contrast is not
  collapsible, so a model that leaves variance out of the linear
  predictor is consistent for the population-averaged contrast rather
  than the conditional one. The section "The offset in the low rungs is
  an estimand" below predicts both offsets to within 2 percent from the
  generator's constants with nothing fitted. No interval width repairs
  that, and none should be offered as though it would.
  `dev/scale-findings.md` had this right before this lane ran and this
  document's first version moved away from it; the details are in that
  section.
- The addition, which was not asked for and is the reason the fifth
  rung was added: **dropping `(1 | id)` instead goes the other way.**
  `... + (1 | id:cond)` with no `(1 | id)` is 1.90 times TOO WIDE, with
  a standardized error of 0.559 and coverage 148 of 148. So "fewer
  random effects, narrower interval" is not the rule. The rule is which
  one: `cond` varies INSIDE a subject, so `(1 | id)` shifts both of a
  subject's conditions together and cancels out of the contrast, while
  `(1 | id:cond)` is the term that carries the contrast's
  between-subject spread. Omit the crossed term and the interval
  collapses; omit the subject term and its variance is forced into the
  crossed term, which inflates it.

## What the replicate count can resolve, stated first

148 complete replicates put a Wilson 95 percent half-width of about
3.5 points on a coverage near 0.95 and about 8 points on one near 0.5.
That separates 0.95 from 0.53 many times over; it would NOT separate
0.95 from 0.92, and no claim below rests on such a difference. The null
arm has 36 replicates, whose half-width near 0.95 is about 11 points,
so it is read as a control on the WIDTH, which is a continuous quantity
measured far more precisely, and not as a coverage measurement.

## The design and the seeds

The Phase 0 design, unchanged: 40 subjects x 2 conditions x 60
frequencies, 8 segments per row, 4,800 rows. The truth puts every
structure in the coherence, on the logit scale: intercept -0.6,
condition contrast 0.5, a Gaussian bump of height 0.8 centered at 0.35,
`sd(id)` 0.35, `sd(id:cond)` 0.20; the two powers and the phase are
constants at 0.3, 0.1 and 0.4.

The simulator is `dev/coh-sim.R`. It draws the SAME data as the Phase 0
tier at the same seed, and `dev/coh-probe.R` checks that with
`identical()` rather than by eye: TRUE at seed 20260908.

| arm | truth for `sd(id:cond)` | seeds | replicates |
| --- | --- | --- | --- |
| main | 0.20 | within 20260911 to 20261105 | 148 |
| null | 0 | 20270911 to 20270946, all | 36 |

The main arm's seeds are not the whole range: six processes each took a
contiguous block of the 200 that were queued and the run was stopped
before the blocks met, so the 148 are the completed prefixes of six
blocks. Each replicate's seed is on its own row in the TSV, and the
seed is `20260910 + the replicate number`.

Scripts: `dev/coh-recovery.R` for one seed block, `dev/coh-run.ps1` to
run several blocks as several processes, `dev/coh-summarize.R`,
`dev/coh-extra.R` and `dev/coh-counts.R` for the tables. The raw rows
are `dev/coh-recovery-main.tsv` and `dev/coh-recovery-null.tsv`, one
line per rung per replicate, each carrying its own seed.

## Phase 0 reproduces, and its 53 percent is its own seed's number

`dev/coh-probe.R`, seed 20260908, one process:

| model | estimate | se | interval | components |
| --- | --- | --- | --- | --- |
| `... + (1 \| id) + (1 \| id:cond)` | 0.4877 | 0.0372 | (0.415, 0.561) | 0.313, 0.146 |
| `... + (1 \| id)` | 0.4858 | 0.0175 | (0.451, 0.520) | 0.328 |

Every figure matches `dev/scale-findings.md` and the plan's item 2.6
row. The width ratio at that seed is 0.4716, that is 52.8 percent
narrower, which is where the recorded "53 percent" comes from. It is
correct at its seed and it UNDERSTATES the typical case: over 148
replicates the mean ratio is 0.374, 62.6 percent narrower, and the
paired ratio of the correct width to the survey's runs from 1.900 to
3.646 with a median of 2.680.

## Coverage and width, 148 replicates

Truth 0.5. `width` is the mean interval width relative to the correct
model's. `sd(z)` is the spread of `(estimate - truth) / se` across
replicates: 1.0 is a calibrated interval, above 1 is too narrow, below
1 too wide. `mean(z)` is the same statistic's center.

| coherence formula | covers | 95 percent | width | sd(z) | mean(z) |
| --- | --- | --- | --- | --- | --- |
| `cond` | 73/148 | 0.493 (0.414, 0.573) | 0.371 | 2.746 | -0.960 |
| `cond + s(freq, by = cond)` | 72/148 | 0.486 (0.407, 0.566) | 0.375 | 2.753 | -0.591 |
| the same `+ (1 \| id)` | 78/148 | 0.527 (0.447, 0.606) | 0.374 | 2.774 | 0.351 |
| the same `+ (1 \| id:cond)` | 148/148 | 1.000 (0.975, 1.000) | 1.897 | 0.559 | 0.032 |
| both, the correct model | 141/148 | 0.953 (0.906, 0.977) | 1.000 | 1.062 | 0.050 |

Coverage predicted from the location and the width alone, under
normality, against what was observed: 0.499 against 0.493, 0.515
against 0.486, 0.519 against 0.527, 1.000 against 1.000, and 0.941
against 0.953. Nothing but where the interval sits and how wide it is
is at work in any rung.

Paired, since every rung of a replicate saw the same data: the
survey's model missed while the correct model covered on **63 of 148**
replicates, and the reverse happened **0** times.

The `mean(z)` column is the one place this table invites a wrong
reading, and the first version of this document took it: the -0.960 and
-0.591 in the two lowest rungs are not a bias to be fixed. The next
section says what they are.

The correct model also recovers what it is asked to: `sd(id)` 0.3412
plus or minus 0.0488 against 0.35, and `sd(id:cond)` 0.1965 plus or
minus 0.0269 against 0.20.

## The offset in the low rungs is an estimand, not an error

Every rung of a replicate saw the same data at the same `n`, so the
rungs difference seed by seed. That is a far sharper instrument than
five means against 0.5, and it also rules out any finite-sample reading
at once: an `n`-driven explanation cannot move two rungs and leave
three alone when all five have the same `n`. From
`dev/coh-marginal.R`, 148 paired differences:

| paired difference | value | se | t |
| --- | --- | --- | --- |
| `cond` minus the correct model | -0.01969 | 0.00050 | -39.70 |
| `smooth` minus the correct model | -0.01328 | 0.00050 | -26.74 |
| `id` minus the correct model | +0.00360 | 0.00078 | +4.63 |
| `idcond` minus the correct model | -0.00000 | 0.00009 | -0.05 |

Rung to rung, which is the question "what does dropping this term
cost": dropping `(1 | id)` moves the contrast by **-0.01687**
(t = -13.92), dropping the smooth by **-0.00641** (t = -49.62), and
dropping `(1 | id:cond)` by **+0.00360** (t = +4.63). The frequency
term is the smallest of the three, not the driver.

**It is the marginal contrast, and the prediction is quantitative.** A
contrast on a logit scale is not collapsible. A model that leaves
variance in the linear predictor unmodeled is consistent for the
population-averaged contrast, which is attenuated by about
`1 / sqrt(1 + 0.346 V)` with `V` the omitted variance. Taking `V` from
the generator's own constants, with nothing fitted: the bump's variance
over the design's frequency grid is 0.07847, so `V` is 0.24097 for
`cond` (bump plus 0.35^2 plus 0.20^2), 0.16250 for `smooth` and 0.04000
for `id`.

| rung | omitted `V` | predicted offset | observed, paired |
| --- | --- | --- | --- |
| `cond` | 0.24097 | -0.01962 | **-0.01969** |
| `smooth` | 0.16250 | -0.01349 | **-0.01328** |
| `idcond` | 0 | 0 | **-0.00000** |

Within 0.4 and 1.6 percent, from constants rather than from a fit.
Those rungs are not estimating 0.5 badly. They are estimating a
different quantity well, and no interval width repairs that. **The two
failures need different remedies, so the documentation has to say
which.**

**This was already in the repository and this lane's first version
moved away from it.** `dev/scale-findings.md`, above its own ladder
table, calls the low rungs' estimate "the attenuated marginal one",
measures the attenuation at "0.016 across the whole ladder, about a
fifth of a standard error", and says it is real, monotone and small.
That record is right and this document adds to it rather than replacing
it: what is new here is 148 paired replicates instead of one seed, and
the closed-form prediction beside them.

**One residual the model does not explain, on the rung the narrative
rests on.** For the `id` rung the attenuation formula predicts -0.00342
and the paired measurement is **+0.00360**, a gap of **+0.00702** at
t = 4.63. It is small in absolute terms, 0.0036 on a truth of 0.5, and
it does not touch the width conclusion, which rests on `sd(z)` = 2.774
and on 63 of 148 paired misses. But it is named here rather than folded
into a bias figure, because it is the one rung of the five the marginal
model gets wrong and nobody has explained it.

Saying it "sits on one rung" is meaningful only because the residuals
are ADDITIVE BY CONSTRUCTION, and that is arithmetic rather than a
second check. Both the paired means and the predictions subtract across
rungs, so the residual of a rung-to-rung difference is exactly the
difference of the two rungs' residuals:
`r(smooth - id) + r(id) - r(smooth)` is 1.73e-18, which is
`all.equal()` to zero. The prediction for dropping `(1 | id)` is
-0.01007 and the measurement is -0.01687; the -0.00681 between them is
not independent evidence, it is `r(smooth)` minus `r(id)`, and the
+0.00021 left over IS the `smooth` rung's own residual to the last
digit. What the identity buys is that the 0.00702 cannot be an artifact
of how the ladder was differenced.

## The control arm, where the term is really absent

Same generator, `sd(id:cond)` set to 0, so `... + (1 | id)` is the
correct model there and the correct model of the main arm is merely
over-parameterized. 36 replicates.

| coherence formula | covers | width | sd(z) | mean(z) |
| --- | --- | --- | --- | --- |
| `cond` | 28/36 | 0.960 | 1.155 | -1.106 |
| `cond + s(freq, by = cond)` | 30/36 | 0.970 | 1.133 | -0.729 |
| the same `+ (1 \| id)` | 33/36 | 0.962 | 1.153 | -0.218 |
| the same `+ (1 \| id:cond)` | 36/36 | 4.198 | 0.263 | -0.043 |
| both | 34/36 | 1.000 | 1.101 | -0.211 |

The two lowest rungs still sit low here, and it is the same estimand
rather than a bias that survived: paired against `id`, which is the
correct model in this arm, `cond` is -0.01570 (se 0.00089) against a
predicted -0.01653 and `smooth` is -0.00917 (se 0.00095) against a
predicted -0.01027, with `V` smaller by the component the truth no
longer has.

**Call that a second `V` the formula SURVIVES, not a confirmation.**
The two arms are not the same grade of agreement: the null arm's
standard error is 1.9 times the main arm's and its residuals are 5.1
times as large, +0.00083 and +0.00110 at 0.9 and 1.2 standard errors,
against -0.00006 and +0.00021 at 0.1 and 0.4. The obvious excuse is
excluded, since an offset reference would show as `full` minus `id` in
this arm and that is -0.00005 with se 0.00010. The right reading is
about the formula rather than the arm: **`0.346` is a standard
approximation, not an identity**, and its demonstrated accuracy over
all four predictions on this generator is about **0.001 on a contrast
of 0.5**. The main arm's agreement in the fourth decimal is better than
the approximation deserves, and nobody should treat it as the standard
to reproduce.

**The width gap closes.** The paired ratio of the two-component model's
width to the survey's model's is 1.000 at the median and 1.264 at its
largest, against 2.680 at the median in the main arm. The estimated
`sd(id:cond)` is 0.0146 plus or minus 0.0175 against a truth of 0, and
it is below 0.01 on 20 of 36 replicates. That is the control the width
comparison needs: adding a term the truth does not have costs
essentially nothing, so the main arm's factor of 2.7 is the missing
component and not the arithmetic of adding a term.

One number is worth more than the rest of this section. In the main
arm the survey's model reports a mean width of **0.0698**. In the null
arm the CORRECT model reports **0.0723**. The two agree to 3.5 percent:
dropping `(1 | id:cond)` makes the fit report the interval belonging to
an experiment that has no subject-by-condition variance at all, on data
that does. It is the right answer to a different question.

## What ships

**`tests/testthat/test-coherence.R`,** a new test,
"a within-condition contrast needs the within-subject term". The
realistic design costs about 50 seconds a fit and this file cannot pay
that, so the test runs a small version: 16 subjects, 2 conditions, 6
replicate rows a cell, no frequency axis, 8 segments a row. It measures

    ratio = width(cond + (1|id) + (1|id:cond)) / width(cond + (1|id))

in a treatment arm with `sd(id:cond) = 0.5` and in a null arm with 0,
at the SAME four seeds, and asserts `min(treat / null) > 1`, a quotient
paired seed by seed. It also asserts `min(treat) > 1`, which is a sign
check and nothing more: adding the term must never narrow the interval.
Sixteen fits, about 25 seconds, and one assertion beside them that
costs no fit at all.

Everything in it is a ratio between two numbers the run itself
measures. There is no absolute width anywhere, which is the rule four
tests in this project have already broken.

Both powers are held at 1 by construction, so `pow2 ~ 1` is the correct
model for them; the section on defects below says why that is not a
detail.

The two arms differ in ONE thing. `rnorm(n, 0, 0)` returns without
drawing, so generating the null arm at `sd = 0` would have moved every
Wishart draw after it and the arms would have differed in their data as
well as in their truth; the generator draws at sd 1 and scales. **That
property is now asserted**, at no fit cost: `cp_cells()` must leave the
random stream in the same state at `sd_idcond` 0 and 0.5. The assertion
was seen to fail against the one spelling it guards, `rnorm(n, 0, sd)`,
by mutating that line of the shipped generator in `dev/coh-absent.R`:
TRUE as shipped, FALSE mutated.

The mutation probe was checked independently in the re-check
(`dev/cohrev2-mutate.R`) and mutates what it claims: exactly one
deparsed line differs, the shipped spelling is absent from the mutated
body and the new one present, and at `sd_idcond = 0` the mutated
generator produces data that is not `identical()` to the shipped one.
The FALSE it prints is earned rather than a silent `gsub` miss.

**What the guard's own numbers are, observation first.** From the 40
paired quotients in `dev/coh-paired.tsv`, re-derived here in
`dev/coh-falsealarm2.R`: **0 of 40 fired**. With no shape assumed that
bounds a seed below **0.072** by Clopper-Pearson at 95 percent, the
rule of three agreeing at 0.075, and so bounds the minimum of four the
test takes below **0.259, about 1 in 3.9**. That is the measurement.

A fit says far rarer, and the fit is doing the work. The log quotient
has mean 0.6210 and standard deviation 0.1685, so 1 sits 3.69 standard
deviations below the mean and a lognormal reads **1 in 2,202**. Forty
points cannot choose the shape: Shapiro-Wilk rejects neither the
lognormal (W 0.9651, p 0.249) nor a plain normal on the raw quotients
(W 0.9586, p 0.150), and the normal reads **1 in 90**. Nor is the
lognormal's own answer resolved, because it is governed by a standard
deviation estimated from 40 points: a nonparametric bootstrap of that
fit gives a median of 1 in 2,631 and a 95 percent range of **1 in 323
to 1 in 114,149**. The review's independent version of the same
exercise adds a gamma at 1 in 568 and a t(5) at 1 in 99, a factor of 25
across four shapes the data cannot separate.

**So the number to carry is 0 of 40 and a bound near 1 in 4.** The
first version of this document quoted "about 1 in 2,200" to two
significant figures with its construction named but its bound missing,
which is manufactured precision on a quantity whose own model spans
three orders of magnitude and whose model-free bound is 500 times
nearer. The plan row and the test's comment now lead with the
observation and demote the fit.

**Where the guard stops discriminating, which is a limit and not a
defect.** With the effect PRESENT but small, `sd(id:cond) = 0.15` at
seeds 2700:2703, the shipped assertion FAILS at
`min(treat / null) = 0.999999956`, because the correct model's
component collapses to zero on two of those four seeds. At 0.5 on those
same fresh seeds it passes at 1.577, so the firing is not a property of
the pinned seeds in either direction. The discrimination floor is
between 0.15 and 0.5; the test's comment says so, so that a failure
there is read for what it is.

**`tests/testthat/test-scale.R`** now records `coh_cond_se`, the
standard error itself and not only the interval it produced. That is
the half of item 2.6 that Phase 0 left out, and it is the column the
rungs differ in.

## How the shipped assertion was calibrated, and what was rejected

The obvious test, comparing the two arms' width ratios at one seed
each, does not work at a small design and the numbers say so.
`dev/coh-calib.R`, 20 seeds in each arm on each of eight designs:

| subjects | rows a cell | `sd(id:cond)` | treatment arm | null arm |
| --- | --- | --- | --- | --- |
| 16 | 6 | 0.35 | 1.094 to 2.129 | 1.000 to 1.383 |
| 24 | 6 | 0.35 | 1.251 to 1.911 | 1.000 to 1.197 |
| 16 | 10 | 0.35 | 1.347 to 2.479 | 1.000 to 1.379 |
| 24 | 10 | 0.35 | 1.356 to 2.559 | 1.000 to 1.267 |
| 16 | 6 | 0.50 | 1.420 to 2.603 | 1.000 to 1.204 |
| 24 | 6 | 0.50 | 1.384 to 2.227 | 1.000 to 1.208 |
| 16 | 10 | 0.50 | 1.491 to 3.075 | 1.000 to 1.370 |
| 24 | 10 | 0.50 | 2.004 to 2.633 | 1.000 to 1.215 |

The ratio rises with every one of the three knobs, and the null arm
sits just above 1 in all eight. The design the test uses is the fifth
row, the cheapest one whose arms are far apart.

The arms overlap one seed at a time, because a null-arm fit that
happens to estimate a nonzero component inflates its own interval. At
60 seeds an arm (`dev/coh-calib3.R`) the treatment arm at the last of
those designs spans 1.265 to 2.603 and the null arm 1.000 to 1.511, so
they overlap there too.

**The first candidate failed open, and the construction that caught it
is the one the standing rules ask for.** It compared the MEAN of four
treatment seeds against the LARGEST of four null seeds, and resampling
blocks of four out of those 120 calibration measurements put its
failure rate at 4 in 20,000 (`dev/coh-falsealarm.R`), which is a
measurement of the wrong thing: it says how often the statistic fails
when the effect is THERE. Run instead on data whose truth has no
subject-by-condition effect at all, two independent blocks of four, it
**passed**: 1.0889 against 1.0772 (`dev/coh-absent.R`). Two blocks of
four from one distribution land either way often enough that a
statistic built on them asserts nothing.

The pairing was there to be used and costs nothing. The generator draws
the deviations at sd 1 and scales, so at one seed the treatment and the
null data differ ONLY by that scaling. With the effect absent from both
arms the two calls therefore draw the SAME data, the widths are
`identical()`, and the quotient is exactly 1: the shipped assertion
fails by construction where the unpaired one flipped a coin. Measured
on the shipped generator over 40 seeds (`dev/coh-paired.R`), the
quotient runs **1.422 to 2.596** with a median of 1.893 and none at or
below 1, with the treatment ratios spanning 1.438 to 2.596 and the null
ratios 1.000 to 1.271. The four seeds the test pins give 1.568, 2.596,
2.057 and 1.709.

The sign check beside it, `min(treat) > 1`, does NOT discriminate:
`dev/coh-absent.R` shows it passing with the effect absent, because
adding a variance component essentially never narrows an interval. It
is in the test for the direction and it is labeled as such.

**Rejected: a model-free reference standard error.** A two-stage
estimator, pooling each subject-by-condition cell into one Hermitian
matrix, differencing the two conditions inside a subject on the logit
scale and reading the standard error off the between-subject spread of
the differences, is a reference that carries the subject-by-condition
variance by construction and needs no model. Over 30 seeds it separates
the arms cleanly (`dev/coh-calib2.R`): the correct model's standard
error is 0.83 to 0.93 of it in BOTH arms, and the survey's model's is
0.48 to 0.58 of it in the treatment arm against 0.79 to 0.87 in the
null arm. Per seed it does not, and cannot at this size: the reference
is a standard deviation of 16 numbers, and in the null arm alone it
ranges from 0.52 to 1.37 of the model's own answer. A test built on it
would have been a pinned-seed test wearing a reference's clothes. The
measurement is kept because it is the one arm of this study that does
not assume the family's own arithmetic.

**Rejected: 200 replicates.** The main arm was launched at 200 and
stopped at 148 when the round's coordinator capped concurrent R
processes at 3 to return memory to the machine. Nothing in this
document needs the difference: at 148 the Wilson half-width near 0.95
is 3.5 points and near 0.5 is 8 points, and the gap being measured is
43 points.

**And the stop did not select.** Each of the 6 dropped replicates is a
strict PREFIX of the fit order `cond, smooth, id, idcond, full`, which
is what truncation looks like and not what attrition on hard data sets
looks like: seeds 20260935, 20260969, 20261003, 20261038, 20261073 and
20261106, with 3, 2, 3, 4, 2 and 3 rungs. Every row of both arms
carries `ok = TRUE`; there are zero failed fits in 757 main rows and
180 null rows, and 148 x 5 + 17 = 757 exactly, so there are no
duplicates either. **The correct model appears in none of the six**, so
141 of 148 cannot move under any handling. Counting every available row
instead of complete replicates only moves `cond` to 0.4935 of 154,
`smooth` to 0.4870 of 154, `id` to 0.5197 of 152 and `idcond` to 1.0000
of 149: the largest move is 0.7 points, on `id`, and it moves DOWN,
which is the opposite of the direction a selection worry points.

## Defects found, and one nuisance that is not one

**No defect was found in the package.** Every rung fits, reports a
positive definite Hessian, and recovers what its formula can carry.

**One was found in the first draft of the shipped test, by its own
help page.** The test file's existing `cp_units()` builds the second
channel as `cm z1 + b z2` with `cm = b exp(eta)`, so that channel's
POWER is `b^2 (1 + exp(2 eta))`: at `lb = 0` the power is tied to the
coherence and spreads 2.6-fold across the rows of this design, which is
exactly the trap `?cross_wishart` documents and which `pow2 ~ 1` is the
wrong model for. It was caught by reading the help page rather than by
a failing run, which is the weak form of the evidence standard, so what
it is worth is stated plainly: with the power free to move with the
coherence the four pinned seeds still separated the arms, 1.315, 1.714,
1.738 and 1.295 against null ratios of 1.097, 1.000, 1.067 and 1.000,
so nothing measured here was wrong in its direction. What was wrong was
that the calibration runs all use a generator with CONSTANT powers, so
the rate being quoted did not belong to the design being shipped.
Choosing `lb` to cancel the factor holds both powers at 1, leaves the
coherence alone because coherence is invariant to scaling a channel,
and raises the same four seeds to 1.663, 2.596, 2.216 and 1.709. The
lesson is the second half: a calibrated rate belongs to the generator
it was calibrated on.

The nuisance: 3 of 148 replicates of the correct model return
convergence code 1, `nlminb`'s "false convergence", with a maximum
absolute gradient of 0.0023, 0.00175 and 0.00313 and a positive
definite Hessian. All three cover. Dropping them moves coverage from
141 of 148 to 138 of 145, that is 0.9527 to 0.9517, so nothing here
turns on how they are treated; both numbers are reported rather than
one. Seeds 20260948, 20260962 and 20261092.

## Process notes

Two things a reader of the raw files should know.

**The replicate count and the process count disagree, and the reason is
Windows.** `dev/coh-run.ps1` starts each seed block with
`bin\Rscript.exe`, which on Windows is a launcher that re-executes
`bin\x64\Rscript.exe` and waits. Six workers therefore show as twelve
processes: six doing the work at 400 MB to 1.2 GB each, and six
launcher stubs holding 0 MB and about 0.1 seconds of CPU. The stubs are
not a cluster, not a per-replicate spawn, and not workers that failed
to connect. Nothing counts processes: a replicate is counted only when
its row reaches the TSV, and `dev/coh-summarize.R` drops any replicate
whose five rungs are not all present. Six partial replicates were
dropped that way when the run was stopped.

**Timings in the TSV are not timings.** The `secs` column was recorded
with up to nine of this machine's cores held by other lanes, and it
ranges over a factor of three for the same fit. The cost figures to use
are `dev/scale-findings.md`'s, which were taken the way that document
describes.

## Files, and what each test file reported

Changed in the package:

| file | what |
| --- | --- |
| `R/cross-wishart.R` | a recovery section on `?cross_wishart` |
| `man/cross_wishart.Rd` | roxygenised from it |
| `tests/testthat/test-coherence.R` | the new test and its generator |
| `tests/testthat/test-scale.R` | records `coh_cond_se` |
| `vignettes/coherence.Rmd` | a subsection beside the other width failure |
| `NEWS.md` | the entry, under a heading with no number in it |
| `dev/scale-row.md` | points at this study |

Changed at the repository root: `dev/extension-gaps-plan.md`, item
2.6's row only, and
`.github/workflows/check-frmtmb-coupling.yaml`, whose per-file
expectation counts this lane moved. The Phase 2 total and the sequence
table in the plan were left alone deliberately, because every lane in
this round would edit the same two lines.

**The workflow counts, and which drift is this lane's.** The comment
there is a record and not an assertion: no step compares it, so the job
could not have failed on it, which is how the previous set went stale.
Punch round 1 called these counts a breakage that "would fail" the job;
the re-check withdrew that and agreed with the reading here, having
confirmed independently that no step in that workflow parses them, that
no other workflow does, and that no test file references the workflow.
The general form is going into the lane rules: **in this repository a
stale CI count is a record defect and not a breakage unless a step
actually compares it.**
Remeasured in one process with `dev/coh-suite-onefile.R`: 431
expectations and 5 skips with `NOT_CRAN`, 408 and 13 without.
Coherence 58 to 61 is this lane. Cross-spectrum 61 to 96,
cross-wishart 203 to 211, and the absence of `test-bracket-access.R`
and `test-scale.R` from the old list are drift this lane did not cause
and refreshed in the same pass.

Nothing under `R/` changed except roxygen comments, and no core API
that this package did not already use was called, so the
`frmtmb (>= 0.55.1)` floor in `DESCRIPTION` still holds. The package
version is not set here: the NEWS heading reads
"frmtmb.coupling (development version)" and the round's session sets
the number.

The suite, one file per R process, `NOT_CRAN=true`, against the private
library:

| file | pass | fail | error | skip | warn | seconds |
| --- | --- | --- | --- | --- | --- | --- |
| `test-bracket-access.R` | 1 | 0 | 0 | 0 | 0 | 3.7 |
| `test-message-uniqueness.R` | 4 | 0 | 0 | 0 | 0 | 3.5 |
| `test-cross-wishart.R` | 211 | 0 | 0 | 0 | 0 | 9.5 |
| `test-scale.R` | 0 | 0 | 0 | 5 | 0 | 1.2 |
| `test-cross-spectrum.R` | 96 | 0 | 0 | 0 | 0 | 23.0 |
| `test-surface.R` | 58 | 0 | 0 | 0 | 1 | 32.4 |
| `test-coherence.R` | 61 | 0 | 0 | 0 | 3 | 63.4 |

`R CMD check --as-cran` on the final tree, once, with pandoc and
TinyTeX on PATH and `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`
(`dev/coh-check.ps1`, log in `dev/coh-check-00check.log`):
**Status: 1 NOTE**, and the NOTE is the runtime of an example this lane
did not touch. `?cross_wishart`'s example measured 4.69 + 0.59 seconds
of CPU under check against a 5 second threshold; the `man` diff is 56
insertions and 0 deletions with none of them inside `\examples`, and in
its own process the same code runs at 3.04 and 3.28 seconds of CPU over
two rounds (`dev/coh-example-time.R`).

**That argument was one-armed and the review supplied the other arm; it
is load, and the statement to carry forward is not "pre-existing" but
LOAD-DEPENDENT.** Three measurements, the first two the reviewer's
(`dev/cohrev-identity2.R`, `dev/cohrev-extime.R`) and the third checked
here directly:

- The work is identical by COUNT and not by clock: all 22 objects in
  this lane's installed namespace deparse identically to the base
  build's, and all 14 example runs return objective 305.7981502108 at
  2 `nlminb` iterations. A change of 54 roxygen lines and 0 lines of
  code cannot make an example slower.
- Seven interleaved fresh-process pairs, lane against base, each
  carrying a fixed 400 x 400 arithmetic block as a CONTROL that must
  report the same time in every process. The control ran from 0.20 to
  0.95 seconds of CPU on identical arithmetic, a factor of **4.7**, and
  example CPU correlates with it at **0.809**. The BASE build itself
  reached **6.75 s**, over the threshold, and this lane's build has the
  faster minimum, 1.94 against 2.33.
- The base commit's own release check, which I read directly:
  `C:/Users/adf44/source/r/checkout-r3/frmtmb.coupling.Rcheck/frmtmb.coupling-Ex.timings`
  records `cross_wishart` at 0.92 + 0.13 = **1.05 s** of CPU, and the
  `frmtmb.coupling-Ex.R` beside it carries this same example verbatim.
  Same code and same threshold, a factor of five from this lane's
  5.28 s.

So an examples-timing NOTE on this machine reports how loaded the box
was. `Status: OK` at round 3 was as much an artifact of a quiet machine
as `1 NOTE` is of a busy one; the harness half of that belongs to the
round's session.

The tarball the check ran on differs from the final tree by one comment
line inside `test-coherence.R`, corrected after the build, and by the
two punch rounds' changes, which came later. **The check was NOT run
again**, and the reasoning a consolidating session needs to accept that
substitution is this, the last four points supplied by the re-check
which ran them independently:

- What the punch rounds changed is documentation, one plan clause and
  one assertion that costs no fit. No R code changed at all, so
  check's code checks and its code-against-documentation checks are
  unreachable by construction.
- `tools::checkRd()` over all five man pages returns **0 messages**,
  and `tools::parse_Rd()` parses the rewritten page
  (`dev/coh-recheck.R`).
- The new Rd content adds **no `\link`**, so there is no
  cross-reference exposure.
- It sits in `\section` blocks, while check's Rd line-width test covers
  only `\usage` and `\examples`, neither of which changed.
- The suite was run ONE FILE PER PROCESS, which is a stronger run than
  the single process inside check.
- The vignette knits in 11 seconds, which is the check step that would
  otherwise rebuild it.

And this round's own first finding argues against re-running it: on
this box a check costs fifteen to twenty-five minutes whose only
variable outcome is an examples-timing NOTE that measures how busy the
machine was, and the consolidating session runs the authoritative check
anyway. The counts below were remeasured on the final tree.

The counts sum `error` as well as `failed`, because a runner that adds
only `failed` prints a clean line for a file that aborted halfway. The
four warnings are pre-existing and belong to two tests this lane did
not touch. `test-scale.R` skips its five rows without
`FRMTMB_SCALE_TESTS`; it was also run gated at
`FRMTMB_SCALE_SMALL=true` on the `coupling-coh-full` row, where it
recorded `coh_cond_se=0.117012` beside an interval of
(-0.0177, 0.4410) on an estimate of 0.2117, which is that estimate plus
and minus 1.96 standard errors and reproduces the figure the test
file's own comment carries for the small design.
