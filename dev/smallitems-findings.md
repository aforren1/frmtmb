# Phase 1 small items: 1.4, 1.5 and 1.6

Lane `smallitems`, 2026-09-08, on branch `wt-smallitems` from 780dec1.
Packages at frmtmb.eam 0.5.1, frmtmb.learn 0.2.1, frmtmb.coupling
0.2.0, against frmtmb 0.55.0 built from this worktree into a private
library.

Two of the three items are delivered in full. The third, item 1.5, is
delivered in half, because the SPELLING the plan asks for is refused by
frmtmb's own parser. The measurement is in section 2, and the half that
IS buildable turned out to close a silent wrong answer that the item's
own check pointed at.

**This file was rewritten after review.** The first version got two
things wrong and this one says so where they were: the counterfactual
guard was declared on five families and needed all eight, and the units
guard's false-alarm cost was measured only at trial counts where it
cannot fail. Both are corrected below with the measurements, and
section 7 lists what changed.

## 1. Item 1.4, a units guard in eam

### What changed

`ddm_check_units()` and the constant `ddm_seconds_ceiling` in
`extensions/frmtmb.eam/R/ddm-shared.R`, called from the `valid_y` of all
five families: `wiener()`, `gddm()`, `lba()`, `rdm()` and
`wiener_gng()`. It warns when the fastest response is above 20, names
milliseconds, and says to divide by 1000. The message carries the class
`frmtmb_eam_units_warning`.

`wiener_gng()` passes the GO responses rather than the whole response
column, for the reason `gng_finalize()` already reads them: a no-go
row's entry is a placeholder the likelihood ignores, and a user who put
a small number there instead of the deadline would hide the data set's
scale behind it.

### Why the check is worth having at all

Every default in the package is an absolute time. The `ndt` link is a
logit scaled onto `(0, min(y))`. `gddm_control(dt = 0.01)` is a time
step, and `t_max = NULL` takes the window from `max(y)`, so a
millisecond record turns a 300-step solve into a 300,000-step one. None
of the five likelihoods refuses milliseconds; the fit converges and
reports a boundary separation three orders of magnitude out. Under the
plan's rule 3 that is the case that goes first.

### Where 20 comes from

Measured on simulated wiener data at the published parameter range
(`dev` script `small-units-measure.R` in the lane scratchpad): drift
0.5, 1.5, 3.0; boundary 0.8, 1.5, 2.5; non-decision time 0.15, 0.3,
0.6; at 200, 400 and 12000 trials. 81 cells, one data set each.

| quantity | measured |
| --- | --- |
| fastest response over the 81 cells | 0.157 to 0.801 s |
| false alarms at a ceiling of 20 s | 0 of 81 |
| caught when the same data is read as ms | 81 of 81 |
| smallest millisecond fastest response | 157 ms |

The statistic is the MINIMUM, not the mean: one slow trial cannot reach
the ceiling, because the whole data set has to be slower than it.

**That sweep does not measure the false-alarm cost, and the first
version of this file claimed it did.** Holding the trial count at 200
and up pins the fastest response just above the non-decision time, so
the guard reduces to "is the non-decision time above 20 seconds", and
the sweep's answer ("a false alarm needs 20 seconds of encoding and
motor time") is a property of the sweep rather than of the guard. The
exposed design is the opposite one: a SHORT session of SLOW decisions,
where the minimum of a few draws sits far above the floor.

Measured with `ddm_simulate()` at a drift of 0.18 and a non-decision
time of 1.5 seconds, 200 replicates a cell, boundary separations of 14
to 40 (`small-r34.R`). The FASTEST response is in brackets beside each
rate, because the fastest is the quantity the guard reads:

| median | 20 trials | 60 trials | 200 trials |
| --- | --- | --- | --- |
| 28 s | 0.000 (9.2 s) | 0.000 (7.3 s) | 0.000 (6.0 s) |
| 39 s | 0.045 (13.4 s) | 0.000 (10.3 s) | 0.000 (8.5 s) |
| 45 s | 0.155 (15.7 s) | 0.000 (12.2 s) | 0.000 (10.1 s) |
| 61 s | 0.690 (22.7 s) | 0.285 (18.2 s) | 0.000 (15.1 s) |
| 100 s | 1.000 (41.7 s) | 0.995 (34.8 s) | 0.995 (29.1 s) |

**The median does not determine the rate, and a table indexed only by
the median implies it does.** Holding it near 50 seconds and reaching
it three ways at 20 trials:

| drift | boundary | median | fastest | rate |
| --- | --- | --- | --- | --- |
| 0.10 | 18.7 | 57.0 s | 15.7 s | 0.185 |
| 0.18 | 21.5 | 49.4 s | 17.5 s | 0.260 |
| 0.35 | 37.6 | 52.4 s | 26.9 s | **0.935** |

A task slow because the boundary is far and the evidence strong has a
tight response time distribution, so the fastest of twenty trials sits
close to the median and the guard fires; a task slow because the
evidence is weak has a long right tail and a fast minimum. At one
median the rate spans 0.185 to 0.935, so the cells above are the LOW
end of a range for their median rather than the range.

The honest one-line cost is therefore the one the guard implements: it
fires when the FASTEST response passes 20 seconds, and how far the
fastest sits below the median depends on the spread as much as on the
center. A 20-trial session with a 40 second median is a deliberation,
insight or matrix-reasoning design, and its non-decision time is 1 to 2
seconds rather than 20.

The gap between the safe grid (boundary 2 to 5) and the exposed one
(14 to 40) was not swept by this lane, which is the shape a
cherry-picked boundary would also have. The review filled it and found
no cliff: the rise is continuous, the first non-zero cell is at a 34.7
second median, and the safe grid tops out at a 10.47 second median with
0 of 192.

The guard is still right, because nothing a standard two-choice task
produces reaches the ceiling. Over 192 cells at boundary separations of
2 to 5, non-decision times of 0.5 to 5 seconds, drifts of 0.3 to 2 and
12 to 80 trials, the fastest response ran from 0.599 to 6.605 seconds
and there were 0 of 192 false alarms.

The miss curve is flat around 20:

| ceiling | clean millisecond cells missed |
| --- | --- |
| 10 | 0 of 81 |
| 20 | 0 of 81 |
| 50 | 0 of 81 |
| 100 | 0 of 81 |
| 200 | 14 of 81 |
| 400 | 49 of 81 |

So 20 sits a factor of five inside the band where nothing is missed and
dominates 10, which would catch nothing more and would false-alarm on a
seconds design whose fastest trial fell between 10 and 20 s. Above 100
the miss rate is real, because an anticipatory response at 100 to 200 ms
is ordinary in millisecond data.

### Is 20 right for every family in the package

Yes, and for one reason: all five are two-choice or n-choice response
time families whose literature is written in seconds, and the threshold
is set by the scale of a response time rather than by anything a family
does with it. The families differ in how BADLY a units error hurts, not
in where the error is. `gddm()` is the worst, because its `dt` is an
absolute step, so the same mistake that costs `wiener()` a wrong
boundary costs `gddm()` a 300,000-step tape. That argues for the guard
being more valuable there, not for a different number.

### What the guard costs a user who legitimately has a slow task

One warning per FIT, and nothing else. Counted over eleven entry points
on a millisecond fit of 300 rows (`small-f2-units.R`): `frm()` raises it
once and `update()` once, because `update()` reassembles the frame and
is a second fit; `predict(type = "link")`, `predict(newdata = )`,
`fitted()`, `simulate()`, `summary()`, `residuals()`, `logLik()` and
`coef()` raise it zero times each. `valid_y` runs once, at frame
assembly (`?frmtmb_family`, "Slot call order"), and no post-fit method
reassembles the frame.

It is a warning rather than a refusal because rule 3 ranks a refusal
above a missing feature but a design with no trial under 20 seconds is a
correct model, and a check that fires on a correct model is a real cost.
The class lets such a design silence this one condition and keep the
convergence warnings.

### What it misses, measured

A millisecond record with fast-guess contamination in it. 30 replicates
per cell, 5 percent of trials replaced by a uniform draw over the
observed range:

| rows | caught | fastest ms response |
| --- | --- | --- |
| 400 | 24 of 30 | 3.24 to 352 ms |
| 12000 | 2 of 30 | 0.095 to 29.5 ms |

At 12000 rows the contaminant reliably reaches below 20 ms and hides the
whole data set.

**The first version of this file argued that the miss is acceptable
because the record is "already broken read as seconds". That is not
true and the argument is withdrawn.** Measured on one such record,
truth `bs = 1.4` (`small-f2-units.R`):

| reading | fastest | warnings | `diagnose()` | boundary |
| --- | --- | --- | --- | --- |
| seconds | 0.001397 | none at all | clean | 2.06, 47 percent out |
| milliseconds | 1.397 | none at all | clean | 65.2, 4558 percent out |

The seconds reading reports "Optimizer convergence code: 0", a maximum
gradient of 1.2e-05 and a positive definite Hessian. So it is a SECOND
silent wrong answer rather than a visible failure, which is the class
rule 3 ranks first, and at this seed the contaminant reached 1.4 ms so
the guard missed the millisecond reading too. That raises the rank of
item 3.5's contaminant mixture rather than excusing the miss.

A median-based arm would close the millisecond half: a millisecond
median sits near 600 against 0.6 in seconds, and contamination cannot
move a median. It was not shipped because the plan's statistic is the
smallest response, the median arm was not measured for false alarms,
and adding an unmeasured second heuristic to a guard is what this
project's own rules forbid. That is the whole justification now; the
"already broken" half of it was wrong.

### Test

`extensions/frmtmb.eam/tests/testthat/test-units.R`, new. Against the
unmodified base exported from HEAD and installed into its own library:
**5 pass, 24 fail, 1 error** (the error is `ddm_seconds_ceiling` not
existing). After: **32 pass**. It covers all five families on a
millisecond fixture and on a seconds fixture, the band's open lower
edge at exactly the ceiling, the warning class, and reachability
through a real `frm()` call.

## 2. Item 1.5, `reward(pay)` with one column

### The spelling the parser refuses, and what that does and does not block

`reward(pay1)` is refused by frmtmb, not by frmtmb.learn:

```
`reward()` takes 2 arguments, not 1
```

That comes from `R/parse.R:276` in core. An addition term's arity is
fixed when it is registered: `frmtmb_register_aterm(arity =)` validates
`length(arity) != 1L || arity < 1 || arity != round(arity)`, stores one
integer, and the parser compares the argument count to it for equality.
`vint()` and `vreal()` are variadic, but they are core terms handled in
their own branch of the parser rather than through the registry.

Registering `reward` at arity 1 instead is not a workaround: it would
refuse `reward(pay1, pay2)`, which every example, vignette, test and
compatibility row in `frmtmb.learn` writes.

**But what that blocks is the SPELLING, not the capability, and the
first version of this file overstated it.** A one-column route needs no
core change at all. A second term name registered at arity 1 keys its
value at `aterms[["reward1"]]`, which is exactly the key the two-column
spelling's first argument already uses, so a family that reads
`reward2` with a fallback takes both spellings. Verified as far as this
lane could take it without changing a shipped family
(`small-f8-seam.R`): after `frmtmb_register_aterm("reward1", arity = 1)`,
`choice | reward1(rec) ~ 1` gets past the parser and lands on the
`reward1` key, and what refuses it is then the FAMILY's own
`required_aterms` declaration, "the density needs `reward2`", which is
this package's to change and not core's.
`reward(pay1, pay2) + reward1(rec)` is refused as a duplicated term.
The review's `rvsm-p2-seam.R` completes the construction with a probe
family that reads the fallback and returns a bitwise identical
log-likelihood.

So the seam is a convenience for users rather than an unblocking
change, and it should be ranked that way:

> `frmtmb_register_aterm(arity =)` accepts a MINIMUM and a MAXIMUM
> argument count, with the NAMING taken from the maximum alone, so a
> term whose maximum arity is above 1 always delivers `name1`, `name2`
> and so on. A `reward()` registered as "1 to 2" then delivers
> `reward1` for both spellings and `reward2` only for the two-column
> one. Needed by frmtmb.learn's `reward(pay)` spelling.

The derivation adds a constraint to that filing, and it is the reason
the min-and-max shape is not merely tidier. Under a plain RANGE arity
with the naming this lane first proposed, the one-column spelling would
arrive as `aterms[["reward"]]` while the two-column one arrives as
`reward1` and `reward2`. `ln_counterfactual_of()` strips digits, so it
would group `c("reward", "reward1")` together;
`ln_check_counterfactual()` would find `cd[["reward1"]]` absent and
skip; and the ONE-COLUMN SPELLING WOULD BE UNGUARDED. That is exactly
the spelling that cannot supply a counterfactual and most needs the
refusal. With a maximum arity above 1 always delivering `name1`,
`name2` and so on, the one-column spelling arrives as `reward1` alone,
the group has one member, and the guard refuses the draw for a
different and equally correct reason. The filing must carry that or the
seam ships a hole.

That shape is the review's and it is better than the range this lane
first proposed. `arity` is read for its NAMING convention at
`R/parse.R:134`, `:167`, `:212` and `:281`, each as `arity > 1L` or
`== 1L`; a plain range makes every one of those a `max()`, while
min-and-max with naming off the max leaves `registered_aterm_of()`,
`aterm_base()` and `aterm_spelling()` untouched. It is `parse.R:276`
plus a `max()` at `:212`. It also removes the open question this lane
left, "how does a family declare which arity it got": the family
declares `required_aterms = "reward1"` unchanged and needs neither an
any-of group nor an `exclusive_aterms` pair.

Nothing was reached with `:::` and nothing under core's `R/` was
touched.

### The half that could be built, and why it was worth more than expected

The item's second sentence is "`simulate()` refuses on it BY NAME,
because it cannot know what the unchosen arm would have paid". That
refusal is buildable, and it closes a silent wrong answer that exists
TODAY under the spelling the docs recommend.

`bandit2arm_delta()`'s help page says, correctly, that a record of the
received outcome alone can be passed twice, because the likelihood
reads only the chosen arm's entry. It then says the second column is
what makes the simulator coherent. What it did not say is that the
simulator ran anyway.

### The fit is EXACTLY right, not right to printed precision

The first version of this file reported the two spellings as agreeing
"to twelve printed digits" and noted that `identical()` on the two
`fixef()` vectors was FALSE. That reads as a near-identity and it
understates the fact.

The unchosen entry is multiplied by a 0/1 chosen indicator before it
reaches anything, and zero times a finite number is exactly zero, so
the OBJECTIVE is bitwise the same function of the parameters whatever
that column holds. Measured by stopping the optimizer before its first
step, so that three fits sit at one parameter vector and the comparison
is of the objective rather than of three optimizer paths
(`small-f1-eight.R`, and the shipped test asserts it): the log-
likelihood is `identical()` across `reward(pay1, pay2)`,
`reward(rec, rec)` and a spelling whose unchosen entries are `N(100,
50)` noise.

The converged values then agree to testthat's relative default and not
always bitwise, and that is gradient accumulation over a tape whose
CONSTANTS changed rather than a different answer. So nothing about any
density needs changing; only the simulator was ever wrong.

### The draw, measured

Measured on the shipped two-armed design, arm 1 paying with probability
0.7 and arm 2 with 0.3. The statistic is the proportion of trials in
the second half that took arm 1:

| family and route | real schedule | duplicated |
| --- | --- | --- |
| bandit2arm_delta, simulate | 0.7444 (0.0024) | 0.5004 (0.0035) |
| prl_fictitious, simulate | 0.8616 (0.0027) | 0.5470 (0.0090) |
| prl_fictitious, task_simulate | 0.8367 (0.0037) | 0.5072 (0.0083) |
| rlddm, task_simulate | 0.1443 (0.0036) | 0.5040 (0.0071) |

40 draws each, 20 for `rlddm()`; the observed proportions in the data
being fitted are 0.7593 and 0.8633. Every duplicated column gives exact
chance, because with both options paying the same on every trial there
is nothing to learn, and each came back as a plausible-looking vector
of option codes with no error and no warning. `ts_par7()` drew 800 rows
the same way.

### The exemption that was wrong, and what it cost

**The first version of this guard was declared on five families and
exempted `prl_fictitious()`, and that exemption was a sixth silent
wrong answer.** The stated reason was that its counterfactual update
reads both columns. It does not. The update is

```r
oc <- c1 * d[["reward1"]] + c2 * d[["reward2"]]
t1 <- c1 * oc - c2 * oc
t2 <- c2 * oc - c1 * oc
```

with `c1` and `c2` the CHOSEN indicators, so `oc` is the chosen
option's payoff and the fictitious update is a sign flip of it. The
unchosen column never enters. Measured: the log-likelihood is bitwise
identical across the two-column spelling, the duplicated one and the
noise one, at -937.55332137794574, and `identical()` on the `fixef()`
vectors is TRUE as well.

The reason the exemption was written is worth recording, because it is
the general lesson. `?prl_fictitious` says, correctly and in its own
words, that the counterfactual outcome "is the NEGATIVE of the realized
one, which is the model's assumption of an anticorrelated task rather
than a reading of the second `reward()` column". The COMPATIBILITY
ROW in `R/zzz.R` says the opposite, and has since 0.1.0. Two documents
in one package said opposite things about one family, and this lane
read the wrong one and did not check it against the code. **An
exemption justified by what code looks like it does is how five became
six.** Nothing but the measurement should have settled it.

`rlddm()` and `ts_par7()` were missed for a different reason: both set
`sim = FALSE`, so `simulate()` refuses them, and this lane stopped
there. But `frm_task_simulate()` is exactly the route `?rlddm`'s own
refusal message tells users to take, and it is the route this lane
guarded second. Both drew from a degenerate schedule with no error.

### What was built

`ln_check_counterfactual()` in `extensions/frmtmb.learn/R/family.R`,
called from two places: the structure's `sim_ctx` slot, which is the
funnel `simulate()`, `frm_simulate()` and
`frmtmb.sample::posterior_predict()` all reach, and
`frm_task_simulate()`, which builds its own data and can be handed a
duplicated design just as easily.

**Derived rather than declared per family.** `ln_counterfactual_of()`
reads the schedule columns off the `aterms` a family already names,
grouped by term and kept where a group has two or more members, and
restricted to `ln_schedule_terms`, which is `reward` and `payoff` and
nothing else. `ln_family(counterfactual = FALSE)` is the explicit
opt-out and no family in this package takes it. That is the fix for the
failure above: a ninth family cannot be added without the question
being asked, because the answer comes from the terms it already names.

`stage2()` is deliberately not a schedule term. Its two columns are the
observed stage-two state and choice, which are data the likelihood
conditions on rather than what each option would have paid, and a
design in which state and choice agree on every trial must still draw.
The test file pins that.

What the derivation produces, read off the eight built family objects
(`small-f1-cover.R`):

| family | guarded columns |
| --- | --- |
| bandit2arm_delta, bandit2arm_dual | reward1, reward2 |
| prl_fictitious, rlddm | reward1, reward2 |
| bandit4arm2_kalman_filter, ts_par7 | payoff1 to payoff4 |
| igt_pvl_delta, igt_orl | payoff1 to payoff4 |

Eight of eight, and `ts_par7()` picks up its four payoff columns and
not its two `stage2()` ones.

**What the derivation is weaker at than a declaration, stated fairly.**
It reads a NAME. A schedule arriving under a name that is not `reward`
or `payoff` returns NULL and nothing says so, and the review built that
family: the same two-armed delta learner with its schedule through an
`outcome()` term draws 0.4300 against 0.7533 on the real schedule, 600
rows, no error, and the same through `vreal()`, which is the route
core's own `?frmtmb_register_aterm` names for a family with no term of
its own. A forgotten DECLARATION at least sits in the family's source
next to `aterms =`; a derivation that does not recognize a name is
invisible.

This is not hypothetical. Plan item 5.1, `bandit_delta(n_option = K)`,
has to register a term name of its own for K outside 2 and 4, because
an addition term's arity is fixed at registration, which is this lane's
own seam finding. If that name is not literally `reward` or `payoff`,
the guard stops applying.

The close is a test rather than a code change, because the only place
the miss is visible is this package's own registration table.
`test-counterfactual.R` now asserts that every `ln_aterms` entry of
arity 2 or more appears either in `ln_schedule_terms` or in a new
`ln_not_schedule_terms`, which names `stage2` and carries the reason it
is excluded. Registering a fourth term without classifying it turns the
suite red at the earliest moment anyone could catch it. The same file
pins the eight families' derived groups and pins that `outcome1`,
`vreal1` and a lone `reward1` are NOT groups, so the limit is written
down as behavior rather than as prose.

**One accident, recorded so nobody reads it as the rule.** The
predicate strips trailing digits, so a term registered as `reward3` at
arity 3 produces `reward31` through `reward33`, which strip to `reward`
and ARE picked up. That is luck. The rule is the name.

The signature is EVERY column of the term identical on every row. Some
columns identical is not the signature, because the data then still
says what at least one option not taken would have paid; the test file
pins that case as accepted on a four-deck family with three of four
columns duplicated.

### The opt-out argument, which had no validation

`ln_family(counterfactual =)` was read as "FALSE means NULL, NULL means
derive, anything else is taken as the groups", and
`ln_check_counterfactual()` then iterates the groups and skips anything
shorter than two. So every value that was neither `NULL` nor a list of
character vectors turned the guard OFF in silence. Measured on a probe
family naming `reward1`/`reward2`, each fed a duplicated schedule
(the review's `rvsm-r2-derive.R`, reproduced here as the new test's
before-run): `NA`, `0`, `"no"`, `list()` and **`TRUE`** all drew.

`counterfactual = TRUE` is the spelling a maintainer reaches for to
mean "yes, guard this one". In the one argument the whole derivation
hangs on, whose justification is that a ninth family cannot be added
"without the question being asked", the affirmative answer to that
question meant no.

`ln_read_counterfactual()` now takes `NULL` to derive,
`identical(FALSE)` to opt out, and the columns themselves as a
character vector or a list of them, and REFUSES everything else at
construction rather than at the draw. A family object that looks built
and has quietly lost its guard is the failure this replaces, so the
refusal has to come before the object exists. `TRUE` gets its own
sentence in the message, because "there is nothing for TRUE to mean" is
the part a reader needs: deriving is already the default, and reading
`TRUE` as an opt-in would make the one spelling that says yes the one
that turns the guard off.

No shipped family passes the argument at all, so this was a hazard
rather than a live defect. It is fixed because the derivation's whole
claim rests on that argument behaving.

### `newdata` is not a way round it, and the message now says so

`simulate(fit, newdata = <a real two-column schedule>)` is refused too,
and correctly: the formula names one column twice, so `reward1` and
`reward2` are re-evaluated from that same column whatever data frame is
supplied. So a user whose record holds the received outcome alone, who
is the exact user item 1.5 exists to serve, can fit and cannot draw
from that fit. Refusing is still right, because a draw at exact chance
presented as a model prediction is worse than an error, and the table
above is what proves it. The refusal now names the route back: refit
with a column per option, or build a schedule with
`frm_task_design()`.

### False alarms

Zero by construction on any real design, and asserted on four: the
shipped `frm_task_design()` schedules for `bandit2arm`, `igt`,
`bandit4arm_restless` and `twostep` all draw without complaint, on
every family that reads them. The only data the check refuses is data
in which every option pays the same on every trial, which is a task
with nothing to learn and nothing to simulate either way.

### The compatibility table

Three rows changed, and one of them is a correction to a RELEASED
claim rather than to this lane's work:

- the shared `simulate` row now says the check runs on all eight
  families and through `frm_task_simulate()` as well;
- `prl_fictitious()` / `reward()` said since 0.1.0 that "the second
  column enters the likelihood and not just the simulator". It does
  not, and the bitwise measurement above settles it. A user choosing a
  family on that table would have chosen this one for a property it
  does not have;
- `rlddm()` / `reward()` and `ts_par7()` / `payoff()` now say which
  route the check runs on, since `simulate()` refuses them anyway.

The shared `simulate` row stays `works` rather than becoming
`conditional`. The vocabulary's `conditional` marks a capability
limited by an OPTION (the way `predict` is, where `type = "response"`
refuses); this is a data guard on a shape that is wrong for a draw. That
is a judgement rather than a measurement and a reviewer may want it the
other way.

### Test

`extensions/frmtmb.learn/tests/testthat/test-counterfactual.R`, new.
Against the unmodified base exported from HEAD and installed into its
own library: **18 pass, 14 fail**. After: **32 pass**. It covers the
bitwise objective identity at one parameter vector for
`bandit2arm_delta()` and for `prl_fictitious()`, the converged fits
agreeing, the refusal from `simulate()`, `frm_simulate()` and
`frm_task_simulate()`, the refusal naming the term and the missing
column and `newdata`, `newdata` really being refused, the four
`reward()` families and the four `payoff()` families each refusing a
duplicated schedule and drawing from a real one, the partial-duplicate
case being accepted, and `stage2()` not being read as a schedule.

The case that asserted the opposite, "a family that READS both columns
is not touched", is gone. Parameter values in the loops are read from
each family's own `init_dpars` rather than written down, so a family
renaming a parameter cannot make the file assert nothing.

## 3. Item 1.6, `frm_cross_spectrum()` ingestion

### Gaps

`NA` in either signal now marks a sample the pair does not have. The
record is cut at it, and segments are laid inside the clean spans, so no
transform crosses a gap. `NaN` counts as `NA`; an infinity is refused by
name, because it is a value the arithmetic cannot use rather than a
value the record is missing. A sample is usable only where BOTH signals
have it, which is asserted: rejecting a span in `x` alone, in `y` alone,
or in both gives one answer.

The contract, in `?frm_cross_spectrum` under "Gaps in the record":

- the segment LENGTH is the usable sample count divided by `segments`;
- each clean span supplies as many whole segments as fit in it, in
  record order, up to `segments` in total;
- `n` is the count the record supplied, not the count asked for.

The cap at `segments` is load-bearing for backward compatibility.
`usable %/% segments` rounds down, so one ungapped span can hold more
than `segments` whole blocks (a 10-sample record cut into 4 holds 5
blocks of 2), and without the cap an ungapped record would come back
with more degrees of freedom than before.

The cap also settles when the shortfall happens, and it is narrower than
it first looks. `seglen = usable %/% segments` gives
`floor(usable / seglen) >= segments` for any single span, so a record
with ONE clean span always yields the full `segments` however much was
rejected from its ends. Only a record broken into several spans loses
anything, because each span drops its own remainder.

One consequence is worth documenting and now is. The segment length
comes from the GLOBAL usable count, so on a multi-span record asking
for FEWER segments can be refused where more would be accepted.
Measured on 2500 samples with spans of 700, 1043 and 699
(`small-f10-gaps.R`): `segments = 2` is refused, because it wants 1221
samples in a row and no span holds them, while 4, 8 and 16 give
`n` of 3, 7 and 14. The refusal's own advice, "Raise one of the three",
is right; what was missing was any statement that raising `segments`
is what works.

In the matrix form each column is its own record, so a column with more
rejected samples gets a shorter segment, its own frequency grid and its
own `n`. Probed on three columns of 1024 with 0, 100 and 300 samples
rejected at `segments = 4`: 127, 115 and 90 rows, `n = 4` in all three.
That is correct, and it is why `freq` is a column of the frame rather
than an attribute of it.

Nothing warns when a gapped record yields fewer segments than asked
for. That was considered and rejected: on a gapped record the shortfall
is the normal case, so a warning would fire on every correct call, which
is what the lane rules forbid. `n` is not a silent channel: it is a
returned column, the density reads it per row, and every standard error
downstream is formed from it.

### A note for item 3.2

Gap splitting partly anticipates Phase 3.2's "list of epochs of unequal
length with a `group` vector". A user can already concatenate epochs
with one `NA` between them and get the epoch behavior, because a
segment never crosses a gap. What 3.2 still adds over that is a `group`
vector on the output and epochs whose LENGTHS differ enough that one
segment length does not serve them, which the NA route does not
address: the segment length here is one number for the whole record.
3.2 should be written on top of `cp_spans()` and `cp_blocks()` rather
than beside them.

### Hann

`window = c("none", "hann")`, applied to each segment before its
transform, scaled by `sqrt(N / sum(w^2))`, which is exactly `sqrt(8/3)`
for the periodic raised cosine. The periodic form and not the symmetric
one, because the three-tap transform the neighbor correlation rests on
is a property of the periodic window.

### The degrees of freedom it costs, which is the number the plan asked for

**None, standing alone.** Segments are disjoint, so a window inside each
one leaves the number of independent complex draws per frequency at
one per segment. Measured the way every other row of the package's
degrees-of-freedom table is measured, `1 / mean(coherence)` at a true
coherence of zero, 3000 replicates on white noise at 512 samples
(`small-hann-probe.R`):

| configuration | nominal | measured (se) |
| --- | --- | --- |
| `segments = 8` | 8 | 8.008 (0.130) |
| `segments = 8, window = "hann"` | 8 | 8.021 (0.127) |

**Everything, in combination with `smooth`:**

| configuration | nominal | no window | Hann |
| --- | --- | --- | --- |
| `segments = 4, smooth = 2` | 8 | 8.188 | **5.660** |
| `segments = 2, smooth = 4` | 8 | 8.023 | **4.983** |
| `segments = 1, smooth = 8` | 8 | 8.047 | **4.678** |

Standard errors 0.067 to 0.132. So `window = "hann"` with `smooth > 1`
is REFUSED by name, exactly as `tapers > 1` with `smooth > 1` already
is, and for the identical reason: both widen the same spectral window.
`window` with `tapers` is refused too, since both are a taper on the
same segment.

### The warning the task gave, checked

The task warned that declaring the equivalent degrees of freedom as
`tapers` does not rescue a smoothed spectrum from core's `whittle()`
refusal at kernel widths of 7 and above
(`dev/reviews/2026-09-08-spectral2.md`: an honest `tapers = 7`
declaration is refused 36.7 percent at 127 ordinates, 88.6 percent at
1023). **This lane's `n` accounting does not rely on that working**, and
that is the design decision the measurement above forces. The
alternative was to admit `window = "hann"` with `smooth` and declare a
reduced `n`. It was rejected for two reasons, one measured here and one
read from the review:

1. There is no one number to declare. Against the nominal 8, the same
   window delivers 71, 62 and 58 percent at widths 2, 4 and 8, so any
   single correction factor is wrong at two of the three.
2. The review established the general form: a wide frequency smooth
   violates the independence the check assumes, whatever number is
   attached to it, because the correlation is a property of the estimate
   rather than of the declaration.

The mechanism is the same in both places, and this lane measured it on
its own code to confirm rather than assume it. Log auto-spectrum
autocorrelation on white noise, 400 replicates at 1024 samples:

| window | lag 1 | lag 2 | lag 3 |
| --- | --- | --- | --- |
| none | -0.004 | -0.004 | -0.000 |
| hann | **0.301** | 0.010 | -0.002 |

That reproduces `dev/freq-findings.md`'s 0.308 at lag 1 and -0.003 at
lag 3 for core's own Hann taper, on a different code path, which is the
construction check the lane rules ask for.

### What `n` does not say, and a pre-existing gap it exposes

`n` counts the draws behind ONE row and the table above shows a window
leaves that alone. What a window changes is the relation BETWEEN rows.

**The first version of this file measured that at a call the function
refuses.** The 0.301 it quoted came from the prototype at
`segments = 1`, and `frm_cross_spectrum(x, y, segments = 1,
window = "hann")` is refused by the function itself, because the
degrees of freedom are then 1. The record is not wrong (core's
`dev/freq-findings.md` has 0.308 for the same construction, a single
raw periodogram) but it is not the frame a user gets. The mechanism is
the log transform: at one segment an ordinate is exponential and its
log has variance `pi^2/6`, which dilutes the shared part, and the
correlation rises toward the raw value as segments are averaged.

Re-measured on the package's own code path, at calls it accepts, 400
replicates of white noise at 1024 samples, as the autocorrelation of
`log(w11)` (`small-f34-spectral.R`). The last column is the variance
inflation factor of a mean over the rows, `1 + 2 * sum` of the positive
correlations, which is what a smooth in frequency actually pays:

| configuration | lag 1 | lag 2 | lag 3 | inflation |
| --- | --- | --- | --- | --- |
| `segments = 8` | -0.010 | -0.015 | -0.021 | 1.00 |
| `segments = 8, window = "hann"` | 0.397 | -0.009 | -0.027 | 1.79 |
| `segments = 4, window = "hann"` | 0.393 | -0.000 | -0.015 | 1.79 |
| `segments = 8, tapers = 2` | 0.525 | 0.104 | -0.035 | 2.24 |
| `segments = 4, tapers = 4` | 0.745 | 0.502 | 0.271 | 4.07 |
| `segments = 1, tapers = 8` | 0.876 | 0.744 | 0.615 | 8.02 |

The inflation column is summed over TWELVE lags. The first version of
this table stopped at four, which is right for every row but the last:
at `tapers = 8` the correlation is still positive past lag 4, and
truncating there reports 6.45 for a row that is really 8.02, which
makes the worst case look better than it is. Every other row is
identical at four lags and at twelve, because the correlation has died
by lag 3 (`small-r34.R`).

Two corrections follow. At the DEFAULT `segments = 8` the Hann lag-1
correlation is 0.397 and not 0.301. And "about one independent
frequency in every three" was the number of TAPS in the window rather
than a measured cost: the inflation factor is 1.79, so it is about one
in 1.8.

**`tapers` is the larger cost, not a footnote to the window's, and that
is what the first version got backwards by leaving it without a
number.** A `tapers = 4` frame carries about one independent frequency
in four, and at `tapers = 8` ordinates are still correlated 0.6 three
bins apart. That is a property of a RELEASED option: users have been
fitting `s(freq)` on tapered frames since 0.1.0 with nothing saying so.
It now has its own NEWS bullet and a row in the help table.

**`n` itself is right in every one of these configurations**, which
this lane raised as a possible correctness problem, and which had to be
settled before the paragraph above could be written: a wrong `n` would
be a correctness bug in a released package. Re-measured on the tighter
estimator, pooling every retained ordinate rather than reading one bin
per replicate, 600 replicates at 512 samples (`small-f4-n.R`):

| configuration | nominal | effective n |
| --- | --- | --- |
| `segments = 8` | 8 | 7.968 (se 0.051) |
| `segments = 8, window = "hann"` | 8 | 7.964 (se 0.051) |
| `segments = 2, tapers = 4` | 8 | 8.033 (se 0.026) |
| `segments = 1, tapers = 8` | 8 | 8.006 (se 0.018) |
| `segments = 1, tapers = 16` | 16 | 16.144 (se 0.039) |

Per-row degrees of freedom and cross-row correlation are different
quantities and `n` claims only the first, so no released version is
wrong and no correctness bullet is owed there.

That table also settles the Hann claim on a better instrument than the
one this lane first used. Reading one frequency bin per replicate, the
way the package's own degrees-of-freedom test does, gives standard
errors of 0.13; pooling every ordinate gives 0.051 for the same work.
On the tighter one the windowed arm sits 0.004 below the untapered arm
with a combined standard error of 0.072, so "a window costs no degrees
of freedom" is a result and not a blind spot. The shipped test still
reads one bin, because there it compares two arms on the SAME draws
with a standard error the run itself computes, which is the paired
comparison rather than the absolute level.

What is still not measured is the coverage loss itself: how far an
`s(freq)` band on a tapered frame falls below nominal. That belongs to
item 2.6, and the inflation factors above set its scale, with
`tapers = 4` the one to measure first.

### The identity against `stats::spec.pgram()`

On the windowed segments of one record, summed the way the accumulator
sums them, agreement is 4.8e-16 of the largest ordinate; the boxcar arm
is 1.5e-16. The scaling constant `N / sum(w^2) = 2.666666667` is pinned
by the same assertion. The power scale does not move when the window is
turned on: on white noise of unit variance, 400 replicates,
`mean(w11 / n)` is 1.0013 (se 0.0022) with no window and 0.9963
(se 0.0032) with Hann.

### Test

`extensions/frmtmb.coupling/tests/testthat/test-cross-spectrum.R`,
extended. Against the unmodified base the file scored **60 pass, 10
errors**; after, **96 pass**. The new assertions are the three-span sum
identity, a step across each gap leaving the answer untouched, either
signal's gap rejecting the pair, `n` reporting what the record supplied,
the two new refusals by name, two windowed configurations added to the
positive-definiteness sweep, the `spec.pgram()` identity, the effective
degrees of freedom against the untapered arm on the same draws, and the
adjacent-frequency correlation at lags 1 and 3.

Every tolerance in the new assertions is either testthat's relative
default or a multiple of a standard error the run itself computes. No
absolute numeric tolerance was written.

## 4. What was decided against

- **A median arm on the eam units guard.** Section 1. It would close the
  contaminated-millisecond miss; it was not measured for false alarms
  and the plan's statistic is the minimum. The reason is that one and
  no longer the "already broken in seconds" argument, which was wrong.
- **A warning when a gapped record yields fewer segments than asked
  for.** Section 3. It would fire on every correct multi-span call.
- **Admitting `window = "hann"` with `smooth` and declaring a reduced
  `n`.** Section 3, with the three measurements that rule it out.
- **Registering `reward` at arity 1.** Section 2. It refuses the
  spelling every existing user writes. A SECOND term name at arity 1
  was not taken either, but for a different reason: it works, and it
  is a second spelling for the same thing, so it belongs with the core
  seam rather than beside it.
- **Moving the learn `simulate` compat row to `conditional`.** Section
  2. A judgement, flagged as one.
- **Exempting any family from the counterfactual guard.** Section 2.
  The one exemption this lane made was wrong, and the guard is now
  derived from the terms a family names so that the next one has to be
  argued for explicitly.
- **Reading `counterfactual = TRUE` as an opt-in.** Section 2. It is
  refused instead, because deriving is already the default and there is
  nothing left for TRUE to mean, so any reading of it would make the
  one spelling that says yes ambiguous.
- **Making `ln_counterfactual_of()` guess at unrecognized term names.**
  Section 2. A predicate that guessed would be wrong in the other
  direction, and the honest close is a test over the package's own
  registration table.

## 5. Defects found and not fixed

1. **The coverage loss of an `s(freq)` band on a tapered or windowed
   frame.** Section 3 measures the variance inflation factor, 1.79 for
   Hann at the default and 4.16 for `tapers = 4`, but not the coverage
   itself. It belongs to item 2.6.
2. **`gd_tri_df` is an unresolvable roxygen link** at
   `extensions/frmtmb.eam/R/gddm.R:169`. Pre-existing, raised on every
   roxygenise, not touched by this lane.
3. **`rlddm()` draws the LOWER-paying arm on a real schedule.** On the
   two-armed design where arm 1 pays with probability 0.7,
   `frm_task_simulate()` takes arm 1 on 0.1443 of late trials at a
   positive drift, where every other family takes it on 0.75 to 0.86.
   Found while measuring section 2 and NOT settled. `logp` and `draw`
   share one expression, so they cannot disagree with each other and a
   shared sign convention would be invisible to the Stan identity test.
   It is pre-existing and out of this lane's scope; one probe in the
   package that owns `rlddm()` would settle it. Recorded because it was
   seen, not because it was diagnosed.
4. **A released compatibility row was false.** The `prl_fictitious()` /
   `reward()` row has said since 0.1.0 that the second column enters
   the likelihood. It does not, and this lane's own exemption was built
   on it. Fixed here, listed as a defect because it shipped.
5. **A schedule under a term name that is not `reward` or `payoff` is
   unguarded.** Section 2. Not fixed in code, because a predicate that
   guessed would be wrong the other way; closed by a test over
   `ln_aterms` that fails the moment a fourth multi-column term is
   registered, which is the earliest anyone could catch it. Item 5.1 is
   the family that will meet it.
6. **The coupling example NOTE is warm-up, not contention alone, and
   not this lane's.** The review timed both examples with
   `tools::Rd2ex`, one fresh process per measurement on a quiet
   machine: `frm_coherence` costs 1.86 to 2.13 s and `cross_wishart`
   1.59 to 2.13 s on FIRST execution, and a second execution in the
   same process costs 0.03 s, so the whole of it is tape build and
   byte-compilation paid once. Under load the same work reached 5.4 to
   6.5 s. So the examples sit a factor of 2.4 to 3 under the 5 second
   threshold on their own and cross it only when the machine is busy,
   on whichever example happens to be running. Across nine coupling
   checks, the review's six quiet ones and this lane's three, it fired
   twice and named two different examples, and this lane's third run
   came back plain OK. Neither is in this diff. If it keeps firing the
   fix belongs to whoever owns those examples, which at 4096 samples
   and 16 segments do more work than showing a coherence fit needs.

## 6. Verification

### Test files, one R process each, `NOT_CRAN=true`

frmtmb.eam, 20 files:

| file | pass | fail |
| --- | --- | --- |
| test-bracket-access.R | 1 | 0 |
| test-brms-parity.R | 13 | 0 |
| test-defects.R | 59 | 0 |
| test-density.R | 132 | 0 |
| test-extension-api.R | 14 | 0 |
| test-family.R | 51 | 0 |
| test-gddm-family.R | 77 | 0 |
| test-gddm-gradients.R | 22 | 0 |
| test-gddm-recovery.R | 28 | 0 |
| test-gddm-reference.R | 166 | 0 |
| test-gddm-solver.R | 94 | 0 |
| test-lba.R | 109 | 0 |
| test-message-uniqueness.R | 4 | 0 |
| test-moments.R | 29 | 0 |
| test-rdm-gng.R | 221 | 0 |
| test-sampling.R | 97 | 0 |
| test-simulate-density.R | 69 | 0 |
| test-surface.R | 50 | 0 |
| test-units.R (new) | 32 | 0 |
| test-variability.R | 140 | 0 |

Total 1408.

frmtmb.learn, 11 files:

| file | pass | fail | skip |
| --- | --- | --- | --- |
| test-bracket-access.R | 1 | 0 | 0 |
| test-clean-session.R | 4 | 0 | 0 |
| test-counterfactual.R (new) | 67 | 0 | 0 |
| test-engine.R | 19 | 0 | 0 |
| test-factorization.R | 35 | 0 | 0 |
| test-families.R | 45 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 |
| test-reference.R | 12 | 0 | 0 |
| test-stan-identity.R | 0 | 0 | 9 |
| test-surface.R | 146 | 0 | 0 |

Total 333 passing, 9 skipped. `test-counterfactual.R` grew twice: 14
in the first version, 32 after the first review when it reached all
eight families, and 67 after the re-check, when it gained the rule that
produces the guard as well as the guard's behavior. Every one of the 36
files was rerun after each round, one R process each.

frmtmb.coupling, 6 files:

| file | pass | fail |
| --- | --- | --- |
| test-bracket-access.R | 1 | 0 |
| test-coherence.R | 58 | 0 |
| test-cross-spectrum.R | 96 | 0 |
| test-cross-wishart.R | 211 | 0 |
| test-message-uniqueness.R | 4 | 0 |
| test-surface.R | 58 | 0 |

Total 428. The `test-stan-identity.R` skips are the gated tier, which
needs `rstan` compilation and was not run.

### The three seen failures, before the change

Re-taken after the review, against the UNMODIFIED BASE and with the
FINAL test files. The base is `git archive HEAD` (780dec1) unpacked
into the scratchpad and installed into a library of its own,
`smallitems-lib-base`, so the count is what a reader running the
shipped file against the shipped package gets:

| package | file | before | after |
| --- | --- | --- | --- |
| frmtmb.eam | test-units.R | 5 pass, 24 fail, 1 error | 32 pass |
| frmtmb.learn | test-counterfactual.R | 28 pass, 24 fail, 5 err | 67 pass |
| frmtmb.coupling | test-cross-spectrum.R | 60 pass, 10 errors | 96 pass |

The learn file was run against the base a second time after the
re-check, with the assertions the re-check added, and the count above
is that run. The re-check's own two additions were also run against the
code they pin, which is this lane's own previous version rather than
the base: a copy of the learn source with `ln_read_counterfactual()`
and `ln_not_schedule_terms` stripped and nothing else changed, in its
own library, scores **53 pass, 10 fail, 1 error**. The error is the
missing classification table and the ten failures are the eight values
that used to disable the guard, the `TRUE`-specific sentence, and the
character-vector spelling that was not accepted before.

The first version of this file reported 6/23/1 and 60 pass with 9
errors, taken before the last edits to those two files. A count that
does not match the shipped file is not a count, so these were re-run
rather than adjusted. The error in the eam row is `ddm_seconds_ceiling`
not existing yet; the learn row now starts from 18 rather than 8
because the file has grown to cover eight families.

### `R CMD check --as-cran`

With `RSTUDIO_PANDOC` and TinyTeX on `PATH`, from a tarball built by
`R CMD build`, against the private library. Two rounds: the first on
the pre-review tree with the CRAN-incoming remote checks both on and
off, and the second on the POST-REVIEW tree with them off, which is the
row that describes what is being merged.

| package | remote on | remote off | after review | after re-check |
| --- | --- | --- | --- | --- |
| frmtmb.eam | 1 WARNING, 1 NOTE | 1 NOTE | 1 NOTE | **1 NOTE** |
| frmtmb.learn | 1 WARNING | OK | OK | **OK** |
| frmtmb.coupling | 1 WARNING | 1 NOTE | 1 NOTE | **OK** |

No ERROR in any run. eam's NOTE is the expected one on the HTML manual,
"Skipping checking math rendering: package 'V8' unavailable", which the
lane rules name as the baseline.

The WARNING is the same in all three pre-review runs and is entirely
the remote half of `checking CRAN incoming feasibility`: "Strong
dependencies not in the CRAN or BioC software repositories: frmtmb"
(and `frmtmb.eam` for learn), which needs the CRAN listing to notice,
plus a 301 on the `https://aforren1.github.io/frmtmb/<pkg>` URL in
`DESCRIPTION`, which wants a trailing slash. This lane's diff changes
no `DESCRIPTION` at all, and the review confirmed the same WARNING on
the unmodified base, so it is structural.

coupling's NOTE is a clock rather than a finding, and the third run
came back plain `OK`. Over the three runs here it fired twice on two
DIFFERENT examples, "Examples with CPU (user + system) or elapsed time
> 5s" for `frm_coherence` at 6.53 s and for `cross_wishart` at 5.39 s,
and the review's six runs on a quiet machine fired zero times. Both
examples are untouched by this diff (`R/coupling.R`,
`R/cross-wishart.R` and their Rd files are all unmodified). Defect 6
carries the timing that settles it.

## 7. What the review changed, and the lesson from each

`dev/reviews/2026-09-08-smallitems.md`. Every number in it that I
re-measured on my own construction reproduced, and two of its findings
were defects in this lane's work rather than in its prose.

**F1, the counterfactual guard named five families and needed eight.**
The one exemption was wrong. Fixed by deriving the guard from the terms
a family names rather than declaring it per family, so all eight are
covered and a ninth cannot be added without the question being asked.
`prl_fictitious()`, `rlddm()` and `ts_par7()` each drew a chance-level
data set through the shipped guard with no error, measured in section
2. **The lesson: an exemption justified by what code looks like it does
is how five became six.** The update I read as "moves the unchosen
value" is a sign flip of the CHOSEN payoff, the family's own help page
said so, and a compatibility row said the opposite; I read the row.
Nothing but the measurement should have settled it, and the measurement
takes one fit.

**F2, the units false-alarm rate was measured where it cannot fail.**
Sweeping only 200, 400 and 12000 trials pins the minimum at the
non-decision-time floor, so the guard reduces to `ndt > 20` and the
answer is a property of the sweep. At 20 trials of a 40 to 70 second
task the rate is 0.07 to 0.74. The guard is unchanged and right; the
sentence describing its cost is rewritten in three places, and the
small-n cells are in the sweep now. **The lesson: a false-alarm rate
measured on one design family is a rate for that family.**

**F3, the Hann correlation was quoted at a call the function refuses.**
0.301 is the value at `segments = 1`, which
`frm_cross_spectrum()` refuses for want of degrees of freedom. At the
default it is 0.397 and the inflation factor is 1.79, so "one
independent frequency in three" was the tap count and not a
measurement.

**F4, `tapers` is the larger cost and had no number.** 0.753 at
`tapers = 4`, an inflation factor of 4.16, against 1.79 for the window.
It now has a NEWS bullet, because it is a released option and users
have been fitting `s(freq)` on tapered frames since 0.1.0.

**F5, "already broken read as seconds" was false.** The seconds reading
of a contaminated record reports nothing at all, diagnoses clean and is
47 percent out. It is a second silent wrong answer, which raises item
3.5's rank rather than excusing the miss.

**F6, the refusal had no escape hatch and did not say so.** `newdata`
cannot rescue it, correctly, because the formula names one column
twice. The message now names the route back.

**F7, the log-likelihood identity is stronger than I claimed.** The
objective is bitwise identical as a function of the parameters, not
identical to printed precision, and the shipped test now pins that at a
fixed parameter vector rather than pinning a tolerance on a converged
one.

**F8, the seam blocks the spelling and not the capability.** A second
term name at arity 1 reaches the same key with no core change. The
filed seam is re-shaped as a minimum and a maximum with naming from the
maximum, which is smaller than the range I proposed and removes the
open question I left.

**F9, two of three before-counts did not match the shipped files.**
Re-taken against the unmodified base with the final files.

**F10, a multi-span record can refuse at low `segments` and pass at
high.** Documented, with the construction reproduced here.

The review also settled a suspicion this lane raised as a possible HIGH:
shipped `n` for `tapers > 1` is NOT wrong. Cross-row correlation and
per-row degrees of freedom are different quantities, and `n` claims
only the second. No released version is affected and no correctness
bullet is owed there.

## 8. What the re-check changed

`dev/reviews/2026-09-08-smallitems.md`, second pass. It confirmed the
derivation covers all eight families and excludes `stage2()`, built a
family naming `reward1`/`reward2` and declaring nothing and found the
guard derived onto it, and reproduced every figure in section 2 to the
digit once it found the construction in `small-f1-eight.R`. Two things
were open.

**R1, `counterfactual = TRUE` silently turned the guard off.** Fixed by
validating the argument at construction. Section 2 carries the table of
values that used to disable it. **The lesson: the argument a guard
hangs on is part of the guard.** I wrote a three-branch `if` with no
`else` that refuses, in the one place where the cost of a wrong value
is the silent loss of the thing being built. The derivation's claim was
"a ninth family cannot be added without the question being asked", and
I had not checked what the affirmative answer did.

**R2, a schedule under a third name is unguarded.** Not fixed in code,
because a predicate that guessed at unrecognized names would be wrong
the other way. Closed by a test over `ln_aterms`, the package's own
registration table, which is the only place the miss can be seen before
a user meets it. Item 5.1 is the family that will meet it, and it is
named in the source so that whoever writes it reads this first.

**R3**, the `tapers = 8` inflation factor was truncated at four lags:
6.45 where twelve lags give 8.02. Every other row is unchanged, because
its correlation has died by lag 3. Corrected in the help table, the
NEWS bullet, the vignette and section 3. **The lesson: a summary
statistic with a cutoff in it has to say where the cutoff is**, and
mine flattered the worst row.

**R4**, the false-alarm rate is not a function of the median response.
At a median near 50 seconds it runs from 0.185 to 0.935 depending on
whether the task is slow from a far boundary or from weak evidence.
The table is re-indexed with the fastest response beside each rate,
which is what the guard actually reads. **The lesson: index a rate by
the quantity the code tests, not by the one that describes the design.**

**The stale before-count** in section 3 said 9 errors where section 6
said 10. Corrected; the count re-taken against the base is 10.

**The coupling example NOTE** is settled and recorded as defect 6: 1.6
to 2.1 s of one-off warm-up per example, 0.03 s on a second run in the
same process, 5.4 to 6.5 s under load, zero NOTEs in six quiet checks.
Nothing to fix and nothing of this lane's in it.
