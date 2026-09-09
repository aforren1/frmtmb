# Review: lane `smallitems`, items 1.4, 1.5 and 1.6

Reviewer, 2026-09-08. Worktree `frmtmb-wt-smallitems` on branch
`wt-smallitems`, uncommitted, based on 780dec1. Reference: the same
commit exported with `git archive` into the scratchpad, so the main tree
was never touched. Two private libraries, `rvsm-lib` (the lane) and
`rvsm-main-lib` (the unmodified base), plus `rvsm-lib-before` for a
learn build with only the `counterfactual =` declarations stripped.
Every script and log is prefixed `rvsm-`.

## Verdict

**Mergeable after named fixes.** The three items are built, the test
suites are green, `R CMD check --as-cran` is clean once the remote
CRAN-incoming half is turned off, and every quantitative claim I
re-measured reproduced. Two fixes are needed before merge and both are
small:

1. **F1**, the counterfactual guard is declared on five families and
   should be declared on all eight. `rlddm()`, `ts_par7()` and
   `prl_fictitious()` read only the chosen option's payoff exactly as
   the five do, and `rlddm()` and `prl_fictitious()` each draw a
   chance-level dataset today, through the guard, with no error. The
   stated reason for exempting `prl_fictitious()` is false as a matter
   of code reading, which is the exact shape the brief warned about.
2. **F2**, the units guard's false-alarm claim is measured only at
   trial counts where it cannot fail. At 20 trials of a 40 to 70 second
   task the rate is 0.07 to 0.74, not zero. The guard is still right;
   the sentence that describes its cost is not.

The other findings are corrections to prose and to numbers in the
findings file, not to code.

**What I changed in the worktree: this file, and nothing else.** F1 is
the only code fix, and it is not the "small and clearly right" kind: it
reverses a decision the lane took deliberately and wrote into four
documents, and it needs a test that currently asserts the opposite
inverted. A reviewer making that change quietly would hide the fact that
it was ever made. It is filed here with the measurement instead.

## What I could not falsify

Recorded first, because most of the lane's numbers held.

| claim | lane | this review |
| --- | --- | --- |
| the two spellings give one log-likelihood | -1445.79870595 twice | bitwise identical, and identical at my own seed too |
| Hann costs no degrees of freedom | 8.021 (0.127) against 8.008 (0.130) | 7.924 (0.051) against 7.992 (0.051), 600 reps, all 31 ordinates |
| `sqrt(8/3)` scaling | exact | `N/sum(w^2)` is bitwise `8/3` |
| Hann against `stats::spec.pgram()` | 4.8e-16 | 3.06e-16 windowed, 1.24e-16 boxcar |
| Hann with `smooth` loses | 5.66 / 4.98 / 4.68 | 5.750 / 5.022 / 4.667 |
| lag-1 log-spectrum correlation, lane's cell | 0.301, lag 3 -0.002 | 0.306 (0.002), lag 3 -0.003 (0.002) |
| `reward(pay1)` is refused by core's parser | quoted message | reproduced exactly, and the arity-1 registration refuses the two-column spelling |
| the gap contract | the sum of the clean pieces | exact to 0.00e+00 against an independent reimplementation, 7 layouts by 4 segment counts |
| ungapped answers unchanged | backward compatible | `identical()` TRUE on 6 configurations, lane against base |
| the guard fires once per fit | 1 at `frm()`, 0 elsewhere | confirmed over 11 entry points |

## Findings

### F1. HIGH. The guard names five families and needs all eight

`ln_family(counterfactual =)` is declared on `bandit2arm_delta()`,
`bandit2arm_dual()`, `bandit4arm2_kalman_filter()`, `igt_pvl_delta()`
and `igt_orl()`. The package has eight families that carry a
`reward()` or `payoff()` term. All three of the others have exactly the
property the guard exists for.

#### `prl_fictitious()`, the deliberate exemption

The lane's reason, in the findings, the NEWS bullet, the help page and
the compatibility note: "its counterfactual update READS both columns,
so identical columns there are a statement about the task that the
likelihood has already used".

Read the update (`R/prl-fictitious.R:86-101`):

```r
oc <- c1 * d[["reward1"]] + c2 * d[["reward2"]]
t1 <- c1 * oc - c2 * oc
t2 <- c2 * oc - c1 * oc
```

`c1` and `c2` are the 0/1 CHOSEN indicators, so `oc` is the chosen
option's payoff and nothing else in the update touches `reward2`. The
fictitious update is a SIGN FLIP of the received outcome, not a reading
of the unchosen column. `prl_fictitious()` reads one column, exactly
like the five that are declared.

Measured (`rvsm-p1c-prl.R`, lane library, 30 subjects by 100 trials,
`alpha = 0.3`, `bias = 0`, `tau = 2`, arm 1 paying 0.7 and arm 2 0.3):

| `prl_fictitious()` | value |
| --- | --- |
| log-likelihood, `reward(pay1, pay2)` | -953.03713558427876 |
| log-likelihood, `reward(rec, rec)` | -953.03713558427876, `identical()` TRUE |
| log-likelihood, unchosen entries replaced by `N(100, 50)` | -953.03713558427876, `identical()` TRUE |
| P(arm 1) over trials 51-100, real schedule, 40 draws | 0.8560 (se 0.0028) |
| P(arm 1), `pay2 <- pay1`, 40 draws | **0.4990 (se 0.0092)** |

So the fit is bitwise unchanged by the second column, which is the
lane's own definition of a family that "carries the column without
reading it", and the draw from a duplicated schedule is exact chance:
`|t|` against 0.5 is 0.10. The drawn subjects do wander (spread across
draws 0.058, range 0.352 to 0.613), because the sign-flip update makes
them perseverate on whichever option they happened to pick, but with no
systematic preference at all.

`prl_fictitious()` has a simulator, so the fitted route is open too.
Same design, `simulate()` on the fitted object, 40 seeds each
(`rvsm-p1d-prlsim.R`), with the two log-likelihoods again bitwise
identical:

| `simulate()` on a fitted `prl_fictitious()` | P(arm 1), trials 51-100 |
| --- | --- |
| observed in the data being fitted | 0.8220 |
| `reward(pay1, pay2)` | 0.8123 (se 0.0042) |
| `reward(rec, rec)` | **0.4862 (se 0.0074)** |

Both routes, `sim_ctx` and `frm_task_simulate()`, are open on this
family. It is the same wrong answer the lane closed for five families,
left open on the sixth by an explicit decision.

**The family's own help page already says this**, and says the
simulator is what needs the column (`R/prl-fictitious.R:15-18`, present
at 780dec1 and unchanged by the lane):

> The counterfactual outcome is the NEGATIVE of the realized one, which
> is the model's assumption of an anticorrelated task rather than a
> reading of the second `reward()` column. The second column is still
> what the simulator pays a counterfactual choice with.

What the lane read instead is the COMPATIBILITY ROW, at
`R/zzz.R:311-315` on the base commit and `:324-328` after the lane's
edit, which says the opposite:

> Both columns are READ here rather than only carried: counterfactual
> updating moves the unchosen option's value too, so the second column
> enters the likelihood and not just the simulator.

That row is false, it is pre-existing, and it has been in the released
table since 0.2.0. The bitwise measurement above settles it: replacing
the unchosen entry with `N(100, 50)` noise leaves the log-likelihood
unchanged to the last bit, so the second column does not enter the
likelihood at all. Two documents in one package have said opposite
things about the same family, and the lane's exemption inherited the
wrong one.

#### `rlddm()` and `ts_par7()`

`rlddm()`'s `update` is character for character the same shape as
`bandit2arm_delta()`'s: `pe2 <- d[["reward2"]] - q2` is multiplied by
the 0/1 indicator `c2`, so the likelihood reads only the chosen arm's
payoff (`R/rlddm.R:193-204`). `ts_par7()` does the same thing over four
stage-two options through its `sel` indicators (`R/ts-par7.R:148-156`).

The lane's own reasoning covers them. Both families set `sim = FALSE`,
so `simulate()` and `frm_simulate()` refuse, but `frm_task_simulate()`
does not, and that is the route the lane guarded second precisely
because it "reads them by those names" and "can duplicate a column in
just as easily" (`R/task.R`). For `rlddm()` it is also the route the
family's own `sim_refusal` message TELLS the user to take: "Use
frmtmb.learn::frm_task_simulate(), which returns whole data frames".

Measured against the LANE library, with the guard installed
(`rvsm-p1-routes.R`). 20 subjects by 100 trials, `alpha = 0.4`,
`drift = 3`, `bs = 1.5`, `ndt = 0.2`, `bias = 0.5`; the statistic is the
proportion of trials 51 to 100 that took arm 1, over 20 draws:

| `rlddm()` through `frm_task_simulate()` | proportion |
| --- | --- |
| the real two-column schedule | 0.1535 (se 0.0048) |
| `pay2 <- pay1` | **0.5034 (se 0.0093)** |

The duplicated schedule draws with no error and no warning, and 0.5034
is exact chance, which is the same signature the lane found for
`bandit2arm_delta()`. `ts_par7()` with all four payoff columns
duplicated drew 800 rows with no complaint on the same run.

The same script confirms `bandit2arm_delta()` refuses, so the guard
works; the list is short.

**Fix.** Declare it on all eight. Better, declare it by DEFAULT from the
`aterms` a family already names, and make the exemption an explicit
opt-out, so that the ninth family cannot be added without the question
being asked. There is no family in this package for which the guard
would be a false alarm: a duplicated column leaves the fit bitwise
unchanged in every one of them, which is the property the guard tests
for, and the draw is then from a task nobody ran in every one of them.

`test-counterfactual.R` should lose its "a family that READS both
columns is not touched" case, which asserts the false statement, and
gain the `prl_fictitious()`, `rlddm()` and `ts_par7()` refusals plus a
bitwise `expect_identical()` on the three `prl_fictitious()`
log-likelihoods, which is what shows the exemption was never earned.

Recorded while there: a separate, pre-existing oddity that is NOT this
lane's. On the real schedule `rlddm()` draws arm 1 on only 15 percent
of late trials, although arm 1 pays with probability 0.7 and the drift
is positive. `logp` and `draw` use the same expression
`drift * (q2 - q1)`, so they cannot disagree with each other, and the
Stan identity test would not see a shared sign convention. It may be
nothing but a convention I read the wrong way round. It is out of this
lane's scope and I did not settle it; it is worth one probe in the
package that owns `rlddm()`.

### F2. MEDIUM. The units false-alarm rate is measured only where it cannot fail

`dev/smallitems-findings.md`, section 1: "So a false alarm needs 20
seconds of encoding and motor time before any evidence accumulates."
That follows from the sweep, but only because the sweep held the trial
count at 200, 400 and 12000. At those counts the smallest response is
pinned just above the non-decision time, so the guard reduces to
`ndt > 20`. A short session of slow decisions is where the minimum sits
far above the floor.

Measured with `frmtmb.eam::ddm_simulate()` (`rvsm-p3b-units.R`),
`mu = 0.15` to `0.2`, `ndt` 1 to 2 seconds, 40 to 100 draws a cell:

| boundary | achieved median response | trials | false-alarm rate |
| --- | --- | --- | --- |
| 14 | 25 s | 20, 40, 200 | 0.000 |
| 18 | 41 s | 20 | 0.070 |
| 18 | 43 s | 60 | 0.000 |
| 20 | 43 s | 20 | 0.175 |
| 26 | 71 s | 20 | **0.740** |
| 26 | 71 s | 60 | **0.420** |
| 40 | 116 s | 20, 60 | **1.000** |

A 20-trial session with a median around 40 seconds is a deliberation,
insight problem solving, Tower of London or matrix-reasoning design, not
a stress test, and the non-decision time in those cells is 1 to 2
seconds rather than 20.

The guard is still worth having. Over 144 cells at the boundary
separations a two-choice task actually produces (2 to 5), with 12 to 80
trials and non-decision times of 0.5 to 5 seconds, the smallest response
ran from 0.580 to 7.588 seconds and there were **0 of 144** false
alarms, so nothing a standard task produces reaches the ceiling. The
exposure starts at median responses above about 30 seconds.

**Fix.** Replace the "20 seconds of encoding and motor time" sentence in
`dev/smallitems-findings.md`, in `extensions/frmtmb.eam/NEWS.md` and in
`ddm_seconds_ceiling`'s `@noRd` block with the measured rate: a
seconds-scale design false-alarms when its median response passes about
40 seconds AND its session is short, at 0.07 at 20 trials and 0.74 at a
70 second median. That is the honest cost, and it is what makes
`frmtmb_eam_units_warning` earn its class.

### F3. MEDIUM. The Hann correlation is measured at a refused call

`?frm_cross_spectrum`, `NEWS.md` and `vignettes/coherence.Rmd` all give
0.301 as the adjacent-frequency correlation of a Hann-windowed frame.
The lane's probe (`small-hann-probe.R`) measures it at `segments = 1`.
`frm_cross_spectrum(x, y, segments = 1, window = "hann")` is refused by
the function itself, because `nseg * tapers * smooth` is then 1. So the
number describes a call the package will not make.

Measured on the package's own code path (`rvsm-p4b-spectral.R`), 400
replicates at N = 1024, `acf` of `log(w11)`. The prototype used for the
`segments = 1` row is bitwise identical to the package where both are
legal (max relative difference 0.00e+00 on `w11`, four configurations).

| configuration | lag 1 | lag 2 | lag 3 | lag 4 | 1 + 2*sum |
| --- | --- | --- | --- | --- | --- |
| segments = 1, hann (the documented cell) | 0.306 (0.002) | 0.015 | -0.003 | -0.003 | 1.63 |
| segments = 2, hann | 0.369 (0.003) | 0.013 | -0.006 | -0.012 | 1.73 |
| segments = 4, hann | 0.398 (0.003) | 0.013 | -0.007 | -0.014 | 1.78 |
| **segments = 8, hann (the default)** | **0.397 (0.005)** | -0.011 | -0.025 | -0.012 | 1.70 |
| segments = 8, none | -0.015 (0.006) | -0.011 | -0.015 | -0.016 | 0.89 |

The mechanism is the log transform. At one segment each ordinate is
exponential, whose log has variance `pi^2/6`, which dilutes the shared
part; at eight segments the log variance is smaller and the correlation
rises toward the raw value of 4/9 that a three-tap window gives.

The unwindowed rows read slightly below zero for the usual reason: the
sample autocorrelation of a series of length `m` has a bias of about
`-1/m`, and `segments = 8` at N = 1024 leaves 63 ordinates, so -0.015
is the bias and not a signal. That is also why the unwindowed inflation
factor is 0.89 rather than 1.00, and why the windowed ones are, if
anything, slight underestimates.

The second half of the same sentence is also wrong in the safe
direction. "About one independent frequency in every three" is the
number of taps in the window, not the measured cost. The variance
inflation factor of a mean over the rows is 1.70, so the frame carries
about one independent frequency in 1.7.

**Fix.** Quote the correlation at `segments = 8`, which is the default
and the shape a user gets, and state the inflation factor rather than
the tap count. Three places: the help section "What a window does to the
rows, which `n` does not say", the coupling NEWS bullet, and the
vignette paragraph.

### F4. MEDIUM. The same sentence understates `tapers` by over 2x

The new help section closes with "The same is true of `tapers` above 1
and always was", with no number. On the same code path and the same
construction:

| configuration | lag 1 | lag 2 | lag 3 | lag 4 | 1 + 2*sum |
| --- | --- | --- | --- | --- | --- |
| segments = 8, tapers = 2 | 0.528 (0.005) | 0.114 | -0.025 | -0.022 | 2.19 |
| segments = 4, tapers = 4 | 0.745 (0.002) | 0.492 | 0.259 | 0.044 | 4.08 |
| segments = 1, tapers = 8 | 0.875 (0.001) | 0.745 | 0.618 | 0.495 | over 6 |

So a `tapers = 4` frame carries about one independent frequency in four,
against one in 1.7 for Hann, and at `tapers = 8` the correlation is
still 0.495 four bins apart. `tapers` is the larger cost, not a
footnote to the window's.

"Was not written down" is not quite right either.
`dev/reviews/2026-09-08-xspec.md` records that tapers "correlate
neighboring bins". What was missing is a user-facing statement, which is
what the lane added, and it should carry these numbers.

**Fix.** One more row in the help section and one NEWS bullet in
`frmtmb.coupling`, which the lane owes and did not write. See
"Consolidation".

### F5. LOW to MEDIUM. "Already broken read as seconds" is not true

`dev/smallitems-findings.md` and the eam NEWS bullet argue that the
contaminated-millisecond miss is acceptable because "That record is
already broken read as seconds, where the same contaminant puts the
`ndt` bound at 9.5e-05, so the miss is a contamination failure rather
than a units failure going unseen."

Measured on one 12000-row record with 5 percent of trials replaced by a
uniform draw over the observed range (`rvsm-p3b-units.R`), truth
`bs = 1.4`:

| reading | smallest response | units warning | other warnings | `diagnose()` | boundary estimate |
| --- | --- | --- | --- | --- | --- |
| milliseconds | 21.578 | fires | large gradient | max grad 0.0222 | 65.72 |
| **seconds** | 0.021578 | silent | **none** | **"No convergence problems detected"** | **2.078** |

The seconds reading is biased by 48 percent and reports nothing at all.
So it is not visibly broken; it is a second silent wrong answer, which
is the class the plan's rule 3 ranks first. The lane's decision stands
on its own merits: the median arm was not measured for false alarms, and
shipping an unmeasured second heuristic inside a guard is what the lane
rules forbid. But the justification should be "the contaminated seconds
fit is itself silent, and item 3.5 is where that is closed", not
"already broken".

**Fix.** One sentence in the findings and in the eam NEWS bullet.

### F6. LOW. The refusal has no escape hatch and its message does not say so

Measured (`rvsm-p1b-escape.R`) on a fit of `reward(rec, rec)`:

| call | result |
| --- | --- |
| `simulate(fit)` | refused |
| `simulate(fit, newdata = <a real two-column schedule>)` | **also refused** |
| `predict(type = "link")`, `summary()`, `frm_value_trace()` | work |
| `fitted()`, `residuals()`, `pp_check()` | refused already, for unrelated reasons |

`newdata` cannot rescue it, and correctly so: the formula names one
column twice, so `reward1` and `reward2` are re-evaluated from the same
column whatever data frame is supplied. The only route left is
`frm_task_simulate()`, which draws from parameters typed by hand rather
than from the fit, and that is what the refusal offers.

So a user whose record holds the received outcome alone, who is the
exact user item 1.5 exists to serve, can fit and can never draw from
that fit. Refusing is still right: a draw at exact chance presented as a
model prediction is worse than an error, and the lane's measurement is
what proves it. But the message should name the route back.

**Fix.** Add one clause to `ln_check_counterfactual()`'s message: once
the schedule is available, refit with the two-column spelling, because
`newdata` cannot change which columns the formula names.

### F7. LOW. The log-likelihood identity is stronger than claimed

`dev/smallitems-findings.md` section 2 says the two fits "agree to
twelve printed digits" and that `identical()` on the two `fixef()`
vectors is FALSE, "which is what a different data vector through the
same tape gives".

Measured (`rvsm-p1-tape.R`, base library, seed 1, 30 by 100). The
OBJECTIVE is bitwise identical at 14 of 14 parameter vectors between
`reward(pay1, pay2)`, `reward(rec, rec)`, and a third spelling whose
unchosen entries are replaced by `N(100, 50)` noise. The GRADIENT
differs at 1 of 14 and 2 of 14, in the last unit in the last place.

At my own seed (`rvsm-p1-silent.R`, 30 by 100, seed 20260908) the fitted
values are bitwise identical as well: `identical(logLik)` TRUE at 0 ulp,
`identical(fixef)` TRUE, and the noise spelling gives the same bits
again. At the lane's seed 1, `alpha` differs by one ulp
(-0.336866175328572326 against -0.336866175328572270) with the
log-likelihood still bitwise identical.

The reason is arithmetic, not data. The update multiplies the unchosen
arm's prediction error by a 0/1 indicator, and `0 * finite` is exactly
0, so the unchosen entry cannot reach the objective at all. The one-ulp
difference is gradient accumulation on a tape whose constants changed,
not a different answer.

So the correct sentence is shorter and stronger: the fit is EXACTLY
right, not right to printed precision, and only the simulator is wrong.
That also settles the remedy question: nothing about the density needs
changing.

**Fix.** Replace the "different data vector through the same tape"
sentence. The noise test is worth adding to `test-counterfactual.R`,
which already has a version of it; making it `expect_identical()` on the
log-likelihood rather than `expect_equal()` would pin the fact rather
than a tolerance.

### F8. LOW. The seam blocks the spelling, not the capability

Verified by construction (`rvsm-p2-seam.R`), which is what the brief
asked for rather than by reading:

- as shipped, `reward(rec)` gives exactly "`reward()` takes 2
  arguments, not 1", from `R/parse.R:276`;
- re-registering `reward` at `arity = 1` refuses `reward(pay1, pay2)`
  with "`reward()` takes 1 argument, not 2", and `reward(rec)` then
  fails one step later on a missing `reward2`, which is exactly what
  `frmtmb.learn/R/zzz.R`'s own comment predicts;
- **but a one-column spelling under a SECOND term name needs no core
  change at all.** `frmtmb_register_aterm("reward1", arity = 1)` keys
  its value at `aterms[["reward1"]]`, which is the same key the
  two-column spelling's first argument already uses, so a family that
  reads `reward2` with a fallback takes both spellings. A probe family
  built in the review fits `choice | reward1(rec) ~ 1` and returns a
  log-likelihood BITWISE identical to `reward(rec, rec)` through the
  shipped family, and `reward(pay1, pay2) + reward1(rec)` is refused as
  a duplicated term.

Two costs of that route the filing should record, neither of which bites
frmtmb.learn today. With both names registered, `aterm_spelling`
returns `reward1(<column>)` for the key `reward1`, so a refusal about
the two-column spelling's FIRST argument would tell the user to write
the wrong thing; and `aterm_base("reward1")` returns `reward1`, so a
family declaring `accepts_aterms = "reward"` would stop accepting it.
frmtmb.learn declares no `accepts_aterms`.

So the seam is worth filing, and what it buys is the SPELLING
`reward(pay)` rather than the capability. The findings should say that,
because "cannot be built in an extension" reads as the stronger claim.

**On the shape of the seam.** The lane proposes a range `arity` in
`R/parse.R`, "about six lines, plus a decision about how a family
declares which arity it got". The declaration half is indeed already
solved by 0.54.0's any-of `required_aterms` and `exclusive_aterms`. The
parser half is larger than one comparison: `arity` is read for its
NAMING convention at `R/parse.R:134` (`registered_aterm_of`), `:167`
(`aterm_spelling`), `:212` (the duplicate check) and `:281` (the key
branch), each as `arity > 1L` or `== 1L`. A range makes every one of
those a `max()`.

A cleaner shape, and the one I would file: store a MINIMUM and a
maximum, and let the naming depend on the maximum alone, so a term whose
maximum arity is above 1 always delivers `name1`, `name2`, and so on. A
`reward()` registered as "1 to 2 arguments" then delivers `reward1` for
both spellings and `reward2` only for the two-column one. The family
declares `required_aterms = "reward1"` unchanged, needs no any-of group
and no `exclusive_aterms` pair, and `registered_aterm_of()`,
`aterm_base()` and `aterm_spelling()` are untouched. That is
`parse.R:276` plus a `max()` at `:212`, and it removes the "how does a
family declare which arity it got" question rather than answering it.

### F9. LOW. The recorded before-counts do not match the shipped test files

Rerun independently against the unmodified base, and for learn against
my own strip of the five `counterfactual =` lines (`rvsm-before.sh`):

| file | lane recorded | measured here |
| --- | --- | --- |
| eam `test-units.R` | 6 pass, 23 fail, 1 error | **5 pass, 24 fail, 1 error** |
| learn `test-counterfactual.R` | 8 pass, 6 fail | 8 pass, 6 fail |
| coupling `test-cross-spectrum.R` | 60 pass, 9 errors | **60 pass, 10 errors** |

All three fail first, which is the claim that matters and the one the
lane rules demand. The two mismatches are consistent with the lane's own
note in section 6 that eam and coupling were edited after their before
runs. The counts should be re-taken so that the record matches the file
a reader will run.

### F10. LOW. A gapped record can refuse at low `segments` and pass at high

The segment length comes from the GLOBAL usable count, so on a record
broken into several spans a small `segments` can produce a length no
span can hold. Measured on 2500 samples with three spans of 700, 1043
and 699 (`rvsm-p4e-gaps.R`):

| `segments` | result |
| --- | --- |
| 2 | REFUSED, "this record supplies 0 segment(s) of 1221 samples" |
| 4 | n = 3 |
| 8 | n = 7 |
| 16 | n = 14 |

The refusal's advice, "Raise one of the three", is correct, and the help
page does say "a span shorter than one segment supplies nothing at all".
What neither says is that asking for FEWER segments can refuse a record
that MORE segments accepts. One sentence in "Gaps in the record" would
close it.

This is a documentation point, not a defect. The contract itself is
exact: see below.

## Verifications in full

### The gap contract is exactly what the help page says

I reimplemented the documented rule from `?frm_cross_spectrum` alone
(usable where both signals have a sample; `seglen = usable %/%
segments`; each span supplies as many whole segments as fit, in record
order, capped at `segments`; mean removed per segment; accumulate and
divide by `seglen`) and compared it to the package over 7 gap layouts by
4 segment counts, including layouts whose spans hold no whole number of
segments (`rvsm-p4e-gaps.R`).

**Maximum relative difference 0.00e+00 in `w11`, `w22`, `w12r` and
`w12i`, and `n` equal, in all 23 cases that were not refused.** The help
page is sufficient to rebuild the function, which is the strongest
statement a contract section can earn.

Backward compatibility is bitwise. On an ungapped record, `identical()`
is TRUE between the lane and the base at
`segments/tapers/smooth` of 8/1/1, 4/2/1, 2/1/4, 1/8/1, 13/1/1 and
5/3/1 (`rvsm-p4d-compat.R`).

### The degrees-of-freedom estimator can see a loss

The brief asked whether `1 / mean(coherence)` at a true coherence of
zero is an instrument that could report a cost if there were one.
Measured over all retained ordinates rather than one bin per replicate,
600 replicates at N = 512 (`rvsm-p4-spectral.R`):

| configuration | nominal | measured (se) |
| --- | --- | --- |
| segments = 2 | 2 | 2.000 (0.004) |
| segments = 4 | 4 | 3.971 (0.016) |
| segments = 8 | 8 | 7.992 (0.051) |
| segments = 16 | 16 | 16.130 (0.158) |
| segments = 8, hann | 8 | 7.924 (0.051) |
| segments = 2, tapers = 4 | 8 | 8.107 (0.026) |
| segments = 1, tapers = 8 | 8 | 7.999 (0.018) |
| segments = 1, tapers = 16 | 16 | 16.106 (0.048) |
| segments = 1, tapers = 32 | 32 | 33.282 (0.103) |

It tracks the nominal count over a factor of eight, and on the same code
it reports 5.750, 5.022 and 4.667 for the three Hann-plus-smooth cells
against a nominal 8, which is a 30 to 40 percent loss. So it is
sensitive, and "Hann alone costs nothing" is a result rather than a
blind spot.

The Hann-plus-smooth reproduction, which is what the refusal rests on
(1200 replicates, N = 512, all ordinates):

| configuration | nominal | none | hann | lane's hann |
| --- | --- | --- | --- | --- |
| segments = 8, smooth = 1 | 8 | 8.027 (0.037) | 7.989 (0.037) | 8.021 |
| segments = 4, smooth = 2 | 8 | 7.972 (0.036) | 5.750 (0.024) | 5.660 |
| segments = 2, smooth = 4 | 8 | 8.049 (0.037) | 5.022 (0.021) | 4.983 |
| segments = 1, smooth = 8 | 8 | 7.990 (0.037) | 4.667 (0.019) | 4.678 |

The refusal is right and I would not trade it for a declared reduced
`n`: there is no single factor, the three cells want 0.72, 0.63 and
0.58, and a user who wants a window and a wider band can raise
`segments`, which is what the message says.

One note on instrument quality. The lane reads one frequency bin per
replicate (`k <- length(w11) %/% 4`), which is the sampling
`?frm_cross_spectrum` already warns about further up the same page for
a different claim. Its standard errors are 0.067 to 0.132 where mine,
pooling every retained ordinate, are 0.019 to 0.037: about three and a
half times wider for the same 3000 against 1200 replicates. Its 8.188 is
1.5 of its own standard errors from 8 and so is not a discrepancy, and
none of its conclusions move. But a claim of the form "this costs
nothing" is worth measuring on the tighter estimator, and on the tighter
one the Hann arm sits 0.068 below the boxcar arm with a combined
standard error of 0.072.

### The guard fires once per fit

Counted on a 300-row millisecond fit (`rvsm-p3-units.R`), by class:

`frm()` 1; `predict(type = "link")` 0; `predict(newdata = )` 0;
`fitted()` 0; `simulate()` 0; `summary()` 0; `residuals()` 0;
`predict(type = "response")` 0; `logLik()` 0; `coef()` 0;
**`update()` 1**, which the lane did not test and which is correct,
because `update()` reassembles the frame and is a second fit.

### Every test file, one R process each, `NOT_CRAN=true`

Lane library, `rvsm-suite.sh`, 36 files:

| package | files | pass | fail | error | skip | warn |
| --- | --- | --- | --- | --- | --- | --- |
| frmtmb.eam | 20 | 1408 | 0 | 0 | 0 | 1 |
| frmtmb.learn | 10 | 280 | 0 | 0 | 9 | 0 |
| frmtmb.coupling | 6 | 428 | 0 | 0 | 0 | 4 |

Every per-file count matches the lane's table. The eam warning is a
pre-existing large-gradient message at `test-rdm-gng.R:733`, not the new
units warning; the learn skips are the gated Stan tier. The lane's table
heads its eam section "19 files" and lists 20.

### House style

No hand-written line in the diff or in the two new test files exceeds 80
columns. The one 100-column line is in the generated
`man/bandit2arm_delta.Rd`, which roxygen wrote. No em dashes, no en
dashes, no emojis, no spaced hyphen standing in for a dash in prose. No
British spelling added, and two removed (`modelling`, `neighbouring`).

No absolute numeric tolerance in any new assertion. Every one is
testthat's relative default, a multiple of a standard error the run
itself computes, a t-statistic threshold, or a structural bound on an
integer count. `just_over[1L] <- ceil * (1 + 8 * .Machine$double.eps)`
in `test-units.R` is relative to the constant it perturbs, which is the
right shape.

## Suspicions chased and disproved

**Shipped `n` is wrong for `tapers > 1`.** The brief asked whether sine
tapers carry the same cross-row correlation as Hann and, if so, whether
released versions ship a wrong `n`. They carry MORE of it (F4), and the
`n` is still right. Per row, at a true coherence of zero,
`segments = 1, tapers = 8` measures 7.999 (se 0.018) against a nominal
8, `segments = 2, tapers = 4` measures 8.107 (0.026), and
`tapers = 16` measures 16.106 (0.048): within about 1.5 percent at every
taper count the package documents, and inside the same band as the pure
`segments` controls. Two earlier measurements agree
(`dev/reviews/2026-09-08-xspec.md`: 7.996, 7.982, 7.982). Cross-row
correlation and per-row degrees of freedom are different quantities, and
only the second is what `n` claims. There is no HIGH finding here and no
release is wrong. What IS missing is the user-facing number, which is
F4.

**The two spellings give materially different coefficients.** The
findings' `identical(fixef) == FALSE` invited that reading. The
objective is bitwise identical as a function of the parameters and the
difference is one ulp of gradient accumulation (F7).

**The refusal could be routed around with `newdata`.** It cannot, and
that is correct behavior, but it makes the capability cost real (F6).

**`conditional_effects(method = "predict")` is a fourth route to a
draw.** It is not: core refuses it for any family without a rowwise
`sim` slot (`R/conditional-effects.R:1766`), which is every learn
family. The three routes through `sim_draw()` are `simulate()`
(`R/predict.R:2889`), `frm_simulate()` (`R/simulate-new.R:663`) and
`frmtmb.sample::posterior_predict()`
(`methods-draws.R:525`), and all three reach `fam_sim_ctx()`. The fourth
route is `frm_task_simulate()`, which the lane already guards. The gap
is the family list, not the route list (F1).

**The 0.301 lag-1 correlation does not reproduce.** My first pass
measured 0.370 at `segments = 2` and I nearly filed the record as
wrong. The lane rules say to find the construction first, and the
construction is in `small-hann-probe.R`: one segment, not eight. At one
segment I get 0.306 (se 0.002), and core's `dev/freq-findings.md`
records 0.308 on a raw periodogram, which is the same construction. The
record is right; what is wrong is that it is quoted as the property of
a frame the function will not produce (F3). The mechanism is the log
transform: at one segment an ordinate is exponential and its log has
variance `pi^2/6`, which dilutes the shared part, and the correlation
rises toward the raw 4/9 as segments are averaged.

**Gap splitting silently drops data or double-counts it.** It does
neither. An independent reimplementation written from the help page
agrees to 0.00e+00 over 22 configurations, `n` never exceeds
`segments`, and an ungapped record is bitwise what the base produced.

**The coupling example timing is a real regression.** It is not: see the
check table.

**Left open, not settled.** `rlddm()`'s drawn choices favor the LOWER
paying arm on a real schedule (0.15 rather than 0.85). Its `logp` and
`draw` share one expression so they cannot disagree with each other,
and a shared sign convention is invisible to the Stan identity test. It
is pre-existing, out of this lane's scope, and I did not chase it. One
probe in the package that owns `rlddm()` would settle it.

## `R CMD check --as-cran`

From a tarball built by `R CMD build`, against the matching private
library, with `RSTUDIO_PANDOC` and TinyTeX on `PATH`. Each package
twice on the lane, once with the CRAN-incoming remote checks on and once
with `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`, and once on the unmodified
base with them on, so the WARNING can be attributed.

| package | lane, remote on | lane, remote off | base, remote on | base, remote off |
| --- | --- | --- | --- | --- |
| frmtmb.coupling | 1 WARNING | **OK** | 1 WARNING | not run |
| frmtmb.learn | 1 WARNING | **OK** | 1 WARNING | not run |
| frmtmb.eam | 1 WARNING, 1 NOTE | 1 NOTE | see below | **1 NOTE** |

eam's base run with the remote checks ON could not complete: twice,
`checking CRAN incoming feasibility` aborted the whole check on a
60 second timeout fetching `PACKAGES.in` and `archive.rds` from CRAN.
That is the network, not the package, and the same step succeeded for
coupling and learn minutes earlier. I ran eam's base with the remote
checks OFF instead, which answers the question the column was there for:
the base raises exactly the same single V8 NOTE that the lane does, so
the lane adds nothing to eam's check.

No ERROR anywhere. The WARNING is the same in all runs and is entirely
the remote half of `checking CRAN incoming feasibility`: "Strong
dependencies not in the CRAN or BioC software repositories: frmtmb"
plus a 301 on the `https://aforren1.github.io/frmtmb/<pkg>` URL in
`DESCRIPTION` that wants a trailing slash. The lane's diff changes **no
DESCRIPTION file at all** (`git diff --stat 780dec1 -- '*DESCRIPTION'`
is empty), and the base runs raise the identical WARNING, so the lane's
attribution is right.

eam's NOTE is "Skipping checking math rendering: package 'V8'
unavailable" on the HTML manual, which is the baseline the lane rules
name.

**The coupling example timing is contention, not a regression.** The
lane reported a NOTE for coupling with the remote checks off,
"Examples with CPU (user + system) or elapsed time > 5s",
`frm_coherence` at 6.53 s, seen in one run and not another. On a quiet
machine, with the test suites finished and nothing else running, both
of my coupling runs report `* checking examples ... OK` and the
remote-off run is a plain **OK** with no NOTE at all. `R/coupling.R`
and `man/frm_coherence.Rd` are untouched by the diff. Settled.

## Consolidation

### Core seam to file

**`frmtmb_register_aterm(arity =)` accepts a range of argument counts.**
Needed by frmtmb.learn's one-column `reward(pay)` spelling (item 1.5).
File it as the lane wrote it, with two amendments from F8:

- what it buys is the SPELLING `reward(pay)`, not the capability. A
  one-column route already exists inside the extension under a second
  term name, verified by construction, so the seam is a convenience for
  users rather than an unblocking change. Rank it accordingly.
- the shape should be a minimum and a maximum with the NAMING taken from
  the maximum alone, so a term whose maximum arity is above 1 always
  delivers `name1`, `name2`, and so on. The family then declares
  `required_aterms = "reward1"` unchanged and needs no any-of group and
  no `exclusive_aterms` pair, and `registered_aterm_of()`,
  `aterm_base()` and `aterm_spelling()` are untouched. Two lines in
  `R/parse.R`, at `:276` and a `max()` at `:212`, and the open question
  the lane left ("how does a family declare which arity it got")
  disappears rather than being answered.

Everything the lane says about the parser's equality check is confirmed
by construction and nothing under core's `R/` was touched.

### NEWS bullets the lane owes and did not write

1. **frmtmb.coupling.** `tapers > 1` correlates adjacent frequencies,
   with the numbers: lag-1 correlation 0.53 at `tapers = 2` and 0.75 at
   `tapers = 4`, against 0.40 for Hann and -0.02 untapered, so a
   `tapers = 4` frame carries about one independent frequency in four.
   The lane lists this as defect 1 in its own section 5 and documents
   the general fact in the help page without a number; a released
   package's users have been fitting `s(freq)` on tapered frames since
   0.1.0 and nothing has told them.
2. **frmtmb.learn**, after F1 is fixed. The existing bullet says
   "`prl_fictitious()` is untouched: its counterfactual update READS
   both columns". That sentence is false and it must not ship: the
   update reads the chosen option's payoff and flips its sign. The
   corrected bullet should say that all eight families carry the column
   without reading it, that a duplicated column leaves every one of
   their fits bitwise unchanged, and that `frm_task_simulate()` is the
   route `rlddm()` and `ts_par7()` reach because `simulate()` already
   refuses them. The same sentence appears in
   `?bandit2arm_delta`-adjacent help text, in `ln_check_counterfactual()`'s
   own `@noRd` block and in the `prl_fictitious` compatibility note in
   `R/zzz.R`, and all four copies need the correction.
3. **frmtmb.learn**, a released-table correction of the kind the lane
   already lists as its own defect 3. The `prl_fictitious` / `reward()`
   compatibility row has said since 0.2.0 that "the second column
   enters the likelihood and not just the simulator". It does not: the
   log-likelihood is bitwise unchanged when the unchosen entry is
   replaced by noise, and `?prl_fictitious` itself says the
   counterfactual outcome is the negative of the realized one "rather
   than a reading of the second `reward()` column". Correct the row and
   say in NEWS that it was wrong, because a user choosing a family on
   that table would have chosen this one for a property it does not
   have.

### Items for later phases, unchanged by this review

- Item 2.6 still owns the coverage measurement for a tapered or
  windowed frame. This review sets its scale: the variance inflation
  factor is 1.70 for Hann at the default and 4.08 for `tapers = 4`, so
  an `s(freq)` band on a `tapers = 4` frame is the one to measure first.
- Item 3.5 still owns the contaminated record. F5 records that the
  seconds reading of such a record is itself silent, which raises its
  rank rather than lowering it.
- Item 3.2 should be written on `cp_spans()` and `cp_blocks()`, as the
  lane says. Both are correct and their contract is exactly documented.

---

# Re-check, 2026-09-08 (second pass)

The lane closed both must-fixes and rewrote its findings file. It did
not close F1 the way this review or the coordinator specified: instead
of declaring the guard on all eight families it DERIVED it, from the
`aterms` a family already names, restricted to a new `ln_schedule_terms`
set of `reward` and `payoff`, with `ln_family(counterfactual = FALSE)`
as an opt-out no family takes. The five per-family declarations are gone
and four family files have left the diff.

Libraries rebuilt from the new worktree into `rvsm-lib`; `rvsm-main-lib`
is still the unmodified 780dec1. Scripts `rvsm-r2-*`.

## Verdict on the re-check

**Mergeable after one small fix, and one test worth adding.** The
derivation is a better shape than a declaration and it holds for
everything that exists today. Both of its failure directions were tested
by construction, and one of them has a hole:

- **R1 (MEDIUM), the opt-out is unvalidated and `counterfactual = TRUE`
  silently turns the guard OFF.** Four lines fix it.
- **R2 (MEDIUM, a test rather than a bug), a schedule under a name that
  is not `reward` or `payoff` is silently unguarded**, and the planned
  item 5.1 is the family that hits it. The fix is a test that fails when
  a new multi-column addition term is registered without being
  classified.

Everything else re-measured and held. Every number the coordinator asked
me to verify reproduces to the last digit.

## What the re-check confirmed

### The derivation covers all eight, and excludes `stage2()`

Read off the eight built family objects (`rvsm-r2-derive.R`):

| family | guarded columns |
| --- | --- |
| bandit2arm_delta, bandit2arm_dual | reward1, reward2 |
| prl_fictitious, rlddm | reward1, reward2 |
| bandit4arm2_kalman_filter, ts_par7 | payoff1 to payoff4 |
| igt_pvl_delta, igt_orl | payoff1 to payoff4 |

Eight of eight, and `ts_par7()` picks up its four `payoff` columns and
not its two `stage2()` ones, exactly as the lane reports.

**The "no one has to remember" claim is true for the names in the set.**
A family built inside the review, naming `reward1`/`reward2` and passing
no `counterfactual` argument at all, has the guard derived onto it and
refuses a duplicated schedule. Nobody edited anything.

### The reproduction, at the lane's own construction

Found the construction in `small-f1-eight.R`: seed 11, 30 subjects by
100 trials, choices drawn from `prl_fictitious()` itself, and
`set.seed(4242)` before the noise columns. My first attempt used a
different seed and got a different log-likelihood; the record is right
and I had not found the construction. Re-run on both libraries
(`rvsm-r2-exact.R`):

| quantity | base (unguarded) | lane (guarded) |
| --- | --- | --- |
| logLik, two-column | -937.55332137794574 | -937.55332137794574 |
| logLik, duplicated | -937.55332137794574 | -937.55332137794574 |
| logLik, unchosen replaced by `N(100, 50)` | -937.55332137794574 | -937.55332137794574 |
| `identical()` on all three | TRUE | TRUE |
| `identical(fixef)` | TRUE | TRUE |
| observed P(arm 1), trials 51-100 | 0.8633 | 0.8633 |
| `simulate()`, two-column | 0.8616 (se 0.0027) | 0.8616 (se 0.0027) |
| `simulate()`, duplicated | **0.5470 (se 0.0090)** | **REFUSED** |
| `frm_task_simulate()`, two-column | 0.8367 (se 0.0037) | 0.8367 (se 0.0037) |
| `frm_task_simulate()`, duplicated | **0.5072 (se 0.0083)** | **REFUSED** |
| `rlddm()`, two-column | 0.1443 (se 0.0036) | 0.1443 (se 0.0036) |
| `rlddm()`, duplicated | **0.5040 (se 0.0071)** | **REFUSED** |
| `ts_par7()`, duplicated | **DREW 800 rows** | **REFUSED** |

Every figure the coordinator listed, to the digit, on the base; and
every one of the four silent draws is closed on the lane.

### The units false-alarm line is drawn honestly

The lane measures a SAFE grid at boundary 2 to 5 and an EXPOSED grid at
boundary 14 to 40, with nothing in between, which is the shape a
cherry-picked boundary would also have. It is not one. I reproduced the
safe grid and then filled the gap with no break (`rvsm-r2-units.R`, 200
draws a cell, 20 trials, drift 0.18, ndt 1.5):

| boundary | median response | fastest response | false-alarm rate |
| --- | --- | --- | --- |
| 5 | 5.9 s | 2.6 s | 0.000 |
| 8 | 11.6 s | 4.1 s | 0.000 |
| 12 | 21.7 s | 7.0 s | 0.000 |
| 14 | 27.2 s | 8.7 s | 0.000 |
| 16 | 34.7 s | 11.0 s | 0.010 |
| 18 | 38.4 s | 12.9 s | 0.035 |
| 20 | 44.2 s | 15.7 s | 0.135 |
| 26 | 61.5 s | 22.6 s | 0.645 |
| 33 | 79.3 s | 31.1 s | 0.935 |
| 40 | 97.0 s | 39.8 s | 0.995 |

The rise is continuous, there is no cliff hidden in the gap, and the
first non-zero cell is at a 34.7 second median. The lane's safe grid,
reproduced at 192 cells, tops out at a median of 10.47 seconds with a
fastest response of 6.197 seconds and **0 of 192** false alarms. So the
two grids do not straddle a concealed failure: the region between them
is genuinely 0.000 to 0.010, and the line is in the right place.

### The other confirmations

**The eleven entry points.** `frm()` 1, `update()` 1, and
`predict(type = "link")`, `predict(newdata = )`, `fitted()`,
`simulate()`, `summary()`, `residuals()`,
`predict(type = "response")`, `logLik()` and `coef()` zero each. Exactly
as the lane now records.

**F5, withdrawn and replaced by a second silent wrong answer, holds.**
On my own contaminated record (12000 rows, 5 percent uniform, truth
`bs = 1.4`, `rvsm-r2-units2.R` and `rvsm-r2-warn.R`):

| reading | fastest | units warning | `diagnose()` | boundary |
| --- | --- | --- | --- | --- |
| seconds | 0.0167 | none | conv 0, max grad 0.00116, pdHess TRUE | 2.087, 49 percent out |
| milliseconds | 16.65 | none, 16.65 is under the ceiling | conv 0, max grad 0.00113, pdHess TRUE | 65.98, 4613 percent out |

One warning is raised in each, and it is the generic "Large maximum
absolute gradient at the optimum (0.00116)" that fires identically on
both readings and says nothing about either problem. `diagnose()`
reports no problem in either. So the lane is right that the seconds
reading is a second silent wrong answer rather than a visible failure,
and right that it raises item 3.5 rather than excusing the miss. It
belongs in the backlog: file it under 3.5 as "a 5 percent contaminant
costs the boundary 47 to 49 percent with every diagnostic clean, in
EITHER unit", because that is the sentence that makes 3.5 a correctness
item rather than a modeling nicety.

**The spectral numbers, my own seeds, 400 replicates at N = 1024**
(`rvsm-r2-spec.R`). The lane's value in brackets:

| configuration | lag 1 | lag 2 | lag 3 | inflation |
| --- | --- | --- | --- | --- |
| `segments = 8` | -0.023 [-0.010] | -0.014 | -0.031 | 1.00 [1.00] |
| `segments = 8, window = "hann"` | 0.399 [0.397] | -0.003 | -0.026 | 1.80 [1.79] |
| `segments = 8, tapers = 2` | 0.529 [0.522] | 0.108 | -0.038 | 2.28 [2.25] |
| `segments = 4, tapers = 4` | 0.749 [0.753] | 0.498 | 0.268 | 4.14 [4.16] |
| `segments = 1, tapers = 8` | 0.875 [0.875] | 0.743 | 0.615 | 6.45 [6.45] |

Every cell inside Monte Carlo error. "One independent frequency in 1.8"
for Hann and "one in four" for `tapers = 4` are both right.

And the re-measured `n`, 600 replicates at N = 512, all retained
ordinates:

| configuration | nominal | mine | lane |
| --- | --- | --- | --- |
| `segments = 8` | 8 | 8.029 (0.052) | 7.968 (0.051) |
| `segments = 8, window = "hann"` | 8 | 7.993 (0.052) | 7.964 (0.051) |
| `segments = 2, tapers = 4` | 8 | 8.010 (0.026) | 8.033 (0.026) |
| `segments = 1, tapers = 8` | 8 | 8.043 (0.018) | 8.006 (0.018) |
| `segments = 1, tapers = 16` | 16 | 15.916 (0.039) | 16.144 (0.039) |

`n` is right in every configuration, on two independent runs. The
suspicion that a released `n` is wrong for `tapers > 1` stays disproved,
and the lane is right that no correctness bullet is owed.

**F9 confirmed.** My independent before-run against the unmodified base
gave eam `test-units.R` 5 pass / 24 fail / 1 error and coupling
`test-cross-spectrum.R` 60 pass / 10 errors, which is what the lane now
records in its section 6. Both mismatches were the lane's, as it says.
One stale copy survives: `dev/smallitems-findings.md` section 3 still
says "60 pass, 9 errors" where section 6 says 10.

## R3. The recurring coupling example NOTE, settled

Neither `frm_coherence` nor `cross_wishart` is touched by the diff, and
the NOTE has now appeared on each of them and on neither. Timing the
example code itself, extracted with `tools::Rd2ex` and sourced, one
fresh `Rscript` per measurement, on a quiet machine
(`rvsm-r2-extime2.R`):

| example | first execution, 6 fresh processes |
| --- | --- |
| `frm_coherence` | 2.13, 2.01, 1.86, 1.92, 1.86, 1.88 |
| `cross_wishart` | 2.13, 1.80, 1.66, 1.93, 1.59, 1.67 |

A SECOND execution of the same example in the same process costs 0.03 s,
so the whole of that 1.6 to 2.1 s is one-off warm-up rather than
arithmetic: both examples fit a `cross_wishart()` model on 4096 samples
cut into 16 segments, and the tape build and the byte-compilation are
paid once. Reversing the order changed nothing, so it is not "the first
example pays for all of them"; each pays its own.

So the answer is neither pure contention nor a genuine borderline:

- The intrinsic cost is **1.59 to 2.13 s**, a factor of 2.4 to 3 under
  the 5 s threshold. The example is not near the line on its own.
- Under load the same work reached 5.91 s in one of my own runs and
  5.39 to 6.53 s in the lane's checks. That is a 3x inflation from
  contention on top of a 2 s floor.
- **It will therefore fire intermittently on a shared CI runner, and on
  whichever example is running when the runner is busiest**, which is
  exactly the pattern observed. Six `R CMD check` runs of coupling on
  this quiet machine, across both rounds of this review, produced the
  NOTE zero times.

Not a blocker and not this lane's: neither example is in the diff. If it
keeps firing, the cheap fix is to shrink the example, which at 4096
samples and 16 segments does more work than it needs to show a coherence
fit.

## Findings

### R1. MEDIUM. The affirmative disables the guard

`ln_family()` reads the argument as

```r
counterfactual <- if (isFALSE(counterfactual)) NULL
  else if (is.null(counterfactual)) ln_counterfactual_of(aterms)
  else counterfactual
```

and `ln_check_counterfactual()` then iterates `for (cols in groups)` and
`next`s on anything shorter than two. So any value that is neither
`NULL` nor a list of character vectors disables the guard silently.
Measured on a probe family that names `reward1`/`reward2`
(`rvsm-r2-derive.R`), each fed a duplicated schedule:

| `counterfactual =` | derived groups | duplicated schedule |
| --- | --- | --- |
| (absent) | reward1, reward2 | REFUSED |
| `FALSE` | NULL | DREW, the intended opt-out |
| `NA` | NA | DREW |
| `0` | 0 | DREW |
| `"no"` | no | DREW |
| `list()` | empty | DREW |
| **`TRUE`** | TRUE | **DREW** |

`counterfactual = TRUE` is the spelling a maintainer would reach for to
mean "yes, guard this one", and it turns the guard off. That is a silent
inversion of the argument's meaning, in the one argument the whole
derivation hangs on: the lane's justification is that a ninth family
cannot be added "without the question being asked", and the affirmative
answer to that question currently means no.

`ln_family()` is internal, so only a package author reaches it, and no
shipped family passes anything. This is a hazard rather than a live
defect, which is why it is MEDIUM. It is also four lines:

```r
if (!is.null(counterfactual) && !isFALSE(counterfactual) &&
      !(is.list(counterfactual) &&
        all(vapply(counterfactual, is.character, TRUE)))) {
  stop(nm, "(): ln_family(counterfactual =) takes NULL to derive the ",
       "guard from `aterms`, FALSE to opt out, or a list of column ",
       "groups", call. = FALSE)
}
```

### R2. MEDIUM. A schedule under a third name is unguarded

This is the axis a derived rule is weakest on, and the answer is that
the predicate is right for every name that exists and has nothing to say
about a name that does not.

Constructed (`rvsm-r2-derive.R`): the same two-armed delta learner,
identical in every respect except that its schedule arrives through an
addition term registered as `outcome()` rather than `reward()`.
`ln_counterfactual_of()` returns NULL and the guard never runs:

| `outcome()` family, 10 subjects by 60 trials | P(arm 1), late trials |
| --- | --- |
| real two-column schedule | 0.7533 |
| `pay2 <- pay1` | **0.4300, drew 600 rows, no error** |

The same holds for a family carrying its schedule through `vreal()`,
which is the route core's own `?frmtmb_register_aterm` names for a
family that does not register a term of its own: NULL, drew 600 rows.

This is not hypothetical. Plan item 5.1 is `bandit_delta(n_option = K)`
covering `banditNarm_*`. An addition term's arity is fixed at
registration, which is this lane's own core seam finding, so a K-armed
family for K outside 2 and 4 must register a new term name. If that name
is not literally `reward` or `payoff`, the guard stops applying and
nothing says so. A declaration that someone forgets at least sits in the
family's source next to `aterms =`; a derivation that returns NULL for
an unrecognized name is invisible.

One accident cuts the other way and is worth knowing: a term registered
as `reward3` at arity 3 produces the keys `reward31`, `reward32`,
`reward33`, which strip to `reward` and ARE picked up. That is luck
rather than design.

**The fix is a test, not a code change.** The package already keeps its
own registration table, `ln_aterms <- c(reward = 2L, payoff = 4L,
stage2 = 2L)` in `R/zzz.R`. A test that every entry with arity 2 or more
is classified, either in `ln_schedule_terms` or in an explicit
`ln_not_schedule_terms` naming `stage2` and saying why, turns the silent
miss into a red suite the moment a fourth term is registered. That is
the property the derivation was chosen for, and it is the only thing
that makes it stick. A second assertion, that every exported family
whose `aterms` hold a repeated multi-column term has a non-NULL
`counterfactual`, closes the same gap from the other side.

Both belong in `test-counterfactual.R`, which currently pins the
behavior of the guard on the eight families that exist and nothing about
the rule that produces it.

### R3. LOW. The `tapers = 8` inflation factor is truncated at four lags

The shipped help table gives 6.45 for `segments = 1, tapers = 8`. That
is `1 + 2 *` the sum of the positive correlations over lags 1 to 4.
Correlations there are still positive past lag 4, and over lags 1 to 8
the factor is **8.02** (`rvsm-r2-spec.R`). For every other row the two
windows agree to 0.01, because the correlation has died by lag 3; only
the `tapers = 8` row moves. Either say the sum is over four lags or
quote 8.02. It changes no conclusion, and it makes the worst row look
better than it is.

### R4. LOW. The false-alarm rate is not a function of the median response alone

The lane's false-alarm table is indexed by median response, and along
its own sweep the numbers are right: I reproduce 0.645 at a 61.5 second
median where it reports 0.560 to 0.710. But the median does not
determine the rate. Holding it near 50 seconds and reaching it along
three different paths (`rvsm-r2-units.R`, 200 draws a cell, 20 trials):

| drift | boundary | median response | fastest response | false-alarm rate |
| --- | --- | --- | --- | --- |
| 0.10 | 18.7 | 53.8 s | 14.7 s | 0.165 |
| 0.18 | 21.5 | 48.8 s | 16.9 s | 0.255 |
| 0.35 | 37.6 | 52.1 s | 26.7 s | **0.940** |

A task that is slow because the boundary is far and the evidence is
strong has a tight response time distribution, so the fastest of twenty
trials sits close to the median and the guard fires; a task that is slow
because the evidence is weak has a long right tail and a fast minimum.
At one median the rate spans 0.165 to 0.940.

The lane's cells are therefore the LOW end of the range for their
median, not the range. The honest one-line summary is the one the guard
actually implements: it fires when the FASTEST response passes 20
seconds, and how far the fastest sits below the median depends on the
response time spread as much as on its center. Adding the fastest
response as a column to the existing table would say that without
another sweep, and the numbers are already in it.

## Test suites and checks, re-run

Lane library, one R process per file, `NOT_CRAN=true` (`rvsm-suite.sh`,
36 files):

| package | files | pass | fail | error | skip | warn |
| --- | --- | --- | --- | --- | --- | --- |
| frmtmb.eam | 20 | 1408 | 0 | 0 | 0 | 1 |
| frmtmb.learn | 10 | 298 | 0 | 0 | 9 | 0 |
| frmtmb.coupling | 6 | 428 | 0 | 0 | 0 | 4 |

`test-counterfactual.R` is 32 assertions, up from 14, and covers all
eight families. The eam warning is still the pre-existing large-gradient
message at `test-rdm-gng.R:733`.

`R CMD check --as-cran`, from a tarball, on a quiet machine:

| package | remote on | remote off |
| --- | --- | --- |
| frmtmb.coupling | 1 WARNING | **OK** |
| frmtmb.learn | 1 WARNING | **OK** |
| frmtmb.eam | 1 WARNING, 1 NOTE | 1 NOTE |

Unchanged from the first pass, and no example NOTE in any of the six
runs. The WARNING is still entirely the remote CRAN-incoming half and no
`DESCRIPTION` is touched by the diff; eam's NOTE is the expected V8 one.

## Consolidation, updated

**The core seam to file is unchanged in substance**, but the derivation
adds a constraint the filing should carry. If `frmtmb_register_aterm()`
gains a range arity and `reward(pay)` becomes legal, the one-column
spelling arrives as `aterms[["reward"]]` under the lane's proposed
naming while the two-column one arrives as `reward1` and `reward2`.
`ln_counterfactual_of()` would then group `c("reward", "reward1")`
together, `ln_check_counterfactual()` would find `cd[["reward1"]]`
absent and skip, and **the one-column spelling would be unguarded**,
which is precisely the spelling that cannot supply a counterfactual and
most needs the refusal. The minimum-and-maximum shape this review
proposed avoids it: a term whose maximum arity is above 1 always
delivers `name1`, `name2`, so the one-column spelling arrives as
`reward1` alone, the group has one member, and the guard refuses the
draw for a different and equally correct reason. Note that in the
filing.

**NEWS bullets.** Both that this review said the lane owed are now
written: the `tapers` correlation bullet in `frmtmb.coupling` and the
corrected `prl_fictitious` compatibility row in `frmtmb.learn`, which
now says the row was wrong from 0.1.0 and what it should have said. The
false claim is gone from all four copies. Nothing further is owed.

**One stale number** in `dev/smallitems-findings.md`: section 3 still
reports the coupling before-count as "60 pass, 9 errors" where section 6
correctly has 10.

## What I changed in the worktree

This file, and nothing else. R1 and R2 are the lane's to apply: R1
changes behavior on inputs no family passes today, and R2 is a test
whose content is a judgement about which terms are schedules.

---

# Final pass, 2026-09-08 (third)

Both re-check findings are closed. **Mergeable.**

## R1, and a hole in the fix I supplied

The lane is right, and it is worse than one hole. Running my own patch
text as written (`rvsm-r3-probe.R`), against the shapes the lane now
refuses:

| value | my patch | groups actually checked |
| --- | --- | --- |
| `TRUE`, `NA`, `0`, `"no"`, `list(1:2)` | refused | -- |
| **`list()`** | **passed** | **0** |
| **`list("reward1")`** | **passed** | **0** |
| **`list(character(0))`** | **passed** | **0** |
| `list(c("reward1", "reward2"))` | passed | 1 |

Three shapes, not one, and all three leave the guard silently off. Two
causes, both mine: `all()` over an empty list is TRUE, which the lane
found; and my predicate tested `is.character` and never the LENGTH, so a
list of one-element vectors passed and then met the `length(cols) < 2L`
skip inside the check. A reviewer's patch is not exempt from the rule
that a guard has to be measured rather than read, and this one was not.

`ln_read_counterfactual()` closes all three with `length(x) > 0L` and
`length(z) >= 2L`, and refuses at family construction rather than at
draw time, which is where a mistake in it can still be seen.

**The widening refuses nothing a real family would pass.** All eight
bad shapes stop the family being built; and:

| passed | result |
| --- | --- |
| absent / `NULL` | built, guards `reward1, reward2` |
| `FALSE` | built, guard off |
| `c("reward1", "reward2")` | built, guards `reward1, reward2` |
| `list(c("reward1", "reward2"))` | built, same |
| `list(reward = c("reward1", "reward2"))` | built, same, name dropped |
| `list(c("reward1","reward2"), paste0("payoff", 1:4))` | built, guards both groups |

All eight shipped families still build. The bare character vector and
the named list are both NEW acceptances relative to my patch, which
would have refused a bare character vector outright, so the lane's
version is more permissive where a family is right and stricter where it
is wrong. Refusing `TRUE` by name, with its own sentence, is the right
call for the reason given: deriving is already the default, so there is
nothing left for an opt-in to mean.

## R2, and whether the test bites

It bites, and it can also be satisfied, which is the part worth checking
because a tripwire that can only fail is not a classification test.
Measured by adding a fourth multi-column term to `ln_aterms` in the
loaded namespace and rerunning the file (`rvsm-r3-bite.R`), using a name
the test file does not itself mention:

| state | test-counterfactual.R |
| --- | --- |
| as shipped | 67 pass, 0 fail |
| `bonus = 3L` registered, classified nowhere | **66 pass, 1 fail** |
| `bonus` added to `ln_schedule_terms` | 67 pass, 0 fail |
| `bonus` added to `ln_not_schedule_terms` with a reason | 67 pass, 0 fail |

So the assertion is on the tables and not on today's agreement between
them, and either classification clears it. (My first probe used
`outcome`, which the file separately pins as an unrecognized name, so
classifying it moved the failure instead of clearing it. That is the
test doing two jobs correctly, not a defect.)

The before-and-after against the unmodified base reproduces:
`test-counterfactual.R` scores **28 pass, 24 fail, 5 error** there.

Closing R2 with a test rather than a predicate is the right call and it
is the one I asked for. A predicate that guessed at unrecognized names
would have to guess wrong in one direction, and the direction it would
guess wrong in is the loud one.

## The seam filing: the reasoning holds, and understates itself

Put through the shipped derivation directly (`rvsm-r3-seam.R`), with a
duplicated schedule in each case:

| naming | family names | one-column data | two-column data |
| --- | --- | --- | --- |
| today | `reward1, reward2` | -- | REFUSED |
| plain range arity | `reward, reward1, reward2` | drew | **drew** |
| min and max | `reward1, reward2` | drew | REFUSED |

Under a plain range arity the family has to name all three keys, the
derivation groups them as one group of three, and whichever spelling the
user writes leaves one key absent, so `ln_check_counterfactual()` finds
a NULL and skips. **That disables the guard for the TWO-column spelling
as well**, which works today. So the seam as first proposed would not
merely leave the new spelling unguarded; it would undo item 1.5's own
fix on the spelling every existing user writes. The min-and-max shape
keeps the two-column case refused.

One caveat the filing should carry, because it is not automatic: the
min-and-max shape does not by itself guard the one-column spelling
either. It arrives as `reward1` alone, the group has one member, and
there is nothing to compare it against. What the shape buys is that the
shortfall becomes a COUNT the family can see, keys present against
`n_option`, rather than an absent key hidden inside a mixed group.
Whoever implements `reward(pay)` owes an explicit refusal for it, and
the derivation will not supply one.

## Bookkeeping, spot-checked

R3 is corrected to 8.02 with the truncation stated in the help page
("summed over twelve lags ... stopping there gives 6.45 instead of
8.02"). R4's table now carries the fastest response beside each rate and
the three-path measurement, at 0.185 / 0.260 / 0.935 against my 0.165 /
0.255 / 0.940 on different seeds. The stale coupling before-count is 10
in both places. The example NOTE is recorded with the timing and the
third quiet run.

## Final state

| | |
| --- | --- |
| frmtmb.eam | 20 files, 1408 pass, 0 fail, 1 pre-existing warning |
| frmtmb.learn | 10 files, **333 pass**, 0 fail, 9 gated skips |
| frmtmb.coupling | 6 files, 428 pass, 0 fail |
| eam `test-units.R` | **32** |
| learn `test-counterfactual.R` | **67** (28 / 24 / 5 against the base) |
| coupling `test-cross-spectrum.R` | **96** |

`R CMD check --as-cran`, remote incoming off: frmtmb.learn **OK**,
frmtmb.coupling **OK**, frmtmb.eam **1 NOTE** (the expected V8 one). With
remote on, the same single CRAN-incoming WARNING as the base, on a
`DESCRIPTION` the diff does not touch. No example NOTE in any of the
nine coupling checks this review has now run on a quiet machine.

Worktree unchanged by me except this file. No commits.
