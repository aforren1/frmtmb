# The non-decision-time bound in `rlddm()`, and the seam it comes from

**Punch rounds 1 and 2 are applied.** Read the retraction below before
anything else: this file first reported a defect in `bf(ndt = )` that
does not exist, and the number it gave was an artifact of my own
reconstruction. `dev/reviews/2026-09-10-rlddm.md` falsified it and I
reproduced the falsification. Round 1's four blockers are closed except
`frmtmb.learn`'s dependency floor, which belongs to consolidation.

**Round 2 changed one number this file had been reporting wrongly in the
other direction**: the `learn-rlddm` row recorded `cor_true = 0`, and
the truth of the block that row estimates is about -0.78. The row now
computes it. Round 2 also documented the one hole the seam cannot close,
measured at 37.5 percent, and put the `sd(ndt)` warning into `?rlddm`
rather than only here.

Item 1.0b of `dev/extension-gaps-plan.md`, in frmtmb.learn and
frmtmb.eam. Written on the `wt-rlddm` worktree off `f9ee297`. R 4.6.1,
Windows 11, 16 logical cores. frmtmb.eam and frmtmb.learn are built
from this worktree into a private library at
`C:/Users/adf44/source/r/rlddm-lib`; frmtmb 0.55.2 and the BEFORE arm of
everything come from the round's shared reference build at
`C:/Users/adf44/source/r/rellib-0552`, which is READ-ONLY and was not
rebuilt.

Every script named below is in `dev/rlddm-scripts/` and carries its seed
in its header. `dev/rlddm-scripts/rlddm-prelude.R` fixes the
`.libPaths()` order the lane rules require, and every script takes the
library as its first argument so that the two arms differ in one thing.

## The claim, and whether it holds

**The `learn-rlddm` scale row converges with a positive definite Hessian
and no NaN standard errors.** It does.

| | 0.3.0, one global bound | this change, `ndt_group(id)` |
|---|---|---|
| convergence code | **1** | **0** |
| maximum absolute gradient | **6.36e+09** | **3.12e-03** |
| positive definite Hessian | **no** | **YES** |
| standard errors that are `NaN` | **4** | **0** |
| log-likelihood | **-8333.1** | **-7781.2** |
| parameters | 14 | 14 |
| `sd(ndt)` on the link | 2.236 | 0.3261 |
| `sd(alpha)` on the link, truth 0.5 | 0.520 | 0.5073 |
| `sd(drift)`, truth 1.0 | 1.138 | 1.094 |
| `sd(bs)` on the log link, truth 0.2 | 0.2119 | 0.1915 |
| `alpha` intercept, truth -0.6190 | -0.5506 (-0.6660, -0.4352) | -0.5186 (-0.6435, -0.3938) |
| learners below their OWN fastest response | (n/a) | 100 of 100 |
| tightest such margin | (n/a) | 20.9 ms |
| per-learner `ndt` RMSE against the drawn truths | (n/a) | 13.996 ms |
| per-learner `ndt` correlation with the truth | (n/a) | 0.9404 |
| `sd(ndt)` natural, truth 0.03930 | (n/a) | 0.03337, and read the caveat below |
| `cor_max_abs`, truth 0.7768 not 0 | 0.3964 | 0.9160 |
| whole `frm(se = TRUE)` call | 1101 s | 574 s |

Construction:
`extensions/frmtmb.learn/tests/testthat/test-scale.R`, row
`learn-rlddm`, seed 20260908, 100 learners by 200 trials, 20,000 rows,
run by `dev/rlddm-scripts/rlddm-scale.R` under
`FRMTMB_SCALE_TESTS=true NOT_CRAN=true`, one row per fresh process. The
records are `dev/rlddm-scripts/rlddm-scale-before.tsv` and
`rlddm-scale-after.tsv`; the logs are beside them.

The BEFORE row is this lane's own run against the shared reference
build, not a figure carried from elsewhere, and it reproduces what
`dev/round-handoff.md` recorded to every digit: `conv=1`,
`maxgrad=6.36e+09`, `pdHess=FALSE`, `nbadse=4`, `logLik=-8333.1`. The
AFTER row was run TWICE, in independent processes, and reproduces
itself to every digit recorded: `-7781.2`, `conv=0`, `maxgrad=0.00312`,
`pdHess=TRUE`, `nbadse=0`, population fraction 0.828719, per-learner
RMSE 13.9959 ms, tightest margin 20.8741 ms. Read the two wall clocks as
upper bounds: both ran beside other R processes, and
`dev/scale-findings.md` records this row at 738 s on a quiet machine.

**What this establishes is that the model is now ESTIMABLE, and the
log-likelihood gap is a symptom rather than a gain.** The review was
right to press on this and the first version of this section overclaimed
it. The before fit has `conv = 1` and a maximum gradient of 6.36e+09, so
-8333.1 is where an optimizer stopped on a model it could not optimize
and not that model's maximum, which is unknown and is at least that
much. Nor are the two arms the same parameter space: both take the bound
from the response, but the per-group bound admits, for 99 of the 100
learners, non-decision times the global bound cannot express at any
parameter value, and a higher likelihood on a larger admissible set is
expected rather than evidence. So the row's claim is the left column
turning into the right one: convergence code 0, maximum gradient
3.12e-03, a positive definite Hessian and 14 finite standard errors
where four were `NaN`. The 551.9 units are reported as the size of the
constraint that was binding.

The `ndt` variance component is the clearest single symptom: on the link
it comes back at 2.236 under the global bound, which is a random effect
running away on a scaled logit against a wall, and at 0.3261 under the
per-learner one.

**Two rows of that table need reading with care, and neither is
asserted.**

`sd(ndt)` natural, 0.03337 against a truth of 0.03930, is the row a
reader is most likely to take as "the variance component is recovered",
and it is the row the FLOORS do better on. This is 1.0a's B5 recurring
in a second family and it gets the caveat `?wiener` had to gain:

| estimator | `sd(ndt)` natural, truth 0.0393041 |
|---|---|
| the fit, with the random effect | 0.0333713 |
| the fit with NO random effect on `ndt` | **0.0380049** |
| 0.812475 times each learner's own floor | 0.0370422 |
| 0.828719 times each learner's own floor | 0.0377828 |

Measured by the review, `dev/rev-rlddm-floors.R` arm `nore` and
`dev/rev-rlddm-sdfloor.R`. The statistics that DO separate the model
from its floors are the per-learner error and the correlation, which are
the two the row asserts; see the next section.

**A caveat in this file was not enough, and punch round 2 said so.** A
shipped `.tsv` field with no comparison beside it is what produced
1.0a's B5 and was producing it a second time. Two things changed rather
than one:

* `test-scale.R` now records `sd_ndt_floor_only` beside
  `sd_ndt_natural`, the spread a COMMON fraction of each learner's own
  floor would give. It costs no second fit, and it puts the comparison
  in the row rather than in a document nobody reads beside it. The
  no-random-effect FIT's 0.0380049 is cited in the same comment, at the
  script that measured it;
* `?rlddm` gains the section `?wiener` gained under 1.0a, with the
  no-random-effect fit's number rather than the arithmetic bound,
  because that is the honest comparison: 0.0380049 against a truth of
  0.0393041 where the full model returns 0.0333713. On this statistic
  the random effect is not neutral, it is worse, and a user reading the
  row could not know that from the row.

`cor_max_abs` goes from 0.396405 to **0.916032**, and the first version
of this file did not mention it at all. Punch round 2 settled it, and it
is not a defect: **the row's stated truth of 0 was wrong.**

`dev/rlddm-scripts/rlddm-cortrue.R`, seed 20260908, the tier's own
design, NO FIT. The design draws four independent deviations, but the
model does not estimate the four it drew. Three it does: `alpha` on the
logit, `drift` on the identity, `bs` on the log. The fourth it does not:
under a per-group bound it estimates `qlogis(ndt_i / floor_i)`, and
`floor_i` is that learner's own fastest response, which is a draw from
its WHOLE parameter vector. Expressing the drawn truths in the
parameterization the row now fits:

| pair | as drawn | in the parameterization now fitted | fitted |
|---|---|---|---|
| bs, ndt | +0.0604 | **-0.7768** | **-0.9160** |
| drift, ndt | +0.1591 | +0.0880 | +0.1163 |
| alpha, ndt | +0.0028 | +0.0085 | +0.0458 |
| alpha, drift | -0.1158 | -0.1158 | -0.2545 |
| drift, bs | +0.0695 | +0.0695 | +0.1034 |
| alpha, bs | +0.0139 | +0.0139 | -0.0224 |
| max abs | 0.1591 | **0.7768** | **0.9160** |

The mechanism, measured in the same script:
`cor(log floor, bs deviation) = +0.4954`, because a learner with a
higher boundary separation is slower and its own fastest response is
later; and `cor(log floor, ndt deviation) = +0.8433`, the floor mostly
being the non-decision time. So a learner with a high `bs` needs a LOWER
fraction to express the same non-decision time, and the two deviations
are negatively dependent by construction, before any estimator sees the
data. My untested hypothesis was right, including which parameter.

**So the row recorded 0.916 beside a truth of 0 and read as though the
fit invented a dependence. It recovered one, and overstated it by
0.14.** `test-scale.R` now COMPUTES `cor_true` per replicate from the
drawn truths and the floors, with no fit, so it moves with the seed
instead of being pinned to a number the design does not support, and it
records which pair carries the maximum. At this seed it is 0.776847.

**This matters most one item away.** Phase 2's item 2.2 is 60
replicates of exactly this design with the correlated block on every
parameter. A diagonal truth there would fail on nearly every replicate
for a reason that has nothing to do with the estimator, and the plan row
below says so.

## The floors do not hand the fit its numbers

The hardest attack on this item is 1.0a's own finding: that the observed
per-subject floors carry the spread, so a recovery statistic can look
good with the model doing no work. The review ran it, on the tier's own
seed and design, and **it does not carry over**
(`dev/rev-rlddm-floors.R`, `dev/rev-rlddm-sdfloor.R`):

| estimator | per-learner RMSE | correlation with the truth |
|---|---|---|
| the fit, `ndt ~ 1 + (1 \| p \| id)` | **13.996 ms** | **0.9404** |
| the best constant fraction of the OWN floor | 21.075 ms | 0.8478 |
| the fit's own fraction, 0.828719, times the floor | 21.692 ms | 0.8478 |
| one number for everybody | 39.107 ms | (undefined) |
| the floor itself, unscaled | 62.921 ms | 0.8478 |
| the fit with the random effect OFF, `ndt ~ 1` | 22.107 ms | 0.8478 |

The correlation of ANY constant multiple of the floor with the truth is
`cor(floor, truth) = 0.847780`, so **no estimator without a random
effect on `ndt` can beat 21.075 ms or 0.8478 on this data.** That is an
analytic bound rather than a fit. The `nore` arm confirms it by
construction: its fitted values are a constant times the floor to
1.7e-13 ms, `cor(fitted, own floor) = 1`, and it costs 97.16
log-likelihood units at four fewer parameters. The two statistics the
row asserts are the fit's own.

The tier still WARNS on the maximum gradient at 0.00312, which is the
same warning the `eam-unbounded` row carries at 0.00293. That is
frmtmb's threshold rather than a property of this fit, and the row
records it.

## What was wrong, and what it cost

`rlddm()` took its diffusion parameterization from frmtmb.eam and, at
0.3.0, took the non-decision time's bound by writing its own copy of the
scaled logit rather than importing one. `?rlddm` said so plainly:

> `frmtmb.eam` builds the same link for the same reason, and this is the
> one thing this family reproduces rather than imports.

The copy is how the family inherited the defect item 1.0a removed in
frmtmb.eam. The bound was `min(rt)` over the WHOLE data set. With 100
learners that is the fastest response of all of them, so a subject
deviation on `ndt` is a deviation on a fraction of somebody else's
floor, and a learner whose true non-decision time is above the global
minimum cannot be expressed at all.

### The behavioral failure, seen

`dev/rlddm-scripts/rlddm-group.R`, seed 909, log at
`rlddm-group-log.txt`. Eight learners in two groups of four, 120 trials
each, true non-decision times 0.18 and 0.40 s, so the slow group's truth
is above the FAST group's fastest response and one bound cannot express
it. That is the construction `dev/ndt-scripts/ndt-before-after.R` used
for frmtmb.eam, on the family this item is about.

| | one global bound | one bound per group |
|---|---|---|
| global `min(rt)` | 0.225935102 | the same data |
| the two groups' own `min(rt)` | 0.225935102, 0.455662732 | |
| fitted `ndt`, truths 0.18 and 0.40 | 0.14740766, **0.22593510** | 0.19269606, 0.40614955 |
| slow group's relative error | **43.5 percent** | 1.5 percent |
| fast group's relative error | 18.1 percent | 7.1 percent |
| log-likelihood | -392.339337 | **-237.703404** |
| convergence code | 0 | 0 |
| maximum absolute gradient | 2.86e-07 | 2.80e-05 |
| positive definite Hessian | YES | YES |
| standard errors that are `NaN` | 0 | 0 |
| margin below the group's own floor | (n/a) | 33.2 ms, 49.5 ms |

154.6 log-likelihood units, at the same parameter count.

**The global arm converges cleanly and is wrong**, which is the part
worth carrying. Its slow group sits at 0.22593510 against a bound of
0.22593510: at the wall to nine digits, with convergence code 0, a
maximum gradient of 2.9e-07, a positive definite Hessian and no bad
standard errors. Nothing in the fit's own diagnostics says anything. At
the scale row's 100 learners the same failure shows up as an exploded
gradient instead, and the two are the same defect at two sizes; a check
that tests either alone is fooled by the other.

## What changed

### In frmtmb.eam: a second export seam

`R/extension-api.R` promised the Wiener density and said a caller
wanting more "should ask for a second seam rather than reach for a
colon". This is that second seam, in `R/ndt-seam.R`, and it is about the
BOUND rather than about any density.

| export | what it is |
|---|---|
| `ndt_bound(y, aterms, max_ndt, what)` | derives the bound from the response, `ndt_group()` and `max_ndt`, with every refusal that goes with it. Returns a `"frmtmb_eam_ndt_bound"` |
| `ndt_bound_pending(max_ndt, what)` | the bound a family carries BEFORE `frm()` has seen a response: a refusing placeholder, or the stated `max_ndt` |
| `ndt_bound_attach(fam, bound)` | puts one on a family: the `ndt` link, the per-row floor as `aterm_data`, the record, and for a grouped bound the starting fraction |
| `ndt_bound_of(fam)` | the settled bound a family already carries, or `NULL`. This is how a `family_finalize` KEEPS a bound instead of re-deriving it |
| `ndt_apply(dpars, aterms, what)` | the non-decision time on the response's own scale, wherever a density would have read `dpars$ndt`, WITH the refusal that goes with it |
| `ndt_bound_key(x)` | the coercion an `ndt_group()` column goes through, for a caller assembling `aterms` by hand |

**`ndt_apply()` is a sixth export added in punch round 1, and the review
was right to ask for it.** The seam first stopped at the bound and
documented the multiplication as three lines a consumer copies. The
arithmetic is three lines; the REFUSAL that goes with it, for a grouped
model reached without its per-row bound, is eleven, and it is the half
that stops a fraction being read as seconds. So the seam was shipping
the safe part as code and the dangerous part as prose, and the review
constructed the consequence: the documented three lines, run with the
floor absent, return 0.5 as a time, silently. `rlddm()`'s `ln_ndt_at()`
is now one line that calls it, and a consumer's copy is zero lines.

The boundary is otherwise where it was, and the review verified the
reason by reading both: `ddm_ndt_install()` stays `@noRd` because it
rewrites five slots BY NAME and `rlddm()`'s likelihood is in none of
them. `ndt_apply()` takes a list in and returns a number, which assumes
nothing about family layout, which is the whole objection to exporting
the installer.

**What is promised.** The bound. How it is derived, the refusals that
derivation owes, the link that makes the constraint structural, the
per-row floor arriving in the addition-term values under the reserved
name `ndt_floor`, and the record a fitted family carries at
`family(fit)$ndt_bound`, which `ndt_time()` and the unread-grouping
frame check both read.

**What is NOT promised.** The four densities' internals. The slot
WRAPPING this package does for its own families: `ddm_ndt_install()`
rewrites `lpdf`, `lcdf`, `lccdf`, `sim` and `post$mean_fn`, which
assumes this package's own family layout, and it stays `@noRd`. The
across-trial range `st` and its doubled bound. The label hash the floor
table is keyed on. The starting value of anything but `ndt`. And
`ndt_apply()` returns a VALUE rather than a modified `dpars`, so a
consumer using frmtmb's robust dpar accessors must drop `.eta_ndt`
itself; this package's own internal wrapper does that for its four
families and the seam says so.

### Two holes in the seam that punch round 1 found

**`ndt_bound_attach()` accepted `gddm()` and broke it.** The guard was
`!is.null(fam[["ndt_raw"]])`, the marker `ddm_ndt_install()` leaves, and
only FOUR of this package's five families carry it. `gddm()` has an
`ndt`, has no `ndt_raw`, and its density reads no `ndt_floor` at all, so
it passed. `dev/rlddm-scripts/rlddm-guard.R`, no seed, exhaustive over
the five:

| family | `ndt_raw` | attach, before the fix | after |
|---|---|---|---|
| `wiener()` | TRUE | refused by name | refused by name |
| `lba(2)` | TRUE | refused by name | refused by name |
| `rdm(2)` | TRUE | refused by name | refused by name |
| `wiener_gng()` | TRUE | refused by name | refused by name |
| `gddm()` | **FALSE** | **ACCEPTED** | **refused by name** |

What the accepted call did to `gddm()`: its `ndt` link went from
`scaled_logit(0, from data)` to plain `logit`, an `ndt_floor` its
density never reads was planted in `aterm_data`, and its `ndt` starting
value moved from 0.155 s to 0.5. That is the silent-wrong-answer shape
this whole item exists to remove, reachable through the item's own new
export, in the item's own package. **My guard tested a slot that
happened to correlate with the property instead of the property**, which
is the same mistake in guard design the round has now seen several
times.

The fixed guard is `!is.null(fam[["ndt_raw"]])` OR the family's name
being in `ddm_accepts`, the table this package's five families are
already enumerated in and that the compatibility rows already derive
from. A sixth family added there is refused without anyone remembering.
`?ndt_bound_attach` now names `gddm()` and says why it is excluded, and
the same script shows a family from ANOTHER package is still accepted.


**The same hole one package further out, DOCUMENTED and not fixed.** A
foreign family is accepted, which is the seam's only reason to exist,
and nothing on this side can check that its density calls
`ndt_apply()`. Punch round 2's review built the silent case: two
`custom_family()` shifted gammas differing in one line, one calling
`ndt_apply()` and one reading `dpars[["ndt"]]` directly, four groups of
300 rows with true shifts 0.60 to 0.90 (`dev/rev-rlddm-foreign2.R`,
seed 707).

| | calls `ndt_apply()` | reads the dpar directly |
|---|---|---|
| attach, fit | accepted, returned | accepted, returned |
| convergence code | 0 | 0 |
| positive definite Hessian | TRUE | TRUE |
| non-finite standard errors | 0 | 0 |
| log-likelihood | 1345.58425 | **1345.58425** |
| worst `ndt_time()` error | 1.3 percent | **37.5 percent** |
| every fitted time below its own floor | TRUE | **TRUE** |

Identical log-likelihood, clean diagnostics, nothing refused, and the
non-decision time up to 37.5 percent too small. It also confirms round
1's point about the floor check: it is TRUE on the broken arm, because
it tests the floor LOOKUP and not the arithmetic.

**Decision: documented, not guarded.** No package can verify another
package's arithmetic, and the review ranks it a nit for that reason. The
review's suggested mitigation is an explicit argument the consumer must
set to declare that its density calls `ndt_apply()`; this lane declines
it, because 1.0a already measured what an affirmative flag gating a
guard does (`ln_family(counterfactual = TRUE)`, the spelling a
maintainer reaches for to mean "yes, guard this one", turned the guard
OFF), and a consumer careless enough to skip a one-line call is careless
enough to set a flag. What the number buys instead is visibility, so it
is now in `?ndt_apply` and pointed at from `?ndt_bound_attach`, together
with the finding that the 0.5 starting fraction is a partial guard only
by luck: on a faster design the broken family died at once with
`NA/NaN gradient evaluation`, and on a slower one it converged silently.

**`ndt_bound()` validated almost nothing, and disagreed with its own
sibling.** `?ndt_bound` invites a caller to assemble `aterms` by hand
and call it outside `frm()`, where no `valid_y()` runs in front of it.
The review's ten probes: `ndt_bound("a")` came back with `ub = "a"`,
`max_ndt = 0` and `max_ndt = -1` were accepted, and `ndt_bound(c(0.3,
NA))` came back with `ub = NA_real_`, which is the exact sentinel
`ddm_scaled_logit()` reads as PENDING, so one missing response would
produce "the bound is not set yet" from a family whose bound WAS set.
Meanwhile `ndt_bound_pending(-1)` refused by name: two entry points,
two contracts, one argument, one of them exported as API.

Both now go through one `max_ndt` check, and `y` must be numeric,
non-empty, finite and positive. `test-ndt-seam.R` pins all eleven cases.

### Why the seam is not the four functions the plan named

The plan named `ddm_ndt_spec()`, `ddm_ndt_scaler()`,
`ddm_ndt_install()` and `ddm_ndt_preinstall()`. Three of the four are
reachable through the seam under other names; `ddm_ndt_install()` is
NOT, and exporting it would have been the wrong promise.

`rlddm()`'s likelihood is not `lpdf`. Its trial contribution depends on
every earlier trial of the same subject, so the family declares
`lpdf = ln_no_rowwise()` and the real likelihood lives in the structured
protocol's `loglik`, `loglik_row` and `loglik_group` slots, which take
`(y, dpars, aterms, weights, block, extra)`. `ddm_ndt_install()` wraps
`(y, dpars, aterms, extra)` slots by name. Exporting it would have
promised frmtmb.eam's family LAYOUT to a package that does not share it,
and it would not have worked for the one consumer that asked for it.

So the seam stops at the bound and the consumer's own density does the
multiplication, which for `rlddm()` is three lines:

```r
ndt <- d[["ndt"]]
fl <- d[["ndt_floor"]]
if (!is.null(fl)) ndt <- ndt * fl
```

That is `ln_ndt_at()` in `R/rlddm.R`, plus a refusal for the case below.

### The reviewer's fifth entry

The 1.0a reviewer asked for "a supported way to obtain a scoreable
family with a stated bound", because `test-extension-api.R` had had to
change to `wiener(max_ndt = 1)`.

**Half of that premise no longer holds, measured.**
`dev/rlddm-scripts/rlddm-smoke.R` scores a BARE `wiener()` outside
`frm()` on the reference build of the base commit:
`wiener()$lpdf(c(0.4, 0.9), list(mu = 1.2, bs = 1.5, ndt = 0.15,
bias = 0.3), list(dec = c(1, 0)))` returns `0.075869859909` and
`-2.586299817283`. It is only the LINK that refuses:
`wiener()$links$ndt$linkinv(0)` errors with "bound is not set yet". The
awkwardness the reviewer saw was round one's, when the bound lived in
the density, and round two removed it: `test-extension-api.R` is back at
plain `wiener()` and passes with 14 assertions.

**The other half is real, and it is the fifth entry.**
`ndt_bound_pending()` is `ddm_ndt_preinstall()` exported under a name
that says what it is for. A family gets one of two honest states before
`frm()`: PENDING, whose link refuses by name, or a stated `max_ndt`,
whose family scores and whose `ndt` is a time.

**The justification I first gave for it was a defect that does not
exist**, and the retraction below replaces it with the two that do: the
out-of-range `bf(ndt = )` case, which reached the objective as `NaN`
and surfaced as an optimizer failure naming nothing, and one contract
across all five families rather than four and an exception.

For frmtmb.eam's own five families the supported spelling stays
`wiener(max_ndt = )`, and `ndt_bound_attach()` refuses all five BY NAME
and says so.

### In frmtmb.learn

`ln_ndt_link()` is gone. `rlddm()` derives its bound through
`ndt_bound()`, attaches it through `ndt_bound_attach()`, keeps a settled
one through `ndt_bound_of()`, and carries `ndt_bound_pending()` from
construction. `ndt_group()` is registered by frmtmb.eam's `.onLoad()`
and `rlddm()` now reads it.

The order in `family_finalize` is the one frmtmb.eam's own families use
and it is not arbitrary: `ndt_bound()` runs FIRST, so a `max_ndt` above
the fastest response is still refused on the second finalize of a family
whose bound is already settled, and only then is the carried bound read.

## RETRACTED: the pinned-`ndt` defect I reported does not exist

**I claimed that `bf(ndt = 0.2)` on `rlddm()` at 0.3.0 was silently
accepted and fitted 0.0435872023 s, an error of 78.2 percent. That is
false. It fitted 0.2, exactly.** Punch round 1's review
(`dev/reviews/2026-09-10-rlddm.md`, B1) falsified it and I reproduced
the falsification before accepting it. The number, the percentage and
the NEWS bullet built on them are withdrawn from NEWS, from the test
comment, from the proposed plan row and from here.

**Why my number was an artifact of my own script.** frmtmb transforms a
constant dpar TWICE and only the second one reaches the parameter:

* `R/parse.R:1266`, `plain_dpar()`, calls `linkfun(constant)` and uses
  the result ONLY to decide whether the constant is in range. The value
  is discarded and the raw constant is stored;
* `R/frame.R:2581` does the transform that counts,
  `par_template[["betad"]][idx] <- lp[["link"]]$linkfun(lp[["constant"]])`,
  and by then `family_finalize()` has run, so `lp[["link"]]` is the
  SETTLED scaled logit.

My first `rlddm-pin.R` never read the template. It pushed 0.2 through
the DECLARED `log` link by hand, got -1.60943791, and read that back
through the fitted link. Nothing uses that number. The rewritten script
reads the template instead, which is where the answer lives.

`dev/rlddm-scripts/rlddm-pin.R`, seed 4242, log at `rlddm-pin-log.txt`,
`min(rt) = 0.261523214`:

| | reference build (0.3.0) | this branch |
|---|---|---|
| `bf(ndt = 0.2)`, no `max_ndt` | ACCEPTED | REFUSED, by name |
| eta in the parameter template | **1.1789028** | (not reached) |
| the time it decodes to | **0.2** | |
| relative error | **0.00 percent** | |
| negative log-likelihood | 88.624345656 | |
| `bf(ndt = 0.2)` with `max_ndt = 0.26` | eta 1.2039728, ndt 0.2 | eta 1.2039728, ndt 0.2 |
| the same, negative log-likelihood | 88.624345656 | 88.624345656 |

1.1789028 is `log(0.2 / (0.261523214 - 0.2))`, the settled link's own
transform of 0.2. And the identical nll at two different bounds, 0.2615
and 0.26, is impossible if either arm had decoded a parse-time eta:
`linkinv(-1.6094)` is 0.04359 under one and 0.04334 under the other. My
own smoke control recorded the reference arm's `pin_02` as 88.62435 and
I did not read what it was telling me.

**The lesson, in the standing rules' own words.** "A printed zero is not
a measured zero" has a sibling I broke: a RECONSTRUCTED number is not a
measured number. I reconstructed a linear predictor from two links I
could see and never checked whether anything read it, on a code path I
had not traced. The check that would have caught it is the one the
rules already name, constructing the case where the guarded thing is
absent: had I asked what the fit did with the constant when the pending
bound was NOT there, I would have had to read the template.

## The defect the pending bound really does catch

Same script, same seed. The parse-time call is a RANGE CHECK, and at
0.3.0 it ran against the placeholder `log` link, which accepts any
positive constant. A constant above the bound the fit would actually use
therefore passed it, and `R/frame.R:2581` has no finiteness check of its
own.

| `bf(ndt = 0.3)` with `rlddm(max_ndt = 0.26)` | reference build | this branch |
|---|---|---|
| eta in the parameter template | **NaN** | (not reached) |
| negative log-likelihood at the start | **NaN** | |
| what a user sees from `frm()` | `NA/NaN gradient evaluation. The likelihood was undefined or unbounded somewhere the optimizer stepped` | `Constant ndt = 0.3 is not in the range of the scaled_logit link` |
| warnings raised on the way | `NaNs produced`, `NA/NaN function evaluation` | `NaNs produced` |

The base build is loud, and it is loud about the OPTIMIZER. The
diagnosis a user needs, that the constant they wrote is above the bound
their data set, is not in it. On this branch it is the whole message and
it arrives at parse.

**So does the refusal survive? Yes, and on two justifications, neither
of which is the one I first gave.**

1. It moves the out-of-range case from an optimizer failure that names
   nothing to a parse-time refusal that names the constant and the
   link.
2. It makes `rlddm()`'s `bf(ndt = )` contract the same as the four
   frmtmb.eam families it borrows its parameterization from. Measured in
   the same script: `wiener()` refuses `bf(ndt = 0.15)` with "the
   non-decision time bound is not set yet", and `wiener(max_ndt = 0.2)`
   accepts it. That refusal is 1.0a's B3 and is deliberate. Leaving
   `rlddm()` as the one family of five that behaves differently would be
   the worse outcome for a seam whose whole purpose is one contract.

**And the price is real and is now stated as a break, not a fix.**
`bf(ndt = c)` on a bare `rlddm()` WORKED for `c` below the fastest
response and now requires `rlddm(max_ndt = )`. NEWS says so under
BREAKING, `?rlddm` gains a section that says so, and the test that pins
it says in its comment that it pins a removal.

**The capability is respelled, not lost, and that is worth a clause in
NEWS because "this is refused now" reads worse than it is.** Same
script, same seed: `bf(ndt = 0.2)` with `rlddm(max_ndt = 0.26)` gives a
log-likelihood of **-86.3374111** on this branch and **-86.3374111** on
the reference build, which is also what `bf(ndt = 0.2)` on a bare 0.3.0
family gave. The user adds one argument and no number moves. The review
measured the same identity on its own design and got -85.7085293 there;
two designs, one conclusion.

**A core defect this exposed, and it is not mine to fix.**
`R/frame.R:2581` transforms a constant with no finiteness check, while
the check at `R/parse.R:1267` is against a link that a
`family_finalize()` may replace. Every family that derives its link from
the data has this hole, and a pending bound is a per-family workaround
for it. It deserves a core item.

## The controls: nothing that did not opt in moved

### frmtmb.learn, a model without `ndt_group()`

`dev/rlddm-scripts/rlddm-smoke.R`, seed 4242, six learners by 60 trials,
run once against `rellib-0552` and once against `rlddm-lib` and compared
by `dev/rlddm-scripts/rlddm-compare.R` under `identical()` with a ulp
count beside it. A printed zero is not a measured zero, so nothing here
is compared with a tolerance.

**Nineteen of 22 quantities, `identical()` TRUE, 0 ulp on every numeric
one:** the simulated response and choice, the log-likelihood
(-85.589871250631), `fixef()`, the optimizer's parameter vector, the
standard errors, the whole `vcov()`, `confint()`, `VarCorr()`,
`predict(dpar = "ndt", type = "response")` over all 360 rows,
`predict(dpar = "drift", type = "link")`, the whole 360-by-5 value
trace, the convergence code, the `max_ndt` arm's log-likelihood
(-85.708399866799) and coefficients, and a bare `wiener()`'s own density
and refusing link.

Three quantities differ and all three are intended: the two
pinned-constant probes now refuse rather than returning a number (which
is the REMOVAL, not a fix; see the retraction above), and the
`max_ndt`-too-high message is the seam's rather than the family's own
copy (both still contain "above the fastest response", which is what
`test-families.R` matches on). Log at `rlddm-smoke-compare-log.txt`.
Rerun after punch round 1 routed `ln_ndt_at()` through the exported
`ndt_apply()`: still 19 of 22, still 0 ulp.

**Regenerated in punch round 2, because the review checked the logs'
modification times and the first pair predated the install they
described by five minutes.** The lane library was installed at 02:30:19
and 02:30:23; the two compare logs are now 02:34:41 and 02:35:08. The
review had already re-run both on its own 02:09 install and got the same
answers, so this is bookkeeping rather than a result, but a log that
predates its install is not evidence and the rule is the lane's own.

**A FOURTH difference the 22 quantities did not cover, found by the
review.** `summary()` on an ungrouped fit prints
`ndt = scaled_logit` where 0.3.0 printed `scaled_logit(0, 0.2615)`. The
link's arithmetic is identical to the last digit and the bound is now on
the fitted family at `family(fit)$ndt_bound$ub`, which 0.3.0 did not
carry, so nothing numeric moved. NEWS says a model without
`ndt_group()` is "identical in every digit", which is true and is not
the same as unchanged, and it now carries a clause saying so.

**Decided NOT to restore the bound in the printed name.** It is
`ddm_scaled_logit()`'s bare name, and 1.0a chose that deliberately with
its reason written in the source: the bound belongs on the fitted family
where `?wiener` sends a reader, not in a link's display string, and two
eam tests pin the bare word. Restoring it for `rlddm()` alone would make
two families print the same link differently, and changing
`ddm_scaled_logit()` would reverse a 1.0a decision from this lane. It is
a display regression, it is recorded in NEWS, and the number it used to
show is now reachable by name.

### frmtmb.eam, which this lane also edits

`dev/rlddm-scripts/rlddm-eam-control.R`, seed 4242, six subjects by 60
trials. Six fits (`wiener`, `wiener(variability = "sv")`,
`wiener(max_ndt = )`, `rdm(2)`, `lba(2)`, and `wiener` WITH
`ndt_group(s)`), each contributing the log-likelihood, `fixef()`, the
parameter vector, `fitted()` and `predict(dpar = "ndt")`; the grouped
fit's floor table, its names, its sizes, its `ub` and `ndt_time()`; and
fixed-parameter probes of three densities, which are the density rather
than an optimizer path.

**Forty-one of 41 quantities `identical()` TRUE, 0 ulp.** Item 1.0a is
not weakened, including its own grouped arm. Log at
`rlddm-eam-compare-log.txt`.

## The defect filed against this item: keying on the label

**Decided NOT to do, and the reason is a measurement rather than a
preference.** `dev/rlddm-scripts/rlddm-labelkey.R`, log at
`rlddm-labelkey-log.txt`.

The recorded obstacle was that "the registry's coercion must return
numbers for the tape and is stateless". The fact underneath it is
stronger: **the label cannot come back OUT of the coercion.** frmtmb
wraps every registered coercion's result in `as.numeric()` at both call
sites, `R/frame.R:1219` and `R/predict.R:770` (the first version of
this file said 771, which is the `eval()` inside that call), and
`as.numeric()` strips every attribute. Those two are the only aterm
coercion sites.

Constructed rather than argued. `ndt_group` was re-registered with a
coercion that returns the code AND attaches the labels to it, and the
value was read back everywhere a family or a prediction can see it:

| where | attributes on the value |
|---|---|
| `fit$frame$aterm_values[["rt"]][["ndt_group"]]` | NONE |
| `aterms_for_newdata()`'s output | NONE |
| after an ordinary `v[1:3]` | NONE |

So recording the fitted labels on the family does not close it through
the coercion's RETURN VALUE: the newdata path has no label to compare
them with.

**My first conclusion, "not closable inside frmtmb.eam", is one word too
strong and the review found the route I missed.** The coercion cannot
return the label, but it is CALLED with it, on BOTH paths: the review's
`dev/rev-rlddm-labelkey.R` counts two coercion calls, one at frame
assembly and one inside `aterms_for_newdata()`, each handed the labels.
A side map written there and read at floor-lookup time closes the
collision residual with NO change to frmtmb.

The accurate sentence is therefore: **not closable through the
coercion's value, and closable in-package only through process-global
mutable state** in a registry this project deliberately made stateless.
The review measured that cost too: after one further unrelated fit the
map holds four entries in one process-wide store with no fit to scope
them to, and nothing to evict them.

**Do I want that trade? No.** The residual is 1.8e-15 per pairing of an
unseen label with a fitted group, structurally zero for every identifier
shape measured, and refused by name wherever both labels are in one
column. Against that, a process-global label map is unbounded state
shared by every fit in a session, is invisible to `refit()` and to
serialization, and would make two fits in one process able to affect
each other through a table neither of them declares. Trading a 1.8e-15
residual for that is the wrong direction, and it would be the kind of
hidden state this project has spent two rounds removing. The right fix
is still the core seam, and the plan row now says why in those terms
rather than saying "cannot".

**The residual that stays open**, quoted with the construction rather
than re-derived: a group label the fit never saw, whose code collides
with a label it did see, is paired with that group's bound instead of
being refused. `dev/reviews/2026-09-09-ndt.md` R2-1 bounds it at 1.8e-15
per pairing, proves that labels of one to three printable ASCII
characters cannot collide at all, and found 0 collisions in 10^8
identifiers of the form `S00000000`. This lane reproduced the shape of
that at 10^6: `ndt_bound_key()` over `S000000` through `S999999` gives
**1,000,000 distinct codes out of 1,000,000**. Where both labels are
visible in one column the fit is refused by name.

What this lane DID check is that the residual is the only hole: a group
the fit never saw is refused on both paths that can reach the floor
lookup from `rlddm()`, `ndt_time(newdata = )` and
`frm()` on the new rows, and the level-order invariance the label key
buys holds on this family too (see below).

## The one assertion I wrote and then took out, with the measurement

The first version of the rlddm scale row asserted the population
non-decision time as a `scale_z()` against the fit's own standard error,
under 4, which is what the eam row asserts. **It came back at 5.28, and
the instrument is what fails, not the fit.**

`dev/rlddm-scripts/rlddm-popndt.R`, log at `rlddm-popndt-log.txt`,
arithmetic on the row's own recorded fields. The reported quantity is
`plogis(b0) * mean(floors)`, where `b0` is the population intercept at
`re.form = NA`. The 9.29 ms gap to the design's stated 0.25 s
decomposes:

| step | value | moves by |
|---|---|---|
| the truth the design states | 0.250000 | |
| `E[ndt] = 0.25 * exp(sd_log^2 / 2)`, sd_log 0.15 | 0.252828 | +2.83 ms |
| the realized mean of the 100 drawn truths | 0.254744 | +1.92 ms |
| the fitted mean of the 100 learners | 0.257234 | **+2.49 ms** |
| `plogis(b0) * mean(floors)` | 0.259290 | +2.06 ms |

**Only the fourth row is the fit**, and it is 0.98 percent of a 254.7 ms
quantity, 0.063 of the truths' own spread. The other three are the
lognormal mean-against-median offset in how the design states its truth,
the realized draw's own deviation from that expectation, and a
conversion that is not the quantity: `mean(frac_i * floor_i)` is not
`mean(frac_i) * mean(floor_i)` when the two are correlated, and they
are, because a learner with a high floor needs a lower fraction for the
same non-decision time. The standard error, 0.00176, covers the
fraction's sampling uncertainty and none of the other three.

**The one-number form, which the review measured and which is stronger
than I wrote it.** The standard error the z divides by covers the
population fraction's sampling uncertainty. The realized-draw term it
does not cover has a standard deviation of its own
(`dev/rev-rlddm-instrument.R`):

| row | subjects | sd of the draw's own mean | the se the row reports | ratio |
|---|---|---|---|---|
| `learn-rlddm` | 100 | 0.003814 | 0.001759 | **2.17** |
| `eam` | 30 | 0.005537 | 0.001964 | **2.82** |

A z whose divisor is two to three times smaller than the standard
deviation of a term in its numerator is not a z.

**A FINDING AGAINST MERGED WORK: the eam tier's `ndt_z < 4` is a
seed-pinned assertion.** I wrote that its z of 1.60 was luck; the review
turned that into a probability and it is sharper than an observation.
The same decomposition on the shipped `dev/ndt-scripts/ndt-scale-r1.tsv`
gives a fitted mean of 0.246146 against a realized truth mean of
0.243452, that is +2.69 ms, essentially the same estimation error as
this row's +2.49 ms. The row passes because its 30 drawn truths averaged
8.35 ms below their own expectation `E[ndt] = 0.251806`, which cancels.
(The first version of this file said "8.35 ms below the stated 0.25";
that is wrong twice over, it is 6.55 ms below 0.25 and 8.35 ms below
`E[ndt]`.) Holding the fit's +2.69 ms and the conversion's +0.72 ms at
what this seed produced,
`extensions/frmtmb.eam/tests/testthat/test-scale.R:207` breaks once the
realized draw mean passes 0.254441, which is **P = 0.317** under the
design's own draw distribution.

A shipped tier assertion with roughly a one-in-three chance of failing
on a different seed is the class `dev/lane-rules.md` names. **It is not
fixed from this lane**, because the eam tier is 1.0a's and this lane
already changes eam; the plan row below proposes the fix, which is the
one this row took: record `ndt_z`, assert the per-subject statistics.

So the rlddm row RECORDS `ndt_pop`, `ndt_pop_se` and `ndt_pop_z` and
asserts the per-learner statistics instead, which is what
`dev/round-handoff.md` says item 2.1 should assert and for the same
reason: every learner below its own floor, the per-learner RMSE below
the spread of the truths it is estimating, the population mean inside
that same spread, and the correlation with the truth. Each is a ratio to
something the run measured.

## What else was measured

### The level order does not matter, on this family too

`dev/rlddm-scripts/rlddm-group.R`, seed 909. The same two rows through
four framings, and the character spelling of the same grouping:

| framing | `ndt_time()` |
|---|---|
| levels as fitted | 0.192696057, 0.406149555 |
| subset then `droplevels()` (slow only) | 0.406149555 |
| levels reversed | 0.192696057, 0.406149555 |
| character column, whole fit | log-likelihood -237.703403967, against -237.703403967 for the factor |

### The refusals reached from `rlddm()`

| construction | outcome |
|---|---|
| `max_ndt` and `ndt_group()` together | refused, "both set the non-decision time's upper bound" |
| a group the fit never saw, `ndt_time(newdata =)` | refused, "was not fitted to" |
| the grouping column absent from newdata | refused, and it names the term |
| a missing value in the grouping | refused, "every row needs one" |
| `ndt_group()` on a softmax family in this package | refused at frame assembly, "no family read it" |
| a grouped model whose density is reached without `ndt_floor` | refused by `ln_ndt_at()` |

### One thing that is NOT refused, and should not be

`predict(fit, newdata =, dpar = "ndt", type = "response")` on a grouped
`rlddm()` fit returns the FRACTION, and it returns it for a group the
fit never saw as well, because a fraction is the inverse link of a
linear predictor and no bound enters it. That is not the eam behavior,
where `predict(type = "response")` goes through a wrapped `mean_fn` and
refuses; `rlddm()` has no fitted mean at all, by design, so there is
nothing for the bound to reach. `ndt_time()` is the call that needs a
bound and it is the call that refuses. `?rlddm` and the compatibility
row both say `predict()` reports the fraction under this
parameterization.

### The bound is kept across a refit on fewer rows

`influence()`, `frm_simulate()` and the prior-predictive path all run
`family_finalize()` again on new rows. Dropping the fastest trial in the
whole data set and re-finalizing, `rlddm-group.R`:

| | floors |
|---|---|
| all rows | 0.455662732, 0.225935102 |
| fastest response dropped | 0.455662732, 0.225935102 |

and the re-finalized family's `aterm_data()` still returns ONE
`ndt_floor` entry rather than two.

### The per-trial factorization is untouched

`sum(log(trace$dens))` is -237.703403967 against `logLik(fit)`
-237.703403967 on the grouped fit, which is the identity
`test-engine.R` asserts for the ungrouped one.

## What was NOT done

**`frm_task_simulate()` does not gain a per-group bound, and needs
none.** It builds its per-row values out of `pars` and the data map, so
no `ndt_floor` and no `ndt_group` is ever in them, and `ndt` arrives as
a time already. `rlddm()`'s `draw` goes through `ln_ndt_at()` anyway,
which is a no-op there today; the reason is written in the source, and
it is that a family reading one parameter two different ways in its two
halves is how a fraction gets drawn as seconds the day a simulator is
turned on.

**frmtmb.eam's own four families were not moved onto the exported
seam.** They keep `ddm_ndt_finalize()` and `ddm_ndt_install()` exactly
as item 1.0a shipped them, which is why the 41-quantity control above
comes back at 0 ulp. The exported functions are thin wrappers over the
SAME internals, `ddm_ndt_spec()`, `ddm_ndt_scaler()`,
`ddm_scaled_logit()` and `ddm_coerce_ndt_group()`, so there is one
implementation of the arithmetic and of the refusals rather than two.

**`ndt_time()` is NOT re-exported from frmtmb.learn.** A user of a
grouped `rlddm()` fit writes `frmtmb.eam::ndt_time(fit)`, because
`library(frmtmb.learn)` does not attach frmtmb.eam. Re-exporting is two
lines and was declined for consistency: this package already imports
`ddm_simulate()` and `wiener_lpdf()` without re-exporting either, and an
export added to frmtmb.learn is a promise frmtmb.learn then owns. `?rlddm`,
its compatibility row and NEWS all spell the call qualified. Worth
revisiting if a second frmtmb.eam function joins the list.

**`gddm()` still does not take `ndt_group()`**, on the scope argument
item 1.0a recorded, and the unnumbered defect that lane filed against it
is untouched here.

**No absolute tolerance went into a test.** The scale row's assertions
are a convergence code, a Hessian flag, a count of non-finite standard
errors, two ratios to a spread this run measured, and a correlation.
`test-rlddm-ndt.R`'s before-and-after test compares the two arms of ONE
run against each other rather than against any constant.

**`dev/scale-findings.md` is not edited here**, on the same rule that
kept item 1.0a out of it: the tier's own `.tsv` is the record and the
lane's findings carry the comparison. What that file's `learn-rlddm` row
should gain when it is next rebuilt is the diagnostic string, which it
does not carry at all today: `conv=0,maxgrad=0.00312,pdHess=TRUE,
nbadse=0` where it used to be `conv=1,maxgrad=6.36e+09,pdHess=FALSE,
nbadse=4`.

## Which version bump this needs

Not chosen here.

**frmtmb.eam.** Six new exports, a new class on a slot that was a plain
list, and one new field on it. Nothing a 0.7.0 model can express
changed, and the 41-quantity control says so at 0 ulp, rerun after the
punch round. That reads MINOR.

**frmtmb.learn.** A new opt-in parameterization, the family's bound now
coming from another package, and one thing that STOPS WORKING:
`bf(ndt = )` on a bare `rlddm()`. A model that does not write
`ndt_group()` is identical in every digit, and `summary()` prints one
link name differently.

The first version of this section called the `bf(ndt = )` refusal a fix
for a silent wrong answer and read the bump as minor on that basis. That
was wrong. It is a REMOVAL of working behavior, bought with a real
refusal on the out-of-range case and with consistency across the five
families. A user who pinned `ndt` has to add `max_ndt`. **That is my
reading of a minor bump with a BREAKING bullet rather than a patch, but
it is the user's call and the NEWS bullet is written so the trade is
visible without reading this file.**

**`frmtmb.learn`'s floor on frmtmb.eam has to rise.** `DESCRIPTION:39`
still says `frmtmb.eam (>= 0.6.0)` and `NAMESPACE` now imports five
functions frmtmb.eam gains in this round, so `library(frmtmb.learn)`
against eam 0.6.0 fails at NAMESPACE LOAD rather than at a call site.
It is not edited here because the version those exports ship in is not
this lane's to choose, and the review agrees the reasoning is right and
the constraint is still a merge blocker: it must rise in the same commit
that lands the exports.

## What the plan should say instead

`dev/extension-gaps-plan.md` is not edited here. Item 1.0b's row should
read:

* the deliverable landed and the acceptance criterion is met: the
  `learn-rlddm` row goes from `conv=1`, a maximum gradient of 6.36e+09,
  a Hessian that is not positive definite and four NaN standard errors
  at -8333.1, to `conv=0`, 3.12e-03, a positive definite Hessian and no
  NaN standard errors at -7781.2. **Say "now estimable" and not "551.9
  log-likelihood units better":** the before fit is not at a maximum of
  anything, and the per-group bound admits, for 99 of the 100 learners,
  times the global bound cannot express at any parameter value, so a
  higher likelihood on a larger admissible set is expected rather than
  evidence;
* the recovery statistics are the fit's OWN and 1.0a's floors finding
  does not carry over: 13.996 ms per-learner RMSE and 0.9404 correlation
  against an ANALYTIC bound of 21.075 ms and 0.847780 for any estimator
  without a random effect on `ndt`, confirmed by a no-random-effect fit
  at 22.107 ms with `cor(fitted, own floor) = 1`. That is the sentence
  the row should carry, and it is the reviewer's measurement;
* **"below their own fastest response" cannot fail** and the acceptance
  criterion should stop reading it as a recovery check: under a
  per-group bound `ndt_i = plogis(eta_i) * floor_i` and `plogis < 1`, so
  100 of 100 is an identity. It checks the floor LOOKUP, which is worth
  having. The tightest margin, 20.9 ms, is the informative number. The
  eam row has the same assertion with the same property;
* the row should carry the `ndt` variance component, which is the
  clearest single symptom: 2.236 on the link under the global bound
  against 0.3261 under the per-learner one;
* **`cor_true = 0` is the WRONG truth for a grouped arm, and item 2.2
  has to know before it runs.** The grouped model does not estimate the
  deviation the design draws; it estimates `qlogis(ndt_i / floor_i)`,
  and the floor is a draw from the learner's whole parameter vector. The
  drawn truths already carry -0.7768 between `bs` and `ndt` in that
  parameterization, where as drawn they carry +0.0604, and the fit
  reports -0.9160. It is computable per replicate from the truths and
  the floors with no fit, and `test-scale.R` now computes it. **Item
  2.2 is 60 replicates of exactly this design with the block on every
  parameter: a diagonal truth there would fail on nearly every
  replicate for a reason that has nothing to do with the estimator;**
* **item 2.2 also needs the `sd(ndt)` warning item 2.1 already carries**,
  and 2.2 carries neither. `sd(ndt)` on the natural scale must not be
  read as the variance component recovering. 1.0a's B5 recurs here:
  0.0334 from the full fit against 0.0380 from the no-random-effect fit
  and 0.0370 from a constant fraction of the floors, truth 0.0393, so
  the full model is WORSE on the statistic than leaving the component
  out. Assert the per-learner error and the log-likelihood;
* the export seam is SIX functions and none of them is
  `ddm_ndt_install()`. That one rewrites five slots BY NAME and
  `rlddm()`'s likelihood is in none of them. The sixth, `ndt_apply()`,
  was added in punch round 1 because the seam had shipped the
  multiplication as three lines of documented arithmetic and its refusal
  as prose, and the documented three lines run with the floor absent
  return a fraction as a time;
* the reviewer's fifth entry, `ndt_bound_pending()`, is justified but
  NOT by the defect this lane first claimed. **That claim is retracted:
  `bf(ndt = 0.2)` on `rlddm()` at 0.3.0 fitted 0.2 exactly.** What the
  pending bound buys is the out-of-range case, `bf(ndt = 0.3)` with
  `max_ndt = 0.26`, which reached the objective as NaN and surfaced as
  an optimizer failure naming nothing, and one `bf(ndt = )` contract
  across all five families. It costs a working spelling: `bf(ndt = )` on
  a bare family now needs `max_ndt`;
* a new CORE item: `R/frame.R:2581` transforms a constant dpar with no
  finiteness check, while the range check at `R/parse.R:1267` runs
  against a link a `family_finalize()` may replace. Every family that
  derives its link from the data has this hole and a pending bound is a
  per-family workaround for it;
* a new item against 0.7.0: `extensions/frmtmb.eam/tests/testthat/`
  `test-scale.R:207` asserts `scale_z(...) < 4` on a quantity whose
  standard error is 2.82 times smaller than the standard deviation of a
  term in its numerator. It passes at this seed by cancellation and
  breaks with **P = 0.317** under the design's own draw. It should
  record `ndt_z` and assert the per-subject statistics, which is what
  the learn row now does;
* **keying on the label is closable in-package only through
  process-global state**, which is the correction the review made to
  this lane's "not closable". frmtmb wraps every registered coercion's
  result in `as.numeric()` at `R/frame.R:1219` and `R/predict.R:770`, so
  the label cannot come back out of the coercion; but the coercion is
  CALLED with it on both paths, so a side map written there would close
  the 1.8e-15 residual with no core change. This lane declines that
  trade, because the map is unbounded state shared by every fit in a
  session with nothing to scope or evict it. The core seam is still the
  right fix and the row should say why in those terms;
* a new item, small: `predict(dpar = "ndt", type = "response")` means a
  FRACTION on a grouped model and SECONDS on an ungrouped one, and
  `ndt_time()` is the call that is the same under both. That asymmetry
  is documented in three places now and is still a trap;
* and a standing limit of the seam, for whoever writes the next
  consumer: a foreign family that attaches a bound and does NOT call
  `ndt_apply()` fits, converges with a positive definite Hessian and no
  bad standard errors, gives an identical log-likelihood, and reports a
  non-decision time measured at 37.5 percent too small. No package can
  check another package's arithmetic, so this is documented in
  `?ndt_apply` rather than guarded, and an affirmative flag was
  considered and declined for the reason 1.0a measured.

## Tests, and the check

### The runner, and the defect it hid

One test file per R process, `NOT_CRAN=true`, the failure cap lifted and
`package =` set so internals resolve, through
`dev/rlddm-scripts/rlddm-testfile.R`. It prints `err` beside `fail` and
lifts testthat's ten-failure cap, because a capped number and a file
that aborted halfway have both been reported as counts in this
repository.

**Its first version attached frmtmb, frmtmb.eam AND frmtmb.learn before
running any file, and that hid a real defect for the whole of this
lane.** A package's own `tests/testthat.R` attaches ONE package. Five
assertions in `test-rlddm-ndt.R` called `ndt_time()` unqualified, which
resolves when frmtmb.eam is on the search path and does not resolve
under `R CMD check`: `Status: 1 ERROR`, `could not find function
"ndt_time"`, five times. `R CMD check` found it, on the final pass,
which is the one place the round's cost rule still runs it.

The runner now attaches only the package under test, and the whole suite
of both packages below was re-run under it. The calls are qualified
`frmtmb.eam::ndt_time()`, which is what a user of `frmtmb.learn` writes
as well: frmtmb.eam is an `Imports` and not a `Depends`, so
`library(frmtmb.learn)` does not put it on the search path, and `?rlddm`
and NEWS both spell it qualified.

Three smaller harness faults were found the same way and are recorded
here because a harness is a guard, and every guard this round built
failed open on its first try.

* The runner's PowerShell wrapper used `2>&1`, which in Windows
  PowerShell 5.1 wraps a native command's stderr in ErrorRecord objects
  and stopped `Select-String` matching the RESULT line on two files that
  pass on their own.
* The joiner against `dev/suite-baseline.tsv` summed such a file as
  ZERO until it was changed to print it as a hole and to count the
  holes.
* **A count read off a STALE log.** In punch round 1 the eam suite's
  summary file was read while the current run was still going, so the
  join reported `test-ndt-seam.R` at 53 over 10 tests when the file on
  disk had 13 tests. The file's modification time is what caught it. The
  counts below are all from logs verified to postdate the install they
  describe, and a one-assertion discrepancy between a standalone run and
  a suite run of `test-rlddm-ndt.R` was chased to the same cause rather
  than being averaged away: 53 is the number, reproduced standalone and
  in-suite.

### Seeing it fail

Both new files, run against the round's shared reference build of the
base commit, `rellib-0552`:

| file | on the base commit | on this branch |
|---|---|---|
| `frmtmb.eam/tests/testthat/test-ndt-seam.R` | `PASS 0 FAIL 1 ERR 13 SKIP 0` over 13 tests | `PASS 83 FAIL 0 ERR 0` |
| `frmtmb.learn/tests/testthat/test-rlddm-ndt.R` | `PASS 5 FAIL 6 ERR 10 SKIP 0` over 12 tests | `PASS 53 FAIL 0 ERR 0` |

The seam file's failures are mostly the weak form, `could not find
function "ndt_bound"`, plus one behavioral: the bound record on a fitted
`wiener()` is a plain list rather than a `"frmtmb_eam_ndt_bound"`.

**The first version of this file said six of the fifteen items were
behavioral. By the standing rule's meaning it was ONE, and the review
was right to correct it. Punch round 1 added a second** by adding the
out-of-range test, which is the case that IS a defect. The honest tally,
over the 16 items the current file reports on the base commit:

| what fails on the base commit | class |
|---|---|
| `ndt_time()` on an `rlddm()` fit refuses by name | **behavioral** |
| `bf(ndt = 0.2)` is accepted where the test now expects a refusal | behavioral, but it pins the REMOVAL, not a bug |
| `bf(ndt = 2 * min(rt))` out of range is accepted, and reaches the objective as `NaN` | **behavioral**, and it is the defect the pending bound really catches |
| `bd[["ub"]]` is `NULL` | the slot does not exist: weak form |
| the link name string differs | vocabulary |
| `frm_compat()` says `untested`, twice | a doc row not written: weak form |
| `ln_ndt_at` not found | weak form, and acknowledged |
| `no family read it`, eight times | the base commit refusing a new term |

So two items are behavioral in the sense the rules mean, and one of the
two pins a removal rather than a fix. That is not fatal and it is worth
saying why: **the thing this item fixes is a new
CAPABILITY, so its before-and-after cannot live in the test file at
all.** A test file can only assert what the branch does; the branch's
`ndt_group()` model has no counterpart on the base commit to fail
against. The before-and-after lives where it can, in
`dev/rlddm-scripts/rlddm-group.R` and in the scale row, and both are
run on both builds.

The eam seam file's tally was honest as written: nine of its ten are
`could not find function` and one is behavioral, the bound record's
class.

### The whole suite of both packages

Joined per file against `dev/suite-baseline.tsv` by
`dev/rlddm-scripts/rlddm-baseline.R`, which prints a file with no
RESULT line as a HOLE rather than summing it as zero.

**frmtmb.learn**, `dev/rlddm-scripts/rlddm-suite-learn.txt`: 12 files,
**386 assertions, FAIL 0, ERR 0**, 11 skips, no non-zero exit, no file
with no RESULT line. Every file is AT its `dev/suite-baseline.tsv`
count except the new one:

| file | this branch | baseline |
|---|---|---|
| `test-rlddm-ndt.R` | 53 | (new) |
| `test-bracket-access.R` | 1 | 1 |
| `test-clean-session.R` | 4 | 4 |
| `test-counterfactual.R` | 67 | 67 |
| `test-engine.R` | 19 | 19 |
| `test-factorization.R` | 35 | 35 |
| `test-families.R` | 45 | 45 |
| `test-message-uniqueness.R` | 4 | 4 |
| `test-reference.R` | 12 | 12 |
| `test-scale.R` | 0, SKIP 2 | 0, SKIP 2 |
| `test-stan-identity.R` | 0, SKIP 9 | 0, SKIP 9 |
| `test-surface.R` | 146 | 146 |

The two zeros are the gated tiers, and both are at their baseline SKIP
counts rather than silently absent.

**frmtmb.eam**, `dev/rlddm-scripts/rlddm-suite-eam.txt`: 23 files,
**1570 assertions, FAIL 0, ERR 0**, 3 skips, no non-zero exit, no file
with no RESULT line. Every file is AT its baseline count except the new
one:

| file | this branch | baseline |
|---|---|---|
| `test-ndt-seam.R` | 83 | (new) |
| `test-ndt-bound.R` | 79 | 79 |
| `test-extension-api.R` | 14 | 14 |
| `test-family.R` | 51 | 51 |
| `test-defects.R` | 59 | 59 |
| `test-brms-parity.R` | 13 | 13 |
| `test-variability.R` | 140 | 140 |
| `test-lba.R` | 109 | 109 |
| `test-rdm-gng.R` | 221 | 221 |
| `test-simulate-density.R` | 69 | 69 |
| `test-density.R` | 132 | 132 |
| `test-gddm-family.R` | 77 | 77 |
| `test-gddm-gradients.R` | 22 | 22 |
| `test-gddm-recovery.R` | 28 | 28 |
| `test-gddm-reference.R` | 166 | 166 |
| `test-gddm-solver.R` | 94 | 94 |
| `test-moments.R` | 29 | 29 |
| `test-sampling.R` | 97 | 97 |
| `test-surface.R` | 50 | 50 |
| `test-units.R` | 32 | 32 |
| `test-message-uniqueness.R` | 4 | 4 |
| `test-bracket-access.R` | 1 | 1 |
| `test-scale.R` | 0, SKIP 3 | 0, SKIP 3 |

`test-sampling.R` at 97 is worth naming: `dev/ndt-findings.md` recorded
it at 94 with one skip on that lane's private library and the reviewer
reran it at 97. It is 97 here.

### The gated tiers

**frmtmb.learn's Stan identity tier**, the one that would see a change
to `rlddm()`'s density: it checks the recursion against an independent
Stan program of the same model at the same estimates.
`dev/rlddm-scripts/rlddm-gated.R`, `FRMTMB_BRMS_FIT_TESTS=true`,
`NOT_CRAN=true`, log at `rlddm-stan-log.txt`:
**FAIL 0, ERR 0, SKIP 0, PASS 55** over 9 tests, which is the count
`dev/round-handoff.md` records for this tier. Two warnings, both on
`bandit2arm_dual` and both pre-existing.

The library was checked before the result was believed, which is what
the lane rules require: **StanHeaders 2.32.10 from
`C:/Users/adf44/source/r/pinlib`**, not the user library's 2.39.1, and
rstan 2.32.7. The runner refuses outright if that is not what it finds.

And the honest limit on what that run proves: the cache held 113 files
before and 113 after, so nothing compiled fresh and this is a check of
the R side against previously compiled programs. That is the right
answer here, because this change alters no Stan program: the tier's
programs are in `helper-stan-programs.R` and are byte-identical to the
base commit's.

The eam scale tier was rerun as well, because this lane changes the
SHAPE of the record `test-scale.R` reads at `family(fit)$ndt_bound`.
`dev/rlddm-scripts/rlddm-scale-eam.tsv`, row `eam`, seed 20260908,
30 subjects by 400 trials, `FAIL 0 | SKIP 2 | PASS 7`. It reproduces
`dev/ndt-scripts/ndt-scale-r1.tsv` to EVERY digit that file records:

| | item 1.0a's record | this branch |
|---|---|---|
| log-likelihood | -7003.01 | -7003.01 |
| convergence code | 0 | 0 |
| maximum absolute gradient | 0.000987 | 0.000987 |
| positive definite Hessian | TRUE | TRUE |
| population `ndt` | 0.246867, se 0.0019639 | the same |
| `sd(ndt)` natural | 0.0254253 | 0.0254253 |
| tightest margin | 27.3347 ms | 27.3347 ms |
| subjects below their own floor | 30 of 30 | 30 of 30 |
| condition effect | 0.910877 (0.858162, 0.963592) | the same |

That is a fourth independent reproduction of 1.0a's headline, on a
fourth build. The review reproduced it a fifth time and compared it
field by field rather than by eye: 34 of the record's 41 fields
identical, and all 7 that differ are clocks or memory.

**That row was run BEFORE punch round 1 and was not rerun after**, on
the round's cost rule. What the punch round changed in frmtmb.eam is
`ndt_bound()`'s input validation, `ndt_bound_attach()`'s guard and the
new `ndt_apply()`, and eam's own four families reach none of the three:
they go through `ddm_ndt_spec()` and `ddm_ndt_install()`, which are
untouched. The 41-quantity eam control WAS rerun after the punch fixes
and comes back 41 of 41 identical at 0 ulp, which is the evidence that
eam did not move.

### `R CMD check --as-cran`

`R CMD build` then `R CMD check --as-cran` on the built tarball, with
pandoc and TinyTeX on PATH and `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE`,
through `dev/rlddm-scripts/rlddm-check.ps1`. No
`--no-build-vignettes`: that flag manufactured two WARNINGs and a NOTE
on all eight packages in the 0.55.2 round.

| package | status | log |
|---|---|---|
| frmtmb.learn | **OK**, no NOTE and no WARNING | `rlddm-check-learn.txt` |
| frmtmb.eam | **1 NOTE**, and it is the expected one, "Skipping checking math rendering: package 'V8' unavailable" | `rlddm-check-eam.txt` |

frmtmb.learn's first check run was `Status: 1 ERROR`, and that is the
one recorded above: five unqualified `ndt_time()` calls that the lane's
own runner resolved and `R CMD check` did not.

Both packages were checked again after punch round 1, because the punch
changed `R/` in both, and the table above is that second pair of runs.


## The scripts, and what each one is for

| script | what it establishes | log |
|---|---|---|
| `rlddm-prelude.R` | the `.libPaths()` order, sourced by `rlddm-probe.R` | |
| `rlddm-probe.R` | the design probes run BEFORE any change: what `ndt_group()` did to `rlddm()` at the base commit | (stdout) |
| `rlddm-smoke.R` | the control: no `ndt_group()` is the 0.3.0 parameterization, 22 quantities under `identical()`. Writes an `.rds` per arm, which is not kept in the tree | `rlddm-smoke-compare-log.txt` |
| `rlddm-eam-control.R` | the other control: six eam fits and three density probes, 41 quantities under `identical()` | `rlddm-eam-compare-log.txt` |
| `rlddm-compare.R` | `identical()` and a ulp count over two such records | the two compare logs |
| `rlddm-popndt.R` | why the population non-decision time is recorded and not asserted | `rlddm-popndt-log.txt` |
| `rlddm-suite.ps1` | one package's whole suite, one file per R process | `rlddm-suite-eam.txt`, `rlddm-suite-learn.txt` |
| `rlddm-gated.R` | one gated tier, refusing to run at all unless StanHeaders is the pinned 2.32.10 | `rlddm-stan-log.txt` |
| `rlddm-baseline.R` | joins a suite run against `dev/suite-baseline.tsv`, per file | (stdout) |
| `rlddm-check.ps1` | `R CMD build` then `R CMD check --as-cran`, once | `rlddm-check-eam.txt`, `rlddm-check-learn.txt` |
| `rlddm-group.R` | the behavioral failure, both parameterizations, and every refusal and invariance reached from `rlddm()` | `rlddm-group-log.txt` |
| `rlddm-pin.R` | what `bf(ndt = )` actually fitted, read off the parameter template, both builds. REWRITTEN in punch round 1: its first version reconstructed an unused eta and produced a retracted number | `rlddm-pin-log.txt` |
| `rlddm-guard.R` | which families `ndt_bound_attach()` refuses, exhaustive over this package's five, with the marker the old guard rested on printed beside the outcome | `rlddm-guard-log.txt` |
| `rlddm-labelkey.R` | what the coercion's VALUE can and cannot carry, and the residual's own numbers. Its conclusion was narrowed in punch round 1: the label cannot come back out of the coercion, but the coercion is called with it, so the in-package route is process-global state and this lane declines it | `rlddm-labelkey-log.txt` |
| `rlddm-cortrue.R` | the truth of the correlated block the grouped row estimates, computed from the drawn truths and the floors with no fit | `rlddm-cortrue-log.txt` |
| `rlddm-scale.R` | one scale row in a fresh process, either arm | `rlddm-scale-before-log.txt`, `rlddm-scale-after-log.txt`, and the two `.tsv` |
| `rlddm-testfile.R` | one test file per process, failure cap lifted, `err` printed beside `fail` | (stdout) |
