# The within-condition contract in `gddm()`

The unnumbered Phase 1 item that item 1.0a filed and did not fix:
`gd_densities()` reads every distributional parameter at the first row
of its condition, so a parameter that varies WITHIN a condition is
silently ignored. Written on the `wt-gddm` worktree off f9ee297.
R 4.6.1, Windows 11, 16 logical cores.

Arms. The BEFORE arm is the round's shared reference build,
`C:/Users/adf44/source/r/rellib-0552`, frmtmb 0.55.2 and frmtmb.eam
0.7.0. The AFTER arm is this worktree's frmtmb.eam installed into
`C:/Users/adf44/source/r/gddm-lib`. `frmtmb` itself is untouched on this
branch, so both arms load the SAME core from `rellib-0552` and no core
difference can contaminate a comparison. Every script takes its arm from
`GDDM_LIB` and its seed from its own header, and lives in
`dev/gddm-scripts/` with its log beside it.

## The claim, and whether it holds

**A dpar that varies within a condition is refused by name, and a dpar
that does not is unaffected.** Both halves hold.

| | before | after |
|---|---|---|
| `mu ~ x` varying inside a condition | fits, no error, no warning | refused, naming `mu`, `x` and the condition |
| `ndt ~ x` | fits | refused, naming `ndt` |
| `bs ~ x` | fits | refused, naming `bs` |
| `bias ~ x` | fits | refused, naming `bias` |
| `mu ~ 1 + (1 \| s)`, `s` crossing the index | fits, `sd(s)` = 7.6e-11 | refused, naming `mu` and `s` |
| 21 correct designs | accepted | accepted, 0 false alarms |
| 26 ways of dropping one variable from those indices | accepted | 26 refusals, each naming the dropped variable |
| gaussian, wiener, lba, rdm, wiener_gng | accepted | accepted |

## The mechanism, reproduced

`dev/gddm-scripts/gddm-firstrow.R`, seed 5, 120 rows in two conditions,
the BEFORE arm. The objective is evaluated at a FITTED parameter vector
after adding 100 to a covariate on 118 of the 120 rows. At the starting
values every regression coefficient is zero, so the linear predictor
does not depend on the covariate and the probe would pass for the wrong
reason; the reviewer of item 1.0a records that trap and this
reproduction keeps it.

| model | fn as drawn | fn, `x + 100` on 118 of 120 rows | `identical()` | fn, `x + 0.01` on the 2 FIRST rows |
|---|---|---|---|---|
| `mu ~ x` | 18.962737109172043 | 18.962737109172043 | TRUE | 18.963062916902754 |
| `ndt ~ x` | 18.810414088214987 | 18.810414088214987 | TRUE | 18.811489838876835 |
| `bs ~ x` | 18.948550650674157 | 18.948550650674157 | TRUE | 18.948948928461569 |
| `bias ~ x` | 17.746051487984953 | 17.746051487984953 | TRUE | 17.746570516609207 |

Bitwise identical on all four, by `identical()` rather than by a printed
zero, and moving the two first rows alone does move it. The fit runs to
convergence with no error and no warning and reports `mu.x = 0.329958`
from two rows out of 120. That is `dev/reviews/2026-09-09-ndt.md`'s F1
reproduced, on a later base commit and a different library, plus a
fourth dpar it did not cover.

**The `vreal()` control fires.** The family's existing check on the
drift covariate refuses a coherence that varies inside a condition:
"gddm: vreal1 is not constant within every condition". That is the
guarded thing PRESENT, and it works.

### The random-effect half, and the instrument it needed

The same defect from the other side: a random-effect grouping that
crosses the condition index. Two probes were wrong before one was right,
and both are recorded because each failed for a plausible-looking
reason.

**The marginal likelihood is the wrong instrument.** With `s` crossing
the index, `dev/gddm-scripts/gddm-firstrow2.R` reports the objective
unchanged when the FIRST rows are relabelled, which looks like evidence
and is not: the subject deviations are integrated out and exchangeable,
so permuting which subject sits on a condition's first row leaves the
Laplace likelihood alone whether or not the other rows are read. The
fitted variance component also collapses to 7.639656e-11, so a probe at
the fitted vector has every `b` at zero for the wrong reason.

**The inner solution is the right one.** A deviation the likelihood
never reads is pinned at exactly zero by its own prior, whatever the
outer parameters are. `dev/gddm-scripts/gddm-firstrow-z.R`, seed 5, 120
rows, 4 subjects, outer vector set away from the optimum:

| condition index | subjects on a condition's first row | inner `b` | bitwise zero |
|---|---|---|---|
| `blk` alone, `s` crosses it | 1, 2 | -1.7661182357184526, -1.8037352708150676, 0, 0 | **2 of 4** |
| `gddm_conditions(blk, s)` | 1, 2, 3, 4 | -1.714, -1.732, -1.681, -1.740 | 0 of 4 |

Two subjects, 60 rows each, never reach the density at all, and the fit
reports a between-subject spread of 7.6e-11 rather than refusing.

## The decision: refuse rather than read, and what that cost

The alternative is to make the density READ the varying parameter. There
is exactly one way to do that, and it is not a solver change: the tape
cannot compare parameter values, so the distinct parameter vectors have
to be found from the DESIGN at frame assembly, and one Fokker-Planck
solve run per distinct vector. That operation is what
`gddm_conditions()` already performs. **"Read it" and "refine the
condition index" are the same thing**, so the only question is who does
it and whether the cost is visible.

The cost, `dev/gddm-scripts/gddm-cost.R`, seed 77, 120 rows, grid
`dt = 0.05, ny = 51, t_max = 2`, arms interleaved in a palindrome with
`gc()` before each and the minimum over 3 rounds:

| ncond | tape build s | one `fn` s | ms per solve | `fn` against ncond = 2 |
|---|---|---|---|---|
| 2 | 0.970 | 0.0002 | 0.081 | 1.0 |
| 4 | 1.500 | 0.0003 | 0.076 | 1.9 |
| 8 | 3.650 | 0.0007 | 0.092 | 4.6 |
| 20 | 8.220 | 0.0017 | 0.083 | 10.4 |
| 60 | 24.440 | 0.0066 | 0.111 | 41.2 |
| 120 | 39.710 | 0.0084 | 0.070 | 51.9 |

One solve per ROW against one per condition is **40.9x the tape build
and about 52x the evaluation**, on a grid five times coarser in time and
four times coarser in space than the shipped default. The exact
quantity needs no clock: solves per likelihood evaluation IS `ncond`, so
the ratio is 60x by construction and the measured 52x is that within
the instrument's own drift.

**The controls.** The tape column carries `ncond = 2` at both ends of
the palindrome at 0.990 and `ncond = 4` at 1.067, so read it as
measured. The `fn` column's same controls come back at 1.348 and 1.328,
a systematic drift toward the end of each round rather than noise, so
read the `fn` column to within about a third. The `ms per solve` column
is flat at 0.070 to 0.111 across a 60-fold range of `ncond`, which is
the structural check the clock cannot give: the cost really is linear in
the number of solves.

So the refusal is the answer, for three reasons that are measurements
rather than preferences.

1. **A silent refinement would hide a 60-fold cost.** A continuous
   covariate on any parameter gives one condition per row. A user who
   wrote `mu ~ rt_prev` would get a fit 52 times slower with no
   statement that they had asked for it, and on the shipped grid that
   is minutes of tape build per evaluation rather than one second.
2. **The contract is already documented.** `?gddm` has stated, in bold,
   that every row sharing a condition must share every parameter value,
   since the family shipped. This is enforcement of a stated contract,
   not a new restriction.
3. **The remedy is one call and it is in the message.** Naming the
   varying variable in `gddm_conditions()` gives exactly the model the
   user wrote, at exactly the cost correctness requires, and the cost is
   then visible as a condition count.

**What the 60x rejects, stated precisely.** It rejects refining the
index WITHOUT SAYING SO. The review is right that this is not the same
as rejecting every middle path, and the distinction matters because the
whole objection was that the cost would be hidden. An opt-in
`gddm_control(refine_conditions = TRUE)` that performs the same
operation `gddm_conditions()` performs, reports the resulting condition
count and the multiple it implies, and lets the user decide, makes the
cost visible and is not what the measurement above argues against. It
is not built here, because the refusal plus a message that hands back
the exact `gddm_conditions()` call is most of the same benefit at none
of the risk, and because an opt-in that silently changes `ncond`
between a fit and a `predict(newdata = )` needs its own study. It is
filed as a plan row rather than closed.

**What a user with a legitimately varying parameter is told.** The
refusal names the parameter, the variable, one condition it varies
inside and how many conditions in all, and then says to name that
variable in the index as well. There is no case where the varying term
was legitimate under the old behavior: it reached nothing.

## What changed

| file | what |
|---|---|
| `extensions/frmtmb.eam/R/gddm.R` | `gd_check_condition_constancy()` and its helpers `gd_ne()`, `gd_col_floor()`, `gd_ne_col()`, `gd_varying_groups()`, `gd_dpar_vars()` and `gd_and()`; the refusal in `gddm_simulate()`; the `?gddm` "What the data must carry" section rewritten, because it said the family cannot check this; `?gddm_conditions` no longer says an under-specified index is undetectable; `?gddm_simulate` states the constraint on `...` |
| `extensions/frmtmb.eam/R/zzz.R` | the registration, and three compat rows corrected (below) |
| `extensions/frmtmb.eam/NEWS.md` | a bullet under a new development heading. No released section touched, no version bumped |
| `extensions/frmtmb.eam/tests/testthat/test-gddm-conditions.R` | NEW, 58 assertions |
| `extensions/frmtmb.eam/man/gddm.Rd`, `man/gddm_conditions.Rd`, `man/gddm_simulate.Rd` | roxygenised |

`NAMESPACE` is unchanged: every new function is `@noRd`. Nothing in
`R/ddm-shared.R` was touched.

### The compatibility rows this change made wrong, and what they say now

A guard that refuses something is a change to what composes, so it is a
change to the table. Three rows in `ddm_compat_rules()`:

- **`gddm x kind:covstruct`, new.** Core's kind-level default says
  every family by covariance-structure pair works, because a covariance
  structure acts on the linear predictor and not on the response. That
  sentence is still true and the conclusion it reaches is now false
  here: a deviation enters the predictor, the predictor is read at its
  condition's first row, so `(1 | g)` is refused unless `g` is in the
  index. One family-side row beats the kind-level default on
  specificity, which is the same route `frmtmb.latent`'s
  `hmm x kind:covstruct` row takes. It reads `conditional`, and it
  carries the 7.6e-11 measurement, so 24 rows that read `works` now
  read `conditional` with a reason. Core's `R/compat.R` is not touched.
- **`gddm x mixture`, `untested` to `refused`.** Exercised by this lane
  and by the review: accepted at frame assembly, refused at objective
  build by the family's own pending-link stop, with and without
  `max_ndt`. The old note worried about an unfloored density inside a
  log-sum-exp, which is moot because the path stops earlier.
- **`gddm x predict`, `works` to `conditional`.** See the defect below.
  The clause the old row carried, that both `vint()` columns are
  mandatory on newdata, was measured and is false in both directions:
  dropping `cond` and dropping `upper` each leave the link prediction
  returning the same values, and the response path that would read them
  fails for every newdata. The row now says that rather than repeating
  it.

### Where the check runs, and why it can see what it needs

`frmtmb_register_frame_check()`, from `.onLoad()`. Three reasons.

- **It is the only seam that sees a design matrix.** The family's own
  `valid_y()` seam, which is where the `vreal()` check lives, is handed
  the response and the addition terms and nothing else; the predictors
  do not exist when it runs. `frmtmb/R/frame.R:2586-2589` says so in a
  comment beside the frame handed to the checks.
- **The precedent is in the repository twice.** `frmtmb.ode`'s
  `check_ode_constancy()` refuses the identical problem through the same
  seam, for a dynamics parameter that varies inside a `frm_ode()` solve
  group, and frmtmb.eam already registers `ddm_check_ndt_group_read()`
  there.
- **The family-side alternative does not exist for this family.**
  `frmtmb_structure(check_frame =)` is the per-family version and runs
  only for a family that declares a structure. `gddm()` is a rowwise
  `frmtmb_family()` with no structure, and giving it one to reach a
  check would be a much larger change than the check.

The check leaves before parsing anything when no response carries a
`gddm()` family, which is measured only as a `vapply` over the response
list; it runs on every frame frmtmb assembles once this package is
loaded.

### What it compares, and why not the design matrix

Each parameter's variables come from its parsed right-hand side, and
what is compared is the model-frame COLUMNS those variables build. Not
`X` and `Z`, which is the shape `check_ode_constancy()` uses, and the
reason is a place that shape fails open: a smooth and a `gp()` put their
basis in `Z`, but `mo()` leaves a column of MULTIPLIERS in `X` with the
ordered level codes beside it in `lp[["mo"]]`, so an `X`-and-`Z` check
passes a `mo()` term that varies. Measured by inspecting the linear
predictors of a gaussian model with each term type:

| term | where its data is |
|---|---|
| `x` | `X` column `x` |
| `s(x, k = 5)` | `X` column `s(x).fx1` plus a 60 x 3 `Z` |
| `(1 \| g)` | 60 x 6 `Z` |
| `offset(w)` | `lp[["offset"]]`, not in `X` at all |
| `mo(m)` | `X` column `mom`, a multiplier; the codes are in `lp[["mo"]]` |
| `gp(x, k = 5)` | 60 x 5 `Z` |

The right-hand side names the variable whatever the term does with it,
so the column route covers all six and any future term type. Under
`nl = TRUE` the primary parameter has no right-hand side at all and its
covariates sit in the nonlinear BODY, so the body is read too;
constructed both ways, `mu ~ a * exp(-b * x)` is refused naming `mu`
and `a ~ x` is refused naming `a`, and both are accepted once `x` is in
the index. It is
sufficient rather than exact, in the direction that refuses: a design
constant within a condition gives a parameter constant within it, which
is the same conservatism `gd_check_response()`'s `vreal()` check has.
It is also the more useful message, because the remedy takes variable
names.

### Where it disagrees with the `vreal()` precedent, and why

`gd_check_response()` compares the `vreal()` columns EXACTLY,
`length(unique(z)) == 1L`. This check does not, and the reason is a
measurement rather than a preference.

`poly()` orthogonalizes over the whole column, so two rows built from
BITWISE IDENTICAL inputs come back apart. On the sweep's own design,
`dev/gddm-scripts/gddm-sweep.R`, `poly(cohn, 2)` differs by **6.03e-15
between rows whose input is the same number**, on a column whose
largest entry is 0.06944: 8.68e-14 relative to that maximum. An exact
comparison therefore refuses `mu ~ poly(coh, 2)` on an index built from
`coh`, which is a correct model, and it did: that was the sweep's only
false alarm before the tolerance was added. The review reproduced
6.030e-15 independently and confirmed that `scale()`, `splines::bs()`,
`splines::ns()`, `log1p()`, `x / sd(x)` and `cut()` are all bitwise
row-deterministic, so `poly()` is the only source anybody has found.

(The ulp figure quoted in an earlier draft, 391, takes one ulp as
`.Machine$double.eps * colmax`. The true spacing at 0.06944 is
`2^-4 * eps = 1.388e-17`, so it is 434 ulp. The relative figure
8.68e-14 is the one the tolerance is set against and it needs no
convention.)

The difference between the two checks is the difference between the two
columns. A `vreal()` column is data as supplied and is never
recomputed, so an exact comparison there has nothing to trip on. A
design column is computed, and how much it moves depends on the term.
`?gddm` now says which of the two comparisons a reader is looking at.

### The tolerance, and the hole the first one had

The tolerance is

    tol = 1e-8 * max(|a|, |b|)  +  1e-11 * max|finite column|

per pair of entries, per column of a matrix. It is not the first form
shipped, and the first form failed open.

**What was wrong.** The first form was `1e-8 * max|finite column|`,
applied to every pair in the column. "Relative" there meant relative to
the column MAXIMUM, not to the two entries being compared, so a column
whose dynamic range exceeds 1e8 carried a blind band wider than its own
small entries. The review built one and it is a shape the field
actually writes: an intertemporal-choice delay of 1, 2 and 3 seconds in
the short blocks beside ten years, 3.15e8 s, in the long ones, with the
index built from the block and `delay` left on the right-hand side.

`dev/gddm-scripts/gddm-band.R`, seed 909, 120 rows in four blocks:

| | old rule, `1e-8 * colmax` | new rule, per pair plus a floor |
|---|---|---|
| tolerance on the pair (1 s, 2 s) | 3.15 s | 3.15e-3 s |
| within-condition spread it is compared against | 2 s | 2 s |
| conditions flagged | **0** | **2** |
| the same column CENTERED, conditions flagged | 2 | 2 |
| the model | **accepted** | refused, naming `mu` and `delay` |

Both rules are evaluated on the same column in one process at the end
of that script, so the flip is measured rather than cited. Two things
it shows. The old rule accepted the model, which is this item's own
defect intact inside its fix: on the reference build the same design
gives `fn(A)` and `fn(B)` bitwise equal at the same parameter vector
across data sets differing on 116 of 120 rows, and reports
`mu.delay = 3.22315e-10` fitted from 4 rows of 120. And the old rule's
answer FLIPPED when the column was centered, so it depended on whether
the user had centered; the new rule gives the same answer either way.

**Why these two constants.** Measured, not chosen. The review
registered a second frame check after this one, so it sees only frames
this check accepted, and recorded the largest relative deviation any
correct design produced from its own condition's first row across the
26 accepted gddm frames of the sweep:

| column | relative deviation |
|---|---|
| `poly(cohn, 2)[2]` | 1.025e-13 |
| `poly(cohn, 2)[1]` | 4.987e-14 |
| every other accepted comparison | exactly 0 |

So the whole requirement is `poly()`'s, and it is 1.03e-13. The first
term of the tolerance is scale-free and carries the ordinary entries;
the second is the floor the entries near zero need, where a per-pair
term goes to zero with them. `1e-11` of the column maximum keeps 97x
over the worst measured requirement by the review's own arithmetic, and
this lane's rerun of the sweep measures the headroom directly, as the
largest ratio of a real deviation to the tolerance it was compared
against:

| | value |
|---|---|
| worst real deviation / its own tolerance | 1.169e-05 |
| headroom, the reciprocal | **8.55e+04 x** |
| smallest within-condition step still refused | 1e-11 on a column of maximum 0.512, that is 2.0e-11 relative |
| largest step accepted | 1e-12, that is 2.0e-12 relative |
| the same two figures under the first shipped tolerance | 1e-08 and 1e-09, a thousand times coarser |

`poly()` is the shape that would produce a false alarm first, because
its entries near zero are where the per-pair term nearly vanishes and
only the floor is left, so it is scanned wider than the sweep's own
design. `dev/gddm-scripts/gddm-polyscan.R`, seed 202609, two arms of
the same family, both reading the floor out of the package rather than
retyping it:

| arm | worst dev / its own tolerance | headroom | cells flagged |
|---|---|---|---|
| 480 rows, integer levels, interleaved | 1.88e-04, at 8 levels degree 4 | 5.33e+03 x | **0 of 15** |
| the review's own: 120 rows, levels on (0, 1), blocked | 4.53e-05, at 6 levels degree 4 | 2.21e+04 x | **0 of 15** |

Both scan 4, 6, 8 and 12 levels by degrees 1 to 4. The second arm
reproduces the review's table to every digit it printed, 2.30e-05,
4.53e-05, 8.29e-06 and 9.94e-06 at 4, 6, 8 and 12 levels. It is carried
because two scans of one family that disagree mean one of them has not
been understood, and the first arm is harsher on three counts: four
times the rows, integer levels, and an interleaved rather than blocked
layout. Read the harsher one, 5.33e+03x, and nothing is flagged either
way.

**The band that is left.** A dpar that genuinely varies by less than
about 1e-8 of the larger of the two entries being compared, and by less
than 1e-11 of its column's largest entry, is still accepted. That is
the measured limit and it is what "sufficient rather than exact" costs.
The per-pair term does not widen with the column's range, so no choice
of units reaches it.

The floor is the absolute half, so the general statement is: **the
guard misses a within-condition difference `s` on a column of maximum
`M` once `M / s` passes about 1e11.** The review measured the flip on
the ten-year delay column at exactly `1e-11 * 3.15e8`, then looked for
a reachable case and could not find one in a covariate: the widest
meaningful ranges in this field are a delay of one second to ten years
(3.15e8), money from a cent to a million dollars (1e8), and an
epoch-millisecond onset against a one-second difference (1.76e9), all
two to three orders short. It reached the band only with a sentinel
1e12 times the data, and behavioral data uses 999, -999, 9999 and `NA`.
Reproduced at `dev/gddm-scripts/gddm-predicate.R` section 3: refused at
sentinels of 1e6, 1e9 and 1e11, accepted from 1e12 up.

The other side of that column is closed by frmtmb rather than by this
guard: a fixed-effect column whose total relative range falls below
about 1e-6 is dropped as rank deficient with a warning naming it, so
the band is five orders narrower than the narrowest range a column can
have and still be fitted at all.

### The form that would close the band, measured and declined

Apply the floor only to a COMPUTED column, one whose model-frame name
parses to a call. The argument is good: the floor exists for the
arithmetic the model frame runs when it evaluates `poly()`, and a
bare-symbol column is data as supplied, which is the same reason a
`vreal()` covariate is compared exactly. On a symbol column the
tolerance would collapse to eight significant digits of the two entries
themselves, which no dynamic range can widen. It is one predicate,
`is.call(str2lang(cn))`, on a parse the check already performs.

It was built, measured and taken out again.
`dev/gddm-scripts/gddm-predicate.R`, seed 4242.

**What it classifies, section 1.** Correctly, for everything reachable
through a formula. `poly()`, `scale()`, `bs()`, `ns()`, `log1p()`,
`cut()`, `I()` and `offset()` all keep the floor. A bare symbol,
`mo(ord)` whose model-frame column is the factor `ord`, and
`s(z, k = 4)` whose column is `z` because the basis lives in `Z`, do
not, and none of those three can need it: a factor is compared exactly,
and a basis is a deterministic function of the column that IS compared.

**What it buys, section 3.** The review's sentinel column, values 1 and
3 in one condition and the sentinel alone in another:

| sentinel | floor | with the floor | without it |
|---|---|---|---|
| 1e+06 | 1e-05 | refused | refused |
| 1e+09 | 1e-02 | refused | refused |
| 1e+11 | 1 | refused | refused |
| 1e+12 | 10 | **accepted** | refused |
| 1e+15 | 1e+04 | **accepted** | refused |
| 1e+18 | 1e+07 | **accepted** | refused |

It does close the band, at any dynamic range.

**What it costs, sections 2b and 2c, and this is what decided it.** A
column computed UPSTREAM and stored under a bare name carries exactly
the noise the floor exists for, and would lose it. Precomputing
`poly()` into a data frame column is one line and a real habit, and its
entries sit at zero whenever a level sits at the mean, which for a
three-level low, medium, high design is the middle one. `poly()`
precomputed at 3, 5, 7 and 9 symmetric levels by degrees 1 to 3, 21
columns:

| | with the floor | without it |
|---|---|---|
| correct designs accepted | **21 of 21** | 18 of 21 |
| false alarms | 0 | **3** |

The worst deviation the floor has to carry there is 2e-14 of the column
maximum, 90 ulp. Section 2b runs one of the three as a model rather
than on the helper: `mu ~ p1 + sc` on precomputed `poly()` and
`scale()` columns is accepted on the shipped rule and refused without
the floor.

**The tension is structural rather than a tuning problem.** The noise a
computed transform carries is scaled to the COLUMN, not to the entry,
so any tolerance that carries it is column-scaled, and any
column-scaled tolerance can be widened by one unrelated large entry.
The predicate resolves that by guessing from the NAME which columns
carry the noise, and the guess is wrong for a column computed upstream.
Splitting the constant instead, a smaller floor for bare columns, buys
about one order of magnitude of band before it meets the same 2e-14
requirement, at the price of a second constant nothing measures.

So the trade is a false alarm on 3 of 21 precomputed orthogonal
polynomials, a shape the field writes, against a silent acceptance
needing 1e12 of dynamic range, a shape the review hunted and could not
find. The lane rule that decides it is the one about false alarms: a
check that fires on a correct model is worse than no check. Rule 3
ranks a silent wrong answer above a refusal, and that orders which
defect to fix first; it does not license trading a reachable refusal
for an unreachable acceptance. The refusal that fires here also gives
WRONG advice: it names `p1` and says to put it in the index, and `p1`
takes as many distinct values as there are rows once the 1e-17 digit is
read, so following the message would ask for one solve per row.

Kept for a later round if anyone finds a principled floor. What such a
floor has to do is carry a 2e-14 relative deviation on a column whose
scale is set by an accumulating computation, without being inflated by
an unrelated entry. Nothing in this lane found one.

It departs from `ode_varying_cols()`'s `1e-8 * max(1, max(abs(M)))` in
three ways, each constructed as a test:

- **no floor at 1**, because a covariate measured in small units would
  otherwise get a tolerance orders of magnitude above its own scale. A
  column whose whole range is 1e-9 is refused here and would be accepted
  under that floor;
- **per column rather than per matrix**, for the same reason one step
  down: an intercept column of ones sets the whole-matrix maximum;
- **per pair rather than per column**, which is the blocker above.

A non-finite entry cannot open the guard: the scale is taken over the
finite entries and an infinite one is compared exactly. Without that,
one `Inf` anywhere in a column makes the tolerance `Inf` and nothing is
ever refused again. That is constructed in the test file rather than
reasoned about, because it is the exact shape of guard failure this
round has seen four times. The review confirms it is defense in depth
rather than a reachable path: frmtmb refuses a non-finite design column
before this check runs.

### A condition of one row

Accepted, always. One row cannot distinguish a constant parameter from a
varying one, and there is nothing to get wrong: the first row IS every
row, so the density reads exactly what the model says. Measured at the
degenerate end, 480 rows and 480 conditions with every parameter varying
between rows: accepted. `?gddm` says so.

## The false-alarm measurement

`dev/gddm-scripts/gddm-sweep.R`, seed 202609, 480 rows, six subjects,
four coherences, two blocks, two cues, run at `dry_run = "frame"` so the
check runs and no solve is taped. Two paired arms over the SAME designs,
because neither is readable alone:

- **ok**: the index is built from every variable in the model, which is
  what `?gddm` asks for. The check must not fire.
- **drop**: one variable is removed from the index and left in the
  model, which is the defect. The check must fire and must name that
  variable.

The 21 designs are the shapes a response-time paper writes: intercept
only; drift by coherence; drift by coherence with the boundary by block;
their interaction; non-decision time by subject; start point by cue; all
four parameters at once; a subject random intercept on the drift; a
subject random slope in coherence; a subject random intercept on the
boundary; a polynomial in coherence; a spline in coherence; a monotonic
ordered coherence; an offset; cell means with no intercept; a numeric
coherence; a collapsing boundary with the rate by block; a uniform start
width by block; a lapse rate by block; the coherence drift nonlinearity
with its `vreal()` covariate; and leaky integration with the leak by
block.

| arm | cases | result |
|---|---|---|
| ok | 21 designs, of which **20** are designs the check COULD fire on | **0 false alarms** |
| drop | 26 one-variable removals | **26 refusals, 0 missed, 0 that failed to name the dropped variable** |
| one row per condition, 480 conditions | 1 | accepted |
| gaussian, wiener, lba, rdm, wiener_gng with a covariate varying across every grouping | 5 | all accepted |

The one design the check cannot fire on is the intercept-only model, and
it is counted separately rather than folded into a rate, because a
refusal that cannot fire is not evidence that it does not false-alarm.

**Rerun after the tolerance was tightened, which is the change that
could have invalidated it.** Tightening a tolerance is exactly how a
guard acquires false alarms, so the whole sweep was rerun on the build
carrying `1e-8 * max(|a|, |b|) + 1e-11 * colmax`: **0 false alarms on
the same 21 designs, 26 of 26 on the drop arm, 0 missed, 0 unnamed.**
Every number in the table above is from that run. The review also
reproduced the previous sweep independently, on its own install rather
than this lane's, and got the same 0 and 26.

**The sweep's first version reported one miss, and the miss was the
sweep.** With `subj` cycled at period 6 and `blk` cycled at period 2,
`blk` is a function of `subj`, so dropping `blk` from an index that
still names `subj` changes nothing and nothing should be refused. The
three crossed factors are drawn rather than cycled now, and the script
prints the distinct block count within each subject so the aliasing
cannot come back unseen.

**Post-fit paths, on a correct model.** `dev/gddm-scripts/gddm-postfit.R`,
seed 11: `fitted()`, `residuals()`, `predict(type = "link")`,
`simulate()`, and `predict(newdata =, type = "link")` all run on both
arms, unchanged.

## The behavioral failure, seen

`tests/testthat/test-gddm-conditions.R` run against the reference build
of the base commit, failure cap lifted, one process:

**`PASS 13 | FAIL 17 | ERROR 5 | WARN 0 | SKIP 0`**
(`dev/gddm-scripts/gddm-testfile-before-log.txt`)

Of the 22, **19 are behavioral** and 3 are the weak form.

- All 17 failures are "Expected ... to throw a error" on a call that
  SUCCEEDS on 0.7.0: 14 `gc_frame()` models, one `frm()` under the
  `dec()` spelling, and two `gddm_simulate()` calls. That is a
  difference in the answer, not in the vocabulary.
- 2 of the 5 errors are `conditionMessage()` applied to `NULL`, which
  is what `expect_error()` hands back when the call did not throw, so
  they have the same behavioral root as the failures above them.
- 3 of the 5 errors are `could not find function`, on
  `gd_varying_groups` twice and `gd_col_floor` once: the weak form,
  from the blocks that construct the fail-open cases on the internal
  helpers because `frm()` cannot reach them.

The same file on this worktree: `PASS 58 | FAIL 0 | ERROR 0`.

An earlier draft recorded `PASS 14 | FAIL 13 | ERROR 2` against a
version of the file that had since grown, and the review caught it. A
count that does not reproduce is not a count; the numbers above are
from the current file.

## What was NOT done, and why

### The solver was not changed

Making `gd_densities()` group by distinct parameter vector rather than
by the given index would be a silent refinement of the condition count,
which is the 60x cost above with no statement to the user. It is also
not obviously stable across `predict(newdata = )`, which reassembles the
index from whatever rows it is handed, so the fitted model and the
predicted one could disagree about how many conditions there are.

### The `vreal()` check was not rewritten to match

It keeps its exact comparison and its own message. Two reasons: a
`vreal()` column is supplied data and does not carry recomputation
noise, so the tolerance would buy nothing there; and rewriting a
shipped refusal to say something new is a change to a message users may
already have in their notes. The two are separate on purpose and this
file records the disagreement rather than hiding it.

### `frm_compat("gddm")`'s `ndt_group()` note was left alone

Item 1.0a's reviewer called it too strong. Reading it now, it already
says the exclusion is on SCOPE rather than on impossibility and prices
it at one solve per subject, so it is accurate. What it does not say is
that the first-row read is now enforced rather than merely documented;
that is a one-line addition somebody should make when the round
consolidates, and it is not made here because the row belongs to
1.0a.

### `mixture(gddm(), ...)` is not covered, and cannot be reached

The check keys on `family[["gddm"]]`, which a mixture family does not
carry: `mixture()` builds a new family with renamed dpars. So a mixture
with a `gddm()` component would keep the defect. It cannot get there:
measured on both arms, with and without `max_ndt`, the model is accepted
at frame assembly and refused at objective build by the family's own
unresolved-link stop, because `mixture()` never runs its components'
`family_finalize`. If that ever changes, this check has to learn about
mixtures at the same time. The review reproduced the refusal by
building the mixture with the dpar names it wants, and
`frm_compat("gddm")`'s `mixture` row is corrected from `untested` to
`refused` with that measurement, because the lane exercised it.

### The one residual gap in "design constant implies parameter constant"

`mi()`. Two missing entries in one condition are two distinct latent
parameters, and the check treats two `NA`s in a column as one datum on
purpose, so a condition whose only variation is between two imputed
values is not seen. The review tried to build it: every mixed
construction is refused, because the observed rows differ from the
missing ones, and the fully unobserved case is degenerate. `mu ~ mi(x)`
with four latent values inside two conditions is refused. Recorded so
the next reader does not have to rediscover it.

## What `gddm_simulate()` did, and what it does now

The review found the same defect in an exported function and it is
fixed here rather than only disclosed.

`gddm_simulate()` recycles every parameter to `n` and then solves once
per distinct value of `coh`, reading each parameter at that value's
first trial. A per-trial parameter vector was accepted and half of it
ignored. `dev/gddm-scripts/gddm-simulate.R`, seed 31, 400 trials,
`mu = c(rep(-2.5, 200), rep(2.5, 200))`, grid `dt = 0.02, ny = 101`:

| | reference build 0.7.0 | this branch |
|---|---|---|
| `coh = 0`, upper rate first half / second half | 0.010 / 0.000 | refused by name |
| error or warning on that call | **none** | the refusal |
| `coh` separating the halves, upper rate | 0.010 / 0.990 | 0.010 / 0.990 |

Both halves of the first call are drawn from `mu = -2.5`. The second
call is the same 400 draws with the same parameters and a `coh` that
separates them, and it is 99 times the upper rate in the half the first
call threw away. The refusal names every offending parameter and says
what to do: give `coh` a distinct value per parameter setting, or call
the function once per setting.

It cannot false-alarm on a correct call, and the four shapes that must
keep working were run on both arms and agree in every digit: scalar
parameters; a length-`n` parameter that is constant within each `coh`;
a coherence drift; and `?gddm_simulate`'s own example. `?gddm_simulate`
now states the constraint on `...` and says `coh` is what separates
parameter settings whether or not a drift term reads it.

## A defect found and not fixed

**`predict(type = "response")` on newdata is broken for `gddm()`, and
so is `conditional_effects()`.** Measured on BOTH arms, so it is
pre-existing and this change neither caused it nor fixed it.
`dev/gddm-scripts/gddm-newdata.R` and `gddm-postfit.R`, seed 11, a
correct four-condition model:

| call | result, both arms |
|---|---|
| `fitted(fit)` | works |
| `predict(fit, type = "link")` | works |
| `predict(fit, newdata =, type = "link")` | works |
| `predict(fit, newdata =, type = "response")` | **error: invalid 'length' argument** |
| `conditional_effects(fit)` | **error: invalid 'length' argument** |

The failing call is `vector("list", 2L * d$ncond)` in `gd_densities()`:
the `.gddm` index the density reads is not rebuilt on the newdata path,
so `d[["ncond"]]` is absent and the length is `integer(0)`. It is a loud
failure rather than a silent wrong answer, so by Rule 3 it ranks below
this item and it is left. `frm_compat("gddm")`'s `predict` row said
"works" and is corrected to `conditional` here, with the measurement in
the note, because a compatibility table that says a broken path works is
itself a wrong answer to the user who consults it. That is the one edit
in this change that is outside the item; revert it if the round would
rather file it separately.

`influence()` on the same fit also fails, differently on the two arms:
"'list' object cannot be coerced to type 'double'" on this build, and on
the reference build the process ends without printing. Not investigated
further.

**The check does not run on the newdata path.** `run_frame_checks()` is
called from `assemble_frame()`, which `frm()`, `influence()`,
`get_prior()` and `frm_simulate()` reach and `predict(newdata = )` does
not. Constructed: a newdata frame with the index collapsed by hand is
accepted at `type = "link"` on both arms. The exposure today is nil,
because a link prediction never touches the density and a response
prediction on newdata fails for the reason above; it becomes real the
moment that defect is fixed, and whoever fixes it should run the check
there too.

## Tests and the check

Every file, one per process, `NOT_CRAN=true`, failure cap lifted,
`package = "frmtmb.eam"`, by `dev/gddm-scripts/gddm-runtest.R`. The
runner counts errors as well as failures, because a runner that sums
`failed` alone prints a clean line for a file that aborted halfway. The
BASELINE column is `dev/suite-baseline.tsv`, joined per file, which is
what that file exists for. Log at `dev/gddm-scripts/gddm-suite-log.txt`.

| file | baseline | this branch | |
|---|---|---|---|
| test-gddm-conditions.R | (new) | 58 | |
| test-bracket-access.R | 1 | 1 | |
| test-brms-parity.R | 13 | 13 | |
| test-defects.R | 59 | 59 | |
| test-density.R | 132 | 132 | |
| test-extension-api.R | 14 | 14 | |
| test-family.R | 51 | 51 | |
| test-gddm-family.R | 77 | 77 | |
| test-gddm-gradients.R | 22 | 22 | |
| test-gddm-recovery.R | 28 | 28 | |
| test-gddm-reference.R | 166 | 166 | |
| test-gddm-solver.R | 94 | 94 | |
| test-lba.R | 109 | 109 | |
| test-message-uniqueness.R | 4 | 4 | |
| test-moments.R | 29 | 29 | |
| test-ndt-bound.R | 79 | 79 | |
| test-rdm-gng.R | 221 | 221 | WARN 1, pre-existing |
| test-sampling.R | 97 | 97 | |
| test-scale.R | 0, SKIP 3 | 0, SKIP 3 | the gated tier |
| test-simulate-density.R | 69 | 69 | |
| test-surface.R | 50 | 50 | |
| test-units.R | 32 | 32 | |
| test-variability.R | 140 | 140 | |

**23 files, FAIL 0, ERROR 0, 1545 passing assertions**, which is the
baseline's 1487 plus the new file's 58. Every pre-existing file is at
its baseline count exactly.

Two totals are both right and a reader needs both. The per-file sum is
**1545 with SKIP 3**. The `R CMD check` whole-suite run, which is one
process rather than 23, reports
**`[ FAIL 0 | WARN 1 | SKIP 5 | PASS 1538 ]`**: two more blocks skip in
that environment, `test-message-uniqueness.R` because it needs the
package SOURCES which a built tarball does not carry, and
`test-sampling.R` because it needs `frmtmb.sample`. Neither is a
failure and neither number contradicts the other. The one WARN is the
pre-existing convergence warning `test-rdm-gng.R` also reports.

`R CMD check --as-cran` on the built tarball, once, on the final
pass, with pandoc and TinyTeX on PATH and
`_R_CHECK_CRAN_INCOMING_REMOTE_` off, by
`dev/gddm-scripts/gddm-check.ps1`: **Status: 1 NOTE**, and the NOTE
is the expected one, "Skipping checking math rendering: package
'V8' unavailable". `checking tests ... [12m] OK`, which is the whole
suite in one process, and `checking examples with --run-donttest`,
`re-building of vignette outputs` and `PDF version of manual` are
all OK. Logs at `dev/gddm-scripts/gddm-check-log.txt` and
`gddm-check/frmtmb.eam.Rcheck/00check.log`.

One check run in this lane is NOT a check result, and it is worth a
line because its shell reported exit 0. Its `00check.log` stopped at
"checking tests ..." with no `Status:` line and its `testthat.Rout`
was 750 bytes, so the test process died about twelve minutes in while
other R processes were being started beside it. A zero from the
wrapper is not a check result; the `Status:` line is. Every check
quoted here comes from a log carrying one.

`test-scale.R` is the gated tier and skips three, as the baseline
records. In the first round's suite run `test-variability.R` came back
`EXIT 127`, a shell "command not found" rather than a test result,
because another R process was starting at that moment; it was rerun on
its own at `PASS 140 | FAIL 0`. The punch-round run recorded above is a
single clean pass over all 23 files with every process exiting 0.

## What the plan should say

`dev/extension-gaps-plan.md` is not edited here. The unnumbered item its
Phase 1 prose describes should become **1.0e**, following the pattern
1.0d already sets: a new small item found by an earlier one, in the same
family of rows.

- **1.0e. DONE.** `gddm()` refuses a distributional parameter that is
  not constant within its condition, naming the parameter, the variable
  and a condition it varies inside. This was a silent wrong answer of
  the widest kind the family has: at a fitted parameter vector, adding
  100 to a covariate on 118 of 120 rows leaves the objective BITWISE
  identical for `mu`, `ndt`, `bs` and `bias` alike, with no error and no
  warning, and a random-effect grouping that crosses the index leaves
  the deviations of every subject not on a condition's first row at
  exactly zero while the fit reports a between-subject spread of
  7.6e-11. The check runs at frame assembly through
  `frmtmb_register_frame_check()`, which is the only seam that sees a
  design; it compares the model-frame columns each parameter is built
  from against the condition's first row, so it is sufficient rather
  than exact. 0 false alarms on 21 designs a response-time paper writes,
  20 of which the check could fire on, against 26 refusals on the 26
  ways of dropping one variable from those same indices. Refusing rather
  than reading was priced: reading means one Fokker-Planck solve per
  distinct parameter vector, which at 120 rows is 40.9x the tape build
  and about 52x the evaluation, and the exact count needs no clock
  because solves per evaluation IS the condition count. The comparison
  carries a measured tolerance,
  `1e-8 * max(|a|, |b|) + 1e-11 * max|column|`, with 8.55e+04x headroom
  over the largest deviation any correct design in the sweep produced
  and a band under it: a parameter varying by less than about 1e-11 of
  its column is still accepted, which is the limit of "sufficient
  rather than exact". The first form of that tolerance scaled by the
  column MAXIMUM and failed open on a column of wide dynamic range, and
  the review built one; the per-pair form closes it at any range and
  no longer depends on whether the user centered. What is LEFT is the
  floor: the guard misses a difference `s` on a column of maximum `M`
  once `M / s` passes about 1e11, which needs a sentinel 1e12 times the
  data to reach and which neither the review nor this lane could
  produce from a covariate the field writes. Charging the floor only to
  a computed column would close that completely and was measured and
  declined: it false-alarms on 3 of 21 precomputed orthogonal
  polynomials, which is a shape the field does write, and its refusal
  there gives advice that would ask for one solve per row. The plan row
  should carry the band as a stated limit, not as a defect.
  `gddm_simulate()` had the same first-row defect in an exported
  function, reading every parameter at the first trial of each `coh`
  level, and it is refused by name too: `mu = c(-2.5 x 200, +2.5 x 200)`
  at `coh = 0` drew both halves from -2.5 with no error and no warning.
  `frmtmb.eam/R/gddm.R`. About half a day, plus half a day of punch.
- The prose paragraph at the end of Phase 1 that files this as "a new
  small item, not yet numbered" should point at 1.0e instead.
- A small Phase 1 or Phase 5 row, not built: an opt-in
  `gddm_control(refine_conditions = TRUE)` that performs the refinement
  itself, reports the resulting condition count and the multiple it
  implies, and lets the user accept the cost. The 60x above rejects
  refining SILENTLY, not refining with the cost stated, and the two
  should not be conflated. What it needs studying for is what happens
  to a `predict(newdata = )` whose refinement differs from the fit's.
- Phase 4, the post-fit surface, should carry a new row:
  `predict(type = "response")` on newdata and `conditional_effects()`
  both fail for `gddm()` with "invalid 'length' argument" at
  `vector("list", 2L * d$ncond)`, because the `.gddm` index is not
  rebuilt on that path. Pre-existing on 0.7.0. Whoever fixes it must
  also run the frame checks on that path, which today they are not:
  `run_frame_checks()` is reached from `assemble_frame()` and
  `predict(newdata = )` does not go through it. Note for that row what
  was measured and what the old compat row claimed: NEITHER `vint()`
  column is mandatory on newdata today, because dropping `cond` and
  dropping `upper` each leave the link prediction returning the same
  values.
- Rule 2's shipping bar mentions the hazard-container guard for a new
  `R/` file. This item added no `R/` file, so `test-bracket-access.R`
  covers the new functions as it stands, and it does: `PASS 1 | FAIL 0`
  with `frm_hazard_reads("frmtmb.eam")` returning `character(0)`.

## The scripts

| script | what it establishes | log |
|---|---|---|
| `gddm-firstrow.R` | the mechanism on four dpars, bitwise, at a fitted vector; and the same five models refused on this build | `gddm-firstrow-before-log.txt`, `gddm-firstrow-after-log.txt` |
| `gddm-firstrow2.R` | the two probes the first script got wrong: an accidentally aligned `vreal()` control, and the marginal-likelihood random-effect probe | `gddm-firstrow2-before-log.txt` |
| `gddm-firstrow3.R` | why the joint parameter vector is not reachable through `obj$fn` | `gddm-firstrow3-before-log.txt` |
| `gddm-firstrow4.R` | the same, through `obj$env$f`, which refuses the length | `gddm-firstrow4-before-log.txt` |
| `gddm-firstrow-z.R` | the random-effect half, through the INNER solution: 2 of 4 deviations bitwise zero | `gddm-firstrow-z-before-log.txt` |
| `gddm-sweep.R` | the paired false-alarm and true-positive arms, and the tolerance at both ends | `gddm-sweep-after-log.txt` |
| `gddm-cost.R` | what one solve costs and what reading the varying dpar would cost | `gddm-cost-log.txt` |
| `gddm-postfit.R` | the post-fit paths on a correct model, both arms | `gddm-postfit-log.txt`, `gddm-postfit-before-log.txt` |
| `gddm-newdata.R` | the newdata defect, and that the check does not run there | `gddm-newdata-log.txt` |
| `gddm-runtest.R` | one test file per process, counting errors as well as failures | per run |
| `gddm-band.R` | punch round 1: the review's wide-range column, both tolerance rules in one process, both arms | `gddm-band-log.txt` |
| `gddm-simulate.R` | punch round 1: `gddm_simulate()` reading one row per coherence, both arms, and the four uses that must keep working | `gddm-simulate-log.txt` |
| `gddm-check.ps1` | the tarball and `R CMD check --as-cran`, once | `gddm-check-log.txt` |
| `gddm-predicate.R` | punch round 2: the computed-column predicate, what it classifies, what it buys and the 3 of 21 false alarms that declined it | `gddm-predicate-log.txt` |
| `gddm-polyscan.R` | punch round 2: `poly()` at 4 to 12 levels by degrees 1 to 4, two arms, the shape a tighter tolerance would false-alarm on first | `gddm-polyscan-log.txt` |
