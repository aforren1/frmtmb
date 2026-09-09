# Scale lane: how Phase 0 was run

Companion to `dev/scale-findings.md`. That file holds the numbers; this
one holds what a reader needs in order to distrust them correctly, and
what the next person to run the tier should know.

## What was added

| file | what |
|---|---|
| `extensions/<pkg>/tests/testthat/helper-scale.R` | the gate, the instrument and the record. BYTE-IDENTICAL in all seven packages |
| `extensions/<pkg>/tests/testthat/test-scale.R` | the designs, one file per package |

No `R/` file in any package was touched, so no new hazard-container
guard is needed: `frm_hazard_reads()` covers code, and this lane added
none.

`helper-scale.R` is copied rather than shared, and the reason is
mechanical rather than stylistic: a testthat helper lives under
`tests/testthat/`, which is NOT installed, so an installed frmtmb.eam
cannot read a helper that lives in frmtmb. Sharing it would mean
putting it in core's `inst/` and reaching it with `system.file()`,
which is a public promise about a test helper for no benefit.

The precedent this file first cited, `helper-sampling.R`'s "a
three-line expectation is cheaper to duplicate than to make into API",
is weaker than it was made to look: that is two copies of three lines,
and the same header two lines above documents the OPPOSITE move,
`sampler_gates_on()` "left rather than being duplicated" once its whole
user base was in one package. This is seven copies of 197 lines.

The drift risk is not hypothetical here. `test-message-uniqueness.R` is
also seven copies, and they are seven DIFFERENT files: seven distinct
hashes at 101, 114, 80, 80, 86, 72 and 114 lines. So the copies carry a
guard: `.github/workflows/scale-helper-identical.yaml` hashes all seven
on any change to one and fails if they differ. It is repository-level
because the comparison is impossible from inside any one package.

That guard asserts the COUNT as well as the hashes, and the count
assertion is not decoration. Without it the job passes when a copy is
DELETED, because one surviving file still hashes to one distinct value.
That is the empty-sweep shape two other lanes hit this round, and it
was in the first version of this workflow. All four cases were run
against the workflow's own shell body before it shipped:

| tree | copies | distinct hashes | result |
|---|---|---|---|
| as shipped | 7 | 1 | pass |
| one copy edited | 7 | 2 | fail, "has drifted" |
| one copy DELETED | 6 | 1 | fail, "expected 7 copies, found 6" |
| one copy added | 8 | 1 | fail, "expected 7 copies, found 8" |

The deletion case passed before the count assertion and fails after it,
which is the only evidence that makes the line worth having.

## Running it

    $env:NOT_CRAN = "true"
    $env:FRMTMB_SCALE_TESTS = "true"
    $env:FRMTMB_SCALE_OUT = "<a file to append to>"
    $env:FRMTMB_SCALE_ROW = "eam"          # one row, in its own process
    Rscript --vanilla -e "testthat::test_file('.../test-scale.R')"

The scratchpad scripts that did this are
`scale-run-row.ps1` (one row, one fresh process, with the process peak
polled from outside) and `scale-queue.ps1` (every row, SERIAL).

`FRMTMB_SCALE_SMALL=true` shrinks every design to a size that exercises
the code path in about a minute per package, ode excepted at about two.
Every row it records carries `small=TRUE`, because nothing measured at
a toy size is a measurement.

**It is a plumbing check and not a smoke test for the result**, and the
difference matters. At the small size the `eam` and `eam-unbounded`
arms return the SAME log-likelihood, -77.9875, and the same estimates,
because at four subjects the `ndt` component COLLAPSES to 5e-05 and
2e-05 on the link instead of running away. The whole finding this tier
exists to record is invisible in small mode, and the assertions differ
there for the same reason: the eam row asserts recovery at the small
size and the defect at the realistic one, and the coupling row's
recovery assertion is switched off at the small size, where 6 subjects
by 15 frequencies gives 0.212 (-0.018, 0.441) for want of data rather
than because anything is wrong.
It exists because the tier's own code has to be debuggable without
paying for the measurement, and it caught five defects before the timed
run:

- a non-conformable profile matrix in the lca row, because
  `lca_profiles()` returns one table per ITEM and the row assumed one
  per class;
- a refused initial distribution in the hmm row, because the default
  stationary one cannot be used when a transition carries a predictor;
- a tape-build subtraction that came out NEGATIVE, because the first
  `frm()` call in a process pays one-time costs that belong to neither
  arm. Two interleaved rounds with the minimum of each fixed it;
- a recovery assertion on three coupling models that do not contain the
  truth, since three rungs of the ladder omit random effects the
  simulator has;
- variance components read by POSITION out of `VarCorr()`. An mgcv
  smooth is a random-effect block too, so position 1 is a different
  block from one coupling rung to the next. Every row now records its
  components by name.

## Three things about the instrument

**`proc.time()` is not usable here and is not used.** It ticks at
10.0 ms on this machine. `Sys.time()` resolves to about 2 us, so every
elapsed time in the tier is a `Sys.time()` difference.

**A gradient is timed over a block that is grown, not over one call.**
`scale_grad()` doubles the batch until the block passes 1.2 s or the
batch reaches 256 calls, then reports the best of three such blocks. A
gradient that is already seconds long stays at one call per block,
which is the right answer for it.

**Every row carries a control.** `scale_control()` runs two blocks of
the SAME work through the same path and reports `max/min`. It must be
near 1. It is in the recorded table for every row, and it is how a
reader can tell a row whose wall clock was disturbed from one that was
not. The rule that produced it is the 0.55.0 round's finding that two
timing claims had measured below the clock's tick.

## What contaminates these numbers, and by how much

The standing rules say other lanes run on this machine at the same
time. The runner records `meta_R_procs_before` and
`meta_R_procs_after` for every row, the count of R processes alive at
the row's start and end with this lane's own included. Nothing there
can see a compiler or a browser, and the first pass counted only
`Rterm` and so MISSED sibling lanes that had launched through
`Rscript`; the second pass counts both. Read the column as a lower
bound on how busy the machine was.

The queue is SERIAL on purpose. Two timed rows at once would make each
row's wall clock a measurement of the other. What the column is
actually for is visible on `coupling-coh-full`, which took 11.5 s in
one pass and 30.6 s in the other and produced identical estimates to
every digit. That is why every number in `dev/scale-findings.md` is a
minimum over passes: contention can only make a measurement longer.

The machine is 16 logical cores and 32 GB, which matters for the
memory column more than for the clock: the rlddm row peaks at 8.0 GB
and the eam rows at about 4.8 GB, so three of those at once would be a
memory problem before it was a scheduling one.

## Memory: two numbers, neither of them the whole story

`mem_mb` in the table is `gc()`'s peak R heap since a reset at the top
of the row. It does NOT see RTMB's tape, which is C++ memory, so it is
a floor.

`meta_peak_ws_mb` is the process peak working set, polled twice a
second from the parent PowerShell process while the child runs. It sees
everything, including R's own startup, the loaded namespaces and the
tape. Neither number is the fit's marginal cost; the difference between
them is roughly what lives outside the R heap.

The runner launches `Rterm.exe` and not `Rscript.exe`, because
`Rscript` is a launcher that starts `Rterm` as a CHILD: polling the
`Rscript` handle reported 7.2 MB for a fit that peaked at 400 MB.

## What was NOT run, and why

The plan's "Realistic scale" table has a second spline line, the
Royston and Parmar design at 2000 subjects with a frailty. The Phase 0
table's spline row names only the curves design, so that is what this
tier runs. The survival design belongs to Phase 2's item 2.5.

## The search that stopped one directory short

The first draft of these notes said the survey's four coupling
benchmarks "could not be recovered". They can. The scripts are gone,
but the constructions are recorded in the survey's own session
transcript under
`C:/Users/adf44/.claude/projects/c--Users-adf44-source-r-frmtmb/`, and
that is where the plan's four numbers came from in the first place.

The search that failed covered the repository, `dev/`, and every
session scratchpad on the machine. It did not cover the transcripts.
The standing rule says "search `dev/` and `dev/reviews/` for the
construction before concluding a record is wrong"; the transcripts are
the third place and they are the place a number quoted in a plan
document is most likely to live, because a plan is written from one.
Add them to the list.

The consequence was not only a wrong sentence. It sent the plan's item
2.6 note to "not reproducible and no claim rests on it", when the
model behind the survey's 0.77 was a specification the ladder did not
carry. It does now, as the `coh-id` rung.

## Verification

**The tier is green under its own gate in BOTH modes.** That was not
true of the first version: the full-size `eam` row failed on purpose,
which is a bad idea because a shipped assertion that fails whenever its
gate is on trains a reader to expect red, and the small-size
`coupling-coh-full` row failed by accident, which went unreported. Both
are fixed, the eam one by asserting the DEFECT so the file turns red on
the day the defect is fixed, the coupling one by gating its recovery
assertion on the design having enough power for it.

| mode | rows | fail | error |
|---|---|---|---|
| `FRMTMB_SCALE_TESTS=true`, full size | 15 | 0 | 0 |
| the same with `FRMTMB_SCALE_SMALL=true` | 15 | 0 | 0 |

**The gate itself.** Every package's `test-scale.R` was run twice with
the tier's variable UNSET, once with `NOT_CRAN=true` and once without.
Fourteen tests skip, none passes, none fails, none errors:

| package | tests | with NOT_CRAN=true | without |
|---|---|---|---|
| frmtmb.eam | 3 | 3 skipped | 3 skipped |
| frmtmb.learn | 2 | 2 skipped | 2 skipped |
| frmtmb.latent | 2 | 2 skipped | 2 skipped |
| frmtmb.ode | 1 | 1 skipped | 1 skipped |
| frmtmb.spline | 1 | 1 skipped | 1 skipped |
| frmtmb.coupling | 5 | 5 skipped | 5 skipped |
| frmtmb.sample | 1 | 1 skipped | 1 skipped |
| total | 15 | 15 skipped | 15 skipped |

Under the gate the same files run every row in this document.

**The public-surface substitution was checked and costs nothing.** The
`unbounded_dpar` measurement in `dev/scale-findings.md` takes its
estimate from `fixef()` and its standard error from the width of
`confint()`'s Wald interval, rather than reaching into core for the
quantity the check itself reads. The reviewer ran both routes side by
side in one process on each eam arm and the difference was exactly 0 on
both the estimate and the standard error, so the table is right to
every printed digit and nothing in it depends on a `:::`.

**`R CMD check --as-cran`**, with the pandoc and TinyTeX directories on
PATH as the standing rules require:

| package | status | tests |
|---|---|---|
| frmtmb.eam | 1 WARNING, 1 NOTE | 324 s, OK |
| frmtmb.learn | 1 WARNING | 68 s, OK |
| frmtmb.latent | 1 WARNING, 1 NOTE | 41 s, OK |
| frmtmb.ode | 1 WARNING | OK |
| frmtmb.spline | 1 WARNING, 1 NOTE | 20 s, OK |
| frmtmb.coupling | 1 WARNING | 14 s, OK |
| frmtmb.sample | 1 WARNING | 21 s, OK |

The NOTE is the expected one: "checking HTML version of manual ...
Skipping checking math rendering: package 'V8' unavailable".

The WARNING is `checking CRAN incoming feasibility`, and it is not
this lane's. It reads "New submission", "Strong dependencies not in the
CRAN or BioC software repositories: frmtmb", and a 301 on the pkgdown
URL in DESCRIPTION. To prove rather than argue that, frmtmb.coupling
was copied WITHOUT the two files this lane added and checked the same
way: `Status: 1 WARNING`, the same one. The lane adds files under
`tests/testthat` and `dev/`, and `dev` is in every package's
`.Rbuildignore`, so nothing it adds reaches the metadata the warning is
about.

**Every test file rerun ONE PER R PROCESS**, `NOT_CRAN=true` and the
scale gate unset, which is the standing rule because a whole-suite run
in one process has repeatedly hidden state leakage here:

| package | files | pass | fail | error | skip | warn |
|---|---|---|---|---|---|---|
| frmtmb.eam | 20 | 1376 | 0 | 0 | 3 | 1 |
| frmtmb.learn | 10 | 266 | 0 | 0 | 11 | 0 |
| frmtmb.latent | 6 | 228 | 0 | 0 | 2 | 0 |
| frmtmb.ode | 6 | 212 | 0 | 0 | 1 | 0 |
| frmtmb.spline | 10 | 333 | 0 | 0 | 1 | 0 |
| frmtmb.coupling | 7 | 393 | 0 | 0 | 5 | 4 |
| frmtmb.sample | 14 | 983 | 0 | 0 | 3 | 0 |
| total | 73 | 3791 | 0 | 0 | 26 | 5 |

Fifteen of the 26 skips are this tier behind its gate; the rest are
each package's own gated tiers. The five warnings are in
`test-rdm-gng.R`, `test-coherence.R` and `test-surface.R` and are not
this lane's: none of those files was touched.

One thing about that runner, because it produced a false alarm first.
`testthat::test_file()` must be given `package = "<pkg>"`, or the test
environment does not descend from the package namespace and every test
that reaches an internal function errors. Without it this runner
reported 20 errors across frmtmb.eam that do not exist: `R CMD check`
was green on the same files at the same moment, which is what said the
runner and not the suite was wrong.

## What the next person should know about the assertions

**Four rows assert interval coverage on ONE seed.** That is a 1-in-20
failure by construction even when the model is correct, and the
coupling row has already shown what it does at a size where the
interval is wide. It is a deliberate trade: this is a cost tier and
Phase 2 is where coverage becomes a rate over 60 replicates. If a row
starts failing intermittently, suspect the seed before the model.

**No row asserts on a status code any more.** The two latent rows used
to assert `identical(fit$opt$convergence, 0L)`. `nlminb` returns 1 on
perfectly good fits, as this tier's own ode row does at a maximum
gradient of 0.0061 with a positive definite Hessian, and the code can
differ with a different BLAS. They now assert a positive definite
Hessian and a gradient scaled by the log-likelihood the run itself
reached, which is a ratio to something measured rather than a code.

**The eam row asserts the DEFECT at the realistic size.** Until item
1.0a lands, the design cannot converge, so the row asserts
non-convergence and NaN intervals rather than recovery. The recovery
assertion is in the file, commented, directly beneath. The file turns
red the day the bound is fixed, which is the day someone should come
and swap them back. A tier that ships a failing assertion trains
readers to expect red.

## One number is not a measurement

Three claims in the first round of this lane came from a single
unreplicated pair, and one of them was wrong by 28 percent: "ten times
the draws cost 12.75 times the time" was `1.59 / 0.1249`, two
single timings, where three interleaved rounds give 9.985x for a draw
ratio of exactly 10.

`scale_interleave()` exists so that the tier cannot make that mistake
again. It takes a named list of zero-argument functions, runs one of
each per round in order, and reports the minimum per arm with that
arm's round-to-round spread beside it. Every post-fit timing in the
hmm, lca, spline and sample rows now goes through it, and the recorded
row carries `rounds` and the spread. The whole-fit timings do not,
because a fit is minutes; there the `passes` column carries the
replication instead.

Where this document or `dev/scale-findings.md` reports a RATIO, it
says how many rounds and what the spread or the control was. Where it
reports a single fit's wall clock, it says how many passes.

## What `scale_interleave()` cannot see

It reports a `first` per ARM, and that is first in its own arm, not
first on the fit. If an earlier arm warms a cache, every later arm's
`first` is a warm number wearing a cold label, and the spread will not
show it either, because after round 1 every round is warm and the
spread is computed over all of them.

That is not hypothetical: it is how this tier first reported
`frm_curve_feature()` at 0.0096 s when its cold cost is 0.398 s, a
factor of 41. The feature is arm 3 and the two band arms had already
paid the memoized joint-precision solve. The spline row now measures
that solve explicitly, on a cold fit, before the interleave starts.

Every other interleaved row was checked for the same shape and none
has it. Measured by calling each arm first on its own fresh fit and
then again:

| row | arm | first on a fresh fit | second call | first when the OTHER arm ran first |
|---|---|---|---|---|
| hmm | `hmm_probs()` | 0.0765 | 0.0783 | 0.0703 |
| hmm | `hmm_viterbi()` | 0.0493 | 0.0592 | 0.0422 |
| lca | `lca_probs()` | 0.0094 | 0.0095 | |
| lca | `lca_profiles()` | 0.0062 | 0.0050 | 0.0030 |

No first-to-second step anywhere, and the order of the arms does not
move the numbers, so neither latent row has the shape.

The sample row was NOT measured this way and the claim for it is
weaker. What is known is that its per-arm spreads move in both
directions across runs rather than falling monotonically after round 1
(the epred ratio came out 10.50 in one run and 11.72 in the next), and
that the reviewer checked every row for a shared cache and found only
the spline one. A cold/warm split there would not change the row's
conclusion, since the whole post-fit surface is 2.2 percent of the
sampling either way, but it has not been ruled out by measurement here.

So the spline row is the only one where this is known to have
mattered, and it is fixed there.

The general rule for anyone reusing the helper: if two arms can share
a cached quantity, measure the cold one OUTSIDE the interleave, first,
on a fit nothing else has touched.

## The minima have not converged

Every timing in `dev/scale-findings.md` is a minimum over one to four
passes, on the argument that contention can only make a wall clock
longer. That argument is sound and the minimum is still the right
statistic, but it should not be read as the true cost.

The reviewer reran three rows on a quiet machine and got numbers 5 to
8 percent BELOW the stated minima. Then the final verification pass of
this lane, one more run of every row, moved several of them again:

| row | minimum before | minimum after | passes |
|---|---|---|---|
| spline | 6.02 s | 3.83 s | 6 |
| eam | 309 s | 269 s | 4 |
| coupling, intercept rung | 0.150 s | 0.110 s | 4 |
| sample | 120 s | 119 s | 4 |
| learn-rlddm | 738 s | 738 s | 3 |

At one to six passes each, on a machine that had sibling lanes on it
for most of the queue, the minima are still descending. Treat every
wall clock in that table as an upper bound, and one that is probably
several percent high, not as a converged figure. The spline row moved
by 36 percent between its fourth and sixth pass.

This is a limit on the whole table rather than a defect in one row,
and it is the reason this lane stopped leaning on `learn-rlddm`'s
clock. That row is 738 s against a 600 s rule, a margin of 1.23x, and
a margin that thin cannot survive a statistic that is still moving by
5 to 8 percent between measurements. Item 1.0b moves with item 1.0a
under Rule 3 instead, which does not depend on a clock at all.
