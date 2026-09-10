# Lane `latent`: items 2.3 and 2.4 of the extension gaps plan

Date: 2026-09-09 into 2026-09-10. Worktree `frmtmb-wt-latent`, branch
`wt-latent`, off `f9ee297` (frmtmb 0.55.2, frmtmb.latent 0.2.2).
R 4.6.1, Windows 11, 16 logical cores, 32 GB. Third-party references:
depmixS4 1.5.4, hmmTMB 1.1.2, poLCA 1.6.0.2.

Libraries, in the order `dev/latent-env.R` sets them: this lane's own
`C:/Users/adf44/source/r/latent-lib` (holds only frmtmb.latent, built
from this worktree), the round's shared reference build
`rellib-0552`, the StanHeaders pin `pinlib`, then the user library.
Nothing in this lane touches Stan, so the pin was checked and not
otherwise exercised.

## The short version

1. **Item 2.4's "none expected" was wrong, and the defect is now
   fixed.** At the plan's own realistic scale, `lca()`'s old
   deterministic start reached a local optimum on **8 of 200** replicate
   data sets, **243 to 284 log-likelihood units** below the one
   `poLCA(nrep = 10)` found on the same data, and seven of the eight
   passed the pair this project's scale tier asserts. The cause was
   diagnosable: the start scored subjects by the MEAN item code, and on
   a design whose classes differ in WHICH items they endorse that score
   explains **0.08 percent** of its own variance between the true
   classes. `lca_init_extras()` now slices on the response PATTERN and
   shrinks the class profiles nine tenths of the way to the pooled
   ones, and reaches poLCA's optimum on **200 of 200 twice**: on the
   200 seeds the shrink weight was tuned on, and on 200 fresh ones it
   was not. Same cost per fit.
2. **The remedy `?lca` prints works, and the one `hmm_starts()` uses
   does not.** Ten refits from PERTURBED DEFAULT STARTS, which is the
   recipe on the help page, reached poLCA's optimum on all eight of the
   old start's failures in 6 to 17 seconds. Ten refits jittered two
   standard errors around the LOCAL OPTIMUM left 72 to 273
   log-likelihood units on the table. That is the single most important
   thing this lane learned about the `frm_multistart()` core seam, and
   it survives the fix: the seam as filed would not have fixed this.
3. **The 8.1-unit probe reproduces exactly** against the shipped
   `hmm()` family: -1096.09575602 cold against -1087.99646521,
   8.09929 units, with `convergence == 0`, a positive definite
   Hessian, `max|grad|` 7e-4 and `diagnose()` silent. `frm_allfit()`
   agrees on the wrong answer across all four optimizers with a
   log-likelihood spread of 7.8e-07.
4. **`hmm_starts()` recovers from it**, on 5 of 5 seeds at `jitter` 1,
   2, 4 and 8 with `n = 8`, and on 1 of 5 at `jitter = 0.5`.
5. **`hmm_starts()` told a correct fit it had found a local optimum, on
   6 of 6 unimodal chains, and now does not.** Found in review. Three
   places asked "is this a different optimum" and answered it three
   ways; two carried no tolerance at all. One definition now serves all
   three, measured at 0 of 6 false alarms with the 8.099-unit detection
   intact at a 7389x margin.
6. **hmmTMB's `initial_state = "estimated"` is a different model from
   frmtmb's `init = "estimated"`**: one initial distribution per
   SEQUENCE against one shared. Left alone the two log-likelihoods sat
   **19.678** apart on 20 sequences, from 38 extra parameters, and
   nothing said why. Pinned on both sides, they agree to 1e-11
   relative. (A first draft of `?hmm` attached 21.231 to that sentence,
   which is the same hmmTMB fit against `init = "uniform"`, a
   40-parameter difference. The 1.553 between them is frmtmb's own two
   initial-distribution parameters. Corrected in punch round 1.)

## Punch round 3: what the third review changed

### R3-BLOCKER. `?lca` pointed at a remedy that does not work

Round two ended the coverage paragraph with "prefer
`confint(method = "profile")` on such a term". The reviewer measured it
and it does nothing:

| | class3:x2 | class4:x2 | pooled |
|---|---|---|---|
| Wald coverage on the subset | 51 / 59 | 52 / 59 | 103 / 118 = 87.3% |
| profile coverage, same 59 | 51 / 59 | 52 / 59 | 103 / 118 = 87.3% |
| median width ratio | 1.001 | 1.001 | 1.001 |

Profile is wider on 118 of 118 intervals and **not one of them changed
a verdict**, at about 4 s each. And it reaches only the 59 replicates of
200 where the scored quantity is a raw parameter: for the other 141 it
is a CONTRAST, `theta_{perm[a]} - theta_{perm[1]}`, which
`confint(method = "profile")` cannot interval at all.

**The diagnosis and the remedy were not the same shape.** A profile
interval widens a Wald one when the profile log-likelihood is not
quadratic, which is a CURVATURE problem. What is wrong here is the
MAGNITUDE of the standard error. Recomputed on this lane's own files
over all 400 replicates (`dev/latent-2p4-pooled.R`), with `z` the error
over the reported standard error:

| coefficient | sd(z) | IQR(z)/1.349 | ratio | kurtosis |
|---|---|---|---|---|
| class2:x2 | 1.0040 | 0.9478 | 1.059 | 2.938 |
| class3:x2 | 1.0600 | 0.9997 | 1.060 | 2.816 |
| class4:x2 | 1.0800 | 1.0332 | 1.045 | 3.100 |

`sd(z)` is 1 when the reported standard error matches the spread. It is
1.06 and 1.08 on the two that under-cover, the robust spread agrees to
within 6 percent, and the kurtosis is below or at 3. No tail, no
curvature: the whole distribution is inflated. Profiling cannot reach
that.

**The pointer is gone.** `?lca` now says the interval is 3 to 4 points
narrow on those two coefficients and that nothing measured here fixes
it, and names a bootstrap interval as the obvious untried candidate
while saying it is untried. An honest gap beats a remedy that does not
remedy.

### The pooled number, which round two did not compute

Round two attributed the under-coverage to the fresh block. That was
the wrong reading. At 200 replicates the Monte Carlo half-width on one
coefficient is about 3 points, so neither block alone can separate a 91
percent row from 95: on the tuning block the same two coefficients are
its two lowest rows, at 93.0 and 92.5, and neither is distinguishable
from nominal on its own. Pooled over all 400
(`dev/latent-2p4-pooled.R`, reading the raw `err_` and `se_` columns of
`dev/latent-2p4-recheck.tsv` and `-recheck-oos.tsv`, nothing refitted):

| coefficient | covered | rate | 95 percent interval | se/sd |
|---|---|---|---|---|
| class2:(Intercept) | 380 / 400 | 95.00% | 92.86 to 97.14 | 0.955 |
| class2:x1 | 378 / 400 | 94.50% | 92.27 to 96.73 | 0.984 |
| class2:x2 | 379 / 400 | 94.75% | 92.56 to 96.94 | 0.994 |
| class3:(Intercept) | 376 / 400 | 94.00% | 91.67 to 96.33 | 0.956 |
| class3:x1 | 384 / 400 | 96.00% | 94.08 to 97.92 | 1.014 |
| **class3:x2** | 367 / 400 | **91.75%** | **89.05 to 94.45** | 0.940 |
| class4:(Intercept) | 376 / 400 | 94.00% | 91.67 to 96.33 | 0.962 |
| class4:x1 | 377 / 400 | 94.25% | 91.97 to 96.53 | 0.970 |
| **class4:x2** | 363 / 400 | **90.75%** | **87.91 to 93.59** | 0.922 |
| **overall** | 3380 / 3600 | **93.89%** | **93.11 to 94.67** | |

Exactly two rows and the overall figure have intervals that exclude 95;
every other row contains it. So it is a property of the DESIGN and not
of one draw. And it is **two of the three `x2` slopes, not all three**:
`class2:x2` pools to 94.75 percent and is fine, which rules out "the
binary covariate" as the whole story and leaves the two classes whose
`x2` coefficient is largest in absolute value (0.9 and 0.4 against
-0.5) as the ones affected.

### R3-NIT. Nothing pinned the shrink weight

The constant that was round two's blocker, and that moved twice,
appeared once, at `R/lca.R:304`, with no test mentioning it. An edit in
either direction left the suite green.

`test-lca.R` gains `the starting rule's shrink weight is pinned`, on
seed **20270476**, which discriminates in BOTH directions. From this
lane's own `dev/latent-2p4-oos.tsv` (poLCA start seed 930076):

| shrink weight | logLik at seed 20270476 |
|---|---|
| 0 | -11763.41 |
| 0.25 | -11749.43 |
| 0.5 | -11764.04 |
| 0.75, **0.9**, 0.95 | **-11507.6035576**, which is poLCA's |
| 0.99 | -11763.41 |
| 1 | -12864.24 |

The block asserts the optimum relative to its own magnitude and a
clearance of more than 200 units over the losing weights, also as a
ratio. One fit, 1.3 s.

**Seen failing, on both edges.** `dev/latent-falsealarm-check.ps1`
gained two more modes, `w05` and `w099`, which move the constant to
each edge of the usable interval on a copy. Both fail the pin with 2
failures each, at a relative gap of 0.0218 and 0.0217 from the optimum.
The full harness now runs four reverts:

| revert | file | result |
|---|---|---|
| `notol` | test-hmm-starts.R | 7 failures |
| `refs` | test-hmm-starts.R | 1 failure, `grad_tol = 0.086` |
| `w05` | test-lca.R | 2 failures |
| `w099` | test-lca.R | 2 failures |

## Punch round 2: what the second review changed

### R2-BLOCKER. The shrink weight was tuned and scored on one set

Round one picked `w = 0.5` as "the smallest weight that loses none" on
the 200 seeds the sweep ran on. That is a constant chosen as the unique
zero of a loss curve on its own tuning set, and it did not survive a
fresh draw. Re-derived here independently, on this lane's own runs, on
**two blocks of 200 that are never pooled**: the tuning block 20260910
to 20261109 (`dev/latent-2p4-shrink.tsv`, `-shrink2.tsv`) and a fresh
block 20270401 to 20270600 with its own `poLCA(nrep = 10)` arm at start
seeds 930001 upward (`dev/latent-2p4-oos.tsv`). Summarised by
`dev/latent-2p4-weight-summarize.R`.

**The reference is adjudicated, not taken on trust.** A seed's
reference is the largest log-likelihood ANYTHING reached on it: poLCA's
best of ten, or any weight in the grid. Taking poLCA on trust would
flatter every start that agrees with it. As it happens no start beat
poLCA on any of the 400, so poLCA is the reference everywhere, and that
is worth knowing rather than assuming.

| w | tuning block, 200 | fresh block, 200 |
|---|---|---|
| the 0.2.2 score cut | 8 lost | not run |
| 0 | 6 | 7 |
| 0.1 | 5 | not run |
| 0.25 | 2 | 1 |
| 0.5, what round one shipped | **0** | **1** |
| 0.75 | **1** | 0 |
| **0.9, what ships now** | **0** | **0** |
| 0.95 | 0 | 0 |
| 0.99 | 0 | **1** |
| 1 | 199 | 200 |

**Only 0.9 and 0.95 lose nothing on either block.** The choice between
them is which has more room on the side that fails catastrophically:
0.9 is 0.09 below 0.99, where a loss reappears, against 0.95's 0.04,
and 0.15 above 0.75. **0.9 ships.**

Three things this measurement adds to the review's:

- **0.99 loses a seed out of sample.** The review had it losing 1 of 17
  hard seeds; on 200 fresh ones it loses 20270476. The upper edge is
  therefore between 0.95 and 0.99 rather than at 0.99 itself, which is
  what makes 0.9 preferable to 0.95 rather than arbitrary between them.
- **The one hard seed is shared by both edges.** 20270476 is lost by
  `w = 0`, 0.25, 0.5 AND 0.99, and survived only by 0.75, 0.9 and 0.95.
  A single data set separating the interior from both boundaries is the
  strongest form the argument for the interior can take here.
- **`w = 1` is the symmetry axis and it is total.** 199 of the tuning
  block and 200 of the fresh one, by up to **1557.9** log-likelihood
  units. Every class profile is the pooled profile there, so the class
  spread in every item logit is exactly 0.

The residual after the change: nothing on either block. The deepest
loss the shipped weight now has anywhere in these 400 data sets is
none; the deepest the 0.2.2 rule had was 283.8 and the deepest `w = 0.5`
had was 256.4, on the seed the review found.

### What the weight change invalidated, re-run

`dev/latent-2p4-recheck.R` at `w = 0.9`, on BOTH blocks, with the poLCA
column read from the sweeps rather than recomputed:

| | tuning block | fresh block |
|---|---|---|
| reaches poLCA's optimum | 200 of 200 | 200 of 200 |
| log-likelihoods agree to, relative | 1.478e-14 | 4.372e-11 |
| Wald coverage, 9 coefficients x 200 | 94.3% (1698/1800) | **93.4%** (1682/1800) |
| its Monte Carlo interval | 93.3 to 95.4 | **92.3 to 94.6** |
| positive definite Hessian | 200 of 200 | 200 of 200 |
| largest bias | 0.012 | 0.032 |
| median fit seconds | 0.330 | 0.311 |

**The per-coefficient table on the tuning block is unchanged from
`w = 0.5`**, to the last digit printed. That is not a coincidence and it
is the right way to read it: both weights reach the same optimum on all
200 of that block, so the estimator is identical there and the weight
can only differ where the optima differ, which on these 400 data sets
is one seed. Median fit time moved 0.341 s to 0.330 s.

**A finding the fresh block turned up, which is not about the start.**
Out of sample the overall Wald rate is 93.4 percent and its own Monte
Carlo interval stops at 94.6. Round two read that as a property of the
fresh block and round three corrected it to a property of the design,
pooled over 400: see "Punch round 3" above for the table and for why
the remedy this paragraph originally recommended was withdrawn.

### R2-NIT 1. One function, three references, so still one band

`hmm_starts_tol()` was called from all three sites but with a different
`ref` at each: the displacement passed the incumbent's log-likelihood
and the merge passed each cluster head. For a negative log-likelihood
the better optimum has the smaller magnitude, so the merge threshold
was always the smaller and a band of `grad_tol` always existed where
the modes table said two optima and the verdict said the incumbent was
best. The review constructed it at `grad_tol = 0.086` on the d4 fit.

**Fixed by carrying a value, not by calling a function.** The threshold
is computed once per call from the original fit's log-likelihood,
stored as `$mode_tol`, and read by all three sites including the print
method. `hmm_starts_modes()` now takes the absolute value rather than
the ingredients, so a local reference cannot be reintroduced there.

**Seen failing, and the first version of the test could not see it.**
`dev/latent-falsealarm-check.ps1` gained a second mode, `refs`, which
puts the three references back. The first invariant test judged the
modes table against `$mode_tol`, the unified value, and so passed on
that build: re-judging the modes by the verdict's own threshold hides
exactly the disagreement it is looking for. The invariant now uses
neither site's threshold. It asks what the TABLE asserts: if the table
names the original as a distinct optimum and names a better one beside
it, the verdict must say so. Against the two reverted builds:

| revert | result |
|---|---|
| `notol` (round 1's fix removed) | 7 failures, `n_false` 6, `n_moved` 6 |
| `refs` (round 2's fix removed) | 1 failure, at `grad_tol = 0.086`: `modes 2, says_local FALSE, better TRUE` |

Against the shipped build, both are green, and the band test sweeps
`grad_tol` over 1e-3, 3e-2, 0.08, 0.085, 0.086, 0.09 and 0.1.

**And the help page now says the two knobs are tied quadratically.**
`$mode_tol` moves as the SQUARE of `grad_tol`, so loosening the
convergence test by 100x loosens the merge threshold by 10 000x; on the
probe the 8.099-unit detection survives every `grad_tol` up to about
0.08 and is gone at 0.086.

### R2-NIT 2. A 0.2.2 fit could not be reproduced

The old rule was gone and `lca_init_extras()` is `@noRd`. `?lca` now
carries it as ten lines of R, in the "Labeling and starting values"
section, with the call that uses it. **Checked rather than asserted**
(`dev/latent-2p4-oldstart.R`): on seeds 20260910, 20260970 and
20261013, two of which the old rule got wrong, the block reproduces the
old starting values with a maximum absolute difference of **exactly 0**
and the old log-likelihood at a relative difference of **exactly 0**.
The block's one limitation is stated on the page: it does not mask
missing responses, which the shipped rule does.

## Punch round 1: what the review changed

`dev/reviews/2026-09-10-latent.md`, with the reviewer's scripts under
`dev/rev-latent-*`. Everything below is a change made in response; the
measurement work the review re-derived is unchanged and is left as it
stood.

### BLOCKER 1. `hmm_starts()` fired on a correct model, 6 of 6

The review's construction: six independent well-behaved two-state
chains (data seeds 4501 to 4506, `hmm_starts` seeds 8801 to 8806,
`n = 8`, `jitter = 2`). On every one, all eight refits reached the same
optimum and `nrow($modes)` was 1, and the summary still printed **"the
original fit found a local optimum"**, on gaps of 5.7e-11 to 1.5e-09.
`$best` was displaced by a refit for a gain of 3e-10.

**The root cause is not a missing constant, it is three copies of a
question.** Whether two log-likelihoods are the same optimum was asked
in three places: `hmm_starts_modes()`, which merged at `grad_tol^2`
relative; the `$best` replacement, which used `>`; and the printed
verdict, which used `gap > 0`. Two of the three carried no tolerance,
and the two thresholds differed by 1.7e05, so the function contradicted
itself inside one printout.

**The fix** adds `hmm_starts_tol(grad_tol, ref)`, one definition of the
threshold with the reason for `grad_tol^2` written beside it, and calls
it from all three. Measured on the review's own construction, with the
lane's build:

| | before | after |
|---|---|---|
| false alarms on the six unimodal fits | 6 of 6 | **0 of 6** |
| `$best` displaced there | 6 of 6 | **0 of 6** |
| the d4 probe still reports its local optimum | yes | **yes** |
| `test-hmm-starts.R` | 57 pass | **75 pass, 0 fail** |

**Seen failing.** `dev/latent-falsealarm-check.ps1` copies the package,
puts both comparisons back the way they shipped
(`dev/latent-falsealarm-revert.R`), installs into a scratch library
that is not this lane's, and runs the test file against it: 3 failures,
reporting `n_false` 6 and `n_moved` 6. The copy and the library are
deleted afterwards.

Two tests are new. `a unimodal fit is not told it found a local
optimum` runs the review's six chains and asserts 0 and 0. `the three
tolerances are one tolerance` asserts the property directly: whatever
the run found, the printed verdict and the modes count cannot
disagree.

### BLOCKER 2. `?hmm`'s 21.2 belonged to a different comparison

The sentence describes hmmTMB's 38 extra initial-distribution
parameters against frmtmb's `init = "estimated"`, and that measures
**19.678** (-3086.49259156 against -3066.81487197, probe 2's own
numbers). 21.231 is the same hmmTMB fit against `init = "uniform"`, a
40-parameter difference, and the 1.553 between them is frmtmb's own two
initial-distribution parameters. Corrected on the help page and here,
with both figures and what separates them.

### The `lca()` decision: option 1, and it works

The coordinator offered improving `lca()`'s own default start or
reclassifying the plan row. **Option 1, measured first.**

**The diagnosis came before the fix.** `lca_init_extras()` scored each
subject by the mean of its item codes and cut the scores into `K`
equal-count slices. On the plan's design every class endorses exactly
three of ten items at 0.85 and seven at 0.2, so every class has the
same expected score. Measured on seed 20260970
(`dev/latent-2p4-startfix.R`): the class means of that score are
0.3949, 0.3919, 0.3956 and 0.3865 against a within-class standard
deviation of 0.127, and the score explains **0.000814** of its own
variance between the true classes. The four slices it produced were
four samples of the same mixture, with essentially the same class
composition in each.

**The replacement.** Slice on the response PATTERN: one-hot the item
codes, seed `K` centers by farthest point, run hard Lloyd iterations to
a fixed point, relabel the clusters by the same mean score, and shrink
each class's smoothed proportions toward the pooled ones (halfway in
round one, nine tenths of the way after round two). It is
deterministic, costs `O(n K J)` and no fitting, and keeps the
documented low-to-high class ordering wherever that score means
anything.

**Measured on all 200 replicates of the sweep, against the poLCA
optima READ from `dev/latent-2p4-lca.tsv` rather than recomputed, so no
EM was re-run:**

| starting rule | reaches poLCA's optimum, on the tuning block |
|---|---|
| the 0.2.2 score cut | 192 of 200 |
| pattern clustering, no shrink | 194 of 200: fixes all 8, breaks 6 others by 79 to 303 units |
| pattern clustering, shrink 0.5 | 200 of 200, but 199 of 200 out of sample |
| pattern clustering, shrink 0.9, what ships | **200 of 200, and 200 of 200 out of sample** |

The shrink is what stops a hard partition from starting more confident
than any fitted profile will be. The weight was measured over
`0, 0.1, 0.25, 0.5, 0.75, 0.9` on all 200 of these seeds
(`dev/latent-2p4-shrink.R`), giving 6, 5, 2, 0, 1, 0, and round one
took `w = 0.5` as "the smallest weight that loses none".

**That basis was wrong and review round two showed why.** It is the
unique zero of a loss curve on the curve's own tuning set, so it moves
when the data does: out of sample 0.5 loses one and 0.75 loses none.
The weight now shipped is **0.9**, chosen from two independent blocks
of 200 and from where both edges of the usable interval are. See
"Punch round 2" above; this section is left as the round-one record.

**The recovery table was re-measured under the start that ships**, at
`w = 0.5` here and again at `w = 0.9` in round two; the round-two
numbers are the ones in "Item 2.4" below. Median fit time is 0.330 s
against the 0.2.2 start's 0.376, so the change costs nothing.

**What did NOT change:** `lca_starts()` was still not built, for the
reason the reviewer agreed with. The `?lca` recipe stays on the help
page, because the surface is still multimodal and a better default
start is not a proof that any particular data set has one mode.

### The surface claim, corrected

`?lca` said "the surface is the cause rather than this
implementation". The review found the counterexample inside the lane's
own range: on seeds 20261010 and 20261013 **all 10 of 10** poLCA single
random starts found the global optimum, while the old deterministic
start landed 274 to 278 units below every one of them, and no poLCA
start ever visited the mode it found. On those two the starting rule
and not the surface is what failed. The help page now says that, and
the fix above is what it implies.

### The smaller ones

- **The alignment check in `hmm_starts_scale()` was a tautology.**
  `confint.frmtmb_fit()` builds its `est` column as `object$opt$par`,
  so the check compared a vector with itself, measured difference
  exactly 0. Replaced by a cross-check of the standard errors against
  `vcov()`, which reaches the same covariance by its own path and names
  its own rows, on the rows the two share. The comment now states
  exactly what is and is not ruled out, including that a permutation
  confined to the `theta_` rows would still not be caught, because
  those have no second path.
- **The `"unit"` fallback is now a construction, not an evidence
  gap.** `test-hmm-starts.R` exercises it on a stub with registered
  `confint` and `vcov` methods: agreeing errors give `"se"`,
  disagreeing errors give `"unit"`, and so does either method being
  unavailable. The reviewer's finding stands and is recorded: no fit
  this family produces can reach the disagreement branch, which is why
  a stub is the right instrument. `"se-partial"` IS reachable from real
  fits (`K = 3` on 8 rows; a response constant within a state).
- **`$spread` read 0 when nothing was measured.** A spread over one
  value is now `NA_real_`, and the summary prints "not measurable, only
  one fit is in this set". A new test constructs the all-errored run
  (`jitter = 1e6`) and asserts both fields are NA.
- **The modes line said "among the converged"** while the set it counts
  includes the incumbent whatever its own gradient did. Relabelled
  "among the original and the converged".
- **`?lca` said "about 5 seconds"** for a recipe measured at 6.10 to
  16.77 s. That whole section is rewritten and the range is now stated.
- **`?hmm` now says what `diagnose()`'s clean verdict is worth here.**
  Its test is `max|grad| < 1e-3` ABSOLUTE and this fit's is 7.00028e-04,
  a margin of 1.43x. The help page says so, and says that what is 8.099
  units wrong is the answer and not the gradient.
- **Two `skip_if()` inside the depmixS4 identity block** could abandon
  it after one or two of five assertions and print a green line. Both
  are now expectations: a reference that did not converge, or a label
  order the two packages do not share, fails loudly.
- **The cost control covers one of three terms, and this is the
  sensitivity.** Only the per-refit 175.7 s is replicated and
  controlled (three blocks, control 1.023, identical gradient counts).
  The plain fit's 183.5 s is context and enters no projection. The
  `sdreport()`'s 37.2 s does enter, as a single reading: at `n = 8` it
  is 2.6 percent of the 24.0 minutes, so a reading wrong by a factor of
  two moves the projection by 2.6 percent and by four moves it by 7.7.
  It is not replicated because `sdreport()` is memoized per fit and a
  second reading needs a second three-minute fit.
- **Two reading traps in the arm tables**, left in the data and named
  here rather than edited after the fact. `identity_ok` is the LABEL
  match against the truth, NOT the third-party match, and it reads 1 on
  all 15 arm B rows including the 9 with no hmmTMB fit; and numeric NA
  was written space-padded, so `read.delim()` returns a character
  column and a denominator taken off the file counts 15 reference fits
  where 6 ran. `dev/latent-2p3-summarize.R` coerces before counting and
  reports the third-party count separately; the scripts that write the
  files are unchanged, because rewriting a measurement output after the
  fact is worse than naming what it holds.
- **"299 assertions" counted 2 skip records as assertions.** testthat's
  `nb` counts every result, including skips and warnings. The counts in
  this document are `nb`; the passing-assertion count is `nb` minus
  those records.

## What changed in the package

| file | what |
|---|---|
| `extensions/frmtmb.latent/R/hmm-starts.R` | NEW. `hmm_starts()`, its print method, and four `@noRd` helpers |
| `extensions/frmtmb.latent/R/hmm.R` | `?hmm`: the "Label switching and local optima" section rewritten with the measured numbers and pointed at `hmm_starts()`; an "Initial distribution" note on the hmmTMB comparison trap; a "Recovery at a realistic scale" section; `hmm_starts()` in `@seealso` |
| `extensions/frmtmb.latent/R/lca.R` | `?lca`: "Labeling and starting values" gains the measured failure rate and the measured verdict on two remedies; a new "Recovery at a realistic scale" section |
| `extensions/frmtmb.latent/tests/testthat/test-hmm-starts.R` | NEW, 13 blocks, 57 assertions |
| `extensions/frmtmb.latent/tests/testthat/test-hmm.R` | two new blocks: a K = 3 depmixS4 identity and a random-effect-on-a-transition hmmTMB identity, both at a size the ordinary suite can afford |
| `extensions/frmtmb.latent/tests/testthat/test-lca.R` | two new blocks: the poLCA identity at the plan's full realistic scale, and the pin on the starting rule's shrink weight |
| `extensions/frmtmb.latent/NAMESPACE`, `man/hmm_starts.Rd`, `man/hmm.Rd`, `man/lca.Rd` | roxygenised |

No version number, DESCRIPTION field or NEWS section was touched. See
"What this needs from whoever consolidates" at the end.

## Item 2.4, lca

### The design and the replicate count

The plan's realistic scale: `K = 4`, `n = 2000`, ten binary items, two
covariates on membership. The simulator is
`dev/latent-lca-sim.R`, lifted from the Phase 0 scale row so the two
rows measure the same construction. Item endorsement probabilities are
0.85 inside a class's own three-item block and 0.2 outside it; the
gating matrix is class 1 referenced with rows
`(-0.4, 0.8, -0.5)`, `(0.2, -0.6, 0.9)` and `(-0.1, 0.3, 0.4)`.

**200 replicates**, seeds 20260910 to 20261109, poLCA start seeds
700001 to 700200. Script `dev/latent-2p4-lca.R`, output
`dev/latent-2p4-lca.tsv`, summary `dev/latent-2p4-summarize.R`. The
count was chosen from a measured 4.2 s per replicate and a wish for a
coverage estimate whose own Monte Carlo interval is about one point
wide; 200 x 9 coefficients gives +/- 1.1 points.

Classes are matched to the truth, and to poLCA, by the item profile
table over all 24 permutations before anything is scored, and the
permutation used is written into every output row.

### The third-party identity

Compared: `lca()` against `poLCA(nrep = 10, tol = 1e-12,
maxiter = 20000)` on the SAME data, both run to convergence, at
`frmtmb_control(rel.tol = 1e-14, x.tol = 1e-14)`. Over the 192
replicates where the two reached the same optimum:

| quantity | worst over 192 | as a ratio to poLCA's own SE |
|---|---|---|
| log-likelihood, `\|d\|/\|ll\|` | 4.27e-14 | |
| item endorsement probabilities, max `\|d\|` | 5.61e-07 | 1.68e-05 |
| posterior class probabilities, max `\|d\|` | 2.12e-06 | |
| gating coefficients, max `\|d\|` | 2.22e-06 | 1.11e-05 |

The tolerance is stated as a ratio in both directions: relative to the
log-likelihood's own magnitude, and relative to the standard error
poLCA reports on the very quantity being compared. Both are quantities
the run measures.

### The recovery table

Two of them, because the starting rule changed in punch round 1 and a
recovery table has to describe the code that ships. Both are over the
SAME 200 seeds and score the same quantities; the only difference is
which start the fit ran from. Gating coefficients are re-referenced to
class 1 and matched to the truth by profile, and the Wald interval is
the fit's own, built from `vcov()` on the contrast `theta_a - theta_b`
that the re-referencing implies.

**The start that ships** (`w = 0.9`), `dev/latent-2p4-recheck.tsv`,
the 200 replicates of the TUNING block, all of which reached poLCA's
optimum. The table is unchanged to the printed digit from the `w = 0.5`
run of round one, because both weights reach the same optimum on all
200 of this block and the estimator is then identical:

| coefficient | truth | bias | RMSE | emp. SD | mean SE | cover95 |
|---|---|---|---|---|---|---|
| class2:(Intercept) | -0.40 | 0.0118 | 0.1517 | 0.1516 | 0.1474 | 95.5% |
| class2:x1 | 0.80 | 0.0002 | 0.1123 | 0.1126 | 0.1114 | 96.0% |
| class2:x2 | -0.50 | -0.0042 | 0.1963 | 0.1967 | 0.2019 | 94.0% |
| class3:(Intercept) | 0.20 | 0.0103 | 0.1111 | 0.1108 | 0.1077 | 94.0% |
| class3:x1 | -0.60 | -0.0070 | 0.0802 | 0.0801 | 0.0807 | 95.5% |
| class3:x2 | 0.90 | 0.0023 | 0.1472 | 0.1476 | 0.1425 | 93.0% |
| class4:(Intercept) | -0.10 | 0.0113 | 0.1159 | 0.1157 | 0.1162 | 95.5% |
| class4:x1 | 0.30 | 0.0044 | 0.0950 | 0.0952 | 0.0886 | 93.0% |
| class4:x2 | 0.40 | -0.0081 | 0.1552 | 0.1554 | 0.1583 | 92.5% |

Overall Wald coverage **94.3 percent, 1698 of 1800**, Monte Carlo
interval 93.3 to 95.4. Positive definite Hessian on 200 of 200;
`max|grad|/|logLik|` between 7.8e-10 and 9.1e-09; median fit 0.330 s.

**And on the 200 seeds the weight was NOT tuned on**, 20270401 to
20270600, `dev/latent-2p4-recheck-oos.tsv`: 200 of 200 reach poLCA's
optimum, the Hessian is positive definite on 200 of 200, the largest
bias is 0.032, and Wald coverage is **93.4 percent, 1682 of 1800**,
Monte Carlo interval 92.3 to 94.6.

That interval does not reach 95, and the right way to read it is the
POOLED one: at 200 replicates a single coefficient's Monte Carlo
half-width is about 3 points, so neither block alone separates a 91
percent row from 95, and the tuning block has the same two coefficients
as its two lowest rows (93.0 and 92.5). Over all 400, `class3:x2` is
91.75 percent (89.05 to 94.45) and `class4:x2` is 90.75 (87.91 to
93.59), both excluding 95, while `class2:x2` is 94.75 and fine and
every other row contains 95. The full table and the mechanism are in
"Punch round 3" above.

| other quantity | min | median | max |
|---|---|---|---|
| item endorsement probabilities, max `\|error\|` over 4 x 10 | 0.034 | 0.058 | 0.098 |
| the same table's RMSE | 0.016 | 0.023 | 0.034 |
| class shares, max `\|error\|` | 0.0027 | 0.011 | 0.031 |
| relative entropy of the classification | 0.784 | 0.807 | 0.830 |
| modal-assignment accuracy | 0.878 | 0.897 | 0.916 |

**The start that shipped before**, `dev/latent-2p4-lca.tsv`, over the
192 replicates where it reached the same optimum poLCA did. It is kept
because it is what the 200-replicate sweep and its identity table were
measured on, and because the two being so close is itself the finding:
the estimator was never the problem.

| coefficient | truth | bias | RMSE | emp. SD | mean SE | cover95 |
|---|---|---|---|---|---|---|
| class2:(Intercept) | -0.40 | 0.0147 | 0.1514 | 0.1511 | 0.1475 | 95.8% |
| class2:x1 | 0.80 | -0.0035 | 0.1087 | 0.1089 | 0.1114 | 96.4% |
| class2:x2 | -0.50 | -0.0071 | 0.1972 | 0.1975 | 0.2021 | 94.3% |
| class3:(Intercept) | 0.20 | 0.0108 | 0.1112 | 0.1109 | 0.1078 | 94.3% |
| class3:x1 | -0.60 | -0.0082 | 0.0794 | 0.0792 | 0.0807 | 96.4% |
| class3:x2 | 0.90 | 0.0025 | 0.1467 | 0.1470 | 0.1426 | 92.7% |
| class4:(Intercept) | -0.10 | 0.0140 | 0.1146 | 0.1141 | 0.1161 | 95.8% |
| class4:x1 | 0.30 | 0.0004 | 0.0938 | 0.0940 | 0.0886 | 93.2% |
| class4:x2 | 0.40 | -0.0084 | 0.1542 | 0.1544 | 0.1583 | 92.2% |

Overall Wald coverage 94.6 percent, 1634 of 1728, Monte Carlo interval
93.5 to 95.6. The empirical spread and the reported standard error
agree on every coefficient in both tables, which is the stronger
statement: the interval is the right WIDTH, not merely the right rate
on average.

**Verdict: the family recovers at the realistic scale and its own
intervals cover at the nominal rate. With the start that now ships it
finds the global optimum on 200 of 200; with the one it replaced, on
192 of 200.**

### The 8 replicates that did not, which is what got fixed

Everything in this section is the OLD start. It is kept because
it is the evidence the fix was built from and because a defect
that was found, measured and removed is a better record than a
defect that was only removed.

Seeds 20260970, 20260990, 20260996, 20260999, 20261010, 20261013,
20261042, 20261043. Script `dev/latent-2p4-modes.R`, output
`dev/latent-2p4-modes.tsv`.

| seed | cold `logLik` | poLCA(nrep=10) | gap | nlminb code | pdHess | `max\|grad\|/\|ll\|` | poLCA nrep=1 also low |
|---|---|---|---|---|---|---|---|
| 20260970 | -11581.784 | -11299.434 | -282.35 | 0 | TRUE | 2.1e-09 | 2 of 10 |
| 20260990 | -11673.907 | -11392.706 | -281.20 | 1 | TRUE | 2.7e-08 | 1 of 10 |
| 20260996 | -11657.291 | -11375.175 | -282.12 | 1 | FALSE | 1.8e-09 | 1 of 10 |
| 20260999 | -11645.779 | -11402.744 | -243.04 | 1 | TRUE | 2.4e-09 | 4 of 10 |
| 20261010 | -11647.382 | -11369.580 | -277.80 | 1 | TRUE | 2.0e-09 | 0 of 10 |
| 20261013 | -11612.796 | -11339.445 | -273.35 | 1 | TRUE | 3.1e-09 | 0 of 10 |
| 20261042 | -11691.789 | -11407.978 | -283.81 | 1 | TRUE | 1.9e-09 | 1 of 10 |
| 20261043 | -11496.205 | -11227.298 | -268.91 | 1 | TRUE | 2.8e-09 | 2 of 10 |

**What `diagnose()` actually says, corrected.** A first draft of this
document claimed it printed "No convergence problems detected" on seven
of the eight. It does not, and the difference matters. It prints that
line on ONE of them, seed 20260970, where `nlminb` returned 0. On six
more it prints three lines ending at "Hessian positive definite: TRUE"
and reports `nlminb` code 1, "false convergence (8)". Code 1 is not
evidence of trouble here: this project's own scale tier says so in a
comment and asserts the pair "positive definite Hessian and a gradient
small relative to the log-likelihood" instead, precisely because
`nlminb` returns 1 on fits that are fine. **Seven of the eight pass
that pair while being 243 to 284 log-likelihood units wrong.** The
eighth, 20260996, has a non positive definite Hessian and one
non-finite standard error, and is the only one a user could catch from
the fit alone.

poLCA's own single-start EM lands low on 0 to 4 of 10 seeds on these
data sets, which says the SURFACE and not this implementation is the
cause. What separates the two packages is that poLCA's `nrep` defaults
users into a multistart and `lca()`'s deterministic start does not.

### Two remedies, measured

`dev/latent-2p4-docrecipe.tsv` and `dev/latent-2p4-modes.tsv`.

| recipe | reached poLCA's optimum | gap left | cost |
|---|---|---|---|
| 10 refits, item parameters of the DEFAULT start perturbed by `rnorm(1)`, which is the recipe `?lca` already prints | 8 of 8 | 3e-11 to 6e-11 | 6.1 to 16.8 s |
| 10 refits, every outer parameter jittered 2 standard errors around the LOCAL OPTIMUM, which is `hmm_starts()`'s recipe | 0 of 8 | 72 to 273 units | comparable |

On one of the eight (20261013) the jittered multistart did not improve
on the cold start at all, in ten refits.

The jitter arm is `hmm_starts()`'s recipe written out by hand in
`dev/latent-2p4-modes.R`, because `hmm_starts()` itself refuses a fit
that is not an `hmm()`. Same construction: the fit's own standard
errors from `confint()`, a normal draw of two of them per outer
parameter, added to the estimates, refit from there.

Between 8 and 10 of the 10 default-start refits beat the cold start on
each data set, so the recipe is not scraping through on one lucky
draw.

**This is the lesson for the core seam.** Jittering around the
incumbent is the right move when the modes are close, and the wrong
move when they are not. On the HMM probe the two modes are 8.1
log-likelihood units apart and one coordinate 3.5 standard errors
apart, and jitter crosses it; on this LCA design they are 280 units
apart and jitter cannot.

### What went into `?lca`

Two sections. "Labeling and starting values" gains the 8-of-200 rate,
the 243-to-284 gap, the fact that `diagnose()` says nothing, the
evidence that the surface rather than the implementation is the cause,
and the measured verdict on both remedies. A new "Recovery at a
realistic scale" section carries the identity table, the coverage and
the recovery numbers above, and ends by pointing back at the starting
values, so a reader cannot take the recovery numbers without the
condition they hold under.

### What went into the suite

`test-lca.R` gains one block, `lca() reproduces poLCA at the realistic
scale (K = 4, n = 2000)`, on the 200-replicate run's own first seed. It
asserts the log-likelihood identity relative to its own magnitude, the
profile and coefficient identities relative to poLCA's own standard
errors on those quantities, and that the model recovers the truth
exactly as well as poLCA does on the same data, which is the only
recovery claim one replicate can carry.

**Seen to fail.** `dev/latent-2p4-testfails.R` runs the same body on
seed 20260970, one of the eight: `|dll|/|ll|` is 0.0244 against a bound
of 1e-10 and the profile gap is 27.85 poLCA standard errors against a
bound of 1e-3. On seed 20260910 the same two are 6.7e-15 and 1.3e-06.

## Item 2.3, hmm

### Reproducing the 8.1 units

`dev/latent-2p3-repro81.R`, and `dev/hmm/probeD4-re-multimodality.R` is
the record it reproduces. The probe predates the shipped `hmm()`
family and used the rung-1 `custom_family()` prototype, so the
question was whether today's family still lands there. It does, to
twelve digits.

| | value |
|---|---|
| probe D4's recorded cold start | -1096.09575602 |
| this family's cold start, same data | **-1096.09575602** |
| probe D4's recorded optimum | -1087.99646521 |
| this family, restarted there | -1087.9964652 |
| gap | **8.09929** |

Data: seed 2026, K = 2, 25 sequences of 30, stationary initial
distribution, `bf(y ~ 1 + (1 | gf), mu2 ~ 1 + (1 | gf))`. At the cold
start `diagnose()` reports `convergence` 0, `max|grad|` 7e-4 (6.4e-07
of the log-likelihood), a positive definite Hessian, no non-finite
standard error, no flat direction, and prints no message.

`frm_allfit()` on the same fit (`dev/latent-2p3-allfit.R`): nlminb,
optim, bobyqa and NLopt L-BFGS all reach -1096.096 with
`convergence == 0` and a log-likelihood spread of **7.84e-07**. Varying
the optimizer is agreement on the wrong answer.

### `hmm_starts(fit, n, jitter, seed, grad_tol, keep)`

What it does: `n` refits, each from the fit's own estimates plus a
normal draw per outer parameter whose standard deviation is `jitter`
times that parameter's own standard error; the best CONVERGED refit is
returned as `$best`, and the spread of the OPTIMA is reported.

Three things the brief asked it to get right.

**The spread is a spread of optima, and non-convergence is not hidden.**
Two spreads are reported, never one: over the original plus the
converged refits, which is the number to quote, and over every refit
that finished at all. A refit counts as converged when it returns
without error and `max|grad| / |logLik| <= grad_tol` (default 1e-3,
the same ratio the Phase 0 scale tier asserts). Refits that do not
converge are counted separately, named in the printed summary, kept out
of the converged spread, kept IN the other one, and can never be
returned as `$best` however high their log-likelihood. Refits that
error are counted, named, and are in neither spread. Three counts and
two spreads is the whole of the contract, and `test-hmm-starts.R`
asserts each branch on a construction rather than waiting for one:
`optCtrl = list(iter.max = 2)` forces non-convergence and
`jitter = 1e4` forces errors.

The INCUMBENT is held to the same test. It stays in both spreads
whatever its own gradient looks like, because it is the point the
spread is measured from and dropping it would leave a reader with a
spread of refits and no anchor, but `$original_converged` records the
verdict and the printed summary says so in as many words when it fails.
A summary that holds the refits to a test it does not hold the
incumbent to is the same guard failing open from the other side.

**The perturbation scale is a measured quantity, not a constant.**
`confint()` reports one row per outer parameter, so the standard errors
come from the fit itself. The row order is an assumption rather than a
contract, so it is CHECKED: `confint()`'s `est` column must equal
`fit$opt$par` to 1e-6 relative or the whole scale vector falls back to
1.0 on the internal scale, and the printed summary says which of the
three cases applied. A parameter whose standard error did not come back
is moved by the median of the ones that did rather than left alone.

**It finds the better mode on the probe that motivated it.**
`dev/latent-2p3-starts-probe.R`, `n = 8`, five `hmm_starts()` seeds per
jitter, on the reproduced probe D4 fit:

| jitter | recovered the -1087.996 optimum | distinct optima found | not converged | errored |
|---|---|---|---|---|
| 0.5 | 1 of 5 | 1 to 2 | 0 | 0 |
| 1 | 5 of 5 | 2 | 0 | 0 |
| 2 | 5 of 5 | 2 | 0 | 0 |
| 4 | 5 of 5 | 2 | 0 | 0 |
| 8 | 5 of 5 | 2 to 4 | 0 | 0 |

At `jitter = 2` a single run of 8 refits split 5 to the global optimum
and 3 to the local one, so the per-refit escape probability on this
surface is about 0.6 and `n = 8` misses with probability about 4e-4. At
`jitter = 8` the refits also find modes 16 to 47 units BELOW the local
one, which is why the summary reports the spread over all finished
refits: a large jitter buys reach and pays in variance, and the reader
should see both.

**Why each refit re-evaluates the call rather than reusing `fit$obj`.**
Reusing the taped object would save the tape build, and it would also
leave RTMB's `last.par.best` at the last refit's parameters. The
deferred `sdreport()` behind `summary()`, `vcov()` and `confint()`
reads exactly that, so `hmm_starts()` would silently move the standard
errors of the fit it was handed. The tape build is 839 ms against 118 s
of optimization at 25 000 rows (`dev/scale-findings.md`), under one
percent, so the safe path is also the cheap one. The cost of that
choice is that the call's data must be visible from where
`hmm_starts()` is called; the frame is assembled once, before any
refit, so that failure is one loud error rather than `n` bad starts.

### What `hmm_starts()` costs

`dev/latent-2p3-starts-cost.R`, run alone on an otherwise quiet
machine, at the plan's realistic design: 50 sequences of 500, K = 3
gaussian, `tr12 ~ (1 | id)`, 25 000 rows, 13 outer parameters, seed
20260910 for the data and 4201 for the starts.

**The instrument first, because two timing claims in an earlier round
did not survive replication.** The same work is run three times, from
the same seed, in one process, and the largest block over the smallest
is the CONTROL. Blocks are minutes long, so no reading is near
`proc.time()`'s 10 ms tick, and `Sys.time()` is what is read anyway.
Beside the clock is a count that does not move with the machine's load:
the same seed gives the same starting values, so the optimizer takes
the same path and the gradient-evaluation counts must be IDENTICAL
across blocks rather than merely close.

| block | seconds for 2 refits | gradient calls |
|---|---|---|
| 1 | 351.40 | 53 |
| 2 | 359.35 | 53 |
| 3 | 354.21 | 53 |

**CONTROL 1.023**, and the gradient counts identical to the unit. The
clock and the count agree, so the number below is a measurement and not
a reading of the machine.

| what | measured |
|---|---|
| one plain `frm()` fit of this model, se deferred | 183.5 s, 16 gradient calls |
| the deferred `sdreport()`, paid ONCE by `hmm_starts()` for the jitter scale | 37.2 s |
| one `hmm_starts()` refit | **175.7 s**, 26.5 gradient calls |
| `hmm_starts(n = 4)` | 12.3 minutes |
| `hmm_starts(n = 8)` | 24.0 minutes |
| `hmm_starts(n = 16)` | 47.5 minutes |

A refit costs what the original fit costs, which is the answer the
design predicts: it is a whole fit, frame and tape and optimization.
It takes MORE gradient evaluations than the original (26.5 against 16)
and slightly less wall clock, because it starts away from the optimum
but on a warm process. The per-refit spread over the six refits was
167 to 187 seconds.

### `n`: what to default to

**The shipped default is `n = 8`, and it is a starting point rather
than a recommendation.** What it is chosen from: at `jitter = 2` on the
probe, one run of 8 refits split 5 to the global optimum and 3 to the
local one, so a single refit escapes with probability about 0.6 on that
surface. At that rate `n = 8` misses with probability about 4e-04 and
`n = 4` with about 2e-02. Both numbers are properties of THAT surface,
not constants, which is why the help page says to raise `n` when the
printed summary shows more than one optimum and to price it from the
`seconds` line when the model is large.

**Would I default to it at the realistic scale? No, and the help page
says so.** 24 minutes is not a thing to do by accident on a model that
took three minutes to fit. What the measured numbers argue for is that
`n` should be chosen against the cost the function itself reports, and
the printed summary carries both the seconds and the load-independent
gradient count so that a second call can be priced from the first. A
default of 8 is right for a model a user is iterating on, which is the
size at which multimodality is cheapest to find, and wrong for a
25 000-row fit, which is why it is documented as a number to change
rather than a number to trust.

### The replicate counts, and how they were chosen

The identity and the recovery table are different questions and cost
different amounts. An identity is a property of the two codes and one
data set establishes it; a recovery table is a property of the
estimator and needs replicates. So the third-party arm runs on the
FIRST few replicates and every replicate contributes to the recovery
table. "First few" rather than "a sample" is deliberate: the first N
replicates are the first N SEEDS, so nothing is selected on its
outcome.

What was measured before any run was launched, on one replicate at full
scale in a quiet process: arm A 22 s (frm 9 s, depmixS4 13 s), arm B
340 s (frm 222 s, hmmTMB 118 s). Both estimates were wrong, and the
reasons are worth writing down because they are the reasons the counts
are what they are.

- **depmixS4's random-start EM varies by more than an order of
  magnitude between replicates.** Its four starts cost 18, 263 and 645
  seconds on the first three replicates, at `tol = 1e-12` and
  `maxit = 5000` on 25 000 rows, where the frmtmb fit beside it cost 12
  to 20 seconds. That variance belongs to the reference, not the model,
  and paying it forty times buys nothing.
- **Two of these processes on one machine is not two processes' worth
  of work.** Running the arms in parallel made each replicate roughly
  five times slower than the single replicate it was priced from. The
  arms were restarted to run SEQUENTIALLY, which is also why the
  seconds in the output tables are labelled indicative: the timing
  CLAIM this row makes comes from `dev/latent-2p3-starts-cost.R`, which
  runs alone and carries its own control.

The counts, all fixed before any result past the third replicate was
read:

- **Arm A identity: 3 replicates**, seeds 20260910 to 20260912, in
  `dev/latent-2p3-hmm-A-id.tsv`, with depmixS4 at `tol = 1e-12`,
  `maxit = 5000`, best of 4 random starts. Three is few, and it is
  three because of the 645 seconds above. What backs it up is a fourth
  data set in the suite (`test-hmm.R`, 20 sequences of 100) and the
  fact that an identity that holds to 3e-14 and 5e-12 relative on
  25 000 rows is not the kind of claim that turns over on the next
  replicate.
- **Arm A recovery: 40 replicates**, seeds 20260910 to 20260949, in
  `dev/latent-2p3-hmm-A.tsv`, no third-party arm. Forty times twelve
  parameters give a coverage estimate whose own Monte Carlo interval is
  about two points wide.
- **Arm B: 15 replicates**, seeds 20260910 to 20260924, with hmmTMB on
  the first 6. Fifteen times thirteen parameters give a coverage
  estimate whose own interval is about three points wide. That is thin,
  and it is stated rather than dressed up: one fit of this model is two
  to five minutes and the budget bought fifteen.

### Arm A: K = 3 gaussian against depmixS4

The plan's design with the random effect removed, since depmixS4 has
none: `bf(y ~ 1)` with `hmm(K = 3, gaussian(), time = t, group = id,
init = "estimated")`, 50 sequences of 500, 25 000 rows. Both sides
estimate ONE initial distribution shared across sequences, which is
depmixS4's `ntimes =` behaviour and frmtmb's `init = "estimated"`.
State means -2, 0, 3 with a common sd of 0.7; the transition matrix is
the plan's, a row-wise multinomial logit with state 1 the reference and
`eta` rows (-1.5, -2.5), (2.0, -1.0), (-1.0, 2.0).

**The identity**, 3 replicates, `dev/latent-2p3-hmm-A-id.tsv`,
depmixS4 at `tol = 1e-12`, `maxit = 5000`, best of 4 random EM starts:

| quantity | over 3 replicates |
|---|---|
| log-likelihood, `\|d\|/\|ll\|` | 3.1e-14 to 4.9e-12 |
| state means and sds, max `\|d\|` | 7.9e-07 to 1.5e-06 |
| transition matrix, max `\|d\|` | 4.0e-07 to 5.9e-07 |

The estimate gaps are stated against the fit's own yardstick: the mean
standard error on a state mean in this design is 0.0081 to 0.0095, so
1.5e-06 is **1.6e-04 of one standard error**.

What the reference itself did, from `dev/latent-2p3-hmm-A-id-log.txt`,
because it is evidence about the surface rather than about either
package:

| replicate | depmixS4's four random starts |
|---|---|
| 1 | -39087.47 three times, -39171.66 once: **84 units** below |
| 2 | -39488.68 three times, -43349.45 once after 3854 EM iterations: **3861 units** below |
| 3 | -38973.33 twice, and two starts that never reported convergence inside `maxit = 5000`, which is where that replicate's 645 seconds went |

So a random EM start on this design lands on a badly wrong optimum
about a quarter of the time, and frmtmb's deterministic quantile-spread
start found the best one on all three. Multimodality here belongs to
the surface, not to either implementation.

**The recovery table**, 40 replicates, `dev/latent-2p3-hmm-A.tsv`,
seeds 20260910 to 20260949. Every parameter is scored on its own
internal scale (a log sd, a transition logit), which is the scale the
fit's interval is built on. The fitted state labels were the true
labels in 40 of 40 replicates, so no permutation enters the table.

| parameter | truth | bias | RMSE | emp. SD | mean SE | cover95 |
|---|---|---|---|---|---|---|
| mu1 | -2.000 | -0.0011 | 0.0092 | 0.0093 | 0.0095 | 95.0% |
| mu2 | 0.000 | -0.0004 | 0.0084 | 0.0085 | 0.0081 | 92.5% |
| mu3 | 3.000 | -0.0003 | 0.0091 | 0.0092 | 0.0091 | 97.5% |
| log sigma1 | -0.357 | 0.0027 | 0.0083 | 0.0079 | 0.0098 | 100.0% |
| log sigma2 | -0.357 | 0.0006 | 0.0076 | 0.0076 | 0.0089 | 97.5% |
| log sigma3 | -0.357 | -0.0022 | 0.0103 | 0.0102 | 0.0094 | 92.5% |
| tr12 | -1.500 | -0.0048 | 0.0364 | 0.0365 | 0.0339 | 97.5% |
| tr13 | -2.500 | -0.0087 | 0.0532 | 0.0531 | 0.0484 | 95.0% |
| tr22 | 2.000 | -0.0066 | 0.0329 | 0.0327 | 0.0357 | 95.0% |
| tr23 | -1.000 | -0.0065 | 0.0615 | 0.0619 | 0.0617 | 97.5% |
| tr32 | -1.000 | 0.0139 | 0.0980 | 0.0982 | 0.0833 | 92.5% |
| tr33 | 2.000 | 0.0098 | 0.0430 | 0.0424 | 0.0418 | 95.0% |

Overall Wald coverage 95.6 percent, 459 of 480, Monte Carlo interval
93.8 to 97.5. The empirical spread and the reported standard error
agree on every parameter, and the largest bias is 0.014 on a
transition logit whose own standard error is 0.083.

On the response scale: the state means come back within 0.021, the
state standard deviations within 0.017, and every cell of the 3 x 3
transition matrix within 0.014, worst case over the 40. Positive
definite Hessian on 40 of 40, `max|grad| / |logLik|` between 4.3e-08
and 4.9e-07, and no replicate needed more than 42 gradient
evaluations.

**Verdict: the family recovers at the realistic scale, and no replicate
of this arm found a local optimum.**

### Arm B: `tr12 ~ (1 | id)` against hmmTMB

The plan's model exactly: `bf(y ~ 1, tr12 ~ 1 + (1 | id))` with
`hmm(K = 3, gaussian(), time = t, group = id, init = "uniform")`, 50
sequences of 500, 25 000 rows, one random intercept per sequence on the
1 -> 2 transition with a true standard deviation of 0.6. `init =
"uniform"` on both sides, and it is the correct model as well as the
allowed one: the simulator draws the first state uniformly, and
`init = "stationary"` is refused when a transition carries a predictor.

**The identity**, 6 replicates, seeds 20260910 to 20260915:

| quantity | over 6 replicates |
|---|---|
| Laplace marginal log-likelihood, `\|d\|/\|ll\|` | 6.1e-12 to 6.9e-11 |
| state means, state sds and `sd(tr12 \| id)`, max `\|d\|` | 5.2e-06 to 1.3e-04 |

The two are the same function of the data once the initial distribution
is pinned on both sides; the three traps that had to be cleared before
that was true are in the next section. The estimate gaps are against the
fit's own yardstick: the mean standard error on a state mean here is
0.008 to 0.009, so 1.3e-04 is **1.5e-02 of one standard error**.

**A gap in the output table that is the HARNESS and not the model, run
down rather than explained away.** `dev/latent-2p3-hmm-B.tsv` has a
`tpm_max_gap` column reading 0.039 to 0.180, next to a log-likelihood
identity of 1e-11. A number that large beside a number that small is
either a real disagreement or a comparison of two different quantities,
and `dev/latent-2p3-tpmgap.R` says which. frmtmb's number there is the
POPULATION transition matrix, built from the `tr{i}{j}` intercepts with
the random effect at zero. hmmTMB's `par()$tpm[, , 1]` is the
transition matrix at DATA ROW 1, which belongs to sequence 1 and
carries sequence 1's own random intercept. On seed 20260910 that
intercept is 0.2828, and:

| comparison | max `\|d\|` |
|---|---|
| frmtmb POPULATION tpm against hmmTMB's `tpm[, , 1]` | 0.0467 |
| frmtmb's tpm with sequence 1's conditional mode added, against the same | **1.8e-06** |

The row-1 fixed effects, which reference state 1 under both
conventions, compare directly and need no such repair: `tr12`
-1.39321375 against -1.39322981 and `tr13` -2.47338704 against
-2.47337785, a largest gap of **3.3e-04 of one standard error** on the
fit's own 0.0486 and 0.0858. The `tpm_max_gap` column is meaningless
for arm B and is left in the file with this note beside it rather than
deleted, since deleting a column a reader can see in the data is how a
number gets re-derived from scratch later.

**The recovery table**, 15 replicates, seeds 20260910 to 20260924. The
fitted labels were the true labels in 15 of 15.

| parameter | truth | bias | RMSE | emp. SD | mean SE | cover95 |
|---|---|---|---|---|---|---|
| mu1 | -2.000 | -0.0019 | 0.0072 | 0.0072 | 0.0093 | 100.0% |
| mu2 | 0.000 | -0.0004 | 0.0065 | 0.0068 | 0.0081 | 100.0% |
| mu3 | 3.000 | -0.0033 | 0.0078 | 0.0073 | 0.0090 | 100.0% |
| log sigma1 | -0.357 | -0.0027 | 0.0088 | 0.0087 | 0.0097 | 100.0% |
| log sigma2 | -0.357 | -0.0027 | 0.0088 | 0.0087 | 0.0090 | 93.3% |
| log sigma3 | -0.357 | 0.0009 | 0.0095 | 0.0098 | 0.0094 | 93.3% |
| tr12 | -1.500 | 0.0150 | 0.1024 | 0.1048 | 0.0941 | 93.3% |
| tr13 | -2.500 | 0.0110 | 0.0500 | 0.0505 | 0.0478 | 86.7% |
| tr22 | 2.000 | -0.0156 | 0.0333 | 0.0304 | 0.0362 | 93.3% |
| tr23 | -1.000 | -0.0034 | 0.0459 | 0.0474 | 0.0615 | 100.0% |
| tr32 | -1.000 | -0.0012 | 0.1144 | 0.1184 | 0.0829 | 86.7% |
| tr33 | 2.000 | 0.0065 | 0.0541 | 0.0556 | 0.0414 | 80.0% |
| log sd(tr12 \| id) | -0.511 | 0.0193 | 0.1070 | 0.1089 | 0.1187 | 100.0% |

Overall Wald coverage 94.4 percent, 184 of 195, Monte Carlo interval
91.1 to 97.6. Fifteen replicates put +/- 10 points on any single row of
that table, so `tr33` at 80 percent and `tr13` and `tr32` at 86.7 are
not distinguishable from 95 at this count; the overall figure is the
one to read, and the thing that would show a real problem, the
empirical spread against the reported standard error, matches on every
row but `tr32` and `tr33`, where the reported error is about 30 percent
narrow.

**The random effect itself recovers.** `sd(tr12 | id)` comes back with
mean 0.615 and standard deviation 0.067 against a truth of 0.6, ranging
0.511 to 0.718 over the fifteen; its own log-scale interval covers the
truth on 15 of 15, which at this count is consistent with 95 percent
and not evidence of anything better.

On the response scale: state means within 0.017, state standard
deviations within 0.018, and every cell of the population transition
matrix within 0.030, worst case over the fifteen. Positive definite
Hessian on 15 of 15, `max|grad| / |logLik|` between 4.5e-08 and
3.5e-07, and no replicate needed more than 31 gradient evaluations.

**No replicate of this arm found a local optimum**, and hmmTMB, which
was started at the truth, reached the same optimum on all six it ran
on.

**Verdict: the family recovers at the realistic scale with a random
effect on a transition, and its own intervals cover at the nominal
rate.**

### The hmmTMB comparison traps

Three, and each cost real time, so each is written down.

1. **`initial_state = "estimated"` is per SEQUENCE.** hmmTMB estimates
   one initial distribution for every time series; frmtmb's
   `init = "estimated"` estimates one shared across sequences. On 20
   sequences of 100 that is 38 extra parameters and a log-likelihood
   **19.678** units higher (-3086.49259156 against -3066.81487197), and
   nothing in either package's output says the two fits are of
   different models. Against `init = "uniform"` on the frmtmb side the
   same hmmTMB fit is 21.231 units higher instead, because uniform
   gives up frmtmb's own two initial-distribution parameters too. The
   first of those is the one that belongs to a sentence about hmmTMB's
   38, and a first draft of `?hmm` used the second.
   `dev/latent-2p3-hmmtmb-probe2.R` is where the gap was found and
   `dev/latent-2p3-hmmtmb-probe4.R` is where it was closed:
   `fixpar = list(delta0 = setNames(rep(NA, 2 * n_seq), rownames(...)))`
   with the default uniform `delta0` pins hmmTMB's initial distribution
   the way `init = "uniform"` pins frmtmb's, and the two log-likelihoods
   then agree to 6.6e-14 relative.
   **A wrong first guess, recorded because it was wrong:** the 21.2
   looked like a missing normalizing constant on hmmTMB's random-effect
   penalty, and `-(k/2) log(2 pi sd^2)` accounts for 5.66 of it. It was
   not that at all. hmmTMB's marginal likelihood is complete; the
   initial distribution was the whole difference.
2. **`ref` cannot be moved once a formula matrix is given.** hmmTMB's
   `MarkovChain$new(ref =)` documents a settable reference cell per
   row, but `check_args()` refuses any formula matrix whose DIAGONAL is
   not `"."`, so `ref = rep(1, K)` is unreachable in practice. Row 1's
   reference is state 1 under both conventions either way, so `tr12`
   and `tr13` compare directly; rows 2 and 3 are compared as a
   transition MATRIX instead.
3. **A column called `state` is read as known states, silently.** This
   is probe D3's finding from the feasibility round and it still holds
   at hmmTMB 1.1.2. The frame handed to hmmTMB drops it.

### What went into `?hmm`

Four changes.

1. **"Label switching and local optima" is rewritten.** It used to say
   multimodality "has been measured" at 8.1 units and send the reader
   to `frm_allfit()`. It now gives both log-likelihoods to the digit,
   names every diagnostic that passes on the wrong one, points at
   `hmm_starts()`, carries the jitter measurement the default of 2 is
   chosen from, and says what `frm_allfit()` actually does on that fit:
   four optimizers, `convergence == 0` each, a log-likelihood spread of
   7.8e-07, all on the wrong answer.
2. **"Initial distribution" gains the hmmTMB trap.** `init =
   "estimated"` fits ONE initial distribution shared across sequences;
   hmmTMB's `initial_state = "estimated"` fits one per sequence. The
   note gives the size of the difference, 38 parameters and 21.2
   log-likelihood units on 20 sequences, and the `fixpar` spelling that
   makes the two comparable.
3. **A new "Recovery at a realistic scale" section**, carrying both
   arms' coverage, bias and identity numbers, and ending on the
   sentence that keeps them honest: no replicate of either arm found a
   local optimum, which is a statement about these designs and not
   about the family.
4. **`hmm_starts()` joins `@seealso`.**

### What went into the suite

`test-hmm.R` gains two blocks, both at a size the ordinary suite can
afford and both with every tolerance a ratio to something the run
measures:

- *a three-state gaussian HMM agrees with depmixS4*, 20 sequences of
  100. The log-likelihood identity is relative to its own magnitude;
  the means, standard deviations and transition matrix are compared
  against the fit's OWN standard error on a state mean. It also
  asserts the relative gradient, because the absolute-gradient warning
  fires on any long chain and suppressing a warning without asserting
  what it was about is how a real one gets lost.
- *a random effect on a transition agrees with hmmTMB*, 15 sequences of
  120, `tr12 ~ (1 | id)`, with the initial distribution pinned on both
  sides. This is new coverage: the existing hmmTMB block is
  fixed-effect only.

`test-hmm-starts.R` is new, 13 blocks and 57 assertions.

## Every test file this lane ran, with its counts

One test file per R process, `NOT_CRAN=true`, through
`dev/latent-runtests.R`, which prints the block, assertion and SKIP
counts before the failure count and exits non-zero on any error. Final
pass, `dev/latent-suite.ps1`:

After punch round 3:

| file | blocks | records | fail | error | skip |
|---|---|---|---|---|---|
| test-bracket-access.R | 1 | 1 | 0 | 0 | 0 |
| test-hmm.R | 25 | 115 | 0 | 0 | 0 |
| test-hmm-starts.R | 18 | 84 | 0 | 0 | 0 |
| test-lca.R | 15 | 105 | 0 | 0 | 0 |
| test-message-uniqueness.R | 1 | 4 | 0 | 0 | 0 |
| test-scale.R | 2 | 2 | 0 | 0 | 2 |
| test-structure-latent.R | 3 | 20 | 0 | 0 | 0 |
| **total** | **65** | **331** | **0** | **0** | **2** |

Before round one: 59 blocks and 299 records; after round one, 63 and
320; after round two, 64 and 329. The column is testthat's `nb`,
which counts every RESULT and not only every assertion: 2 of the 320
are the gated Phase 0 scale tier's skip records, so the passing
assertion count is 329. Those two skips need `FRMTMB_SCALE_TESTS=true`
and are named in the output rather than counted silently.

`R CMD check --as-cran` on the built tarball, with pandoc and TinyTeX
on PATH and no `--no-manual` or `--no-build-vignettes`:
**Status: 1 NOTE**, and the NOTE is the environmental one the lane
rules predict, "checking HTML version of manual ... Skipping checking
math rendering: package 'V8' unavailable". Tests inside the check ran
in 74 s and the vignette rebuilt in 15 s.

## The guards, seen failing

Six, because a guard that has only been seen passing is not evidence,
and one of these was seen NOT failing first, which is why it changed.

0. **`hmm_starts()`'s guard tests go red on the code they pin, both
   of them.** `dev/latent-falsealarm-check.ps1` copies the package,
   reverts ONE fix at a time (`dev/latent-falsealarm-revert.R`),
   installs into a scratch library that is not this lane's, and runs
   the test file. Mode `notol` puts back round one's untoleranced
   comparisons: 7 failures, `n_false` 6 and `n_moved` 6 against an
   expected 0 each. Mode `refs` puts back round two's three reference
   values: 1 failure, at `grad_tol = 0.086`, reporting `modes 2,
   says_local FALSE, better TRUE`. Against the shipped build the file
   is 84 of 84. Copies and scratch library are deleted afterwards.

   **The first version of the round-two test could not see it.** It
   judged the modes table against `$mode_tol`, the unified value, and
   so passed on the `refs` build: re-judging the modes by the verdict's
   own threshold hides exactly the disagreement it is looking for. That
   is a guard failing open, found by running it against the absent case
   rather than by review. The invariant now uses neither site's
   threshold and asks what the TABLE asserts.

   Round three added two more modes to the same harness, `w05` and
   `w099`, which move the starting rule's shrink weight to each edge of
   its usable interval. Both fail `test-lca.R`'s new pin with 2
   failures, at relative gaps of 0.0218 and 0.0217 from the optimum, so
   the constant is now held from both sides.

1. **The hazard scanner reaches the new `R/` file.**
   `dev/latent-hazard-check.R` plants `est$b` in `hmm_starts_scale()`,
   which is a `$` read on one of frmtmb's hazard containers with a name
   that partial-matches `beta`. `test-bracket-access.R` then fails with
   `hits` = `"hmm_starts_scale: est$b"`. Taking it out again returns the
   file to green. An empty scan is only evidence once the scan has been
   shown to find something in that file.
2. **The lca identity assertions can go red.**
   `dev/latent-2p4-testfails.R` runs the new `test-lca.R` block's body
   on seed 20260970, one of the eight where the deterministic start
   reaches a local optimum: `|dll| / |ll|` is 0.0244 against a bound of
   1e-10, and the profile gap is 27.85 poLCA standard errors against a
   bound of 1e-3. On the seed the test actually uses the same two are
   6.7e-15 and 1.3e-06, so the margins are about 1e05 and 750.
3. **`hmm_starts()`'s non-convergence and error branches are reached by
   construction, not by waiting.** `optCtrl = list(iter.max = 2)` makes
   every refit stop far from any optimum, and `jitter = 1e4` puts a
   start far enough out that the objective is not finite there. Both
   are asserted in `test-hmm-starts.R` rather than hoped for, and the
   error block skips loudly if no refit errors rather than asserting
   nothing.

And one that was NOT a guard failing, but a harness failing open, found
the same way: `dev/latent-install.ps1` gated its `R CMD INSTALL` on
PowerShell's `$?`, which goes false the moment a native command writes
anything to stderr. `roxygenise()` writes "Loading required package:
frmtmb" there on every run, so the install was silently skipped every
time and two probes ran against a stale library. It gates on
`$LASTEXITCODE` now and prints the exit code.

## What I decided NOT to do, and why

- **`frm_multistart()` in core.** Out of scope for this lane by
  instruction, and the lca measurement above says the core version
  needs a design decision this lane cannot make on its own: see "what
  the core seam would need" below.
- **A multistart for `lca()`.** Building `lca_starts()` beside
  `hmm_starts()` would be two copies of the thing that is filed as one
  core seam, and the reviewer agreed. What punch round 1 changed is
  that the alternative is no longer a help-page paragraph: the DEFAULT
  START now carries the eight, 200 of 200, deterministically and at the
  same cost per fit. A better default is not a seam and not a second
  copy of one. The `?lca` recipe stays on the page anyway, because a
  start that works on 200 replicates of one design is not a proof about
  anybody's data.
- **A test pinning the lca local optimum.** It is a defect that is not
  fixed in this lane, so a test asserting the current behaviour would
  pin the wrong answer, and a test asserting the right one would fail.
  The measurement is on the help page and in this document instead.
- **Asserting `sd(tr12 | id)` recovery on one replicate in the suite.**
  A variance component on 50 groups has a wide sampling distribution
  and a single-seed assertion on it is the kind of check the 1.0a note
  warns about. The replicate table carries it; the suite asserts the
  identity against hmmTMB, which is a property of the code and not of
  the draw.
- **`hmm_starts()` scattering the family's DEFAULT start as well as the
  incumbent.** That is what the lca measurement says a general
  multistart needs, and it is a feature the plan's row does not
  describe. It is written up below instead of built.

## Defects found and not fixed

1. **FIXED in punch round 1.** `lca()`'s old deterministic start
   reached a local optimum on 4 percent of realistic-scale data sets,
   silently. The replacement reaches poLCA's optimum on 200 of 200. The
   defect is recorded rather than deleted because the shape of it
   generalizes: a starting rule built from a scalar summary of the
   response is blind to any latent structure that summary does not
   carry, and `mixture()`'s response quantiles are the same shape of
   rule.
2. **`diagnose()` has nothing to say about multimodality.** It is a
   local check by construction and this is not a criticism of it, but
   the pair "convergence 0 and pdHess TRUE" reads to a user as "this is
   the answer", and on both families in this package it is not. A
   `diagnose()` line that says which checks were NOT performed would
   help; that is a core change and belongs with the seam.
3. **hmmTMB's `ref =` is documented but unusable with a formula
   matrix.** Upstream, not ours. Worth an upstream issue alongside the
   `state` column one that `dev/hmm-feasibility.md` already filed.
4. **RESOLVED in punch round 1.** `hmm_starts_scale()`'s `"unit"`
   fallback was recorded as an evidence gap. It is now a construction:
   `test-hmm-starts.R` reaches it on a stub with registered `confint`
   and `vcov` methods, in three ways (the two disagreeing, `confint()`
   unavailable, `vcov()` unavailable). The finding underneath stands
   and is the reason a stub is the right instrument: no fit this family
   produces can make the two paths disagree. `"se-partial"` IS
   reachable from real fits, on `K = 3` fitted to 8 rows and on a
   response that is constant within a state.

## What the core `frm_multistart()` seam would need

The plan files `frm_multistart(fit, n, jitter)` as the general form of
`hmm_starts()`. Five things this lane measured that it would have to
carry, and the first is the one that is not obvious.

1. **Two families of starting values, not one.** Jitter around the
   incumbent crosses an 8.1-unit gap on the HMM probe and does not
   cross a 280-unit one on the LCA design. A general multistart needs
   both "perturb the answer" and "re-draw the family's own default
   start", and probably a third, "draw at random the way the reference
   implementations do". `hmm_starts()` implements only the first, and
   the `?lca` recipe is the second written out by hand. Punch round 1
   sharpened this rather than removing it: fixing `lca()`'s default
   start removed one instance of the problem and left the general
   statement exactly where it was, because the seam as filed would
   still have recovered 0 of those 8.
2. **A per-family hook for what to perturb.** The `?lca` recipe
   perturbs the ITEM PARAMETERS specifically, which are family extra
   parameters, and leaves the gating coefficients alone. A core
   function jittering every outer parameter uniformly would not
   reproduce it. `frmtmb_structure()` is where a family would declare
   "these are the parameters a restart should scatter".
3. **The scale question, answered from the fit.** Standard errors from
   `confint()` are the unit a user can reason about, and the row-order
   assumption behind reading them has to be checked rather than
   assumed; core owns `outer_par_names()` and would not need the check
   that `hmm_starts()` does.
4. **The three counts and the two spreads.** Converged, not converged,
   errored; the spread over converged optima and the spread over
   everything that finished. This is not decoration: a multistart that
   reports one spread over the runs it liked is a guard that fails
   open.
5. **A load-independent cost, reported.** `$table` carries the
   optimizer's own function and gradient evaluation counts beside the
   seconds. The same seed gives the same counts exactly, which is a
   check the clock cannot make.
6. **One definition of the threshold, not one per call site.**
   `hmm_starts()` shipped this lane with three copies of "is this a
   different optimum", two of them carrying no tolerance at all, and
   the result was a guard firing on a correct model 6 of 6 times. The
   copies are what made it possible. A core version will have at least
   as many call sites.

## What this needs from whoever consolidates

- **A minor version bump for frmtmb.latent**, 0.2.2 to 0.3.0. Two
  reasons now, not one: `hmm_starts()` is a new export, and
  `lca()`'s DEFAULT STARTING VALUES changed. The second is a behaviour
  change on a shipped family and it belongs in the user's hands: a
  given data set can now end on a different optimum, and will end on a
  better one where the two differ. The class LABELING convention is
  preserved as far as the old scoring statistic can preserve it (the
  clusters are ordered by that same score), but a data set where the
  score is uninformative can come back with the classes in a different
  order than 0.2.2 gave. I have not changed DESCRIPTION or NEWS.
- **Two NEWS bullets**, saying in substance:
  - `hmm_starts()` refits an `hmm()` model from jittered starting
    values and reports the spread of the optima, because the default
    cold start has been measured converging 8.1 log-likelihood units
    below the optimum with every diagnostic clean;
  - **`lca()`'s starting values changed, and a fit can move.** The old
    rule cut subjects on the mean of their item codes, which is blind
    to any design whose classes differ in which items they endorse, and
    it reached a local optimum 243 to 284 log-likelihood units below
    `poLCA(nrep = 10)` on 8 of 200 replicates of the realistic-scale
    design, with a positive definite Hessian on seven of the eight. The
    rule now clusters on the response pattern and shrinks the class
    profiles nine tenths of the way toward the pooled ones, and reaches
    poLCA's optimum on 200 of 200 on the seeds its one constant was
    tuned on and on 200 fresh ones. Class 1 is still the low-score end
    wherever that score means anything; where it does not, the labeling
    is deterministic but not the same one 0.2.2 gave. **A 0.2.2 fit can
    be reproduced**: `?lca`'s "Labeling and starting values" section
    carries the old rule as ten lines of R, checked to reproduce the
    old starting values and the old log-likelihood at a relative
    difference of exactly 0.
- **`dev/extension-gaps-plan.md` should be edited** (I did not edit
  it):
  - row 2.4's "expected trouble" cell should stop saying "none
    expected" and say instead: *WRONG, and FIXED. 8 of 200 replicates
    at this design reached a local optimum 243 to 284 log-likelihood
    units below poLCA's under the 0.2.2 starting rule, and seven of the
    eight passed the pair the scale tier asserts. The cause was the
    start, not the surface, on at least 2 of the 8: there all 10 of 10
    poLCA single random starts found the global optimum. `lca()`'s
    default start now clusters on the response pattern and reaches
    poLCA's optimum on 200 of 200, with coverage 94.3 percent and the
    log-likelihoods agreeing to 1.8e-14 relative. Measured, DONE*;
  - row 2.3 should be marked DONE, with: *`hmm_starts(fit, n, jitter)`
    shipped. The 8.1-unit probe reproduces to twelve digits against the
    shipped family and `hmm_starts(n = 8, jitter = 2)` recovers from it
    on 5 of 5 seeds. Recovery: coverage 95.6 percent over 40 replicates
    against depmixS4's model and 94.4 percent over 15 with
    `tr12 ~ (1 | id)`; identities 3.1e-14 to 4.9e-12 relative against
    depmixS4 and 6.1e-12 to 6.9e-11 against hmmTMB. No replicate of
    either arm found a local optimum. One refit costs 175.7 s at this
    scale, so `hmm_starts(n = 8)` is 24 minutes here*;
  - the Core seams row for `frm_multistart(fit, n, jitter)` should
    gain: *`hmm_starts()` shipped in frmtmb.latent 0.3.0 is the
    jitter-around-the-incumbent half only. Item 2.4 measured that half
    recovering 0 of 8 on the LCA design, where the modes are 280
    log-likelihood units apart; the core version needs a second family
    of starts that re-draws the family's own default, and a per-family
    declaration of WHICH parameters a restart should scatter. It also
    needs ONE definition of "is this a different optimum": the
    extension shipped three, two of them untoleranced, and the guard
    then fired on a correct model 6 of 6 times.*

## Every script this lane wrote, and what each measures

| script | what |
|---|---|
| `dev/latent-env.R` | the library stack every script here sources |
| `dev/latent-install.ps1`, `dev/latent-roxy.R` | roxygenise and install into the lane library |
| `dev/latent-runtests.R` | one test file per process, with the block, assertion and skip counts printed before the failure count |
| `dev/latent-hazard-check.R` | plants and removes a `$` hazard read, to prove `test-bracket-access.R` reaches the new file |
| `dev/latent-lca-sim.R` | item 2.4's simulator and the label-matching helpers |
| `dev/latent-lca-explore.R` | one replicate, to learn the shapes |
| `dev/latent-2p4-lca.R` | the 200-replicate run |
| `dev/latent-2p4-split.R`, `dev/latent-2p4-summarize.R` | the tables above |
| `dev/latent-2p4-modes.R` | the eight disagreeing replicates, one at a time |
| `dev/latent-2p4-docrecipe.R` | whether the recipe `?lca` prints actually works |
| `dev/latent-2p4-testfails.R` | the new lca identity assertions, seen failing |
| `dev/latent-2p4-startfix.R`, `dev/latent-2p4-startfix-fns.R` | why the old lca start failed, and the candidate that replaces it, on all 200 seeds |
| `dev/latent-2p4-shrink.R` | the shrink weight, swept over six values on all 200 seeds |
| `dev/latent-2p4-recheck.R`, `dev/latent-2p4-recheck-summarize.R` | the lca recovery table re-measured under the start that ships |
| `dev/latent-falsealarm-check.ps1`, `dev/latent-falsealarm-revert.R` | both guard tests, seen failing on copies with one fix removed at a time |
| `dev/latent-2p4-oos.R`, `dev/latent-2p4-shrink2.R`, `dev/latent-2p4-weight.ps1` | the shrink grid on 200 FRESH seeds with their own poLCA arm, and the tuning block's right-hand end |
| `dev/latent-2p4-weight-summarize.R` | the grid in sample and out, never pooled, against an adjudicated reference |
| `dev/latent-2p4-oldstart.R` | whether the recipe `?lca` prints for reproducing a 0.2.2 fit actually reproduces it |
| `dev/latent-2p4-pooled.R` | the coverage pooled over both blocks, per coefficient, and whether the shortfall is a tail or the whole distribution |
| `dev/latent-hmm-sim.R` | item 2.3's simulator, with and without the random effect |
| `dev/latent-2p3-repro81.R` | the 8.099-unit cold start, reproduced |
| `dev/latent-2p3-allfit.R` | `frm_allfit()` agreeing on the wrong answer |
| `dev/latent-2p3-scale-probe.R` | what scale a jittered start should use |
| `dev/latent-2p3-starts-probe.R` | `hmm_starts()` across five jitters and five seeds |
| `dev/latent-2p3-starts-cost.R`, `dev/latent-2p3-cost.ps1` | what `hmm_starts()` costs, with its control; the .ps1 waits for the arms so the measurement has the machine to itself |
| `dev/latent-2p3-hmmtmb-probe.R` ... `-probe4.R` | pinning the hmmTMB conventions, in the order they were found |
| `dev/latent-2p3-hmm.R`, `dev/latent-2p3-run.ps1`, `dev/latent-2p3-summarize.R` | the two hmm arms, run sequentially |
| `dev/latent-2p3-tpmgap.R` | why arm B's transition-matrix column reads 0.18 when its log-likelihood identity reads 1e-11 |
| `dev/latent-suite.ps1`, `dev/latent-check.ps1` | the final whole-suite run and the single `R CMD check` |
