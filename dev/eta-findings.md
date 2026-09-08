# eta lane: a public API for the linear-predictor scale

Round 0.55. Worktree `frmtmb-wt-eta`, branch `wt-eta`, based on
`a2deb65`.

## The problem, restated with the measurement

`build_objective()` stores each distributional parameter's linear
predictor beside the parameter under `.eta_<dpar>`. Five accessors read
it, and all five were `@noRd`. A family defined in another package can
SEE the entry in `dpars` and had no sanctioned way to use it, so
`frmtmb.coupling` wrote the arithmetic again in `cw_complement()`.

Why it must be read at all, measured on the `cross_wishart()` log
density at `n = 16`, one process:

| eta | `1 - plogis(eta)` | naive log density | accessor log density |
| --- | --- | --- | --- |
| 10 | 4.54e-05 | -38741.9441435 | -38741.9441435 |
| 30 | 9.36e-14 | -18873016199287 | -18853765065131 |
| 36 | 2.32e-16 | -7945539820084449 | -7606151693197552 |
| 40 | 0 | NaN | -4.1528181132e+17 |
| 100 | 0 | NaN | -4.74254897384e+43 |

`plogis(eta)` is exactly 1 from `eta = 36.7368005696771` up. The naive
form is not merely inaccurate before it dies: at `eta = 30` it is wrong
by 1.0e-3 relative, and at `eta = 36` by 4.5e-2.

## 1. What is exported, and what is not

Exported: `dpar_log()`, `dpar_log1m()`, `dpar_log_complement()`,
`dpar_complement()`. Not exported: `robust_logit()`, `robust_logmu()`.

Four, not three. The first draft exported only the pair-returning
`dpar_log_complement()` for the unit-interval case, and that made
`frmtmb.coupling` compute a `log(C)` it never uses. The conclusion
survived review; **the measurement behind it did not, and the number
reached three shipped documents before it was caught.** See
"I measured the wrong quantity" below.

The cost that repeats is the TAPED GRADIENT SWEEP. Measured at 1000
rows, six interleaved blocks, medians:

| what | base `a2deb65` | one-sided (shipped) | pair |
| --- | --- | --- | --- |
| full `cross_wishart()` sweep, 2000 each | 462 us | 480 us | 633 us |
| accessor alone, 3000 each | | 205 us | 387 us |
| AD nodes recorded, accessor alone | | 10002 | 14002 |

The pair is 1.32 times the one-sided on the full density and 1.89
times on the accessor alone, and the gap clears the run-to-run spread:
the pair's fastest block, 620 us, is above the one-sided's slowest,
540 us. The one-sided form is level with the hand-written arithmetic it
replaces, 480 against 462, which is the number that matters for
`frmtmb.coupling`: adopting the public API costs it nothing.

The node count is the mechanism and is worth more than the timing. The
pair records 4 extra AD nodes per row for the half the density
discards, R's `list()` forces both components so laziness cannot skip
it, and RTMB exposes no simplify or optimize method on a `Tape`, so the
dead half is replayed on every sweep for the life of the fit.

The line is TOTALITY, not usefulness. An exported accessor always
returns a usable value, on the tape and off it, and folds the phase
branch in. `robust_logit()` and `robust_logmu()` return `NULL` off the
tape by design, which forces the caller to write the two-branch density
by hand. That two-branch density is the exact trap this task names: the
`NULL` branch appears only off the tape, so a family author who tests
only there ships the wrong on-tape expression and never sees it.

For each one not exported:

* `robust_logmu()`. A family author gets the same number from
  `dpar_log()`, which adds only the `log(dpars[[dpar]])` fallback the
  author would otherwise write. The one extra thing the `NULL` carries
  is whether the value is EXACT, and the only use core makes of that is
  to switch between two whole density expressions
  (`RTMB::dnbinom_robust()` against `RTMB::dnbinom2()`, and a closed
  form for `log P(0)` against a `logspace_add` form). That switch is a
  core-internal optimization between two spellings of the same
  quantity. Writing the robust spelling unconditionally is correct in
  both regimes, because `dpar_log()`'s fallback puts the plain log
  there.
* `robust_logit()`. Same argument. Every use a density has for a log
  odds is `log(p)`, `log(1 - p)` or the pair `(p, 1 - p)`, and those
  are the two exported accessors. Handing back a raw log odds that is
  `NULL` half the time buys nothing a family cannot get totally.

Keeping them internal also keeps the door open: if a family outside
frmtmb ever needs the exact-or-not distinction, exporting one of them
later is additive. Withdrawing a `NULL` from a public return is not.

`dpar_log()` also gained a branch it did not have as `log_dpar()`. A
unit-interval link carries no `log_eta`, only a `logit_eta`, so the
internal version fell through to `log(dpars[[dpar]])` for such a dpar.
That made one public function answer two different ways: exact at
`eta = -800` on a log link, and `-Inf` on a logit, where
`stats::plogis(-800)` is exactly 0. It now takes `log(x)` off the log
odds when that is what the link offers. No core call site changes,
because `dpar_log()` is called only for `shape` and `phi`, whose links
are positive-support ones.

## 2. The names, and why

Chosen: `dpar_log()`, `dpar_log1m()`, `dpar_log_complement()`,
`dpar_complement()`. Was: `log_dpar()`, `gate_logs()`, `mu_pair()`, and
nothing at all for the one-sided complement.

The package's exported extension vocabulary is unprefixed noun phrases
that open with the object they read: `dpar_linpred()`, `frame_block_of()`,
`fit_extras()`, `response_mean()`, `eval_dpars()`, `single_response()`.
`frmtmb_*` is reserved for constructors and registries (`frmtmb_family()`,
`frmtmb_structure()`, `frmtmb_register_*()`), and `frm_*` for verbs on a
fitted model (`frm_lp_basis()`, `frm_joint_cov()`). These three read a
dpar, so they take the `dpar_` opening that `dpar_linpred()` already
established.

Point by point:

* `log_dpar` to `dpar_log`. Same words, ordered the way the rest of the
  surface orders them: object first, quantity second. `log_dpar` read as
  a verb on a dpar, which is what made it look private.
* `gate_logs` to `dpar_log_complement`. "Gate" is jargon for `zi` and
  `hu` and the function was never limited to those; its own internal
  documentation said "or any other dpar on the unit interval". The
  quantity is a value and its complement, on the log scale, and
  `frmtmb.coupling` had already named its private copy `cw_complement()`
  without having read this code, which is evidence the word is the
  natural one.
* `mu_pair` to `dpar_complement`. `mu_pair` named a dpar (`mu`) in a
  function that takes the dpar name as an argument, so it was wrong for
  every call that passed anything else. The pair is the value and its
  complement; naming it after the complement says which half you called
  for.
* `dpar_log_complement()` is `dpar_complement()` on the log scale, and
  the names say so. The two are the pair a density needs, and the `log_`
  is the only difference between them.
* `dpar_log1m()` is the second half of `dpar_log()`, and `1m` for "one
  minus" is the package's own existing convention: `log1m_inv_logit()`
  is the internal it is built on, and `l1m` is the field name it
  matches inside `dpar_log_complement()`. So the two one-sided
  functions and the two fields of the pair carry the same two names.

Rejected: `*_probs` (`dpar_probs()`), because the package already has a
`*_probs` family (`latent_probs()`, `mixture_probs()`, `lca_probs()`,
`hmm_probs()`) which all return a matrix of class probabilities from a
FIT. A two-element list from a dpar list under the same suffix would
read as one of those.

The argument order was made uniform: `(dpars, dpar, link)` for all
three. It was `(dpars, name, link)` for two and
`(dpars, link, name = "mu")` for the third, which is not a signature to
put in public.

Field names were NOT renamed. `l` and `l1m` follow the package's
existing `log1m_inv_logit()` convention, `p` and `q` are standard for a
probability and its complement, and the two pairs use different letters
so a reader can never mistake which accessor's result is in hand.
Renaming them would have touched about 30 call sites for a documentation
gain the help page already delivers.

### I measured the wrong quantity, and how

My first justification for four exports timed the R-level `cw_lpdf()`
call: "1.9 ms against 1.2 ms per call at 1000 rows". That number is
real and it is irrelevant, because **the R density runs once per fit.**
Counted on a fitted `cross_wishart()`, the lpdf is entered once by
`frm()` to build the tape and once more by `summary()`, which tapes the
report; `fitted()`, `predict()`, `residuals()`, `logLik()`, `confint()`,
`vcov()`, `AIC()` and both coherence extractors enter it zero times.
A cost paid twice in the life of a model is not an argument for
anything. What repeats thousands of times is `$jacobian()` on the tape,
and the R code is not on the tape: only the arithmetic it recorded is.

Three things led me there, all avoidable:

1. **I profiled the thing I had edited.** I changed R code, so I timed
   R code. The edit's whole purpose was to change what gets RECORDED,
   and recording and replaying are different budgets.
2. **I never counted the calls.** One `frm()` fit with a counter in the
   density would have shown "1" and stopped the R-call measurement
   dead. It takes a minute and I did not spend it.
3. **The number looked decisive, so I stopped.** 1.9 against 1.2 is a
   large ratio and it pointed at the conclusion I already believed. A
   measurement that agrees with you deserves the same scrutiny as one
   that does not, and it did not get it.

The check that generalizes: before quoting a cost, say what repeats it.
If the answer is "once per fit", it is a setup cost and belongs in
prose, not in a performance table. For anything inside a density the
unit is the gradient sweep.

A second lesson, cheaper: my first sweep measurement ran each variant
to completion in turn and put base at 1110 us against one-sided at
550 us, for two tapes that record the same nodes. That gap was drift,
not arithmetic. Interleaving the blocks removed it. Wherever the
difference under test is tens of percent, interleave.

### The link is an argument, and it is not optional

`link` takes the dpar's OWN link, as the name the family wrote in
`links = ` or as a link object. `get_link()` is not exported, so a family
outside frmtmb usually holds only the name; refusing a name would have
made the API unreachable from the place that needs it most.

Cost of accepting a name, 50000 replicates, one process:

| | us/call |
| --- | --- |
| `get_link("logit")` | 0.800 |
| `dpar_log_complement()` at n = 1, link object | 78.0 |
| `dpar_log_complement()` at n = 1, link name | 76.0 |
| `dpar_log_complement()` at n = 500, link object | 572.0 |
| `dpar_log_complement()` at n = 500, link name | 538.0 |

The name path costs one registry lookup, 0.8 us, which is below the
run-to-run spread of the accessor itself: the measured name-versus-object
difference came out NEGATIVE at all three sizes. A link object still
passes straight through with no work, which is the path every family
inside frmtmb takes.

(Unrelated observation, not a regression: the accessor itself costs 78 us
at one row, dominated by two `RTMB::logspace_add()` calls. That is the
pre-existing cost of the internal accessors, unchanged by this lane.)

## 3. Is `.eta_<dpar>` API?

No. The accessors are the only sanctioned route, and the reference topic
says so in a section of its own, "The `.eta_` entries are not the API".

Two reasons, both concrete:

* The entry is on the LINK scale, so its meaning depends on the dpar's
  link. `-logspace_add(0, e)` is `log(1 - x)` on a logit and is nothing
  at all on an identity link; `e` is `log(x)` on a log link and is not
  on a softplus. This is not hypothetical: core's own `log_dpar()`
  carried a comment recording that it once assumed the log link
  unconditionally and would have "quietly returned the wrong density"
  once `link_shape` became an argument, and the repository has two tests
  (`test-families.R:200` and `:221`) pinning exactly that bug for a
  positive dpar and for a gate.
* Whether the entry is PRESENT is a property of the phase, not of the
  model. Reading it directly means writing the two-branch density, and
  the branch that only runs on the tape is the one nothing tests.

I found and fixed one live instance of the first hazard inside this
repository. `inst/bcm/binomial-extras.R`, `bcm_contaminant()`, read
`dpars[[".eta_phi"]]` and applied `-logspace_add(0, eta)` to it while
taking `link_phi` as a constructor argument. At the default
(`link_phi = "logit"`) it is right; at `link_phi = "identity"`, which the
unit-interval link set allows, it would have read a probability as a log
odds and fitted a different density with no message. It now calls
`dpar_log_complement(dpars, "phi", link_phi)`.

`vignettes/inputs.Rmd` already described `.eta_` as a reserved prefix
that a data column cannot collide with; that stays, and the new topic is
the other half of the same statement.

## 4. Both paths, in the documentation and in the tests

The reference topic has a section "On the tape and off it" that says the
entries exist only while the objective is taped, that the fallback is
correct there because nothing off the tape is differentiated, and that a
family author must test both, with the recipe for building the on-tape
list by hand.

I also found the reason a family author gets this wrong.
`check_custom_family()`, the package's own family-checking entry point,
tapes the lpdf with respect to the DPAR VALUES and supplies no linear
predictors at all. Every family checked through it exercises only the
fallback branch, which is the plain arithmetic, which is the thing that
saturates. Its help page now says so and points at the new topic. I did
not extend the function itself to sweep the eta scale: that needs the
links and a per-dpar eta grid, it is a change to a public function's
behavior, and it belongs to whoever owns the family-checking surface.

Tests, in `tests/testthat/test-dpar-eta-api.R`, cover both paths for
every accessor:

* value on the tape and off it for all three accessors, at
  `eta` in `{-700, -40, -20, -3, 0, 3, 40, 700}`;
* the off-tape path asserted to be NOT finite at `eta = 40`, which is
  the trap itself under test;
* the wrong-link hazard for a positive dpar (softplus) and for a gate
  (identity);
* `dpar_log()` on a unit-interval dpar, exact where the plain form is
  `-Inf`, which is the branch that function gained;
* `dpar_log1m()` against the closed form at
  `eta` in `{-700, -40, 0, 40, 700}`, asserted IDENTICAL to
  `dpar_log_complement()$l1m` so the two spellings cannot drift;
* a link given as a name string, equal to the same link given as an
  object, and an unknown name refused;
* value AND gradient through `RTMB::MakeTape()` at `eta` in
  `{-40, 0, 40}`, against the closed-form derivative;
* the fallback branch reached ON the tape (identity link), so that it is
  differentiated too;
* end to end through `frm()`: a custom family that records whether it
  saw `.eta_mu`, fitted against `bernoulli()` for the same log
  likelihood and coefficients.

### The test failing before the change

Run against the base commit's installed package, `test-dpar-eta-api.R`
gave 10 failures and errors, the first being

```
Error: 'dpar_log_complement' is not an exported object from
  'namespace:frmtmb'
```

**Weak evidence, and I should have said so.** Every one of those ten is
the same fact, a missing symbol, counted ten times. A file full of
calls to functions that do not exist yet fails whatever the functions
would have done, so this shows the API is new and nothing about
whether it is right. The behavior comparison that carries the weight is
the base-versus-lane table of `log_dpar()` against `dpar_log()` and
`gate_logs()` against `dpar_log_complement()` at the saturating ends,
which the review supplied and which the test file already encodes
inline (`expect_identical(log(stats::plogis(-800)), -Inf)`,
`expect_identical(log(1 - 1e-17), 0)`). Next time: when a failing-first
run is all missing symbols, say the count is one fact and go find the
value comparison.

After the change: 123 passing assertions, 0 failures, which includes
the two blocks the review added.

## The equivalence measurement

`extensions/frmtmb.coupling/R/cross-wishart.R`, `cw_complement()`, now
reads `dpar_log1m(dpars, "coh", "logit")` instead of
`-RTMB::logspace_add(0 * eta, eta)` off `dpars[[".eta_coh"]]`. It is
two lines and the whole helper.

Compared on the FULL `cross_wishart()` log density, not only on the
helper, with the old body kept verbatim beside the new one. 961 values
of the coherence linear predictor: 10 to 100 in steps of 0.25, then 101
to 700 in steps of 1. `n = 16`, `y = 10.5`, `w22 = 10`, `w12 = 9.9 +
0.3i`, `S11 = 1.2`, `S22 = 1.1`, `phase = 0.4`.

**Maximum difference: 0. The two are bit-identical at all 961 points,
and every value is finite.**

| eta | old | new |
| --- | --- | --- |
| 10 | -3.874194e+04 | -3.874194e+04 |
| 30 | -1.885377e+13 | -1.885377e+13 |
| 36.75 | -1.610222e+16 | -1.610222e+16 |
| 40 | -4.152818e+17 | -4.152818e+17 |
| 100 | -4.742549e+43 | -4.742549e+43 |
| 400 | -9.212052e+173 | -9.212052e+173 |
| 700 | -1.789373e+304 | -1.789373e+304 |

Off the tape, where both fall back, they agree exactly at
`C` in `{0.05, 0.25, 0.5, 0.9}` and differ by 2.27e-13 at `C = 0.999` on
a log density of -1693.41653, which is 1.3e-16 relative, under one unit
in the last place. The difference is in `1 / (1 - C)` alone: the old code
spelled it `1 / (1 - C)` and the accessor's caller now spells it
`exp(-log(1 - C))`. I measured the alternative before accepting it:

| spelling | on tape, max rel diff in `inv` | off tape, max rel diff |
| --- | --- | --- |
| `exp(-l1m)` (shipped) | 0 | 9.09e-16 |
| `1 / dpar_complement()$q` | 2.2e-16 | 0 |

Neither is exact on both, so the choice is which path to be exact on,
and the two paths are not comparable in weight. I first justified the
tape spelling as "the path a fit runs on", which is true and too weak.
The sharper reason, counted rather than assumed: I instrumented
`cw_complement()` in the installed namespace and drove fourteen
post-fit entry points against a fitted model.

| path | calls |
| --- | --- |
| `frm()`, on the tape | 1 (then every gradient sweep, taped) |
| `summary()`, on the tape | 1 |
| `residuals(type = "deviance")`, off the tape | 1 |
| `fitted()`, `predict()`, `residuals(type = "response")`, `residuals(type = "pearson")`, `logLik()`, `confint()`, `vcov()`, `AIC()`, `frm_coherence()`, `frm_phase()` | 0 |

`simulate()` never arrives either; it hits the family's own refusal
first. So the off-tape surface is ONE method, and the quantity it
computes is a diagnostic residual, evaluated once. The tape spelling is
exact on the arithmetic the optimizer walks thousands of times and
costs 4 ulp on a number that is reported once and rounded before anyone
reads it. That is the trade, and stated that way it is not close.

### The fallback now uses log1p, and my first reason for it was wrong

`frmtmb.coupling` spelled its off-tape fallback `log1p(-C)` where
`gate_logs()` spelled it `log(1 - p)`, so I changed the accessor to
`log1p(-p)` and wrote in the commit that the subtraction loses digits as
`p` approaches 1. **That reason is false, and the test I wrote to pin it
failed.** `1 - p` is EXACT for `p` in `[0.5, 1]` by Sterbenz's lemma, so
the two spellings are bit-identical at the saturating end. Measured, 12
values:

| p | `log(1 - p)` | `log1p(-p)` | identical |
| --- | --- | --- | --- |
| 1e-18 | 0 | -1.0000000000000001e-18 | no |
| 1e-16 | -1.1102230246251565e-16 | -9.9999999999999998e-17 | no |
| 1e-08 | -1.0000000100247594e-08 | -1.0000000050000001e-08 | no |
| 1e-04 | -0.00010000500033334732 | -0.00010000500033335834 | no |
| 0.25 | -0.2876820724517809 | -0.2876820724517809 | yes |
| 0.999 | -6.9077552789821359 | -6.9077552789821359 | yes |
| 1 - 1e-16 | -36.736800569677101 | -36.736800569677101 | yes |

Over 800 points spread logarithmically toward both ends of `(0, 1)`,
they differ at 397, and every one of those has `p < 0.1`. The largest
relative difference is 1, at `p = 1e-18`, where the subtraction rounds
`1 - p` to 1 and `log(1 - p)` returns exactly 0 for a value of -1e-17.

So `log1p(-p)` stays, for the OTHER boundary: a gate the optimizer has
pushed toward 0, which happens as readily as one pushed toward 1, and
which the fallback branch does reach on the tape whenever the gate's
link carries no `logit_eta` (`identity` is in the allowed set for `zi`,
`hu` and `quantile`). It is never worse and is sometimes better. The
test now pins that, in both directions.

`log1p` has an RTMB tape method (checked: value and gradient both
correct through `MakeTape` at `x = 0.3`), so the fallback is still
differentiable where it is reached on the tape. There is a test for
that too.

The 2.27e-13 off-tape difference in the coupling comparison is NOT from
this term. It is entirely from `1 / (1 - C)`, measured separately in the
table above.

## The frmtmb version floor on frmtmb.coupling

**`extensions/frmtmb.coupling/DESCRIPTION` still says
`Depends: frmtmb (>= 0.53.0)` and that is now WRONG.** The floor has to
rise to the frmtmb release that exports `dpar_log1m()`. I did
not edit it, because naming that release means inventing a core version
number and consolidation sets versions. This is the one thing in this
lane that is deliberately left undone; it is a one-line edit once the
core version for this round is fixed.

A build against an older frmtmb fails at INSTALL time, not at run time,
because the name is in `importFrom` in `NAMESPACE`. That is the right
failure mode and it is stated in the coupling NEWS entry.

`RTMB` also leaves `frmtmb.coupling`'s `Imports:`. The
`-RTMB::logspace_add(0 * eta, eta)` call the accessor replaces was its
only use anywhere in the package, and an `Imports:` entry nothing
imports from is an `R CMD check` NOTE.

## Verification

Test files, one per R process, `NOT_CRAN=true`, private library.

Core (`frmtmb`):

| file | pass | fail | error | warn | skip |
| --- | --- | --- | --- | --- | --- |
| test-dpar-eta-api.R (new) | 113 | 0 | 0 | 0 | 0 |
| test-families.R | 216 | 0 | 0 | 0 | 0 |
| test-custom-family.R | 99 | 0 | 0 | 0 | 0 |
| test-numerical-robustness.R | 679 | 0 | 0 | 0 | 0 |
| test-interop.R | 44 | 0 | 0 | 0 | 0 |
| test-bcm-psychophysics.R | 7 | 0 | 0 | 0 | 2 |

`frmtmb.coupling`, whole suite, one file per process:

| file | pass | fail | error | warn | skip |
| --- | --- | --- | --- | --- | --- |
| test-cross-wishart.R | 211 | 0 | 0 | 0 | 0 |
| test-coherence.R | 58 | 0 | 0 | 3 | 0 |
| test-cross-spectrum.R | 61 | 0 | 0 | 0 | 0 |
| test-message-uniqueness.R | 4 | 0 | 0 | 0 | 0 |
| test-surface.R | 58 | 0 | 0 | 1 | 0 |

Both suites also ran inside `R CMD check`, in one process each, as that
command runs them: core `Running 'testthat.R' [243s] OK`, coupling
`Running 'testthat.R' [23s] OK`.

### R CMD check, final tree after review

Rebuilt and rechecked on the tree as it now stands, with pandoc and
TinyTeX on the PATH.

* **Core: Status 2 NOTEs.** CRAN incoming "New submission"; HTML manual
  "Skipping checking math rendering: package 'V8' unavailable".
  `checking tests ... Running 'testthat.R' [240s] OK`,
  `checking examples ... [33s] OK`,
  `checking examples with --run-donttest ... [36s] OK`,
  `checking re-building of vignette outputs ... [272s] OK`,
  `checking PDF version of manual ... OK`.
* **`frmtmb.coupling`: Status 1 WARNING**, the CRAN-incoming one, which
  the base commit produces identically. Everything else OK, including
  the PDF and HTML manuals and `Running 'testthat.R' [18s] OK`.

Both match the pre-review runs and the base-commit baselines. The
`profile.frmtmb_fit` example timer that appeared once did not return on
either of the two later runs, nor on the reviewer's independent one;
three OK results against one NOTE, so it was contention.

### R CMD check, first round

`roxygen2::roxygenise()` runs clean on both packages.

`R CMD check --as-cran` on core, run twice on this tree:
**Status: 2 NOTEs** the first time and **3 NOTEs** the second, the
extra one being a timer.

* `checking CRAN incoming feasibility ... NOTE`: "New submission".
  frmtmb is not on CRAN. Another lane's baseline check of the same
  package carries the identical NOTE.
* `checking HTML version of manual ... NOTE`: "Skipping checking math
  rendering: package 'V8' unavailable". This is the one the round was
  told to expect.
* `checking examples ... NOTE` in the SECOND run only, naming
  `profile.frmtmb_fit` at 10.69s elapsed against the 5s threshold. It
  is a load-dependent timer and not mine: the same check on the same
  tree reported `checking examples ... [42s] OK` twelve minutes
  earlier, and profile likelihood has nothing to do with these
  accessors. Six lanes were checking packages on this machine at once.

Everything else is OK in both runs, including `checking tests ...
Running 'testthat.R'` at [302s] and [278s] OK, `checking re-building of
vignette outputs` at [359s] and [248s] OK, and `checking PDF version of
manual` OK.

`R CMD check --as-cran` on `frmtmb.coupling`: **Status: 1 WARNING**,
and the SAME package at the base commit `a2deb65`, checked in the same
environment, gives **Status: 1 WARNING** with the identical text. The
warning is `checking CRAN incoming feasibility`: "New submission",
"Strong dependencies not in the CRAN or BioC software repositories:
frmtmb", and a 301 on the pkgdown URL in `DESCRIPTION`. All three are
properties of an extension package living in this repository and none
is introduced here. The changed package and the base package now check
to the same status, line for line.

An earlier run of the changed package added a `checking examples`
timing NOTE naming `frm_coherence` at 7.54s and `cross_wishart` at
6.41s. It did not reproduce, here or on the baseline, which is what a
load-dependent timer does on a machine running six lanes.

An earlier run of the same check reported `1 ERROR, 1 WARNING, 3 NOTEs`
and the error read "pdflatex is not available". That was my shell, not
the package: TinyTeX lives at
`C:/Users/adf44/AppData/Roaming/TinyTeX/bin/windows` and was not on the
PATH the lane's scripts built. With it on the PATH the manual builds,
and `R CMD Rd2pdf` on the source tree confirms it separately: exit 0,
156 pages. Recording it because the same trap will catch the next lane
that builds its own environment.

The four testthat warnings are pre-existing and were shown to be,
rather than assumed. The base commit `a2deb65` was extracted read-only,
installed into a second private library of mine, and the same two files
run against it: identical pass counts (58 and 58), identical warning
counts (3 and 1), identical text, including the same
`Large maximum absolute gradient at the optimum (0.00117)`. They belong
to two deliberately hard models in that suite.

## Round two, after review

The review (`dev/reviews/2026-09-08-eta.md`) applied four fixes itself
and left me three. All three are done, and I reproduced every number
rather than copying it.

* **The published timing measured the wrong quantity.** Corrected in
  all three places it had reached, `?"frmtmb-robust-dpars"`, `NEWS.md`
  and the coupling's `cross_wishart()` Rd, to quote the per-gradient-
  sweep cost and say what it is a cost of. Post-mortem above under
  "I measured the wrong quantity, and how". My sweep numbers (462 base,
  480 one-sided, 633 pair; accessor alone 205 against 387, ratio 1.89)
  differ in level from the review's (325, 340, 460; 240 against 413,
  ratio 1.72) because the machine load differs, and agree on the ratio
  and the ordering, which is the part the argument rests on. I also
  counted the recorded tape, 10002 nodes against 14002, which is the
  mechanism and does not depend on load at all.
* **The off-tape surface was overstated.** Reproduced by instrumenting
  `cw_complement()` and driving fourteen post-fit entry points:
  `residuals(type = "deviance")` is the only one that reaches it, and
  `frm()` and `summary()` are the only on-tape callers, once each. The
  reasoning for the tape-exact spelling is rewritten around that.
* **`bcm_contaminant()` has its own NEWS bullet.** Reproduced the whole
  table: bit-identical at the default logit across 15 linear predictors
  from -700 to 700 (max difference 0), and wrong elsewhere by +0.0303
  nats at cauchit, -0.0965 at probit, -0.1174 at cloglog, all at
  `eta = 2`. I also measured the case the review did not: `identity`,
  which is the one non-logit link the framework's own unit-interval set
  admits, is the worst at **+1.5716 nats per observation** at
  `phi = 0.05`. All four links construct: I called
  `bcm_contaminant(link_phi = )` for each and got a family back.

Two lower findings closed while there: `frmtmb-extension-api` now
points at the new topic and warns that `dpar_linpred()` is a different
kind of `dpar_`, and the new topic says the raw log odds is
reconstructible as `dpar_log() - dpar_log1m()` (worst absolute error
5.0e-17, reproduced), which is why the surface is complete at four.

Not acted on: the NEWS heading spelling, where four of five worktrees
agree with mine and consolidation will settle it.

## Files touched

Core:

* `R/families.R`. The three accessors, exported and documented as the
  new topic `frmtmb-robust-dpars`; `robust_logit()` and
  `robust_logmu()` kept `@noRd` with the reason written down; a new
  internal `robust_link()`; every call site renamed; `frmtmb_family()`'s
  `lpdf` parameter and `@seealso` pointed at the topic.
* `R/objective.R`. The comment beside the `.eta_` write now names the
  public accessors and says the entry is reserved.
* `R/interop.R`. `check_custom_family()` says it exercises the fallback
  path only.
* `inst/bcm/binomial-extras.R`. `bcm_contaminant()` calls the accessor
  instead of reading `.eta_phi`, which also fixes its assumed-logit bug.
* `vignettes/frmtmb.Rmd`, "Custom families": a paragraph, a runnable
  two-line demonstration of the trap and the fix, and the warning about
  `check_custom_family()`.
* `NAMESPACE`, `man/frmtmb-robust-dpars.Rd` (new),
  `man/frmtmb_family.Rd`, `man/check_custom_family.Rd`,
  `_pkgdown.yml` (the topic added under Families), `NEWS.md` (a new
  development heading).
* `tests/testthat/test-dpar-eta-api.R` (new).
  `tests/testthat/test-families.R` and
  `tests/testthat/test-numerical-robustness.R`: comments only, renaming
  the functions they name.

`frmtmb.coupling`:

* `R/cross-wishart.R`. `cw_complement()` is two lines over
  `dpar_log1m()`; the family's "Why the links are the constraint"
  section and the helper's own notes rewritten, including why the
  one-sided accessor and not the pair.
* `R/frmtmb.coupling-package.R`. The "What core supplies and what it
  does not" section no longer says the seam is unpromised, because it
  now is; `importFrom(RTMB, logspace_add)` replaced by
  `importFrom(frmtmb, dpar_log1m)`.
* `DESCRIPTION`: `RTMB` out of `Imports:`. The `frmtmb` floor is
  deliberately NOT touched, see above.
* `vignettes/coherence.Rmd`: the one paragraph that named
  `-logspace_add(0, eta)`.
* `NAMESPACE`, `man/cross_wishart.Rd`,
  `man/frmtmb.coupling-package.Rd`, `NEWS.md` (a new development
  heading), `tests/testthat/test-cross-wishart.R` (one new structural
  test).

## What I did not do

* Did not extend `check_custom_family()` to sweep the eta scale. Stated
  above; it is a behavior change to a public function and it needs a
  design decision about where the etas come from.
* Did not rename the returned field names `l`, `l1m`, `p`, `q`. Stated
  above.
* Did not export `robust_logit()` or `robust_logmu()`. Stated above.
* Did not touch `extensions/frmtmb.coupling/DESCRIPTION`'s frmtmb floor
  or any `Version:` field.
* Did not need to audit the other six extensions for hand-rolled
  copies. I grepped: `.eta_` appeared outside core only in
  `frmtmb.coupling` and `inst/bcm/binomial-extras.R`, and both now go
  through the API. What is left of the string in `frmtmb.coupling` is
  in its TESTS, where a hand-built on-tape `dpars` list is exactly the
  use the new topic sanctions.
* Did not relitigate `bcm_contaminant()`'s link choice, though my
  change weakens its stated reason. Its comment says the contamination
  rate takes a logit "where the book puts it on a probit" and that "the
  reason is numerical". The probit link also carries `logit_eta` (exact
  to `|eta| = 38.2`, where `pnorm` itself underflows), so with the
  accessor a probit gate would be exact too. The link choice also moves
  the group prior, and the file says the port's Stan program encodes
  the logit, so this is a modelling decision with a test behind it and
  not mine to change. Someone should revisit the sentence.
* Did not change `dpar_complement()`'s off-tape `q` to anything other
  than `1 - x`. There is no better spelling on the natural scale, and
  the measurement above says the two candidate spellings of `1/(1-C)`
  trade one ulp for another rather than one dominating.
