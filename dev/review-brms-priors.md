# Review: the `wt-brms-priors` measurement lane

Reviewer's independent reproduction of the eight claims in
`dev/brms-priors-findings.md`. Written incrementally; the verdict is at
the end.

Environment reproduced from scratch: private library
`.../scratchpad/rp-lib`, frmtmb 0.50.0 installed from the worktree
(`R CMD INSTALL --library=`), brms 2.23.0, rstan 2.32.7,
StanHeaders 2.39.1, RTMB 1.9, testthat 3.3.2, R 4.6.1. Stan cache under
`.../scratchpad/rp-stan-cache`, empty at the start. `R_MAKEVARS_USER`
points at `.../scratchpad/rp-makevars`, whose `CXX17FLAGS` ends in
`-std=gnu++17`; without it rstan 2.32.7 with StanHeaders 2.39.1
compiles nothing.

## Claim 8, scope (checked first, since it bounds everything else)

CONFIRMED.

- `git -C <worktree> rev-parse HEAD` is `564e185`, the branch point, so
  every change is uncommitted working-tree state.
- `git diff --name-only 564e185`: `.gitignore`,
  `dev/brms-likelihood-tests.md`. Untracked:
  `dev/brms-priors-findings.md`,
  `tests/testthat/helper-brms-priors.R`,
  `tests/testthat/test-brms-priors.R`.
- No file under `R/` appears in either list. Confirmed: the lane
  changed no package code.
- `.gitignore` gains exactly `dev/stan-cache-priors/` and
  `dev/bp-makevars`, both scratch paths the lane names in its
  environment section.
- Main checkout `C:/Users/adf44/source/r/frmtmb` is at `5dfdd84` and
  `git status --porcelain` is empty, both at the start of this review
  and (re-checked) at the end. I did not touch it.

## Claim 2, default rows apply nothing

CONFIRMED on all seven shapes, and the finding is LARGER than the lane
states. Reproduced with my own script, not the lane's helper.

Handing `brms::get_prior()` output to `frm(prior = )`, per shape: every
live row has `source == "default"`, `fit$prior` is `NULL`, the drop
message fires, and `logLik()` equals the no-prior fit to 1e-10.

| shape | live rows | all source=="default" | fit$prior NULL | logLik == unpenalized |
| --- | --- | --- | --- | --- |
| S1 | 3 | TRUE | TRUE | TRUE |
| S2 | 2 | TRUE | TRUE | TRUE |
| S3 | 4 | TRUE | TRUE | TRUE |
| S4 | 1 | TRUE | TRUE | TRUE |
| S5 | 1 | TRUE | TRUE | TRUE |
| S6 | **6** | TRUE | TRUE | TRUE |
| S7 | 2 | TRUE | TRUE | TRUE |

### Correction to the lane's table

S6 has **six** live rows, not five. The findings summary table reads
`S6 | 5 | 2 | 3`; the correct row is `6 | 2 | 4`. The lane's own test
file pins a six-element status vector for S6
(`test-brms-priors.R`, "every default row's fate is one of five"), and
its own S6 section lists six rows, so the summary table contradicts the
rest of the lane. Minor, but it is a number in a table the maintainer
will read. **Punch-list item.**

### What the README sentence is true of, and what it is not

`README.md:61` says "a prior object brms itself built is translated".
Measured, on S1:

| object | `source` | entries reaching `fit$prior` | logLik |
| --- | --- | --- | --- |
| `brms::set_prior("normal(0, 10)", class = "b")` | `user` | 1 | -900.805166644 |
| `brms::prior(normal(0, 10), class = "b")` | `user` | 1 | -900.805166644 |
| `brms::get_prior(...)` whole table | `default` | 0 | -897.039321503 |
| no prior at all | - | 0 | -897.039321503 |

So the sentence is true of a `brms::set_prior()`- or `brms::prior()`-
built object, which carries `source == "user"`. It is false of the
object `get_prior()` returns, which is the one a porting user most
often has.

### The part the lane missed: an EDITED row is dropped too

This is the sharpest form of the finding and the lane does not state
it. `source` is set when brms BUILDS the row and is **not updated when
the user edits the `prior` cell in place**, which is the ordinary brms
workflow (`gp <- get_prior(...); gp$prior[i] <- "normal(0, 20)"; brm(prior = gp)`).

Measured on S1, editing the `sd` row to `normal(0, 20)`:

- the edited row's `source` is still `"default"` in brms's own frame;
- **brms honors it**: its `lprior` block changes from
  `student_t_lpdf(sd_1 | 3, 0, 59.3)` to `normal_lpdf(sd_1 | 0, 20)`;
- **frmtmb drops it**: `fit$prior` is `NULL` and logLik is unchanged at
  -897.039321503, against -898.925661682 for the same prior respelled
  through frmtmb's own `set_prior("normal(0, 20)", class = "sd")`.

And the message the user is shown is not merely uninformative, it is
wrong about the row it dropped:

> Translating a brms prior: dropped 3 row(s) brms had filled in as its
> own defaults.

One of those three rows was filled in by the user, not by brms. So the
user-facing consequence of D1 is not "brms's defaults are ignored"; it
is **"`frm(prior = )` drops rows by the `source` column, and `source`
does not track authorship after `get_prior()` returns"**. A user who
deliberately writes a prior into a `get_prior()` table loses it
silently and is told brms wrote it. That is a correctness bug, not an
ergonomics gap, and it raises D1's priority above where the lane put
it.

Also confirmed the asymmetry that makes the drop invisible: `brm()` and
`brm(prior = get_prior(...))` generate **byte-identical** Stan code, so
in brms passing the default table back is genuinely a no-op. That is
exactly why a porting user has no reason to suspect the row matters,
and why frmtmb's silent divergence is hard to notice.

## Claim 5, refusal misfires

CONFIRMED, all three, and the lane's classification of them is right.

**(a) `theta2` -- a BUG.** `check_brms_prior_class()` special-cases
`identical(cls, "theta")` (`R/priors.R:509`), and brms spells a mixture
proportion `theta2`, so the generic branch fires and emits:

> ... its frmtmb spelling is class = "Intercept", dpar = "theta2", and
> that density sits on the LINK scale where brms puts it on theta2
> itself.

A user who follows that advice verbatim gets
`Prior target not found (class=Intercept, dpar=theta2)`. Reproduced
directly. The bare class `"theta"` does get the correct special-cased
hint, which confirms the condition is simply too narrow. One line.

**(b) ordinal thresholds -- a documented limit, reported late.** The
class gate accepts `"Intercept"` (returns without error), and the
failure comes later from `resolve_priorlist()` as
`Prior target not found (class=Intercept)`. frmtmb genuinely has no
threshold slot, so this is a feature gap rather than a mistranslation;
the defect is only that the diagnosis arrives from the resolver with no
mention of thresholds. Also confirmed the lane's undocumented route:
`frm(prior = list(tau_raw = prior_normal(0, 5)))` resolves to **2**
entries on component `tau_raw` at `scale = "internal"` and fits.

**(c) `logistic(0, 1)` -- a documented limit, correctly reported.**
`parse_prior_dist()` refuses with
`Unsupported prior distribution 'logistic' (supported: normal,
student_t, cauchy, exponential, lkj)`. Swept the parser: `normal`,
`student_t`, `cauchy`, `exponential`, `lkj` parse; `logistic`, `gamma`,
`dirichlet`, `beta`, `uniform` do not. The message names the supported
set, so this one is honest and needs no fix beyond deciding whether to
add the density.

## Claim 3, the documentation half

CONFIRMED: frmtmb's documentation claims centering **nowhere**.
`man/set_prior.Rd:57` reads in full:

> `"Intercept"`: the intercept of `dpar`. Link scale.

Every `center`-adjacent string under `R/priors.R` is about a NONLINEAR
sub-formula not being centered (`:146`, `:755`, `:1065`), which is a
different subject. So the translation makes no false claim; it is
silent, which is why the gap is invisible. `man/set_prior.Rd:59` does
document class `"sd"` as "on the NATURAL sd scale with the
log-Jacobian applied", which supports the lane's D0.

## Claim 4, second half: the D2a-prime strong variant

The lane says the stronger variant -- flipping `class = "Intercept"`/
`"b"` carrying any `dpar` to the natural scale, one line at
`R/priors.R:1191` -- must be refused because it produces "Two hard
failures". **Reproduced on a scratch copy. The mechanism is confirmed;
the cost claim is half wrong, and the truth is worse.**

Method: `tar`-copied the worktree to `.../scratchpad/rp-strong`
(the worktree itself untouched), applied exactly the one-line change
the lane names:

```r
sc <- if (isTRUE(s$natural) || nzchar(s$dpar)) "sd" else "internal"
```

installed to `.../scratchpad/rp-lib-strong`, and ran each named test
file both ways, one file per process.

### `tests/testthat/test-simulate-ergonomics.R` -- CONFIRMED, hard failure

| build | result |
| --- | --- |
| unpatched worktree | passed 41, failed 0, error 0 |
| D2a-prime | passed 38, **error 1** |

```
Error ('test-simulate-ergonomics.R:120:3'):
  The prior on sigma_Intercept did not produce a draw inside [0, Inf]
  in 1000 tries
```

Exactly the lane's predicted mechanism: the `set_prior(normal(log(0.6),
1e-9), class = "Intercept", dpar = "sigma")` at lines 118-119 is
deliberately on the link scale, and on the natural scale every draw is
negative. The lane cites 118-119 for the `set_prior()` call; the error
surfaces at the `frm_simulate()` call on line 120. Substantively
correct.

### `extensions/frmtmb.ddm/.../test-surface.R` -- REFUTED as a failure

| build | result |
| --- | --- |
| unpatched worktree | passed 46, failed 0, error 0 |
| D2a-prime | passed 46, **failed 0, error 0** |

The test does **not** fail. But this is not a reprieve for D2a-prime,
because the lane's mechanism claim is right and the test simply misses
it. Probed directly on the same model:

| build | `R/priors.R:1191` reaches it? | entry scale | unpenalized `mu.cond` | penalized `mu.cond` | assertion `abs(pen) < abs(unpen)` |
| --- | --- | --- | --- | --- | --- |
| unpatched | n/a | `internal` | 0.9165381216 | **+0.2654012893** | TRUE |
| D2a-prime | yes | **`sd`** | 0.9165381216 | **-0.3301564534** | TRUE |

The patch does reach the case: the resolved entry's scale flips from
`internal` to `sd`. And the estimate behaves exactly as the lane
predicted -- the density on `exp(coef)` pushes the coefficient
negative, from a correct shrink-toward-zero at +0.265 to **-0.330**,
a sign flip. The test survives only because it asserts a magnitude
comparison, `abs(-0.330) < abs(0.917)`, which is still true.

So the correct statement for the maintainer is:

- D2a-prime causes **one** hard failure (`test-simulate-ergonomics.R`),
  not two;
- on the ddm surface it causes a **silent sign flip in a fitted
  coefficient that the existing test does not catch**.

That strengthens the case against D2a-prime rather than weakening it: a
change whose damage is invisible to the suite is worse than one that
stops it. It also flags a second, independent defect that is worth
fixing whatever is decided about D2 -- **`test-surface.R`'s dpar-prior
assertion is too weak.** Asserting only that the magnitude shrank lets
a sign flip through; it should assert the shrunk estimate stays between
0 and the unpenalized value. **Punch-list item.**

The lane's D2a-prime recommendation ("do not") stands, and my
recommendation is the same. Only its evidence line needs correcting.

## Claim 4, D2a's "no existing frmtmb spelling changes meaning"

CONFIRMED, and the lane's one specific warning about how to implement
it is confirmed too, by building both shapes.

I implemented D2a on two scratch copies of the worktree:

- **D2a-A** (the lane's shape): `set_prior()` gains `natural = FALSE`,
  and the field is written on the spec **only when TRUE**;
  `as_priorlist()` routes a brms distributional class to
  `class = "Intercept"`, `dpar = <class>`, `natural = TRUE`;
  `check_brms_prior_class()` keeps refusing only `theta*`.
- **D2a-B**: identical except the field is **always** written
  (`natural = isTRUE(natural)`), so a default spec carries
  `natural = FALSE` rather than no field.

Every file run one per process, against control (the unpatched
worktree) and each variant:

| test file | control | D2a-A | D2a-B |
| --- | --- | --- | --- |
| `test-setprior.R` | 27 pass | **27 pass** | - |
| `test-prior-compat.R` | 105 pass | 104 pass, **1 fail** | - |
| `test-priors-autocor-classes.R` | 63 pass | **63 pass** | - |
| `test-priors-bounds-grcov.R` | 49 pass | **49 pass** | - |
| `test-simulate-ergonomics.R` | 41 pass | **41 pass** | - |
| `frmtmb.sample/.../test-sample-direct.R` | 135 pass | **135 pass** | 134 pass, **1 fail** |

Findings:

1. **D2a-A's only casualty is the one the lane names.** The single
   `test-prior-compat.R` failure is at **line 424**, the
   `expect_error(as_priorlist(brms::prior(student_t(3, 0, 10), class =
   "sigma")), "LINK scale")` assertion. The lane predicted exactly this
   ("`test-prior-compat.R:424-426` greps `"LINK scale"` out of the
   refusal and goes"). Nothing else in the prior suite moves.
2. **The claim of no meaning change holds.** `test-setprior.R`,
   the autocor and bounds prior files, and `test-simulate-ergonomics.R`
   are all bit-identical in outcome. Note that
   `test-simulate-ergonomics.R` is the file D2a-prime breaks, and D2a
   leaves it alone -- which is the concrete difference between the two.
3. **The lane's implementation warning is real and I reproduced it.**
   D2a-B fails at `test-sample-direct.R:220:3` with
   `Expected `sg$natural` to be NULL.` The lane wrote that this test
   "breaks if the new argument defaults to `FALSE` rather than absent".
   Correct, and worth flagging to whoever implements D2a: write the
   field only when TRUE, or update that assertion deliberately.

Also confirmed by reading, as the lane states: `print.frmtmb_priorlist()`
already renders `scale=natural` (`R/priors.R:618`, and
`test-sample-direct.R:282` already greps for it);
`resolve_priorlist()` already reads the flag (`R/priors.R:1191`); and
`natural_dpar_prior()` (`frmtmb.sample/R/sample.R:524-528`) does by
hand what D2a would expose, so it does become redundant.

## Claim 1, the stale premise and the `sd` placement

CONFIRMED, to machine precision, by a route independent of the lane.

`R/priors.R:1307-1312` reads exactly as the lane quotes: for
`scale == "sd"` it sets `jac <- x`, replaces `x` with `exp(x)`, and
returns `prior_base_logdens(x, dist) + jac`. So frmtmb evaluates the
density at the natural sd and adds `log(sd)`.

What brms writes, read off `make_stancode()` output (text generation,
no compile):

```
S1  lprior += student_t_lpdf(Intercept | 3, 288.7, 59.3);
    lprior += student_t_lpdf(sigma | 3, 0, 59.3)     [- lccdf]
    lprior += student_t_lpdf(sd_1 | 3, 0, 59.3)      [- lccdf]
S3  ... plus lkj_corr_cholesky_lpdf(L_1 | 1);
S7  Intercept + sigma only
S5  sigma only,  1 lccdf line
```

The lower bound is written as the density minus one
`student_t_lccdf(0 | ...)` per element, and `lccdf(0)` for a
symmetric t at 0 is `log(1/2)`, so brms's statement is
`st(sd) + k log 2`.

Verified numerically, in both directions:

| quantity | value |
| --- | --- |
| `prior_logdens(log(36.67), student_t(3,0,59.3), "sd")` | -1.7214827919 |
| `st(36.67) + log(36.67)` | -1.7214827919 (equal, 1e-12) |
| `prior_logdens(log(36.67), ..., "internal")` | -5.0859563146 |
| `st(log(36.67))` | -5.0859563146 (the link spelling, no Jacobian) |
| brms half-t at `sd = (24.9060881, 5.9852937)` | -8.9017551850 |
| frmtmb at the same sd | -5.2836318600 |
| **frm - brms** | **+3.6181233250** |
| closed form `sum(log sd) - k log 2` | +3.6181233250 (equal, 1e-12) |

The lane's S3 worked example reports `+3.61812333`. I get
`+3.6181233250` by a separate derivation. **Match to nine figures.**

So the plan's premise really was stale: `class = "sd"` is already
brms's placement, differing only by the constant `k log 2`, which
moves no mode. D0 "nothing to decide" is correct.

## Claim 3, the intercept centering

CONFIRMED, mechanism and magnitude, again without rstan.

**Mechanism, from brms's own generated Stan code for `Reaction ~ Days`:**

```
matrix[N, Kc] Xc;                 // centered version of X without an intercept
vector[Kc] means_X;               // column means of X before centering
means_X[i - 1] = mean(X[, i]);
Xc[, i - 1] = X[, i] - means_X[i - 1];
real Intercept;                   // temporary intercept for centered predictors
lprior += student_t_lpdf(Intercept | 3, 288.7, 59.3);
real b_Intercept = Intercept - dot_product(means_X, b);
```

The prior is on the temporary intercept, i.e. the intercept at the
predictor means. `mean(Days) = 4.5`, so on `sleepstudy` the gap is
structural, not a sampling artifact.

**Magnitude.** I wrote brms's penalized objective by hand in R directly
from that Stan code and maximized it with `optim`, using no Stan and
none of the lane's helper:

| quantity | MLE | brms mode | frmtmb, same prior string |
| --- | --- | --- | --- |
| slope `Days` | 10.46728620 | 10.46728551 | **10.38298713** |
| raw intercept | 251.40509669 | 251.35914902 | 251.93899630 |
| centered intercept | - | 298.46193382 | 298.66243837 |
| sigma | 47.44890396 | 47.48793321 | 47.44976717 |

- **frmtmb - brms slope = -0.0842984.** The lane reports -0.08444059.
  Agreement to about 0.2 percent, which is optimizer tolerance: my
  hand objective's frmtmb slope is 10.38298713 against `frm()`'s own
  10.38302871, a relative difference of 4e-06.
- **brms - MLE slope = -6.9e-07**, i.e. nothing. The lane reports
  3.5e-05, also nothing; mine is tighter-converged. The qualitative
  claim -- brms's Intercept prior does not move the slope because the
  intercept it constrains is orthogonal to it, frmtmb's does -- is
  reproduced exactly.
- `SE(Days) = 1.231295522` from frmtmb's own `summary()`, matching the
  lane's 1.2313, so the shift is **0.0685 SE**. The lane says 0.069.
  CONFIRMED.
- The centered intercepts, 298.46193382 vs 298.66243837, match the
  lane's 298.461943 and 298.6627025.

Cross-checked: `frm()` with the same prior really does give the shifted
slope (10.38302871 against 10.46724855 unpenalized), so this is
frmtmb's behavior and not an artifact of my objective.

## Claim 4, first half: dpar placement on S5

CONFIRMED, and the lane's figures reproduce to seven significant
figures by the same Stan-free route.

| spelling | my sigma | lane's sigma |
| --- | --- | --- |
| MLE (no prior) | 0.1349066835 | 0.134906316 |
| frmtmb LINK | 0.1351048982 | 0.135104527 |
| frmtmb NATURAL | 0.1354694064 | 0.135469736 |
| brms mode (AT=TRUE) | 0.1354694061 | 0.135469737 |

- **natural - brms = 2.5e-10.** The lane says the natural placement
  "reproduces brms's mode to nine figures". CONFIRMED -- and it must,
  because the two objectives differ by the constant `log(2)`.
- **link - brms = -3.645e-04**, against the lane's 3.65e-04.
  CONFIRMED.
- **The link spelling captures 35.22 percent** of the shift; the lane
  says 35.2 percent. CONFIRMED.
- `a` and `b` agree to nine figures across all four columns
  (2.509173..., 0.786629...), so the placement moves the dispersion
  parameter and nothing else, exactly as claimed.

## Claim 7, part 1: do the pins actually watch the package?

**PARTLY REFUTED, and this is my most consequential finding about the
test file.**

The lane's stated justification for the whole file is:

> The new test pins today's behavior, so **whichever way each decision
> goes, `test-brms-priors.R` is where it fails first and says what
> changed.** Under D1 the test at its line 54 inverts; under D2a the
> status vectors and the `link` branch of the row-5 test change.

The D1 half is true. **The D2a half is not.**

`bp_classify_rows()` (`helper-brms-priors.R`) decides `"refused: class"`
and `"refused: distribution"` from its OWN hardcoded lists:

```r
bp_classes <- function() c("b", "Intercept", "sd", "cor")
bp_dists <- function() c("normal", "student_t", "cauchy", "exponential", "lkj")
...
g$status[[i]] <- if (!g$class[[i]] %in% bp_classes()) "refused: class"
                 else if (!g$kind[[i]] %in% bp_dists()) "refused: distribution"
                 else "honored"
```

It never calls `frmtmb:::check_brms_prior_class()` or
`frmtmb:::parse_prior_dist()`. It is a **parallel re-implementation of
the gate**, so the status-vector assertions watch the helper, not the
package.

Demonstrated on S1 against my D2a-A build:

| | control | D2a-A |
| --- | --- | --- |
| helper status vector | `honored, honored, refused: class` | **`honored, honored, refused: class`** (unchanged) |
| `frmtmb:::check_brms_prior_class("sigma", ...)` | refused | **ACCEPTED** |
| `as_priorlist()` on a user-sourced table | ERROR ("LINK scale") | **3 entries: `Intercept`, `sd`, `Intercept/sigma[natural]`** |

frmtmb's behavior changed completely; the pin did not move.

The consequence propagates through the whole file, because `bp_shape()`
derives its `honored` and `dpar` row sets from `bp_classify_rows()`
too, and every other test re-spells rows through `set_prior()` via
`bp_frm_prior()`, which bypasses the gate entirely. So under D2a I
expect the entire file to pass unchanged. (Verified below against the
warm cache.)

**This does not invalidate the measurement.** Every number the lane
reports is still reproducible -- I reproduced the important ones by
independent routes. What it invalidates is the lane's claim about what
the test file will *protect*. As written, `test-brms-priors.R` is a
faithful record of today's numbers and a genuine regression guard for
D1 and D3, but it is **not** a guard for D2.

Fix, and it is small: have `bp_classify_rows()` derive the class and
distribution verdicts by CALLING frmtmb, e.g.

```r
cls_ok <- !inherits(try(frmtmb:::check_brms_prior_class(g$class[[i]],
                                                       g$prior[[i]]),
                        silent = TRUE), "try-error")
```

and likewise `parse_prior_dist()` for the density branch. Then the
status vector is an observation of the package and the lane's stated
guarantee becomes true. **Punch-list item, and the one I would insist
on before merge.**

### Two coverage gaps

1. **The `theta2` message defect is pinned nowhere.** `theta2` appears
   in the test file only inside comments (lines 121, 168); no
   assertion reads the refusal text. So the D4 one-line fix the lane
   recommends would break no test and be signalled by nothing. One
   `expect_error(..., 'dpar = "theta2"')` would close it, and it also
   documents the defect where a maintainer will meet it.
2. **S7 is not in the test file at all.** The lane's own headline
   finding -- the centering biasing a regression slope by 0.069 SE,
   "the largest effect this lane measured" -- lives only in the
   findings document. The centering is pinned on S1 and S3 as a
   density difference inside the check-A decomposition, so the
   attribution is guarded; the user-facing consequence is not. Adding
   S7 costs three Stan programs, two of which the file already
   compiles for other rows.

`S2b` and `S5b` are likewise absent, but the lane says explicitly that
they exist to show the identities do not move with n, so omitting them
from the suite is a defensible choice rather than a gap.

## Claim 6, the residual decomposition

CONFIRMED, exactly, on S1, S3, S5 and S7. Reproduced with neither
rstan nor the lane's helper: frmtmb's log prior from its own two-tape
difference (`-obj$env$f(par)` with the prior minus the same without,
at the penalized fit's parameters), brms's `lprior` evaluated
analytically from the statements its generated Stan code actually
contains.

Every figure below matches the lane's table to all eight decimals it
prints.

| shape | quantity | mine | lane |
| --- | --- | --- | --- |
| S1 | frm log prior | **-7.03864998** | -7.03864998 |
| S1 | stan log prior | **-9.73573525** | -9.73573525 |
| S1 | frm - stan | **+2.69708527** | +2.69708527 |
| S1 | sd Jacobian - log 2 | +2.90884658 | - |
| S1 | centering | -0.21176131 | -0.2118 |
| S1 | **unattributed** | **4.4e-16** | 0.00000000 |
| S3 | frm log prior | **-11.30284344** | -11.30284344 |
| S3 | stan log prior | **-14.69839383** | -14.69839383 |
| S3 | frm - stan | **+3.39555038** | +3.39555038 |
| S3 | sd Jacobian - 2 log 2 | **+3.61812333** | +3.61812333 |
| S3 | centering | **-0.22050516** | -0.22050516 |
| S3 | LKJ coordinate | **-0.00206779** | -0.00206779 |
| S3 | **unattributed** | **0.0e+00** | 0.00000000 |
| S7 | frm log prior | **-5.32456812** | -5.32456812 |
| S7 | stan log prior | **-5.10222613** | -5.10222613 |
| S7 | frm - stan (pure centering) | **-0.22234198** | -0.22234198 |

Supporting estimates also match: S3's `sd = (24.90608813,
5.98529373)` against the lane's `24.906, 5.985`, and
`rho = 0.03711572` against `0.0371157`.

S5's honored set carries no frmtmb prior at all, so its residual is 0
by construction; its `nat` residual is `log(2)` because the natural
objective and brms's differ by exactly that constant, which is also
why the two share a mode (verified separately: `natural - brms` sigma
is 2.5e-10).

**The lane's caveat is correct and worth preserving.** The residual is
NOT the plain Jacobian sum. frmtmb carries a Jacobian only where a
prior asked for one; Stan's `jac` runs over every constrained
parameter whether priored or not. On S1's honored set the priored part
is `log(sd) = 3.6019897` while Stan's `jac` is `7.0324881`, and the
difference is `log(sigma)` -- sigma being constrained in Stan and
unpriored on both sides. That distinction is stated clearly in the
findings and is the kind of thing a later reader would otherwise get
wrong.

## D1's dependency on D2, checked

The lane writes that **"D1a is only safe after D2"**, because a default
table carries a `sigma` row that `check_brms_prior_class()` still
refuses, so honoring the table would turn a silent no-op into an error.

CONFIRMED, and it is broader than stated. I simulated D1a the way the
change would behave -- flipping `source` to `"user"` so the rows are
honored -- and handed each shape's whole `get_prior()` table to
`as_priorlist()`:

| shape | result today, if default rows were honored |
| --- | --- |
| S1 | **ERROR**: `class = "sigma" ... has no faithful frmtmb spelling` |
| S3 | **ERROR**: same |
| S7 | **ERROR**: same |
| S5 | **ERROR**: same (`student_t(3, 0, 2.5)`) |
| S4 | translates (1 entry), then **fails at resolve** with `Prior target not found (class=Intercept)` (shown earlier) |

So D1a shipped on its own would turn `frm(prior = get_prior(...))` from
a silent no-op into a **hard error on every shape I tested** -- four at
translation, the fifth at resolution. The ordering constraint the lane
identifies is real and is the single most important implementation note
in the findings document.

## Claim 7, part 2: running the file

The file runs green and the assertion count is exactly as claimed.

| run | cache state | wall | blocks | passed | failed | error | skipped |
| --- | --- | --- | --- | --- | --- | --- | --- |
| **cold** | **empty**, compiled all 15 | **1925.2 s** | 9 | **49** | 0 | 0 | 0 |
| completing run | 13 of 15 programs cached | 263 s | 9 | **49** | 0 | 0 | 0 |
| **warm** | 15 of 15 cached | **39.9 s** | 9 | **49** | 0 | 0 | 0 |
| **under D2a-A** | warm | 27.9 s | 9 | **49** | 0 | 0 | 0 |

- **49 assertions**, matching the lane's count exactly, and 49 is also
  the static count of `expect_*` calls in the file, so nothing is
  looped or conditional.
- The file compiles **15 Stan programs**, exactly as the lane says.
- Warm wall 39.9 s against the lane's 52 s; same order, and both
  figures are contention-dependent on this machine.
- Nothing skips once `FRMTMB_BRMS_FIT_TESTS=true` is set, and
  `skip_unless_brms_fit()` correctly skips the whole file when it is
  not.

**Cold from a genuinely empty cache: 1925.2 s**, compiling all 15
programs, still 49 passed and 0 failed. The lane reports ~3000 s. Both
are the same measurement under different load: my first cold attempt
ran while this machine carried up to ten competing R processes and
three `cc1plus` compiles from other lanes, and a single mixture program
took over an hour of wall time under that load; the run that produced
1925 s had the machine largely to itself. **The cold number is a
property of the contended machine rather than of the test file**, so I
would not hold the lane to 3000 s either way. What matters, and what is
confirmed, is that the file is green cold as well as warm, and that the
warm cost of 40 s is what CI pays on every run after the first.

### The decisive check: D2a

Running the file against my fully working D2a-A build -- where
`check_brms_prior_class()` accepts `sigma` and `as_priorlist()` emits
`Intercept/sigma[natural]` -- gives **49 passed, 0 failed: completely
unchanged.**

This confirms the prediction from the helper analysis above and
**refutes** the lane's claim that "under D2a the status vectors and the
`link` branch of the row-5 test change", and with it the broader claim
that the file is "where it fails first and says what changed"
whichever way each decision goes. For D1 and D3 that claim holds; for
**D2 the file is blind.**

Two causes, both in `helper-brms-priors.R`:

1. `bp_classify_rows()` re-implements the class and distribution gates
   from its own hardcoded `bp_classes()` / `bp_dists()` lists instead
   of calling frmtmb's, so the status vectors cannot move.
2. `bp_shape()` derives its `honored` and `dpar` row sets from those
   same statuses, and every other test re-spells rows through
   `set_prior()` via `bp_frm_prior()`, bypassing the gate entirely.

So the `link` and `nat` branches keep measuring exactly what they
measure today, whatever `frm(prior = )` now does with the same table.

### Does any pin hide a defect?

No pin asserts something false. Two pins assert defects **on purpose**,
and both say so in a comment:

- lines 54-80 pin the default-drop (`expect_null(fit$prior)` plus
  logLik equality). This is the D1 defect, deliberately recorded. It
  will invert under D1a, which is the intent.
- the S6 status vector pins `"refused: class"` for `theta2`, which is
  the misfire. But it pins only the CLASSIFICATION, not the misleading
  hint text, so the actual defect -- the message naming a spelling that
  does not exist -- is invisible to the suite.

The honest summary is that the file hides no defect, but it under-pins
two of them (the `theta2` message, and everything D2 would change) and
omits the shape that carries the largest measured effect (S7).

### Two lint-style files, for scope

| file | result | what it actually checks |
| --- | --- | --- |
| `test-message-uniqueness.R` | 6 passed, 0 failed | scans `R/` only |
| `test-bracket-access.R` | 7 passed, 0 failed | scans `R/` and `extensions/*/R` only |

Both scan package source, never `tests/`. Their passing confirms the
lane touched no `R/` code -- useful corroboration of claim 8 -- but
says nothing about the quality of the new test or helper.

## The plan document's follow-on section

Accurate against everything I reproduced. It records plainly that its
own prediction was wrong and why, which is the most useful thing it
could do, and its four numbered findings match my measurements. The
methodological note (that the mode-distance question is not answerable
on a random-effect shape, because Stan's joint optimum runs the
standard deviations away, `sd_1 = (110.2, 66.1)` against frmtmb's
marginal `(24.91, 5.99)`) is the right caveat and is consistent with
the flat-prior tier's own reasoning.

## Citation accuracy

I spot-checked every file:line the findings cite. They are precise:

| citation | verified |
| --- | --- |
| `R/priors.R:460` default drop | yes, exact line |
| `R/priors.R:507`, `:509` class gate and theta hint | yes |
| `R/priors.R:618` print renders `scale=natural` | yes |
| `R/priors.R:1191` the `sc <-` line | yes |
| `R/priors.R:1307` `prior_logdens()` | yes |
| `man/set_prior.Rd:59` natural sd scale documented | yes |
| `README.md:59-62` the migration sentence | yes |
| `test-prior-compat.R:424-426` LINK scale | yes, fails at 424 under D2a |
| `test-prior-compat.R:437-441` default-drop message | yes, at 440 |
| `frmtmb.ddm/.../test-surface.R:133` | yes, exact line |
| `test-simulate-ergonomics.R:118-119` | yes, error surfaces at 120 |
| `frmtmb.sample/.../test-sample-direct.R:220` | yes |
| `frmtmb.sample/R/sample.R:524-528` `natural_dpar_prior()` | yes |

This is unusually careful work and it made the review much faster.

---

# VERDICT: GO WITH FIXES

The measurement is sound. Every number I checked reproduced, most of
them exactly, and several by routes deliberately independent of the
lane's own helper: hand-written penalized objectives, analytic
evaluation of brms's `lprior` from its generated Stan text, and
`make_stancode` inspection. The premise correction is right, the four
decisions are the right four, and the costings are concrete and mostly
accurate. The maintainer can act on these numbers.

Three things must be fixed before it lands, one of them substantive.

## Punch list

**Must fix**

1. **`bp_classify_rows()` must ask frmtmb, not its own copy of the
   rules.** `helper-brms-priors.R` re-implements the class and
   distribution gates as hardcoded lists (`bp_classes()`,
   `bp_dists()`), so `test-brms-priors.R` passes **49/49 unchanged**
   against a working D2a build. The file's stated purpose, that it is
   where any of these decisions fails first, is false for D2. Route
   both branches through `frmtmb:::check_brms_prior_class()` and
   `frmtmb:::parse_prior_dist()`.
2. **Correct the two claims that rest on that.** In the findings,
   "under D2a the status vectors and the `link` branch of the row-5
   test change" is wrong as written, and so is the general guarantee in
   "What this lane did not do".
3. **Fix the S6 row count.** The summary table reads `S6 | 5 | 2 | 3`;
   it is `6 | 2 | 4`. The lane's own test file and its own S6 section
   both say six.

**Should fix**

4. **Correct the D2a-prime evidence.** It is **one** hard failure
   (`test-simulate-ergonomics.R`), not two. On the ddm surface the
   change is a **silent sign flip** (`mu.cond` from +0.265 to -0.330)
   that `test-surface.R` does not catch. That makes the case against
   D2a-prime stronger, not weaker, and it should be stated that way.
5. **Report the edited-row case under D1.** The drop is keyed on
   `source`, and brms does not update `source` when a user edits the
   `prior` cell of a `get_prior()` table in place. brms honors such an
   edit; frmtmb silently drops it and reports it as a row "brms had
   filled in as its own defaults", which is false for that row. This
   is the sharpest and most damaging form of D1 and it is missing from
   the findings.
6. **Pin the theta2 message.** The defect is reported but no assertion
   reads the refusal text, so the recommended one-line fix would be
   signalled by nothing.

**Nice to have**

7. **Add S7 to the test file.** The lane's own largest measured effect,
   0.069 SE of slope bias, is pinned nowhere.
8. **Strengthen `test-surface.R:133-135`** to bound the sign, not only
   the magnitude, of a dpar-prior shrink. As written it admits a sign
   flip.

## My recommendation on the four decisions

**D1 (default rows): take D1a, but only after D2, and treat it as a bug
fix rather than an ergonomics improvement.** The lane costs this as a
surprise; my measurements say it is worse. `source` records who *built*
a row, not who *wrote* it, so the ordinary brms workflow of editing a
`get_prior()` table in place loses the user's own prior silently, under
a message that misattributes it to brms. Meanwhile `brm()` and
`brm(prior = get_prior(...))` generate byte-identical Stan code, so
nothing in the user's brms experience warns them. The cost is one
condition in `as_priorlist()` plus threading the calling path through
six known call sites, and two tests invert by design. The sequencing
constraint is real and I verified it: D1a shipped alone turns today's
silent no-op into a **hard error on every shape I tested**, four at
translation and the fifth at resolution. If D2 is not wanted, D1b
(message to warning) is a poor substitute, because the warning would
still misattribute the row; at minimum the message text must stop
claiming brms wrote it.

**D2 (dpar placement): take D2a, refuse D2a-prime.** The natural
placement already exists, is already what `frm_sample()`'s own defaults
use, is already rendered by `print()`, and reproduces brms's mode to
nine figures on S5, while the link spelling captures about a third of
the intended shift on S5 and moves sigma the wrong way by a
ten-thousandth on `sleepstudy`. I verified the "no meaning change"
claim by building D2a: the prior suite is unchanged except the single
assertion that greps the obsolete refusal text. One implementation
note, confirmed by building both shapes: write the `natural` field only
when TRUE; a version that always writes `natural = FALSE` fails
`test-sample-direct.R:220`. D2a-prime must be refused. The lane's
mechanism is right even though its evidence line is half wrong, and the
ddm sign flip being *silent* is the strongest argument against it.

**D3 (Intercept centering): take D3b now, and do not let the number
disappear into a vignette.** This is the only one of the differences
that biases a regression coefficient, at 0.0685 SE on 180 rows, and
frmtmb's documentation currently says only "the intercept of `dpar`.
Link scale.", which is not wrong but tells a porting user nothing. D3a
is genuinely structural, a density on a linear combination of parameter
slots with its own gradient, and I agree it is out of scope for this
lane. But the documentation should carry the measured slope shift, and
the migration vignette should say outright that a brms `Intercept`
prior does not mean the same thing after porting. Add S7 to the suite
at the same time so the number cannot rot.

**D4 (refusal misfires): take the theta2 fix now; treat the rest as
documented limits.** The theta2 hint is an outright bug, naming a
spelling (`class = "Intercept", dpar = "theta2"`) that then fails with
`Prior target not found`. One line, worth taking whatever happens to
D2, and it should arrive with the assertion that is currently missing.
The `logistic`/`gamma`/`dirichlet` gap is honestly reported by the
parser and is a feature request rather than a defect; `logistic` and
`gamma` are one `switch` arm each if wanted. The ordinal threshold gap
is a real feature gap whose only defect is that the diagnosis arrives
from the resolver as a bare `Prior target not found (class=Intercept)`
with no mention of thresholds. The `frm(prior = list(tau_raw = ...))`
route works today (I confirmed two entries at `scale = "internal"`) and
is documented nowhere, which is worth a line in `set_prior()`'s help
even if nothing else changes.

## Edits I made to the worktree

**One file, added; nothing modified, nothing deleted:**

- `dev/review-brms-priors.md` (this file).

I changed no other file in the worktree, nothing under `R/`, and
nothing in the main checkout, which is still at `5dfdd84` with a clean
`git status`. Every experiment ran on `tar` copies under
`.../scratchpad/rp-strong`, `rp-d2a-A` and `rp-d2a-B`, against the
private library `.../scratchpad/rp-lib`. I committed nothing.
---

# Punch re-check, 2026-09-05

Second pass over the punch round, verified against the diff and the
files rather than the worker's prose, reusing my own Stan cache
(`.../scratchpad/rp-stan-cache`) and my own scratch builds. Main is
still `5dfdd84` and clean; worktree HEAD is still `564e185`; `R/` is
still untouched, and the only new tracked change is
`extensions/frmtmb.ddm/tests/testthat/test-surface.R` (item 8).

**Verdict: GO.** All eight items land. Three residual items below, one
of them substantive and new.

## Baseline, measured before anything else

| run | wall | blocks | passed | failed | error | skipped |
| --- | --- | --- | --- | --- | --- | --- |
| worktree, first run | 393.3 s | 11 | **61** | 0 | 0 | 0 |
| worktree, warm | **33.6 s** | 11 | **61** | 0 | 0 | 0 |
| `test-message-uniqueness.R` | - | - | 6 | 0 | 0 | 0 |
| `test-bracket-access.R` | - | - | 7 | 0 | 0 | 0 |
| ddm `test-surface.R`, control | - | - | **47** | 0 | 0 | 0 |
| ddm `test-surface.R`, D2a-prime | - | - | 46 | **1** | 0 | 0 |
| `test-simulate-ergonomics.R`, control | - | - | **41** | 0 | 0 | 0 |
| `test-simulate-ergonomics.R`, D2a-prime | - | - | 38 | 0 | **1** | 0 |

All one file per process.

## Item 1: the helper now asks the package. VERIFIED, and the proof
## reproduces independently

Structurally, in `tests/testthat/helper-brms-priors.R`:

- `bp_class_ok()` (`:28-31`) calls `frmtmb:::check_brms_prior_class()`;
- `bp_dist_ok()` (`:33-36`) calls `frmtmb:::parse_prior_dist()`;
- `bp_row_priorlist()` (`:44-49`) calls `frmtmb:::as_priorlist()` on the
  single row with `source` flipped to `"user"`, and
  `bp_classify_rows()` (`:86-91`) uses it for the target question;
- `bp_classes()` and `bp_dists()` are **gone** (no hits in the file).

`bp_frm_prior()` (`:105-116`) additionally falls back to
`bp_row_priorlist()` when `set_prior()` has no name for the class,
which is the right delegation: it asks the translator for the spelling
instead of guessing one.

**Is the `source`-flipped `as_priorlist()` route the same path
`frm(prior = )` takes, or a third path?** The same path. Traced:

- `as_priorlist()` dispatches on `inherits(x, "brmsprior")`
  (`R/priors.R:434`), and `bp_row_priorlist()` builds
  `structure(row, class = c("brmsprior", "data.frame"))`, so it enters
  the brms-translation branch and no other.
- `frm()` calls `as_priorlist(prior)` at the argument boundary
  (`R/fit.R:440`) and then `resolve_prior_input(...)` at
  (`R/fit.R:681`). The helper calls the same two functions in the same
  order.

The only difference is the value of the `source` cell, which is
precisely the field D1 governs, and the helper says so in its comment.
Not a third path.

**The proof, rebuilt from my own D2a shape.** I did not use the
worker's `psp-d2a-lib`. I re-copied the *current* worktree to
`.../scratchpad/rp-d2a-A2` and re-applied my own review-round patch
(`set_prior()` gains `natural = FALSE` written only when TRUE;
`check_brms_prior_class()` refuses only `theta*`; `as_priorlist()`
routes a distributional class to `Intercept` + `dpar` + `natural`):

| build | blocks | passed | failed | error |
| --- | --- | --- | --- | --- |
| worktree | 11 | 61 | 0 | 0 |
| **my D2a-A2** | 11 | **35** | **15** | **3** |

**35 / 15 / 3, matching the worker's reported numbers exactly.** The
file that passed 49 of 49 unchanged against D2a before the punch round
now moves 18 assertions. Item 1 is real.

### Are the 18 the ones that SHOULD move?

Yes. All 18, enumerated with `SummaryReporter$new(max_reports = 100)`:

| # | line | what moved | D2 consequence |
| --- | --- | --- | --- |
| 1 | 141 | S1 status vector, `sigma` no longer `refused: class` | direct |
| 2 | 163 | S6 status vector, `sigma1`/`sigma2` now honored | direct |
| 3 | 202 | S5 status vector | direct |
| 4 | 203 | `abs(r$hon$dF) < 1e-8` | honored set gained sigma, so the honored fit is now penalized |
| 5 | 204 | `r$hon$gF < 1e-3` | same |
| 6 | 207 | `r$hon$gT == 1` | same |
| 7 | 218 | `r$link$frm_prior` equals the link density | the `link` variant no longer exists |
| 8 | 223 | ERROR on `expect_gt(r$link$gF, 0.01)` | `r$link` is NULL |
| 9 | 256 | ERROR, C++ `std::runtime_error` | `r$nat$pars` is NULL |
| 10 | 321 | S1 check-A decomposition | honored set grew |
| 11 | 334 | ERROR on `expect_lt(r$nat$gFz, 1e-8)` | `r$nat` is NULL |
| 12 | 379 | S3 check-A with the LKJ term | honored set grew |
| 13 | 390 | `bp_half_t_count(r$code$hon)` vs `length(sd1)` | the honored Stan program gained sigma's lccdf line |
| 14 | 452 | S7 status vector | direct |
| 15 | 494 | S7 check-A decomposition | honored set grew |
| 16 | 553 | mixture `ent$comp == c("beta","beta")` | `betad` entries appear |
| 17 | 554 | mixture `unique(ent$scale) == "internal"` | `"sd"` appears: the placement itself |
| 18 | 562 | mixture `betad` sum | the nat set changed |

**None moves for an unrelated reason.** I chased the only ambiguous
one, the bare C++ error at `:256`, to its cause rather than assuming
it. Probed on S5 under both builds:

| | control | D2a-A2 |
| --- | --- | --- |
| row status | `refused: class` | **`honored`** |
| `honored` idx | none (n=0) | 1 |
| `dpar` idx | 1 (n=1) | **none (n=0)** |
| `r$link` / `r$nat` | present | **NULL** |

`bp_shape()` builds the `link` and `nat` variants only
`if (length(dpi))`. Under D2a the class gate honors the row, so `dpi`
is empty, so `r$nat` is never built, so `:256`'s
`rstan::unconstrain_pars(r$sf$full, r$nat$pars)` gets `NULL` and rstan
raises `variable does not exist ... variable name=b_a`. That is a D2
consequence reached by a legitimate chain. Same for `:223` and `:334`.

Presentation nit, not a defect: those three surface as a raw C++ error
and two `Result of comparison must be TRUE, FALSE, or NA` messages
rather than as a sentence naming the cause. A reader of a future CI
failure would have to do the tracing I just did. A guard such as
`expect_false(is.null(r$nat))` ahead of them would make the intent
legible for one line each.

## Items 2-8, checked against the files

| item | claim | verdict |
| --- | --- | --- |
| 2 | two findings claims corrected, old wording quoted | **VERIFIED** at `dev/brms-priors-findings.md:688-697`; the quoted old wording matches what I reported in the review round, and the mechanism (hardcoded lists) is stated correctly. One number wrong, below. |
| 3 | S6 is `6 \| 2 \| 4` in three places | **VERIFIED**: `findings:98` (table), `findings:481` ("six live defaults"), `test-brms-priors.R:536` (comment) |
| 4 | D2a-prime re-measured | **VERIFIED**: `test-simulate-ergonomics.R` 41 control vs 38 + 1 error at `:120`; ddm `mu.cond` +0.2654 to -0.3302 with the entry scale flipping `internal` to `sd` (reproduced in the review round and again here) |
| 5 | edited-row case pinned | **VERIFIED** at `test-brms-priors.R:82-119`: asserts `source` stays `"default"`, the `"dropped 3 row"` misattribution, `expect_null(fit$prior)`, logLik equality with the unpriored fit, and that the respelled prior moves the objective by more than 1 nat. It pins today's wrong behavior and says so in the comment. |
| 6 | theta2 refusal text pinned | **VERIFIED** at `test-brms-priors.R:175-176`: `expect_error(frmtmb:::check_brms_prior_class("theta2", "logistic(0, 1)"), 'dpar = "theta2"', fixed = TRUE)` |
| 7 | S7 added | **VERIFIED** at `test-brms-priors.R:439-498`; pins `abs(frm_days - mle_days) / se_days == 0.0684` (tol 0.02) and the centering identity `raw + mean(Days) * slope == centered`. My own independent figure last round was 0.0685 SE. "Zero new compiles" is wrong, below. |
| 8 | ddm assertion bounds the sign | **VERIFIED** in the diff: `test-surface.R:134-142` now computes `shrink <- pen / unpen` and asserts `expect_gt(shrink, 0)` and `expect_lt(shrink, 1)`. Control 47/0; under D2a-prime 46/1, failing at `:141:3` with `Expected shrink > 0` (shrink = -0.3302/0.9165 = -0.36). This is the fix I asked for and it works. |

## Residual items

**R1. `dev/brms-priors-findings.md:684` says "6 of the 11 blocks"; it is
7.** Measured per block under my D2a-A2 build, the blocks carrying a
failure or an error are:

| # | block | failed | error |
| --- | --- | --- | --- |
| 3 | every default row's fate is one of five, by shape | 2 | - |
| 4 | row 5: a dpar prior is the whole placement question | 5 | yes |
| 5 | row 5: the natural placement reproduces brms's mode | 0 | yes |
| 6 | row C: class sd is brms's placement, up to log(2) | 1 | yes |
| 7 | row C: class cor is the same LKJ in another coordinate | 2 | - |
| 9 | S7: the centering is what biases a regression slope | 2 | - |
| 11 | row 17: a mixture keeps its intercepts and loses the rest | 3 | - |

Seven blocks; 15 failures plus 3 errors is the 18. Blocks 5 and 6 are
easy to miss because block 5 has zero failures and only an error. One
word.

**R2. Item 7's "zero new compiles" is true only against a cache that
already holds the lane's own measurement programs.** I started this
re-check with the 15 programs the pre-punch file needed. The first run
of the new file compiled **three more** (cache 15 to 18, wall 393.3 s,
the three `.rds` written at 10:19, 10:21 and 10:23); the second run was
33.6 s with no new compiles. S7 genuinely adds three Stan programs,
because `Reaction ~ Days` with no random effect shares no program with
S1 or S3. The lane's own findings record "S7 (3 programs) | 723 s" from
the measurement phase, which is where the worker's cache got them.
`dev/brms-priors-findings.md` (the timings paragraph near `:792`)
should say the file now needs **18** programs rather than "no new
compiles", because 18 is what a cold CI pays.

**R3 (new, substantive, and a cost D2 has not been charged).
Naive D2a routing passes brms's `lb` through onto the link-scale
parameter and pins it.** Found while checking why my D2a build warned
about convergence on S5. Measured on the nonlinear shape:

| spec | sigma | max abs gradient |
| --- | --- | --- |
| hand-built `natural = TRUE` (the lane's measurement) | **0.135469736** | 3.5e-05 |
| D2a-routed from the brms row | **1.000000000** | **117** |
| D2a-routed with `lb`/`ub` dropped | **0.135469736** | - |

brms's `sigma` row carries `lb = "0"`, a bound on **sigma itself**.
Routing the row to `class = "Intercept", dpar = "sigma",
natural = TRUE` keeps that bound, but frmtmb applies a bound to the
**internal** parameter, which here is `log sigma`. So `lb = 0` becomes
`log(sigma) >= 0`, i.e. `sigma >= 1`, and the fit pins at exactly
1.000000000 instead of 0.135469736.

This is not a defect in the lane's work: D2a is a recommendation, not
an implementation, and the lane's own measurements used a hand-built
spec that carries no bound. But it is a real constraint on
implementing D2a, it is invisible in the current costing, and it hits
**every** brms dispersion default, since they all carry `lb = 0`.
Whoever implements D2a must transform the bound (`log(lb)` for a log
link) or drop it, and should pin that with a test. Worth one line in
D2's cost column.

It also slightly qualifies the D2a proof: a few of the 18 moving
assertions compare numbers from a fit whose sigma is pinned at 1. The
*conclusion* is unaffected, because the six status and scale
assertions (`:141`, `:163`, `:202`, `:452`, `:553`, `:554`) are pure
gate observations and move regardless of what the fit does.

## Standing items from the review round

Both remain open by design and neither blocks the merge:

- **S2b and S5b** are still absent from the suite. The lane's stated
  reason (the identities do not move with n) still holds.
- **D3a** is still not costed, which the findings say plainly.

## Edits I made in this round

**None.** I appended this section to
`dev/review-brms-priors.md` and changed no other file. All builds were
`tar` copies under `.../scratchpad/rp-d2a-A2`, `rp-strong2`, against
the private library `.../scratchpad/rp-lib`. Nothing committed; the
main checkout was not touched and is still `5dfdd84`, clean.
