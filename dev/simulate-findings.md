# The simulator seam: what was measured, and what it changed

Working notes from the round that added
`tests/testthat/test-simulate-density.R`. Everything below is a
measurement or a decision, not a plan.

## What the gap actually was

The round started from a claim that "about nine of the 45 built-in
family constructors declare no `sim` slot". Measured at the entry
point, by building every family in `frmtmb:::family_registry` and
calling `frm_simulate()` on it, the number is **three**:

| Family | Status before |
| --- | --- |
| `tweedie` | no simulator |
| `compois` | no simulator |
| `hurdle_poisson` | no simulator |
| `cox` | no simulator, and says why (deliberate) |
| `categorical` | **not a gap** |

`categorical()` reads as a gap if you inspect the registry, because
the registry entry is a DEFERRED placeholder: the category count is a
property of the data, so `categorical()` carries a `defer` slot and no
simulator. Every entry point calls `resolve_deferred_families()`
before asking whether the family can simulate, and the resolved family
has always had a `sim`. For `categorical()` the registry is the wrong
measurement; calling `frm_simulate()` is the right one.

The nine is a STALE number, not a misreading of the registry. At
df9d6bc the registry reads five by every slicing: 40 entries, 32
exported zero-argument constructors, and the same five without a `sim`
(`tweedie`, `compois`, `hurdle_poisson`, `cox`, `categorical`). Nine is
the state of an earlier tree. At the 2026-08-31 commit b492cf0
("mixtures, mi, cs, rr") `R/families.R` carries 35 registry entries and
26 `sim` slots, which is nine short. The number was carried forward
from that round rather than re-measured.

`cox()` is the only remaining refusal. It states its reason: inverting
a cumulative baseline hazard identified only on the observed time
window has no answer past the last event.

## The defect the tier found

The ddm extension was described as simulating "only through its own
`ddm_simulate()` and `gddm_simulate()`, outside the family seam".
That is not what the code says. `wiener()` passes `sim` in its
constructor, `lba()` does too, and `gddm()` installs `sim` in
`gd_finalize()`. All three were wired.

What was broken was in **core**. `frm_simulate()` read
`spec$responses[[1]]$family` straight from `parse_spec()`, before
`assemble_frame()` had run the family's `family_finalize()` slot. A
family that derives itself from the response is not itself until
assembly, so:

- `gddm()` was refused by `frm_simulate()` for having no simulator;
- `simulate()` on a **fit** of the same model worked, because a fit
  carries the finalized family;
- `frm()` had the same problem once and already fixes it, at
  `R/fit.R:507-510`, by carrying the finalized responses over from the
  frame.

The fix in `R/simulate-new.R` is the same three lines. Nothing in
`extensions/frmtmb.ddm` needed changing. The regression is pinned in
two places: a core test using a `custom_family()` that installs its own
`sim` in `family_finalize()`, and
`extensions/frmtmb.ddm/tests/testthat/test-simulate-density.R`.

## The three new simulators

Each states a generative process. None inverts the density, which is
the whole point: a simulator written from the density cannot disagree
with it.

**`tweedie`.** The `power12` link guarantees `1 < p < 2`, and on that
range the Tweedie IS a compound Poisson sum of gamma jumps. Draw
a Poisson count of jumps, then a gamma of that total shape. The
reparameterization that puts the mean at `mu` and the variance at
`phi * mu^p` is

```
lambda = mu^(2-p) / (phi (2-p))     alpha = (2-p)/(p-1)
scale  = phi (p-1) mu^(p-1)
```

The point mass at zero falls out rather than being added: a row that
draws no jumps has gamma shape exactly zero, and `rgamma()` returns an
exact zero for shape zero (checked). Verified against the density at
`mu = 2.4, phi = 1.2, p = 1.5` on 2e5 draws: mean 2.397 against 2.4,
variance 4.470 against 4.462, and `P(Y = 0)` 0.07614 against both
`exp(-lambda)` and `exp(lpdf(0))` = 0.075623.

**`hurdle_poisson`.** Clear the hurdle with probability `1 - hu`, then
draw a count that cannot be zero. The positive part is drawn by
inverse transform on the Poisson CDF above its own zero,
`qpois(p0 + u (1 - p0), mu)`, not by drawing Poisson variates until one
is positive. Both are correct; only this one has bounded cost, and the
rejection loop degenerates exactly where hurdle models are used, at a
small `mu` where nearly every plain draw is the zero the hurdle already
accounts for.

**`compois`.** The COM-Poisson has no constructive generative process
the way a Poisson or a gamma does; every published sampler works from
the unnormalized weights `lambda^y / (y!)^nu`. This one does the same,
and solves for the rate that puts the distribution's MEAN where it was
asked for, by bisection on an expanding bracket. That solve is the part
worth doing independently: it is what confirms that
`RTMB::dcompois2()`'s second argument is the mean rather than the mode.
Checked at `(mu, nu)` of `(2, 1)`, `(2, 1.5)`, `(4, 0.6)` and `(1, 3)`:
maximum absolute probability difference against the family's own
density is between 1.4e-16 and 2.1e-13.

Rows sharing a `(mu, nu)` pair share one solve, because the solve is
the expensive part and a design has far fewer distinct pairs than rows.

## The tier

`tests/testthat/test-simulate-density.R`. 435 assertions, about 11
seconds.

Design, and why each piece is there:

- **The reference is the family's own `lpdf`, evaluated numerically.**
  RTMB's `d*` functions are numeric on numeric input, so the density
  the objective tapes can be run in plain doubles. Nothing is compared
  against a closed form written beside the simulator, because a
  convention error the simulator and the density SHARE is exactly what
  such a comparison cannot see.
- **Two cells.** The design is `y ~ 0 + g` with a two-level factor, so
  a simulator that reads `dpars[["mu"]][1]` and recycles it fails
  instead of passing on the first cell's parameters.
- **One representation for every support.** The density becomes a
  discrete measure: exact pmf on an integer support, trapezoid weights
  on a fine grid for a continuous one. Cells, moments and the mass
  check all read off that one object.
- **Cells of roughly equal probability**, cut from the measure's own
  cumulative sum, with integer supports cut at half-integers so a draw
  can never land on a boundary.
- **An atom gets its own statistic.** A zero-inflated or hurdle family
  is tested as `P(Y = 0)` against `exp(lpdf(0))` by a binomial
  deviate, plus a goodness of fit on the nonzero draws against the
  continuous part. Folding the atom into a cell dilutes exactly the
  error a reversed gate produces.
- **Every continuous family asserts its density integrates to one**
  over the bracket used. That is what makes a bracket clipping real
  mass a failure rather than a silent bias.

### Tolerances

One false-alarm rate, `SIM_ALPHA = 1e-6`, drives the chi-square
critical value and the two-sided normal deviate for the moments.

The argument is about maintenance, not statistics. The seed is fixed,
so a false alarm would not be flaky: it would be permanent, and
somebody would have to work out whether a correct simulator or a real
defect was on the other end. Buying that costs almost no detection. On
eight cells the critical value at this alpha is about 41 on 7 degrees
of freedom, while a simulator drawing from the wrong distribution
misses the cell probabilities by whole percentage points and lands a
statistic in the hundreds at `SIM_N = 4000` draws per cell. The margin
between 41 and "hundreds" is what the small alpha buys, and it is why
the number is one constant rather than a value tuned per family.

The variance check uses the delta-method standard error
`sqrt((mu4 - sigma^4) / n)`, with `mu4` taken from the same numeric
measure. `SIM_N` is 4000 per cell, so the mean tolerance is about
`0.077` sample standard deviations.

### It has teeth

Mutation-tested, in a scratch harness, by patching a family's `sim` in
the namespace and re-running the tier. Nine plausible convention slips,
**all nine caught**:

| Mutant | Failures raised |
| --- | --- |
| gaussian: sd/variance swap | 19 |
| gaussian: recycles cell 1's mean over every row | 19 |
| Gamma: shape/scale swap | 4 |
| zero_inflated_poisson: gate read the wrong way round | 6 |
| hurdle_poisson: plain Poisson instead of zero-truncated | 2 |
| tweedie: dispersion left out of the jump rate | 8 |
| compois: `nu` ignored (draws a Poisson) | 4 |
| binomial: `trials()` ignored, draws a bernoulli | 6 |
| cumulative: category order reversed | 4 |

The unpatched tier raises zero. Worth re-running that harness after any
change to the machinery: a tier that cannot fail is worse than no tier,
because it reads as coverage.

### The extension tier

`extensions/frmtmb.ddm/tests/testthat/test-simulate-density.R`, 41
assertions. The shape differs because a choice-RT density is
DEFECTIVE: `lpdf` gives the joint density of (time, boundary), so it
integrates to the probability of that boundary. The `sim` slot draws a
time conditional on the boundary the row already observed. So the
reference is the defective density renormalized on the row's own
boundary, and the two boundaries' masses are asserted to sum to one,
which is the check that no probability has gone missing. Measured for
`wiener()` at `mu = 1.1, bs = 1.5, ndt = 0.2, bias = 0.45`: 0.196881
lower and 0.803119 upper, summing to 1. An earlier draft of this note
gave 0.157095 and 0.842905 at that point. Those are the masses of
`mu = 1.2, bs = 1.4, ndt = 0.25, bias = 0.5`, the parameters of the
simulate-on-a-fit test at line 141 of the same tier file, so the pair
was transcribed from the wrong test.

`gddm()`'s density reads a `.gddm` aterm that frame assembly builds
from the response. A density evaluated on a grid of candidate
responses has to build it the same way, or it scores the wrong rows;
the tier calls the family's `aterm_data()` slot on its own grid. Its
modeled time window is also fixed at finalize time and the density is
undefined outside it, so the tier names `t_max` through
`gddm_control()` rather than inheriting whatever the largest simulated
response happened to be.

## The sweep

### Vignettes

Measured by knitting each vignette five times before and five times
after, in fresh processes, interleaved. Wall-clock timings on this
machine are noisy by a factor of two or more, so the MINIMUM of five is
the honest statistic; a single before/after pair suggested a large
speedup that was entirely cache warming.

| Vignette | Before (min of 5) | After (min of 5) |
| --- | --- | --- |
| `frmtmb.Rmd` | 10.5 s | 10.3 s |
| `inputs.Rmd` | 5.5 s | 5.6 s |

Timing-neutral. The conversion is about the source reading as one
statement of the model, not about speed.

Converted:

- `vignettes/frmtmb.Rmd`, the mixture chunk. This is the conversion
  worth having: `groups = ~g` means the class belongs to the group, and
  `frm_simulate()` draws it that way, which the hand-written
  `rbinom(40, ...)` had to spell out separately from the model.
- `vignettes/inputs.Rmd`, the `dry_run` chunk.

Left alone, with reasons:

| Site | Why |
| --- | --- |
| all of `case-studies.Rmd` | Every one of its sites is a validation site: the chunk fits and then compares against a closed form, `metafor`, `mgcv`, `nlme::gls`, `MASS::polr` or a hand-rolled ML `optim`, and eight end in a literal `stopifnot()`. Independent generation is the point. |
| `frmtmb.Rmd` nonlinear chunk | Prints `fixef()` against the known 2.5 and 0.7. Validation-flavored. |
| `frmtmb.Rmd` custom-family chunk | The response feeds `check_custom_family()` on a family that has no `sim` slot yet, so `frm_simulate()` would refuse it by name. |
| `frmtmb.Rmd` smooth/GP chunk | The generative mean is `sin(x)`, which is not expressible as the fitted model's linear predictor. The point of the example is that a smooth recovers it. |
| `frmtmb.Rmd` `car()` chunk | The hand-written field carries a deliberate spatial GRADIENT that the CAR prior would not produce, and `car()` has no `from_natural` map, so the call would need the internal `theta` spelling. |
| `frmtmb.Rmd` `mo()` chunk | `mo()`'s simplex (`zeta1`) has no natural-scale name. `frm_simulate()` refuses the natural spelling here, and the internal one would put raw simplex parameters in a vignette. |
| `diagnostics.Rmd` | Owned by another lane this round. Both its sites are scaffolding and are good conversion candidates. |
| `brms-migration.Rmd`, `compatibility.Rmd`, `habit.Rmd` | No hand generation to convert. |

### Tests

Converted, all behavior tests whose data is scaffolding:

`test-dry-run.R`, `test-par-template.R`, `test-naming-collisions.R`,
`test-get-prior-route.R`, `test-input-validation.R`,
`test-api-spellings.R`, `test-tabular-inputs.R`, `test-portability.R`.

**The trap, and it is worth knowing about.** The first conversion
passed the fixture's own seed to `frm_simulate()`. That restarts the
same random stream that made the covariates, so the residuals came
back exactly equal to `x` and the fixtures fitted a noiseless line.
It showed up as 43 new convergence warnings in `test-portability.R`
and two outright failures. Every converted fixture now gives the draw
its own seed, with a comment saying why.

Left alone as validation, deliberately:

- `helper-reference.R`. Its `sim_ar1_data()` and `sim_pois_glmm()`
  feed ten reference-comparison files, and the file carries the
  package's `@srrstats` G5.4/G5.5/G5.6 evidence, which claims exactly
  that the data is simulated independently from known parameters.
  Converting it would invalidate the standards text as written.
- `helper-fuzz.R`. Its section 3 simulates from the sampled model and
  its section 5 asserts against that truth, including a coverage
  invariant on the generating slope `FUZZ_BETA_X = 0.4`. Its own
  comment already names the failure mode: "Comparing simulate() with
  predict() cannot see an aterm that BOTH of them ignore."
- Every file the survey classified as validation: the reference
  agreement files against glmmTMB, lme4, mgcv, MASS, survival, nnet,
  mclust, nlme, GLMMadaptive and `mice`, plus the TMB example
  reproductions.
- `test-brms-*.R` and `helper-brms*.R`, owned by other lanes.

Not converted, but candidates once the owning lanes land:
`test-data2.R` (its `gr(cov = A)` field is drawable, since `gr_cov`
has a `from_natural` map) and `test-prior-compat.R` (a nonlinear
loss-reserving curve; drawable through the `nl` formula, but the
generative curve is the subject of the example).

## Gaps worth closing next

1. **`mo()` has no natural-scale name for `zeta1`.** This blocks
   `frm_simulate()` in the natural spelling for every monotonic-effect
   model, which is what kept the `mo()` vignette chunk and its example
   pages hand-written.
2. **Structures without a `from_natural` map** (`car`, `spde`, `ar1`,
   `cs`, `toep`, the GPs) need the internal `theta` spelling. Legal,
   but opaque wherever it appears in documentation.
3. **The roxygen example sweep**, listed in `dev/simulate-examples-todo.md`.

## Punch round, 2026-09-05

Worked from `dev/review-simulate.md`, whose 15-item punch list is the
contract. Every number below was measured in a private library
(`scratchpad/ssp-lib`: this worktree's core, then `frmtmb.ddm` from
`extensions/frmtmb.ddm` on top), R 4.6.1. Nothing was committed; the
main checkout is untouched at 5dfdd84 and clean, and this worktree is
still at df9d6bc on `wt-simulate`.

### The three blockers

**1. The seed collision, fixed in both places it appears.**
`vignettes/inputs.Rmd:161` now draws with `seed = 1001` against the
chunk's `set.seed(1)`, and carries the same three-line comment the
eight converted fixtures use. Reproduced the defect first, to be sure
the fix is the fix: recovering the generating residual exactly (the
draw spends its first six normals on the group effects and the next
thirty on the residuals) gives

| seed | residual sd | max abs cor with the seed-1 stream | at offset |
| --- | --- | --- | --- |
| 1 (as shipped) | 0.91723 | **1.00000** | 6 |
| 1001 (fixed) | 1.16733 | 0.25016 | 3 |

and a refit of the generating model moves the slope from **0.20020**
back to **0.37872** against the generating 0.4. The correlation at
offset 6 falls from 1.00000 to 0.03185.

Then swept the whole tree, by parsing each `frm_simulate()` call to its
closing parenthesis and comparing its literal `seed =` against the
NEAREST PRECEDING `set.seed()` in the same file: 367 files, 34 calls
with a literal seed, **1 further collision**, at
`vignettes/frmtmb.Rmd:279` (pre-existing, not this lane's, and the
review says so). Fixed the same way. The sweep now reports **0
collisions**. The review's larger count of 44 counts calls whose seed
is an expression (`seed + 1000L`); those cannot collide by
construction.

**2. The wiener boundary masses were transcribed from the wrong test,
and the wrong pair IS reproducible.** `dev/simulate-findings.md:202`
now reads 0.196881 lower and 0.803119 upper, which is what the tier's
own `ddm_cells()` path gives at
`mu = 1.1, bs = 1.5, ndt = 0.2, bias = 0.45` (sum 1.000000000). The
review could not place the old pair and said it was unreachable by
varying any single parameter. It is reachable by varying all four:
**0.157095 / 0.842905 are the masses at
`mu = 1.2, bs = 1.4, ndt = 0.25, bias = 0.5`**, which is the
simulate-on-a-fit test at line 141 of the same tier file. The note now
says so. For the record, no single-parameter move reaches the old
lower mass at the stated point except by coincidence: it would need
`mu = 1.293163`, `bs = 1.763404` or `bias = 0.506162`.

**3. The `R/families.R` line numbers, regenerated from the tree.**
Rather than trust the review's five, I audited ALL 65 table rows
against disk by asking whether each named line is an `#' @examples`
tag: **60 of 65** held, and the five that did not were exactly the
`R/families.R` rows the review names. Its replacements are right: I
regenerated them independently and got the same. The shift at that
point was **117 lines**, not the review's 123: the diff is +120/-3,
and 2714 - 2597 = 117.

The renumber then had to be done TWICE, which is worth recording
because it is the same trap the review caught the lane in. Items 10,
11 and 12 below add lines to `R/families.R` and change the length of
`R/fit.R`, `R/priors.R` and `R/conditional-effects.R`, so the numbers
written earlier in this round went stale inside the round. The audit
caught nine rows the second time, five in `R/families.R` and four
elsewhere. The final, verified numbers are:

| Block | pristine | after the lane | after this round |
| --- | --- | --- | --- |
| `R/families.R` `mixture` | 2597 | 2714 | **2724** |
| `R/families.R` `mixture_probs` | 2956 | 3073 | **3083** |
| `R/families.R` `mixture_mvn` | 3190 | 3307 | **3317** |
| `R/families.R` `cox_baseline` | 3910 | 4027 | **4057** |
| `R/families.R` `frmtmb-families` | 4228 | 4345 | **4375** |
| `R/conditional-effects.R` `plot.frmtmb_fit` | 1561 | 1561 | **1570** |
| `R/conditional-effects.R` `pp_check` | 1619 | 1619 | **1628** |
| `R/priors.R` `frmtmb-priors` | 1360 | 1360 | **1366** |
| `R/fit.R` `frmtmb_control` | 1360 | 1360 | **1354** |

Each was resolved by NAME, by mapping every `@examples` tag to the
object it documents, not by taking the nearest tag. `R/families.R:113`
is correct throughout, because it sits above every insertion.
Re-audited last, after all other edits: **65 of 65**. The independent
sweep of `R/` also still reads 79 `@examples` tags, 66 with a hand RNG
call, 1 of those already using `frm_simulate()`, which is the review's
count exactly.

### Should fix

**4. `R/parse.R:491` moved to the convert list**
(`dev/simulate-examples-todo.md:212`). Drew the response myself in the
natural spelling to check the review's claim:
`sd_mmschool1school2__Intercept` on `(1 | mm(school1, school2))`,
sd(y) = 0.9294, and a refit recovers Intercept 0.915 and slope 0.509
against the generating 1 and 0.5. The row now records the true reason
and one thing the review did not mention: the page's WEIGHTED term has
its own natural name,
`sd_mmschool1school2weightscbindshare1share2__Intercept`, so the draw
must use the unweighted spelling. That costs nothing, because the page
draws `y` once and fits it twice.

**5. The "nine" was stale, not a registry misreading**
(`dev/simulate-findings.md:30`). Confirmed both halves by measurement
rather than by reasoning. On this tree the registry reads 40 entries,
32 exported zero-argument constructors, and 2 without a `sim`
(`cox`, `categorical`); the lane's diff adds `sim` to exactly three
families, so pristine df9d6bc reads five by every slicing. And the
nine has a source: at the 2026-08-31 commit b492cf0
("mixtures, mi, cs, rr") `R/families.R` carries **35 registry entries
and 26 `sim` slots**, which is nine short. The paragraph now says the
number was carried forward from that round.

**6. `CONTRIBUTING.md:95`** gains the rule that makes the
`frm_simulate()` policy safe: give the draw a seed distinct from the
covariates', because `seed =` restarts the stream, with the repository
convention (covariate seed plus 1000) named.

**7. The ungrammatical tweedie clause**, at `R/families.R:1158` and
`dev/simulate-findings.md:76`. Both now read "the `power12` link
guarantees `1 < p < 2`, and on that range...".

**8. Spaced hyphens used as em-dash stand-ins**, all 12, restructured
rather than substituted: `tests/testthat/test-simulate-density.R` 8,
9, 32, 33, 78, 79, 601; the ddm tier's 19, 89, 193; and the tweedie
and hurdle_poisson `sim` comments in `R/families.R`. Verified by
re-sweeping every line this lane ADDED: the five remaining ` - ` hits
in the diff are all arithmetic (`lw - max(lw)` and friends), and the
two new test files have none left in prose.

**9. `NEWS.md:25`** no longer names `posterior_predict()` as core; it
now reads "`simulate()` and `frm_simulate()` work for them, and so
does `posterior_predict()` through frmtmb.sample". While there, the
registry sentence at line 38 gained the one qualifier section 7 of the
review asked for, naming the deferred entries.

### Follow-ups, applied rather than deferred

**10. `R/families.R:1430`, the hurdle_poisson floor.** The positive
part is clamped with `pmax(..., 1L)`, which is exact because its
support starts at 1, and the comment now attributes the loss to
`p0 + u (1 - p0)` rather than to "a mu near the double epsilon".
Measured, at `hu = 0.3` over 2e5 draws:

| mu | P(Y = 0) before | P(Y = 0) after |
| --- | --- | --- |
| 1e-8 | 0.2973 | 0.2995 |
| 1e-14 | 0.4293 | 0.2995 |
| 1e-16 | 1.0000 | 0.2995 |
| 1e-300 | 1.0000 | 0.2995 |

and the positive part's minimum is 1 everywhere. The clamp moves
nothing in normal use: at `mu` of 0.4, 2.5 and 8 the draws are
BIT-IDENTICAL to the unclamped expression, same seed, zero differing
rows. That is why the tier's fixed-seed assertions do not move.

**11. `R/families.R:1264`, the compois cache comment.** Qualified
rather than re-engineered, which is the cheaper half of the review's
"either/or". Measured here: at n = 1000 with a continuous predictor
`frm_simulate()` takes **1.02 s** (1.02 ms per row) against
effectively zero for a two-level factor design and for `poisson()` on
the same design. The comment now says the cache pays on a factor
design, never hits on a continuous one, and names the grid-and-
interpolate remedy.

**12. The pre-finalize reads, all four sites through one helper.**
Added `carry_finalized_responses()` at `R/families.R:3797` and used it
at `R/fit.R:504`, `R/simulate-new.R:488`, and the two the review
flagged: `R/par-template.R:155` and `R/priors.R:719`. Verified
behaviour-neutral today on `gddm()`, the only in-repo family whose
finalize matters: `get_prior()` returns 4 rows with identical
class/coef/dpar/group keys on the formula and fit routes, and
`par_template()` gives `beta 1, betad 3` with and without `prior =`.

**13. `R/conditional-effects.R:1192`, and a correction to the
review.** The review calls this latent because no registry family
simulates only through `sim_ctx`. It is NOT latent:
**`mixture_mvn(2, 2)` has `sim_ctx` and no `sim`**, so
`conditional_effects(method = "predict")` refuses a model
`simulate()` accepts today. It escapes the registry sweep the way
`multinomial` does, by needing arguments.

A blind swap to `sim_can(fam)` would have made things worse, not
better: `sim_response()` reads `fam[["sim"]]` directly, and the band
code below takes per-row quantiles of a `replicate()`, so a
whole-draw family returning an n by D matrix would have failed
downstream instead of being refused. So the gate is now two refusals:
`!sim_can(fam)` raises the same message every other entry point
raises, now carrying the family's own `sim_note()` (which is how
`cox()` gets to explain itself here), and a second check refuses a
family that draws whole, saying that rather than claiming it has no
simulator.

**14. `dev/simulate-examples-todo.md:67`** now says 65 pages listed, 64
needing work, and names the 65th.

**15. `tests/testthat/test-simulate-density.R:54`** records the count
the alpha argument rests on. Counted it rather than copying it, by
instrumenting `testthat::expect_lt()` while the tier runs: **251
statistical assertions, 87 chi-square, 77 mean, 77 variance, 10
P(Y = 0)**, which is the review's split exactly. The comment gives the
0.0003 per-reseed false-alarm probability and tells the next reader to
re-check it if the tier grows.

`SIM_ALPHA = 1e-6` is unchanged, as instructed.

### One thing left alone deliberately

The review notes that `test-input-validation.R:14` and
`test-portability.R:8` say the response "only has to come from the
model that is fitted" when the draw is really from a submodel
(`bf(y ~ x)`, fitted as `y ~ x + (1 | g)`). It calls this "worth a
word but not a fix" and does not put it on the punch list, so the
comments stand.

### Verification

One process per test file, counts audited by name against the files on
disk (no duplicates, none missing), against the lane's install in
`ssp-lib`.

```
test-*.R files on disk         : 105
files reported                 : 105   (reported set == on-disk set)
TOTAL pass 5858   fail 0   error 0   skip 21
summed wall clock              : 12.9 min
```

| File | pass | fail | secs |
| --- | --- | --- | --- |
| test-simulate-density.R | 435 | 0 | 4.4 |
| test-dry-run.R | 23 | 0 | 1.8 |
| test-par-template.R | 51 | 0 | 2.8 |
| test-naming-collisions.R | 31 | 0 | 2.7 |
| test-get-prior-route.R | 36 | 0 | 2.0 |
| test-input-validation.R | 43 | 0 | 1.8 |
| test-api-spellings.R | 22 | 0 | 2.3 |
| test-tabular-inputs.R | 9 | 0 | 5.8 |
| test-portability.R | 91 | 0 | 9.0 |
| test-message-uniqueness.R | 6 | 0 | 3.3 |
| test-bracket-access.R | 7 | 0 | 1.4 |

**5858, 0 failures, 105 of 105.** Every count the review published
reproduces after these changes, which is the point: nothing here moved
a number.

The ddm suite, one process, against the `ssp-lib` core:

```
core from: scratchpad/ssp-lib
ddm  from: scratchpad/ssp-lib
DDM SUITE  files: 15  pass: 925  fail: 0  error: 0  skip: 0  secs: 858.7
```

with the extension tier at **41 pass, 0 fail** inside it.

Both converted vignettes render, with pandoc on PATH rather than by
`knit()` alone: `frmtmb.Rmd` 5.4 s, `inputs.Rmd` 1.3 s, **0 errors, 0
warnings**.

`R CMD build` then `R CMD check --as-cran --no-manual`, with
`_R_CHECK_CRAN_INCOMING_=false` and `NOT_CRAN=true` so both agreement
tiers ran inside the check:

```
* checking examples ... [40s] OK
* checking examples with --run-donttest ... [39s] OK
* checking tests ... [581s] OK
* checking re-building of vignette outputs ... [156s] OK
* DONE
Status: OK
```

**0 NOTEs, 0 WARNINGs, 0 ERRORs**; `grep -cE "NOTE|WARNING|ERROR"` over
`00check.log` returns 0.

`roxygen2::roxygenise()` was run after the roxygen edit and is
idempotent: a second run leaves `git status` byte-identical, and `man/`
is untouched (the new helper is `@noRd`).

Two measurement notes, so the numbers above can be trusted later. A
first pass reported spurious errors because test files were run
without the package as their context: tests reach package internals
and use `local_mocked_bindings()`, both of which need
`env = new.env(parent = asNamespace(pkg))` AND `package = pkg` on
`test_file()`, which is what `test_check()` supplies under `R CMD
check`. And `test-perf.R` failed once with 1 of 3 when the suite ran
concurrently with the as-cran check; it passes 3 of 3 unloaded and
passed inside as-cran's own test stage. Its assertion is wall clock,
and the file says so.

### Scope

`git diff --name-only df9d6bc` is now **19 tracked files**, the lane's
15 plus the four that items 12 and 13 touch (`R/fit.R`,
`R/par-template.R`, `R/priors.R`, `R/conditional-effects.R`), plus the
same **5 untracked** files. Nothing was committed. This worktree is
still at df9d6bc on `wt-simulate`; the main checkout is at 5dfdd84 and
clean, and I neither wrote to it nor created a worktree from it.
