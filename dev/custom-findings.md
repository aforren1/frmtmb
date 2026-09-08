# Custom-family seams: `se()` and `cens()`

Lane: custom. Worktree `frmtmb-wt-custom`, branch `wt-custom`, off
`c18253e` (frmtmb 0.53.0, main's head merged with `wt-api`).

Two gates made a family written in plain R a second-class citizen.
`dev/bcm-findings.md` records both with reproductions, as seams 2 and 3
of four. This file records what they were, what replaced them, and what
stays refused.

## Baseline, before any edit

The gated Stan identity tier, run one file per process against the
private library. `options(frmtmb.bcm_report = TRUE)` prints the
residual of each identity: `measured` is
`log_prob(theta_hat) - frmtmb_log_density(theta_hat)` and the claim is
that it equals the stated constant, which is 0 for every model below.

| file | model | ours | measured - stated | max grad |
|---|---|---|---|---|
| test-bcm-esp.R | OptionalStopping | -30.37952411 | 3.552713679e-15 | 2.24e-07 |
| test-bcm-esp.R | Extraversion | -1.747339859 | -2.853717262e-12 | 1.62e-10 |
| test-bcm-data-analysis.R | ChaSaSoon | -118.1411391 | -3.694822226e-13 | 9.07e-08 |

```
RESULT test-bcm-esp.R           OK pass=17 fail=0 skip=0 warn=0  5.9s
RESULT test-bcm-binomial.R      OK pass=37 fail=0 skip=0 warn=0  8.2s
RESULT test-bcm-data-analysis.R OK pass=45 fail=0 skip=0 warn=0  9.8s
```

`test-bcm-binomial.R` holds six identities and none of them touches
either seam; it is in the table of runs because the brief names it, and
its residuals are listed in full further down.

## What the discrete `cens()` gate was protecting against

Not nothing. The censoring algebra in `row_lpdf()` (R/objective.R) is
written for a continuous response, where the endpoints of an interval
carry no probability:

- right censoring at `k` scores `1 - F(k)`, which is `P(Y > k)`;
- an interval `[k, k2]` scores `F(k2) - F(k)`, which is `P(k < Y <= k2)`.

For a count that silently drops the point mass at the recorded value. A
family that carries a CDF but declares `type = "continuous"` slips past
the gate today and reaches exactly that algebra, so the loss is
measurable on the CURRENT build rather than argued. 200 poisson draws
at lambda 4, right censored at 6, 30 rows censored
(`scratchpad/cf-probe-gate.R`):

```
frmtmb logLik        -384.2856427254
hand, Y >  k         -384.2856427254     <- what the core computes
hand, Y >= k         -363.4564099828     <- what "6 or more" means
mu hat (Y > k)         3.8590634324
mu hat (Y >= k)        3.7227288111
```

Twenty log units and a 3.7 percent shift in the estimate, with no
warning. Deleting the gate and keeping the algebra would have shipped
that. So the gate was a placeholder for a decision nobody had made, and
the decision is below.

## The convention, decided

**A censoring bound on a discrete response NAMES a value the response
can take, and is included in the event.**

| code | written | means | scored as |
|---|---|---|---|
| 0 | `"none"` | `Y == y` | `f(y)` |
| -1 | `"left"` | `Y <= y` | `F(y)` |
| 1 | `"right"` | `Y >= y` | `1 - F(y - 1)` |
| 2 | `"interval"` | `y <= Y <= y2` | `F(y2) - F(y - 1)` |

Equivalently, and this is how it is implemented: every LOWER edge
enters the CDF as `F(edge - 1)`, and upper edges are unchanged because
`F` already includes its argument. One line of arithmetic in
`row_lpdf()` and the rest of the block reads as it did.

Three reasons, in the order they decided it.

1. It is the rule the package already follows. `trunc(lb = k)` on a
   discrete response has always meant `Y >= k` and has always been
   scored with `F(k - 1)` (R/objective.R, citing brms#1903). Had right
   censoring at `k` meant `Y > k`, the same number would have meant two
   different things in two addition terms on the same response.
2. It is what a recorded count means. A survey band "5 or more" is
   `Y >= 5`. The exclusive reading throws away the point mass at the
   recorded value, which is the 20 log units measured above.
3. The port needs it. ChaSaSoon's band is `F(hi) - F(lo - 1)` in the
   book's own Stan program, which is
   `binomial_cdf(25 | n, theta) - binomial_cdf(14 | n, theta)`:
   inclusive at both ends.

It DIFFERS from brms for right and interval censoring of a discrete
response, where brms emits `poisson_lccdf(y | mu)` = `P(Y > y)`. That
is the same off-by-one brms fixed for truncation and not for censoring,
and it is stated in the migration vignette rather than left to be
discovered.

A one-point discrete interval (`y2 == y`) is now legal, because under
this convention it is exactly `P(Y = y)`. It stays refused on a
continuous response, where it is an event of probability zero.

## The convention, measured

240 poisson draws at lambda 5, one data set carrying all four codes
(`scratchpad/cf-probe-cens.R`), against a hand-written likelihood:

```
rows  observed 60  right 60  left 60  interval 60
frmtmb   mu 5.282906261974  logLik -272.538062061560
hand     mu 5.282906279932  logLik -272.538062061560
|dlogLik| 5.684e-14   |dmu| 1.796e-08
exclusive convention at the same mu: -317.470925675037 (differs by 44.932864)

right     alone: frmtmb -154.729571329804  hand -154.729571329805  diff 1.137e-13
left      alone: frmtmb -166.024584966086  hand -166.024584966086  diff 5.684e-14
interval  alone: frmtmb -196.220258306889  hand -196.220258306889  diff 1.705e-13
degenerate interval: -516.174510202242 vs plain -516.174510202242  diff 1.137e-13
cens+trunc: frmtmb -148.770872057292  hand -148.770872057293  diff 2.274e-13
```

## `se()`, decided

The name test becomes a DECLARATION test, in the vocabulary that
already exists: a family gets `se()` if it names `"se"` in
`frmtmb_family(accepts_aterms =)` or in `required_aterms`.

`se()` is the one core addition term whose entire effect is inside the
density. The core does not multiply it in (`weights()`), reshape
anything with it (`cens()`, `trunc()`) or create parameters from it
(`mi()`); it hands `aterms[["se"]]` to the family and maps out the now
redundant `sigma`. So "does this family read it" is the only question,
and only the family can answer.

An UNDECLARED family (`accepts_aterms = NULL`) is refused, even though
NULL otherwise means "every term". For an allow-list, "did not say"
safely means "accepts"; for a capability whose whole effect is a read,
"did not say" has to mean no, or a family that ignores the term gets a
silently wrong likelihood. It is also exactly the behavior every custom
family had before, so nothing that fitted before fits differently now.

The built-in gaussian and student already declared `"se"`; no other
built-in family does. The name test and the declaration test therefore
agree on every built-in family but one, and that one is the mixture.

## What `mixture()` gave back, and why it is refused

A mixture's allow-list is the UNION of its components', so
`mixture(gaussian, gaussian)` declared `"se"` and the declaration test
admitted it where the name test had refused it. Measured on 80 draws
from a two-component normal, `v | se(s) ~ 1`:

```
NO ERROR  logLik -228.2941603162
sigma1 (Intercept)  0.9794869      <- the shared starting value
sigma2 (Intercept)  0.9794869      <- unmoved
every standard error: NaN
```

Reading the term is only half of what `se()` means. The other half is
that the residual scale it replaces is mapped out, and that step names
the dpar `sigma`, which `sigma1` and `sigma2` are not. So the union
subtracts `"se"` (R/families.R), and `mixture() + se()` is refused by
name with the message that says which families read it.

## The port, simplified

Both workarounds collapse.

| was | is |
|---|---|
| `bf(xs \| vreal(sx) ~ ...)` + a family reading `aterms[["vreal1"]]` | `bf(xs \| se(sx) ~ ...)` + a family declaring `required_aterms = "se"` |
| `bcm_binomial_band()`, band written into the density, `vint(hi)` | `bcm_binomial_cdf()`, four lines of `lcdf`, band written as `cens(cc, hi)` |

`bcm_binomial_band()` is renamed `bcm_binomial_cdf()` because that is
what is left of it: the core binomial plus the CDF the first gate asks
for. Its `lpdf` is now `RTMB::dbinom` rather than a hand-rolled CDF
difference, and 30 lines of density and validation went with the band.

### Identity residuals, before and after

Same models, same estimates, same claim (`const` 0 everywhere):

| model | ours before | ours after | residual before | residual after | max grad |
|---|---|---|---|---|---|
| Extraversion | -1.747339859 | -1.747339859 | -2.853717262e-12 | -2.853717262e-12 | 1.62e-10 |
| ChaSaSoon | -118.1411391 | -118.1411391 | -3.694822226e-13 | -1.421085472e-14 | 9.07e-08 |
| OptionalStopping | -30.37952411 | -30.37952411 | 3.552713679e-15 | 3.552713679e-15 | 2.24e-07 |

Extraversion's residual is bit-identical: `se()` and `vreal()` carry the
same column into the same `dnorm`. ChaSaSoon's improves by an order of
magnitude, because the CDF difference is now formed once by the core in
`row_lpdf()` rather than inside the family, and the log of a difference
rounds differently from a difference of logs. The estimate and the
gradient are unchanged.

`test-bcm-binomial.R` touches neither seam and its six residuals are
bit-identical before and after: -6.217248938e-15, -1.243449788e-14,
-1.33226763e-14, -4.440892099e-16, 0, -7.283063042e-14.

## What stays refused, and why

Not everything opened. Each refusal below names itself.

1. **A family that does not declare `"se"`.** The message says how to
   opt in and which built-in families do. This is the old behavior of
   every custom family, kept: a term whose whole effect is a read
   cannot default to yes.
2. **`mixture()` with `se()`.** Above: half a capability is worse than
   none.
3. **`cens()` without a CDF.** Unchanged, and the gate that was always
   right. The core binomial is still refused for exactly this, and
   `bcm_binomial_cdf()` still exists to pass it.
4. **A non-integer censoring bound on a discrete family.** The
   inclusive convention reads a lower edge as `F(k - 1)`, which assumes
   the unit integer lattice. A family whose support is not that would
   be shifted onto a point it has no mass at, silently. Refused by
   name, with the convention in the message.
5. **`residuals(type = "osa")` on a censored discrete fit.** The
   one-step window is `[lo, hi]`, which is right where the endpoints
   carry no mass. Under inclusive bounds an uncensored count's support
   is `[lo + 1, hi - 1]` and its PIT renormalizes on `F(hi - 1) -
   F(lo)`. Shifting the window is one line; verifying it against a
   reference is not, and an unverified PIT that looks calibrated is
   the failure mode this whole lane is about. Every other residual type
   works and `dharma_residuals()` covers the same ground.
6. **An `lccdf` for `poisson()`.** Right censoring would use it, and
   `RTMB::ppois(lower.tail = FALSE, log.p = TRUE)` is exact in R
   (-2773.28 at q = 700, where `log(1 - F)` is -Inf) and does not TAPE:
   inside `MakeADFun` it reaches `stats::ppois` and errors with
   "Non-numeric argument to mathematical function". So a censored
   poisson still scores a right-censored row as `log(1 - F)` and meets
   the representable-tail limit the gaussian escaped. Writing an exact
   discrete log survivor is the follow-up.

## Boundary cases, measured

`scratchpad/cf-probe-zero.R`, and both are now tests:

- Right censoring at 0 evaluates the CDF at -1. `P(Y >= 0) = 1`, so
  the row contributes exactly nothing, and the fit equals the
  uncensored-rows-only likelihood to 1.7e-13.
- An interval `[0, k]` is `F(k) - F(-1) = F(k)`, with no `NaN` from
  below the support.

## Touched files

Core, five files:

| file | what |
|---|---|
| `R/frame.R` | the `se()` capability test (~1435); the discrete `cens()` lattice check replacing the blanket refusal (~1371); the discrete interval-bound rule (~1415) |
| `R/objective.R` | `row_lpdf()`: the lower-edge shift for a discrete family (~99), and the two CDF calls that read it |
| `R/families.R` | `family_declares_aterm()`; `accepts_aterms` on `poisson()`; `"se"` removed from the `mixture()` union; the `lcdf`, `accepts_aterms` and lccdf docs; the new "Censoring a discrete response" section |
| `R/predict.R` | `osa_cens_domain(discrete =)`, the refusal and its reason |
| `R/compat.R` | four rows: `cens()`x`poisson`, `cens()`x`group:discrete`, `se()`x`kind:family`, `mixture`x`se()`, plus the `residuals_osa` note and the `cdf_continuous` group comment |

Port and docs: `inst/bcm/binomial-extras.R` (both families),
`tests/testthat/test-bcm-esp.R`, `tests/testthat/test-bcm-data-analysis.R`,
`vignettes/bayesian-cognitive-modeling.Rmd`,
`vignettes/brms-migration.Rmd`, `vignettes/compatibility.Rmd`,
`NEWS.md`, `man/frmtmb_family.Rd`, `dev/bcm-findings.md`.

Tests: `test-cens-trunc.R` (+8 tests), `test-custom-family.R` (+5),
and three corrected in place: `test-review-v25.R` (the discrete
composed reference now shifts the right-censored lower edge),
`test-v14.R` (asserts the declaration message, not the family names),
`test-cens-lccdf.R` (the comment for why poisson has no `lccdf` is now
the taping reason, not the refusal).

## Verification

### The core suite, one file per process, audited by name

131 test files on disk, 131 output files, 131 result lines, 131
distinct names: no file missing and none run twice. Run four shards
wide, each file in its own process writing its own output file, with
the gated tiers ON (`FRMTMB_BRMS_FIT_TESTS=true`, `NOT_CRAN=true`,
`FRMTMB_STAN_CACHE` on this lane's own copy of the cache).

```
pass=9103 fail=0 skip=1 warn=2
```

- The one SKIP is `test-fuzz.R`, which is gated on `FRMTMB_FUZZ=true`.
- The two WARNs are pre-existing and unrelated:
  `test-prior-compat.R` "coef and group narrow the classes that read
  them", and `test-brms-likelihood.R` row 22 "the brms link roster on
  positive-mean responses".

One file failed on the first pass and is the only test this lane broke
by changing behavior rather than by breaking code: `test-compat.R`
asserted `frm_compat("cens()", "poisson")$status == "refused"` in two
places. That is now `"conditional"`, which is the change. Both
assertions were updated and the file re-run: 266 pass, 0 fail.

### The gated BCM chapter files

All pass, and every identity residual is in the table above.

```
test-bcm-esp.R           pass=17 fail=0
test-bcm-binomial.R      pass=37 fail=0
test-bcm-data-analysis.R pass=46 fail=0   (45 before; +1 hand-likelihood test)
```

### R CMD check --as-cran, with the manual

Pandoc from RStudio's quarto tools, TinyTeX prepended to PATH, gated
tiers OFF so the list is the one a submission would get.

```
Status: 3 NOTEs        (0 ERRORs, 0 WARNINGs)
  tests            'testthat.R' [417s] OK
  vignettes        re-building [338s] OK
  PDF manual       [29s] OK
```

Each NOTE, with its cause:

1. **CRAN incoming feasibility [283s]: "New submission".** The package
   is not on CRAN. Not actionable and not this lane's.
2. **Examples with elapsed time > 5s.** `profile.frmtmb_fit` 7.94s and
   `residuals.frmtmb_fit` 5.25s. Both pre-existing; neither example is
   touched here. (Run 1 measured the same two at 5.23s and 6.56s, so
   they sit near the threshold and swap places with machine load.)
3. **HTML manual: "Skipping checking math rendering: package 'V8'
   unavailable".** V8 is not installed, which is the expected NOTE.

The first run of the check, on a tarball built before the
`test-compat.R` fix, was `1 ERROR, 3 NOTEs`; the ERROR was those two
assertions and nothing else. Same three NOTEs both times.

Roxygen is idempotent: the second `roxygenise()` writes nothing, and
`man/frmtmb_family.Rd` is the only generated file that changed.

---

# Punch round, 2026-09-08

Review at `dev/reviews/2026-09-08-custom.md`, verdict PUNCH, three
blocking items, none touching the likelihood. The correctness of the
convention, the capability gate, the OSA refusal and the identity
residuals were all reproduced independently and held.

## The brms divergence: upheld, with the argument I should have made

The reviewer checked the premises rather than accepting them, and found
the fact that settles it. **brms is internally inconsistent here.** Its
discrete truncation emits `poisson_lccdf(lb - 1 | mu)`, an INCLUSIVE
lower bound; its censoring emits `poisson_lccdf(y | mu)`, exclusive. So
in brms `trunc(lb = 6)` means `Y >= 6` while a right-censored row
recorded at 6 means `Y > 6`: one number, two meanings, one response.

frmtmb's `trunc()` reproduces brms bit for bit (poisson on
`y = 3,4,5,6` at `b0 = log 4`: `trunc(lb = 2)` gives 6.9990718955 under
both, against 6.2954805379 exclusive). So frmtmb's `cens()` could match
brms's censoring or match frmtmb's own truncation, and not both. It
matched its own. That turns "we differ from brms" into a reason.

Two more facts from the review, both now in the docs:

- **Left censoring already agrees with brms** (`P(Y <= y)` in both), so
  the divergence is RIGHT and INTERVAL only. My write-up said "right
  and interval" but did not say left agrees; every statement of it now
  does.
- **The inclusive reading is the only one consistent with the censored
  simulator**, which predates this lane. `simulate(censored = TRUE)`
  caps a draw with `pmin(pmax(y, lo), hi)` (`R/predict.R:2727`), which
  records `k` exactly when the latent draw is `>= k`. Measured on 4000
  draws from a fit censored at 7: simulated mass at the point 0.13250,
  against `P(Y >= 7) = 0.12190` inclusive and `P(Y > 7) = 0.05763`
  exclusive. brms's reading would silently decouple the likelihood from
  the simulator `dharma_residuals()` rests on.

Written into `?frmtmb_family` ("Censoring a discrete response"),
`vignettes/brms-migration.Rmd`, `NEWS.md`, and the `frm_compat()` note
for `trunc()` x `poisson`, which now states that the truncation bound
is inclusive because that is the premise the whole argument rests on.

## Item 1 - the refusal pointed at a function nobody can call

`R/frame.R:1451`. The message told family authors to see `resid_sd()`,
which is `@noRd`, absent from NAMESPACE and in no man page; the
reviewer hit "could not find function" while writing a conforming
family. It also named only `accepts_aterms`, though
`required_aterms = "se"` opens the gate too.

The message now names both declarations and states the convention
inline instead of by reference: read `aterms[["se"]]` as the standard
deviation, and honor `se(x, sigma = TRUE)` by reading
`aterms[["se_sigma"]]` and using
`sqrt(dpars[["sigma"]]^2 + aterms[["se"]]^2)`. A test asserts the
message mentions both declarations, states the `se_sigma` convention,
and does NOT contain the string `resid_sd`.

## Item 2 - the half-capability my own change opened

`R/parse.R:1582`. The other half of `se()` - mapping out the scale it
replaces - was keyed on the dpar being literally named `sigma`. Opening
`se()` to custom families opened a second route to the exact pathology
`mixture()` is refused for. Reproduced before the fix, dpars
`c("mu", "tau")`:

```
FITTED (the pathology)
$tau (Intercept)  -0.1372695     <- its starting value, never moved
standard errors all NaN: TRUE
warning: The likelihood is FLAT in 1 direction: tau_(Intercept)
```

Neither of the review's two suggestions works as stated. Keying the
map-out on "the family's scale dpar" needs a slot that does not exist,
and refusing every declaring family with no dpar named `sigma` would
refuse `bcm_gaussian_probit()`, whose `dpars` is `"mu"` alone and which
is the model that motivated the whole lane. The rule that separates
them is what the core can actually map out:

> `se()` without `sigma = TRUE` replaces the residual scale. The core
> maps out the dpar the convention names, `sigma`. A declaring family
> with no `sigma` and no other free dpar has nothing to map out and is
> fine. A declaring family with no `sigma` that carries another free
> dpar is REFUSED, because that dpar would be left free and unread.

A family carrying `sigma` is trusted about its other dpars, which is
why `student()` keeps `nu` estimated beside a known `se()` exactly as
it always has. Three ways out, all named in the message and all
measured in the tests:

| way out | result |
|---|---|
| name the scale `sigma` | fits, all SEs finite |
| pin it in the formula, `bf(y \| se(s) ~ x, tau = 1)` | fits, same logLik to 1e-8 |
| `se(x, sigma = TRUE)` | fits; `tau` 0.905871 against a true 0.9, all SEs finite, and bit-identical to the same family whose scale IS named `sigma` |

The third needed data carrying variance beyond `se` to measure at all:
on data whose whole spread is the known `se`, the extra scale sits on
its boundary at zero and is flat for a real reason rather than a
structural one. The test says so.

## Item 3 - the refusal was missing from `?residuals`

`R/predict.R:2291`, rendered to `man/residuals.frmtmb_fit.Rd`. The page
enumerates which `osa` x `cens()` combinations are refused and did not
list the new discrete one. Added with its reason.

## Also done, from the non-blocking list

- `R/families.R:307` - the poisson upper-tail figure `-2773.28` now
  states its lambda (5).
- The `frm_compat()` note for `trunc()` x `poisson` now says the lower
  bound is INCLUSIVE, `P(Y >= k)` scored with `F(k - 1)`.

## Left for the owner

- **`docs/`** still says `cens()` "additionally refuses discrete
  responses" (`docs/articles/compatibility.md:83` and the rendered
  HTML). This project rebuilds the site in its own commit; that is the
  release step, not this lane's.
- **A one-time runtime message** on the first discrete censored fit,
  which the review floats as the only thing that would close the silent
  porting hole documentation cannot. It is a real gap and a real design
  decision: it needs a policy for once-per-session notices that the
  package does not have, and it would speak to every user of discrete
  censoring to warn the subset porting from brms. Recorded rather than
  guessed at.

## Verification after the punch round

Suite re-run, four shards, one file per process, audited by name: 131
files on disk, 131 output files, 131 result lines, 131 distinct names.
None missing, none twice.

```
pass=9124 fail=0 skip=1 warn=2      (9103 before; +21 from the new tests)
```

Same one skip (`test-fuzz.R`, gated on `FRMTMB_FUZZ`) and same two
pre-existing warns (`test-prior-compat.R`, `test-brms-likelihood.R`).

Every BCM identity residual is bit-identical to the pre-punch run:
Extraversion -2.853717262e-12, ChaSaSoon -1.421085472e-14,
OptionalStopping 3.552713679e-15, and all six of
`test-bcm-binomial.R`. None of the three items touches the likelihood,
and the numbers say so.

`R CMD check --as-cran`, manual and vignettes built:

```
Status: 3 NOTEs        (0 ERRORs, 0 WARNINGs)
  tests            'testthat.R' [340s] OK
  vignettes        re-building [330s] OK
  PDF manual       [10s] OK
```

1. CRAN incoming feasibility: "New submission". The package is not on
   CRAN.
2. Examples over 5s: `residuals.frmtmb_fit` 7.08s and
   `dharma_residuals` 6.46s. Pre-existing and load-dependent - across
   the three checks this round the pair over the threshold has been
   `profile`/`residuals` twice and `residuals`/`dharma_residuals` once,
   all of them sitting near 5s. No example touched by this lane.
3. HTML manual: "Skipping checking math rendering: package 'V8'
   unavailable" - the expected NOTE.

Roxygen idempotent: the second `roxygenise()` writes nothing.
`man/frmtmb_family.Rd` and `man/residuals.frmtmb_fit.Rd` are the two
generated files that changed.
