# Frequency-domain lane: what was built and what was measured

Date: 2026-09-07. Branch `wt-freq` at `c18253e` (main plus the merged
API lane), frmtmb 0.53.0, R 4.6.1, Windows 11. Every number here comes
from a private library; nothing was installed into the shared one.

Specification: `dev/frequency-domain-todo.md`, items 1 to 5.

**Amended 2026-09-08 by the `spectral2` lane.** The smoothness refusal
now compares ordinates THREE apart rather than neighbors. Every number
here that describes that check was re-measured; where a figure changed,
the new one is in place and the old one is named beside it. Two claims
in this file were wrong and are corrected where they stand: the Hann
false-refusal rate, which is 0 rather than 1 to 2 percent, and the cost
of attaching a `frmtmb_structure()`, which is nothing for a structure
carrying no `loglik`. The decisions are in `dev/spectral2-findings.md`.

## What landed

| item | where |
| --- | --- |
| `frm_periodogram()` | `R/spectral.R:157` |
| its per-series core | `R/spectral.R:52` |
| taper weights | `R/spectral.R:16` |
| `whittle()` | `R/spectral.R:360` |
| its two refusals | `R/spectral.R:258` |
| `frm_series_draw()` | `R/spectral.R:460` |
| family-fixed dpar seam | `R/parse.R:1622` |
| tests (94 checks) | `tests/testthat/test-spectral.R` |
| vignette | `vignettes/spectral.Rmd` |

## The identity, re-verified through the shipped code

`exponential()` is mean-parameterized: the intercept of an
intercept-only fit is `log(mean(y))` to the last digit. So a Whittle fit
is an ordinary `frm()` call.

AR(1), n = 1024, phi = 0.7, `frm(bf(pgram ~ ls - log(1 - 2 ph cos w +`
`ph^2), ...), family = whittle())`:

| quantity | frmtmb | `arima(method = "ML")` |
| --- | --- | --- |
| phi | 0.708617 | 0.707178 |
| se(phi) | 0.02218 | 0.02204 |
| innovation variance | 0.98967 | 0.98799 |

## Item 1: the periodogram constructor

`frm_periodogram()` was checked against R's own estimate rather than
against arithmetic written twice. With `taper = 0, detrend = FALSE,`
`fast = FALSE`, `stats::spec.pgram()` and `frm_periodogram()` agree on
every retained ordinate to a ratio of exactly 1, and on every frequency
to 0. That single check pins the scaling AND the frequency axis:
ordinate `|sum_t h_t x_t e^{-2i pi f t / fs}|^2 / (fs sum_t h_t^2)`,
frequency `j fs / m` in Hz.

Measured properties, each with a test that fails without it:

- **Zero and Nyquist are dropped.** A series of length m returns
  `floor((m - 1) / 2)` rows: 3 at m = 8, 4 at m = 9, 7 at m = 16, 8 at
  m = 17. `spec.pgram()` returns 512 rows for n = 1024 where this
  returns 511.
- **Scaling.** White noise: `mean(pgram)` 8.73672 against `var(y)/fs`
  8.73470 at fs = 1, and 0.03413 against 0.03412 at fs = 256.
- **Mean removal is a no-op on what is returned.** Adding 1000 to the
  series and asking for `detrend = "none"` reproduces the mean-removed
  ordinates to 1.9e-12 relative. This is not a coincidence to be
  documented away: the retained frequencies are exactly the ones where
  the DFT of a constant vanishes, so the test doubles as a guard on the
  grid. `detrend = "linear"` does change the ordinates.
- **Tapers preserve `E I`.** `mean(pgram)/var(y)` is 1.0002 (none),
  0.9915 (Hann), 0.9956 (split cosine), because the divisor is
  `sum(h^2)` rather than `m`.
- **Segment averaging.** `segments = k` gives `floor((m/k - 1)/2)` rows
  and a coefficient of variation of 0.981 (k = 1), 0.492 (k = 4), 0.336
  (k = 8) against the model's 1, 0.5, 0.354.
- **Multi-series.** Matrix columns and a `group` vector both work, group
  lengths may differ (600 and 424 give 299 + 211 = 510 rows), and asking
  for both at once is refused.

## Item 2: the family wrapper

`whittle()` is `exponential(log)` at `tapers = 1` and `Gamma(log)` with
the shape fixed at `tapers` above it. The fixing goes through the
existing constant-dpar mechanism (`plain_dpar()` ->
`betad_fixed_idx` -> `map$betad`), reached by a new family slot
`fixed_dpars` that mirrors the existing `default_forms` slot, six lines
in `R/parse.R:1622`. Verified: `fixef(fit)$shape` is `log(4)` exactly,
`map$betad` is `NA`, the fit has 3 free parameters against 4 for a free
shape, and `predict(type = "disp")` returns 4.

### The refusals, and the one that cannot be written

The response is the only thing `valid_y()` sees. A frequency-grid check
would need the assembled model frame, and the only seam that reaches it
from a family is `frmtmb_structure(check_frame =)`.

**Corrected 2026-09-08.** This section used to say that attaching a
structure costs the family REML, quadrature, profiling,
`conditional_effects()` and cluster-robust standard errors, because
`structure_allows()` reads an absent `supports` flag as a refusal. That
is true only for a structure that carries a `loglik`. For one that does
not, `frmtmb_structure()` (`R/structure.R:390`) defaults `supports` to
all THIRTEEN flags TRUE, because such a structure is a capability
declaration on a family whose likelihood is already rowwise. Measured
by construction: with a `check_frame`-only structure attached to
`whittle()`, REML, quadrature, `profile()`, `conditional_effects()` and
`vcov_cluster()` all work, and behave exactly as they do with no
structure at all; the same structure given a `loglik` refuses all
thirteen. So the route is open and free, and the reason the grid is
still guaranteed at construction by `frm_periodogram()` is only that
`check_frame` reaches the frame, which does not carry a taper either.

The family refuses what the response alone can prove:

1. **Not positive.** Log power and power in dB go negative; this is the
   most common wrong response and it is caught exactly.
2. **Too smooth to be raw.** `var(diff(log(I), lag = 3))` is
   `2 trigamma(k)` for ordinates of shape k, whatever the spectrum is,
   because the log ratio of two ordinates with the same mean has that
   variance. A smooth spectrum and shuffled rows both only ADD to it,
   so a value far below the model's is one-sided evidence that the
   ordinates are not independent exponentials: either already averaged,
   or leakage rather than signal. Calibration, power and the leakage
   case are in the two sections below. The lag was 1 in 0.54.0 and is 3
   from 0.55.0, for the reason in "The Hann false refusal, closed".

The literal reading of "a length that does not match the frequency grid"
is vacuous: `m = floor((n-1)/2)` has a solution `n` for every `m`, and a
user fitting one frequency band legitimately passes a subset. The
dispersion check is what replaces it.

## Item 3: post-processing, and why `simulate()` is refused

Measured first. On a whittle fit, `simulate(fit, nsim = 2)` returns a
data frame of columns `sim_1`, `sim_2` with 511 rows for a series of
1024 points. Nothing in the object says the numbers are spectral
ordinates; the only tell is that the row count is not the series length,
and it is a tell only to a user who already knows.

So `whittle()` ships no simulator and states why (`sim_refusal`), which
is the mechanism `cox()` already uses. All five entry points were then
CALLED, not just read: `simulate()`, `pp_check()`, `frm_simulate()`,
`frm_bootstrap()` and `conditional_effects(method = "predict")` each
refuse and each prints the reason. `pp_check()` was missing from this
list, and from `?whittle` and the vignette, until the review found it;
it is the verb a user reaches for first, so the omission mattered.
`fitted()`, `predict()`, `residuals()` and
`conditional_effects(method = "epred")` all still work on the same fit.
The refusal hands over the two replacements:
`rexp(nobs(fit), 1 / fitted(fit))` for ordinates in the open, and
`frm_series_draw(fit)` for a series whose spectrum is the fitted one.
Everything the refusal costs is one line away; the trap it closes is
silent.

What keeps working is about the periodogram, and is documented in one
place: the family's own help, under a section named for the confusion.

`frm_series_draw()` gives each Fourier frequency an independent complex
Gaussian amplitude with the fitted variance and transforms back. The
length is `2 nf + 1`, odd on purpose: no Nyquist ordinate has to be
invented and the series has exactly zero mean. The sampling rate follows
from the grid spacing, and the result is a `ts`, so
`frm_periodogram(frm_series_draw(fit))` lands back on the same grid.

## WITHDRAWN defect: bf(shape = k) is fine. Do not re-file

`bf(y ~ x, shape = 4)` with `family = Gamma()` appeared to converge to
the wrong answer: the slope came back with the wrong sign against
`glm(family = Gamma(link = "log"))`, the objective was `NaN` at the
origin, and `fixef()` looked inconsistent with `fitted()`.

All three are one thing. `stats::Gamma()` defaults to `link =`
`"inverse"`, so that fit was an inverse-link model compared against a
log-link reference. `mu = 1 / eta` is where the `NaN` at `eta = 0` comes
from, and the reported coefficients reproduce the reported fitted values
exactly once inverted. Re-checked against `glm()` with the matching
link, the constant-dpar route is correct.

The gotcha is worth keeping: `whittle()` is `Gamma(link = "log")` with
the shape pinned, and any reference fit written beside it must say
`link = "log"` out loud, because the `family = Gamma()` spelling does
not mean what the log-link context suggests it means.

**Cause, for the record: `stats::Gamma()` defaults to `link =`
`"inverse"`.** There is no bug in the constant-dpar route, nothing to
fix, and nothing to re-file. `tests/testthat/test-spectral.R` pins the
correct comparison (`whittle(tapers = 4)` against
`glm(family = Gamma(link = "log"))`, agreeing to 1e-4) so that the same
question does not have to be asked again.


## Item 2, continued: the dispersion refusal, rebuilt in the punch round

### What it was, and why that was wrong

The first rule was a flat threshold at half the model value with a hard
floor: below 200 ordinates the check returned early and did nothing.
The floor was calibrated honestly (under the flat spectrum, which
MINIMIZES the statistic, 1 sample in 20000 of 200 ordinates reached the
threshold and none of 500 or more did) but the consequence was not
measured, and the review measured it:

| untapered exponent-3 power law | ordinates | refused, old rule |
| --- | --- | --- |
| n = 256 (a one-second epoch at 256 Hz) | 127 | **0%** |
| n = 512 | 255 | 93% |
| n = 1024 | 511 | 97% |

The check was inert at exactly the length this feature is aimed at, and
worst of all in the row of the tapering table with the largest bias
(-1.04). For the error it was built for it was worse still, because
averaging spends the ordinates the check needs: `segments = 4` needed a
1604-point series and `segments = 8` needed 3208 before it could say
anything at all.

### What it is now

`whittle_smooth_frac()` (`R/spectral.R:281`) replaces the floor with a
threshold that shrinks as the row count falls, so a short response buys
its safety by needing more evidence rather than by not being looked at.
Calibrated at the flat spectrum, 120000 samples per cell, four shapes
(k = 1, 2, 4, 8), taking the smallest 1e-4 ratio across k with a 15%
margin, monotone in nf and capped at the old half:

| ordinates | 1e-4 ratio measured | shipped fraction |
| --- | --- | --- |
| 24 | 0.144 | 0.123 |
| 32 | 0.201 | 0.171 |
| 48 | 0.278 | 0.236 |
| 64 | 0.329 | 0.279 |
| 96 | 0.423 | 0.359 |
| 128 | 0.470 | 0.400 |
| 192 | 0.536 | 0.455 |
| 256 | 0.581 | 0.494 |
| 384, 512 and up | 0.651, 0.679 | 0.500 (capped) |

Below 24 ordinates it still does not run: no threshold separates the
two hypotheses there, and a series that short has no spectrum worth
fitting.

### What that bought, measured against the review's own cells

| case | ordinates | old | new |
| --- | --- | --- | --- |
| exponent-3 leakage, n = 256 | 127 | 0% | **89%** |
| exponent-3 leakage, n = 512 | 255 | 93% | 95% |
| exponent-3 leakage, n = 2048 | 1023 | 96% | 100% |
| k = 4 declared raw, n = 512 | 63 | 0% | **98%** |
| k = 8 declared raw, n = 512 | 31 | 0% | **98%** |
| k = 4 declared raw, n = 1024 | 127 | 0% | **100%** |
| k = 8 declared raw, n = 2048 | 127 | 0% | **100%** |
| k = 2 declared raw, n = 512 | 127 | 0% | 46% |

The review's exact reproduction line,
`frm_periodogram(rnorm(1024), segments = 4)` fitted with `whittle()`,
was accepted silently and is now refused.

Nothing was paid for it. Under the flat null, 20000 samples per cell,
the realized false-alarm rate is **0.000%** at 24, 32, 64, 127, 255 and
511 ordinates for k = 1 and k = 4. On legitimate spectra, 200
replicates each: AR(1) at phi 0, 0.5, 0.9 all 0% at 127 ordinates,
AR(2) resonances at r = 0.9, 0.98, 0.995 all 0% at 1023 ordinates, a
1/f background with an alpha peak of amplitude 1, 3 and 10 all 0% at
1023, ten stacked 64-point series 0% at 310. **The lengths are part of
these claims**, added 2026-09-08 because leaving them out let a later
lane read the AR(2) entry as a claim about a short record and report
it as irreproducible. It reproduces at 1023 ordinates and does not
hold at 127: see the AR(2) note in "The Hann false refusal, closed".

Two honest limits, both now in `?whittle`:

- **k = 2 below about 200 ordinates** (46% at 127, re-measured at 54%
  over 2000 replicates). Gamma(2) and Gamma(1) are too close to
  separate on that many rows.
- **A near-unit-root AR(1) at a short length** fires occasionally:
  0.5% at phi = 0.99 with 127 ordinates and 8.5% at phi = 0.999 with
  255. The review saw 10% in the second case independently. Re-measured
  over 2000 replicates with a 5000-sample burn-in, those are 0.2% and
  6.3%. Those spectra fall like `f^-2` across most of the band, so an
  untapered periodogram of one really is leakage-dominated; it is close
  to a true positive rather than a false one, but it is a refusal a user
  can meet on a legitimate model and it is documented as such.

### A third limit the 0.54.0 table did not reach: the Hann taper

A Hann-tapered periodogram is legitimately raw, and at lag 1 it was
refused 1.9% of the time at 255 ordinates. That is item 7 of the
0.54.0 follow-ups and it is closed in the section below.

### The one-sidedness argument, corrected

The old text claimed the two extra terms "can only ADD". A random
reordering does add (a shuffled raw periodogram gives 3.46 against a
model 3.29, and is accepted). A **monotone** one subtracts, heavily: a
raw periodogram sorted by its own power gives 0.002 and is refused
every time. So the statistic assumes frequency order, which is what
`frm_periodogram()` returns, and that is now what the code comment,
`?whittle` and a test all say.

## The Hann false refusal, closed (spectral2 lane, 2026-09-08)

Follow-up item 7. The statistic now compares ordinates THREE apart.
Everything below is measured on the shipped code, from a private
library, at frmtmb 0.54.0 plus this change.

### Why the taper deflated the statistic, and why three

A periodic Hann window's transform has exactly three non-zero taps, so
a tapered ordinate at bin j is built from bins j-1, j and j+1 alone.
Neighbors share two of those bins; ordinates three apart share none.
Log-periodogram autocorrelation of white noise, n = 1024, 2000
replicates:

| window | lag 1 | lag 2 | lag 3 | lag 4 |
| --- | --- | --- | --- | --- |
| none | -0.001 | -0.004 | -0.003 | -0.003 |
| hann | **0.308** | 0.013 | -0.003 | -0.005 |
| hamming | 0.265 | 0.006 | -0.003 | -0.005 |
| split cosine | 0.011 | 0.005 | 0.003 | 0.001 |
| blackman | **0.414** | 0.059 | -0.001 | -0.005 |

Correlation SUBTRACTS from `var(diff(log I))`: under Hann the mean fell
from 3.290 to 2.274 at lag 1, and is 3.292 at lag 3. Lag 2 already
works (3.237) but only nearly; lag 3 is exact for a three-tap window,
and it also takes Blackman, which lag 2 does not.

### The rates, paired on the same responses

Every cell below is one set of simulated responses scored twice, once
with the lag-1 rule and once with the lag-3 rule, so the two columns
are paired and the difference is not a seed difference. `se` is the
binomial Monte Carlo standard error of the cell.

| Hann-tapered, legitimately raw | reps | lag 1 (before) | lag 3 (after) |
| --- | --- | --- | --- |
| white noise, n = 256 (nf 127) | 2000 | 0.8% (0.20) | **0.0%** (0.00) |
| white noise, n = 512 (nf 255) | 2000 | **1.9%** (0.31) | **0.0%** (0.00) |
| white noise, n = 1024 (nf 511) | 2000 | 0.1% (0.09) | 0.0% (0.00) |
| white noise, n = 2048 | 2000 | 0.0% (0.00) | 0.0% (0.00) |
| white noise, n = 4096 | 2000 | 0.0% (0.00) | 0.0% (0.00) |
| AR(1) 0.6, n = 512 | 2000 | 1.9% (0.31) | 0.0% (0.00) |
| AR(1) 0.6, n = 1024 | 2000 | 0.0% (0.05) | 0.0% (0.00) |
| exponent 3, n = 1023 | 2000 | 0.0% (0.05) | 0.0% (0.00) |
| exponent 3, n = 256 | 2000 | 0.4% (0.15) | 0.0% (0.00) |
| hamming, n = 512 | 2000 | 1.0% (0.22) | 0.0% (0.00) |
| hamming, n = 1024 | 2000 | 0.0% (0.05) | 0.0% (0.00) |
| blackman, n = 512 | 2000 | **18.1%** (0.86) | 0.0% (0.00) |
| blackman, n = 1024 | 2000 | **10.8%** (0.70) | 0.0% (0.00) |
| split cosine, n = 1024 | 2000 | 0.0% (0.00) | 0.0% (0.00) |

Every "after" cell is 0 out of 2000, so the one-sided 95% upper bound
on each is 0.15%. The worst cell went from 1.9% to under 0.15%, which
is the debt paid.

Detection is unchanged. Paired, 2000 replicates a cell, with the
McNemar standard error of the difference:

| case | lag 1 | lag 3 | difference |
| --- | --- | --- | --- |
| k = 4 declared raw, n = 512 (nf 63) | 98.25% | 98.60% | +0.35 (se 0.32) |
| k = 8 declared raw, n = 512 (nf 31) | 99.10% | 99.25% | +0.15 (se 0.26) |
| k = 4 declared raw, n = 1024 | 100.0% | 100.0% | 0 |
| k = 8 declared raw, n = 2048 | 100.0% | 100.0% | 0 |
| k = 2 declared raw, n = 512 | 55.65% | 56.45% | +0.80 (se 1.05) |
| k = 2 declared raw, n = 1024 | 96.2% | 96.2% | 0 |
| leakage exponent 3, n = 256 | 92.65% | 91.65% | **-1.00** (se 0.27) |
| leakage exponent 3, n = 512 | 95.40% | 95.05% | -0.35 (se 0.17) |
| leakage exponent 3, n = 1024 | 95.85% | 95.65% | -0.20 (se 0.12) |
| leakage exponent 3, n = 2048 | 96.10% | 95.95% | -0.15 (se 0.09) |
| leakage exponent 2.5, n = 1024 | 76.95% | 76.60% | -0.35 (se 0.21) |
| leakage exponent 2.0, n = 1024 | 12.65% | 11.80% | -0.85 (se 0.28) |

**What the wider step costs is one point of leakage detection at 127
ordinates** and nothing separable from zero above that. It is a real
loss, three to four standard errors at the two shortest cells, and it
is the price of the 1.9% that was being charged to correct users. The
averaged-periodogram rows move the other way and none of those moves
is significant.

At 6000 replicates a cell, all four lengths cut from the same long
draw, exponent-3 leakage is caught 91.5%, 94.6%, 95.7% and 96.3% at
n = 256, 512, 1024 and 2048 (se under 0.4). Exponent 2.5 at nf = 1023
is 76.7% and exponent 2 is 9.8%: the vignette's earlier 81% and 14%
came from 100 replicates and were noise around these.

False alarms on legitimate untapered spectra, 2000 replicates a cell,
AR series with a 5000-sample burn-in:

| case | lag 1 | lag 3 |
| --- | --- | --- |
| AR(1) phi 0, 0.5, 0.9 at nf 127 and 255 | 0.0% | 0.0% |
| AR(1) phi 0.99, nf 127 | 0.2% | 0.2% |
| AR(1) phi 0.999, nf 127 | 5.3% | 4.3% |
| AR(1) phi 0.999, nf 255 | 6.6% | 6.3% |
| AR(2) resonance r = 0.9, nf 127 | 0.0% | 0.0% |
| AR(2) resonance r = 0.98, nf 127 | 0.4% | 0.1% |
| AR(2) resonance r = 0.995, nf 127 | **16.6%** | **5.4%** |
| 1/f + alpha peak, amplitude 1, 3, 10 | 0.0% | 0.0% |
| 1/f alone, nf 511 | 0.0% | 0.0% |
| ten stacked 64-point series | 0.0% | 0.0% |

**The r = 0.995 row is a different cell from the one the 0.54.0 round
reported, not a contradiction of it.** Corrected 2026-09-08 by the
review of this lane. The 0% in the 0.54.0 record was measured at 1023
ordinates, which `dev/reviews/2026-09-08-freq.md` states in its own
table; the row above is 127 ordinates. Re-measured on an unmodified
0.54.0 build, 1000 replicates, poles at radius 0.995 and normalized
frequency 0.1: 16.2% at 127 ordinates, 5.9% at 255, 1.1% at 511, 0.1%
at 1023 and 0 in 1000 at 2047. The rate falls with length because a
longer record resolves the pole, so the periodogram stops being
leakage-dominated. The released claim is reproducible at the length it
was taken; it is length-conditional and the record now says so.

The near-unit-root AR(2) resonance improves by a factor of three,
because a wider step admits more of the spectrum's own drift, and drift
only adds. That was not the point of the change and it is the largest
single improvement in the table.

Under the flat null, 20000 samples per cell at nf 24, 32, 64, 127, 255
and 511 for k = 1 and k = 4, the rate is 0.000% before and after.

### The threshold curve did not have to move, and why

The curve is calibrated at a FLAT spectrum, where the ordinates are iid
whatever the step is, so the null distribution of the statistic does
not depend on the lag. Re-measured at lags 1, 2 and 3 over 120000
samples a cell and four shapes, the 1e-4 quantile of the statistic as a
ratio to `2 trigamma(k)` (minimum over k = 1, 2, 4, 8):

| ordinates | lag 1 | lag 2 | lag 3 | shipped fraction |
| --- | --- | --- | --- | --- |
| 24 | 0.160 | 0.154 | 0.136 | 0.123 |
| 32 | 0.213 | 0.190 | 0.198 | 0.171 |
| 48 | 0.284 | 0.284 | 0.277 | 0.236 |
| 64 | 0.327 | 0.332 | 0.323 | 0.279 |
| 96 | 0.406 | 0.404 | 0.412 | 0.359 |
| 128 | 0.459 | 0.442 | 0.456 | 0.400 |
| 192 | 0.530 | 0.540 | 0.539 | 0.455 |
| 256 | 0.588 | 0.588 | 0.589 | 0.494 |
| 384 | 0.636 | 0.640 | 0.643 | 0.500 |
| 512 | 0.685 | 0.678 | 0.680 | 0.500 |

The three columns agree to the spread of a 1e-4 quantile estimated
from 120000 samples, so `whittle_smooth_frac()` ships unchanged. The
one place the margin is thin is 24 ordinates, where the shipped 0.123
sits above the lag-3 estimate times 0.85. Measured directly rather than
through a quantile, 500000 flat-null samples at the shipped threshold:

| ordinates | k = 1 | k = 2 | k = 4 | k = 8 |
| --- | --- | --- | --- | --- |
| 24, lag 1 | 9 | 6 | 1 | 0 |
| 24, lag 3 | 17 | 3 | 6 | 7 |
| 32, lag 1 | 12 | 2 | 1 | 2 |
| 32, lag 3 | 16 | 5 | 5 | 0 |
| 64, lag 3 | 9 | 1 | 0 | 0 |
| 128, lag 3 | 4 | 0 | 0 | 0 |

Counts, not rates: 17 in 500000 is 3.4e-5, a third of the 1e-4 the
curve is calibrated to and about twice the lag-1 count at the same
cell. The cap at 0.5 for 384 ordinates and up is left alone
deliberately. The null tail would allow 0.65, but the fraction is what
holds the near-unit-root cases at 5% rather than the null at 1e-4, so
raising it would buy detection with false refusals on legitimate steep
spectra.

### Item 8: the assertion that was balanced on its seed

`tests/testthat/test-spectral.R` fitted a Hann-tapered exponent-3
periodogram and needed the refusal not to fire. On the pinned seed 21
the statistic sat 1.289 times the trigger; over 2000 fresh seeds the
rule refused 0.25% of them and the smallest margin seen was 0.907, so
the test failed on about one seed in 400. At lag 3 the same seed sits
at 2.113 times the trigger, the refusal rate over the same 2000 seeds
is 0, and the smallest margin is 1.438.

The test is now two things. A rate over 1000 seeds, read off the
shipped refusal and compared against the untapered response the
threshold was calibrated on, which fails 30 times out of 30 on the
lag-1 rule (10 to 26 refusals per 1000 against a control of 0) and
holds 30 times out of 30 on the shipped one. And a deterministic
margin assertion on the pinned seed, which needs no second copy of the
statistic: the statistic is a variance of log differences, so `y^a`
scales it by `a^2` exactly, and asking the shipped refusal to accept
`y^(1/sqrt(1.5))` asserts that the response clears the trigger by half
again.

## Item 4: the small-sample bias, measured

What is compared: two estimators of the same parameter, on the same
simulated AR(1) series with innovation variance 1. **Whittle** is the
profile likelihood over phi with the innovation variance concentrated
out, evaluated on `frm_periodogram()`s own output (mean removed by the
dropped zero ordinate, no taper). **exact** is
`arima(x, order = c(1, 0, 0), method = "ML")`, exact Gaussian maximum
likelihood on the series itself, which estimates a mean as Whittle
implicitly does. Lengths 64, 128, 256, 512, 1024, 2048, 4096; phi 0,
0.3, 0.5, 0.7, 0.9, 0.95; 400 replicates per cell to n = 1024, 200 at
2048, 100 at 4096. Monte Carlo standard error on each Whittle cell runs
0.0042 to 0.0064 at n = 64, 0.0014 to 0.0032 at n = 256, and 0.0005 to
0.0015 at n = 4096, so the -0.041 against -0.061 comparison at n = 64 is
roughly three standard errors and the sub-0.003 entries at n >= 1024 are
within noise of zero. Bias in phi:

| n | phi = 0 | 0.3 | 0.5 | 0.7 | 0.9 | 0.95 |
| --- | --- | --- | --- | --- | --- | --- |
| **Whittle** | | | | | | |
| 64 | -0.006 | -0.016 | -0.027 | -0.031 | -0.041 | -0.044 |
| 128 | 0.004 | -0.010 | -0.014 | -0.018 | -0.020 | -0.023 |
| 256 | -0.006 | 0.000 | -0.004 | -0.010 | -0.012 | -0.013 |
| 512 | -0.002 | 0.001 | -0.001 | -0.006 | -0.006 | -0.007 |
| 1024 | -0.003 | -0.001 | -0.001 | -0.003 | -0.002 | -0.002 |
| 2048 | 0.000 | -0.004 | 0.000 | -0.001 | -0.001 | -0.002 |
| 4096 | 0.000 | -0.002 | -0.003 | 0.000 | 0.000 | 0.000 |
| **exact ML** | | | | | | |
| 64 | -0.022 | -0.031 | -0.039 | -0.045 | -0.061 | -0.060 |
| 128 | -0.005 | -0.017 | -0.021 | -0.025 | -0.029 | -0.026 |
| 256 | -0.011 | -0.004 | -0.008 | -0.014 | -0.016 | -0.007 |
| 512 | -0.004 | -0.001 | -0.003 | -0.008 | -0.008 | 0.004 |
| 1024 | -0.003 | -0.001 | -0.002 | -0.004 | -0.003 | 0.018 |
| 2048 | -0.001 | -0.004 | 0.000 | -0.001 | -0.002 | 0.022 |
| 4096 | 0.000 | -0.002 | -0.003 | 0.000 | 0.000 | 0.031 |
| **Hann taper** | | | | | | |
| 64 | -0.014 | -0.031 | -0.036 | -0.042 | -0.060 | -0.062 |
| 256 | -0.005 | 0.001 | -0.006 | -0.017 | -0.020 | -0.019 |
| 1024 | -0.003 | 0.001 | -0.002 | -0.005 | -0.003 | -0.004 |
| **debiased** | | | | | | |
| 64 | -0.006 | -0.012 | -0.019 | -0.020 | -0.031 | -0.037 |
| 256 | -0.006 | 0.001 | -0.002 | -0.007 | -0.008 | -0.009 |
| 1024 | -0.003 | 0.000 | 0.000 | -0.002 | -0.001 | -0.001 |

Three conclusions, none of them the expected one.

1. **Whittle is not the worse estimator here.** Its absolute bias is
   the smaller one in 40 of the 42 cells **on this seed**, and the
   count is seed-dependent: the review reran it and reproduced the
   n = 256, phi = 0.95 exception almost exactly (-0.0128 against
   -0.0062, here -0.0129 against -0.0071) but not the n = 512 one,
   where it measured -0.0055 against +0.0085, about two standard
   errors from this run's -0.0072 against +0.0043. So the durable
   statement is "Whittle is not worse anywhere in this grid, and is
   better in nearly every cell", not a count. Both exceptions sit at
   phi = 0.95, for the reason conclusion 2 gives: the reference has
   started drifting upward there, which partly cancels its own downward
   small-sample bias.
   Most of what the table shows is the ordinary small-sample bias of an
   autoregressive estimate, which both methods pay, rather than anything
   Whittle adds. The claim is deliberately narrow: one model (AR(1)),
   one reference implementation, a mean estimated on both sides, and
   this grid. It says Whittle costs nothing here; it does not say
   Whittle beats exact maximum likelihood in general.
2. **At phi = 0.95 the reference is what fails, rather than Whittle
   succeeding.** `arima()` drifts toward the unit root there: median
   0.978 at n = 1024 and 0.990 at n = 4096, against a Whittle median of
   0.949 and 0.950 on the same series. The positive "exact" numbers in
   the last column, and the whole RMSE advantage Whittle appears to hold
   in that column (a ratio of 0.119 at n = 4096), are the behaviour of
   that optimizer near the boundary. That column is evidence about
   `arima()`, not about the two likelihoods.
3. **The debiased variant is worth about a third of the remaining
   bias**, and Hann tapering makes it WORSE for an AR(1) at every
   length, because that spectrum is not steep enough for leakage to be
   the problem.

Safe: n of 1024 or more, where nothing here exceeds 0.003. Marginal: 256
to 512, a one to two second epoch at 256 Hz, where one to three percent
of a strongly autocorrelated parameter is worth stating in a paper.

### Where tapering does pay: a steep spectrum

A power law `f^-chi` synthesized long and then CUT, so the segment is
not periodic. Bias in the estimated exponent, 300 replicates:

| chi | n | raw | Hann | split cosine |
| --- | --- | --- | --- | --- |
| 1 | 256 | 0.000 | -0.024 | -0.003 |
| 1 | 1024 | 0.005 | -0.002 | 0.002 |
| 2 | 256 | -0.084 | 0.017 | 0.061 |
| 2 | 1024 | -0.082 | 0.005 | 0.027 |
| 3 | 256 | **-1.033** | 0.087 | 0.331 |
| 3 | 1024 | **-1.094** | 0.029 | 0.148 |

At chi = 3 the raw periodogram does not estimate the exponent at all: it
returns 1.967 at n = 256 and 1.906 at n = 1024, because the leakage
floor of an untapered transform falls like `f^-2` and anything steeper
disappears under it. The chi = 2 row lands at 1.916 and 1.918, so true
exponents of 2 and 3 both come back between 1.91 and 1.97: the estimate
saturates near the floor instead of following the spectrum, which is
what "the floor is what reaches the model" means. At chi = 1, which is
flatter than the floor, the estimate is unbiased.

**A longer recording does not help; it makes it worse.** This
is the practical result of the lane for the aperiodic-exponent audience,
where exponents of 2 to 3 are ordinary, and the taper argument is what
answers it.

### The dispersion refusal catches most of this, unplanned

Leakage-dominated ordinates are not independent exponentials either, so
they leave the same fingerprint as an already-averaged periodogram and
the same one-sided check sees them. Measured at n = 1023, 100
replicates, raw periodogram:

| true exponent | var(diff(log I)) in one draw | refused |
| --- | --- | --- |
| 1.0 | 3.818 | 0% |
| 2.0 | 0.917 | 14% |
| 2.5 | 3.091 | 81% |
| 3.0 | 0.014 | 97% |

At exponent 3 the ordinates are nearly deterministic (0.014 against a
model value of 3.29): there is no signal left in them at all. The
message therefore names both causes and both remedies, and the vignette
shows the refusal firing rather than a wrong number coming back.

**The gap is at exponent 2**, caught only 14% of the time, where the
untapered fit returns a quiet -0.08 with small standard errors and
well-behaved residuals, because the model is describing the leakage
correctly. The vignette says plainly that the refusal is a backstop for
the extreme case and not a test to lean on: taper on the shape of the
spectrum, not on whether an error appeared.

### The debiased variant is not expressible here, and what it would take

The mean at one frequency is a triangle-weighted sum over lags of the
model autocovariance, so it is not a function of that row's predictors.
Writing it needs a term whose value at row j is a linear functional of a
vector computed once per likelihood evaluation: an O(n log n) transform
of a parameter-dependent sequence, on the tape. RTMB 1.9 has the complex
AD for it and the grammar has no term type that can take it. That is the
same seam the complex cross-spectrum needs, so it belongs with item 5's
extension rather than with core.

## The AR(1) spectrum has two maxima, and the vignette says so

`ls - log(1 - 2 ph cos w + ph^2)` with an unconstrained `ph` returns
1.1288861 as readily as 0.8858289 on the same data, and
1 / 1.1288861 = 0.8858290. The two log spectra differ by a constant to
1.3e-15: a spectrum cannot tell phi from 1 / phi, and only the
stationary root is the answer. `tanh` is what the vignette and the test
use.

With that parameterization the shipped fit reproduces the direct profile
optimizer at every cell checked:

| n | phi | frm | direct | difference |
| --- | --- | --- | --- | --- |
| 128 | 0.3 | 0.3612968 | 0.3612982 | -1.4e-06 |
| 128 | 0.9 | 0.8858288 | 0.8858289 | -1.5e-07 |
| 512 | 0.3 | 0.3292284 | 0.3292284 | -9.8e-10 |
| 512 | 0.9 | 0.8990039 | 0.8990038 | 4.0e-08 |
| 2048 | 0.3 | 0.2824265 | 0.2824265 | -1.7e-09 |
| 2048 | 0.9 | 0.8931927 | 0.8931928 | -3.8e-08 |

So the bias table describes the package rather than a second
implementation of it.

## Item 5: the vignette

`vignettes/spectral.Rmd`, tutorial-shaped, renders in 8.4 s. It fits a
parametric spectrum against `arima()`, a nonparametric one whose
smoothing parameter comes from the same Laplace marginal likelihood, a
hierarchical one in which the per-subject aperiodic exponent is a random
SLOPE on `log f` (correlation 0.997 with the truth over 12 subjects),
and a peak on an 8-segment averaged periodogram with
`whittle(tapers = 8)`, recovering centre 10.298, exponent 1.412, width
1.233 and height 0.994 against truths 10, 1.4, 1.5 and 1. It then starts
the same peak's centre at 30 Hz and shows what the note predicted: the
centre drifts to 28.57 rather than to the peak, its standard error is
NaN, and the log-likelihood is 135.70 against 211.80 for the good fit.

Between the hierarchical fit and the peak it shows the leakage failure
live: an untapered periodogram of an exponent-3 spectrum is REFUSED by
`whittle()`, and `taper = "hann"` on the same series returns 3.018.

It carries the bias table, the tapering rule, and the three things this
framework does not do: phase-amplitude coupling (the bispectrum has no
exponential marginal), a non-stationary epoch (needs segmentation), and
two signals (complex Wishart, a family rather than a formula).

## Verification, punch round (2026-09-08)

Everything below was re-run after the punch fixes, on the frozen source.

### Test suite, one file per process

| | measured |
| --- | --- |
| test files on disk | 132 |
| RESULT lines | 132 |
| distinct names | 132 |
| duplicates | **0** |
| passing | **7549** |
| failing | **0** |
| errors | 0 |
| skipped | 133 |
| `test-spectral.R` | 109 checks, all passing |

Fully clean. The previous round's single failure was `test-perf.R`, and
this run reproduces the review's diagnosis rather than my earlier one:
the file is byte-identical to main, its two assertions are pure wall
clock, and the reviewer got pass / fail / pass on three ISOLATED runs,
so it is flaky in itself on a loaded machine, not merely contended. It
passed here because this run had the machine mostly to itself. Nothing
about it is attributable to this lane, and the flakiness is worth a
ticket from whoever owns it.

The pass count rose 7533 -> 7549 for two reasons: `test-perf.R`
contributes 3 rather than 2, and `test-spectral.R` grew from 94 to 109
with the three tests the punch round added (short-response coverage,
the shrinking threshold, and frequency-order sensitivity).

### Roxygen

Idempotent: the second `roxygenise()` wrote nothing.

### R CMD check --as-cran, with the manual

Two environment fixes from the review, both adopted: PATH set from
PowerShell with Windows semicolons (a `C:/...` entry on a
colon-separated PATH from Git Bash is read as a separator, which is why
pandoc did not resolve last round), and `R_BIOC_VERSION=3.21` pinned.

**Exit 0. Status: 2 NOTEs, 0 WARNINGs, 0 ERRORs**, down from 4 NOTEs.

1. **`checking CRAN incoming feasibility ... [23s] NOTE`** -
   "Maintainer: ... / New submission". The package is not on CRAN.
   Environmental. The pin cut this from 292s to 23s and removed six
   `cannot open URL` warnings.
2. **`checking HTML version of manual ... [21s] NOTE`** - "Skipping
   checking math rendering: package 'V8' unavailable". The expected
   one.

The two that are gone:

- **`checking top-level files`** is now **OK**: with
  `quarto/bin/tools` on PATH, `Sys.which("pandoc")` resolves and the
  NEWS.md check runs. Exactly the behaviour the review predicted by
  removing it.
- **`checking examples`** is now **[46s] OK** with no NOTE at all. Last
  round it named `profile.frmtmb_fit`, and the review saw four topics
  (`residuals.frmtmb_fit`, `profile.frmtmb_fit`, `dharma_residuals`,
  `VarCorr`) on a loaded machine. On a quiet one every example is under
  the threshold, so that NOTE was contention, not content. The review's
  observation that `residuals.frmtmb_fit` used 14.85s of CPU stands as
  a pre-existing cost worth someone's attention; it is not this lane's.

Everything this lane could have broken is OK, including
`checking R code for possible problems ... OK`, `checking Rd files ...
OK`, `checking Rd cross-references ... OK`, `checking for missing
documentation entries ... OK`, `checking examples with --run-donttest
... [33s] OK`, `checking tests ... Running 'testthat.R' [236s] OK` (the
whole suite in one process on the tarball), `checking re-building of
vignette outputs ... [273s] OK` (the vignette builds, refusal blocks and
all), and `checking PDF version of manual ... OK` under TinyTeX.

### The docs index

`_pkgdown.yml` indexes the three exports and the vignette; the YAML
re-parses to 8 reference sections and 10 articles. `R CMD check` does
not read this file, so a lane that only ran the check would have handed
the docs CI a failure.

## Seams recorded for other lanes

1. ~~**A rowwise family cannot reach the model frame.**~~ **WITHDRAWN
   2026-09-08. Do not act on this.** The route is
   `frmtmb_structure(check_frame =)` and it costs nothing: a structure
   with `loglik = NULL` already defaults every one of the thirteen
   `supports` flags to TRUE (`R/structure.R:390`), which is exactly the
   "check_frame-only structure that keeps the rowwise defaults" this
   entry asked someone to write. Verified by construction rather than
   by reading: with such a structure attached to `whittle()`, REML,
   quadrature, `profile()`, `conditional_effects()` and
   `vcov_cluster()` all behave as they do with no structure at all,
   while the same structure given a `loglik` refuses all thirteen. The
   `check_frame` receives the assembled frame, including its
   `data_frame`. It does not help the taper question, because the frame
   holds no record of the taper either.
2. **`sim_refusal` prints after "has no simulator yet"**
   (`R/predict.R:2832`, `R/simulate-new.R:578`, `R/compat.R:1343`). For
   a deliberate refusal "yet" is wrong. `cox()` already lives with it,
   so the wording is a shared decision rather than this lane's to make.
3. **Stopping a background suite driver does not stop the suite.**
   `TaskStop` on the task that launched `powershell.exe -File suite.ps1`
   killed the launcher and left the driver running. Two drivers then
   appended to the same log, which showed 127 result lines for 131
   files with 30 of them duplicated, and the interleaving was visible
   only because the file order stopped being alphabetical. Any lane
   auditing a count from a shared log should check
   `grep -o "test-[a-z0-9-]*\.R" log | sort | uniq -d` before trusting
   it, and confirm by PID that the driver is gone after a stop.
4. **`residuals(type = "osa")`** on an exponential-family fit errors with
   `Not compatible with requested type: [type=S4; target=double]`, which
   names neither the family nor the cause. Reproduced on the untouched
   install, so it predates this lane; not investigated further.
