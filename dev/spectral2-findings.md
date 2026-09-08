# spectral2 lane: the Hann false refusal, and a test that was
# balanced on its seed

Date: 2026-09-08. Branch `wt-spectral2` on `a2deb65` (frmtmb 0.54.0
plus the docs rebuild), R 4.6.1, Windows 11. Every number here comes
from a private library; nothing was installed into the shared one.

Scope: follow-up items 7 and 8 of `dev/feature-gaps.md`. The
measurement tables are in `dev/freq-findings.md`, under "The Hann false
refusal, closed", so that the frequency-domain record stays in one
place. This file is the decisions.

## The change, in one line

`whittle()`'s smoothness statistic is
`var(diff(log(y), lag = 3))` rather than `var(diff(log(y)))`. Nothing
else moved: the threshold curve, the 24-ordinate floor, the positivity
refusal and the family's arithmetic are untouched.

## Why that is the right change and not a weakened refusal

The refusal exists to catch two errors. Neither travels through the
differencing lag in the form this lane measured, but one of them has a
second form that does; see the correction below the two bullets.

- **A SEGMENT-averaged periodogram declared raw** is deflated through
  its MARGINAL shape: ordinates of shape k have
  `var(diff(log I, lag)) = 2 trigamma(k)` at every lag, because Welch
  and Bartlett averaging leaves them independent across frequency. No
  lag can see around that and no lag helps it. Measured: the detection
  rate moves by at most 0.8 points and not once significantly, at
  k = 2, 4 and 8 and at 31, 63, 127 and 255 ordinates. This is the
  averaging `frm_periodogram(segments =)` performs.
- **A leakage floor** is smooth across MANY bins, so it deflates the
  statistic at every small lag. Measured: exponent-3 leakage detection
  falls by 1.00 point at 127 ordinates (se 0.27) and by 0.15 to 0.35
  points above that. The review re-measured this on the construction
  `test-spectral.R` uses, 2000 replicates a cell, and adds the
  shallower exponents at the short end, which this lane measured only
  at 1023 ordinates: at 127 ordinates the loss is 1.25 points at
  exponent 3 (se 0.30), 2.30 at exponent 2.5 (se 0.40) and 2.25 at
  exponent 2 (se 0.38); at 255 it is 0.25, 0.95 and 0.85; at 511 it is
  0.15, 0.45 and 0.85; at 1023 it is 0.05, 0.10 and 0.35. So the cost
  at the short end is two points rather than one for a spectrum that
  is only just steep enough to leak.

**Corrected 2026-09-08 by the review of this lane.** The first bullet
above originally said "an averaged periodogram", without the word
SEGMENT, and that generalization is false. An estimate smoothed ACROSS
FREQUENCY, which is the other way to get non-raw ordinates, correlates
them, the correlation was doing part of the detection at lag 1, and
lag 3 gives three bins of it up. Measured on the shipped threshold,
1000 replicates a cell, a plain three-bin Daniell smooth declared raw:
95.4% caught at lag 1 against 17.7% at lag 3 at 32 ordinates, 99.7%
against 37.7% at 48, 100% against 58.2% at 64, 100% against 92.2% at
100, and the two agree at 100% from 128 ordinates up. `spec.pgram(x,
spans = 3)` is that case and is the commonest smoothed-periodogram
idiom in R. A three-taper sine multitaper behaves the same way (85.4%
against 16.5% at 24 ordinates). The trade is still worth taking,
because at lag 1 the same correlation refused the HONEST declaration
of those estimates: a three-bin Daniell declared as `tapers = 3` was
refused 68.6% of the time at 127 ordinates and 97.7% at 255, and a
single Slepian taper 26.9% and 56.1%, all of which go to 0 at lag 3.
But it is a real loss on a real case and it belongs in the record.
Numbers in `dev/reviews/2026-09-08-spectral2.md`.

The Hann taper is the one thing that lives at exactly lag 1. A periodic
Hann window's transform has three non-zero taps, so ordinates one apart
share two bins (log correlation 0.308) and ordinates three apart share
none (-0.003). That is a property of the window, not of the data, so
moving the step past it removes the false refusal without removing any
of the evidence the check reads.

The threshold curve did not have to move, and that is the strongest
single piece of evidence that this is not a weakening. The curve is
calibrated at a FLAT spectrum where the ordinates are iid at any lag,
so the null distribution is lag-invariant. Re-measured over 120000
samples a cell at lags 1, 2 and 3 and four shapes, the 1e-4 quantile
of the statistic agrees across the three lags to within the spread of
the estimate: above 32 ordinates no lag differs from another by more
than 0.017 in the ratio, and at 24 and 32 the gap reaches 0.024, where
a 1e-4 quantile read off 120000 samples is 12 order statistics and is
not worth that much precision. Those two knots were settled by direct
measurement instead, 500000 flat-null samples at the shipped
threshold: 17 refusals against the lag-1 rule's 9 at 24 ordinates and
k = 1, both a third or less of the 1e-4 the curve is calibrated to.
The realized false-alarm rate under the flat null is 0 in 20000 per
cell before and after, at six row counts and two shapes.

## What it cost, stated plainly

One point of leakage detection at 127 ordinates: 92.65% to 91.65%,
paired, 2000 replicates, McNemar se 0.27. Three to four standard
errors, so real. Above 511 ordinates nothing is separable from zero.

What it bought: the worst Hann cell went from 1.9% (se 0.31) to 0 in
2000, a one-sided 95% upper bound of 0.15%. Nine Hann cells, two
Hamming cells, two Blackman cells and one split-cosine cell are every
one of them 0 in 2000. The near-unit-root AR(2) resonance at r = 0.995
fell from 16.6% to 5.4%, which was not the point of the change and is
the largest improvement in the table.

The trade is one point of a backstop against two points of refusals
charged to correct users, and the check is documented as a backstop.

## Lag 2 or lag 3

Lag 2 also clears every window measured: Hann, Hamming and Blackman are
each 0 in 2000 at lag 2, and the mean Hann statistic there is 3.237
against a model 3.290. Lag 3 was chosen on the residual correlation
rather than on a rate. A three-tap window (Hann, Hamming) leaves lag 3
exactly uncorrelated by construction, and measurement agrees (-0.003
and -0.003); at lag 2 the residual is 0.013 for Hann and 0.059 for
Blackman, which is small but is a number that would have to be
re-measured for any window a user brings from elsewhere. The cost of
the extra step over lag 2 is inside the Monte Carlo error of every
detection cell.

## Two routes examined and rejected, with the measurement

### An attribute on the response

`frm_periodogram()` knows the taper, so the obvious cheap fix is
`attr(out$pgram, "frm_taper") <- taper` read back by `valid_y()`. It
does not work, for two independent reasons, both measured rather than
assumed.

1. **The attribute reaches frmtmb's model frame and is then stripped
   inside `extract_y()`** (`R/frame.R:247`, `y <- as.numeric(as.vector(y))`).
   Traced: present at `extract_y()` entry, absent at its exit, absent
   at `valid_y()`. `stats::model.frame()` itself keeps it, so this is
   frmtmb's own coercion, not R's.
2. **Five of eight ordinary user operations drop it anyway.** Kept by
   `data.frame()`, `rbind()` and `transform()`; dropped by `d[i, ]`,
   `subset()`, `na.omit()`, reordering with `[`, and `merge()`. Two of
   those are inside model-frame assembly for a fit with missing rows.

So an attribute would make the refusal depend on how the user handled
the data frame between the constructor and the fit. A refusal that
fires on some subsets of a response and not others is worse than one
that fires 1.9% of the time, because the user cannot reason about it.
Relaxing the coercion in `extract_y()` would change what every family
in the package sees and is not this lane's to do for one taper flag.

### `frmtmb_structure(check_frame =)`

`dev/freq-findings.md` ruled this out on the grounds that any structure
flips every `supports` flag to refused. **That is no longer true and
the entry is withdrawn.** `frmtmb_structure()` (`R/structure.R:390`)
defaults `supports` to all thirteen flags TRUE when `loglik` is NULL,
because such a structure is a capability declaration on a family whose
likelihood is already rowwise. Verified by construction: with a
`check_frame`-only structure attached to `whittle()`, REML, quadrature,
`profile()`, `conditional_effects()` and `vcov_cluster()` all work and
behave as with no structure; the same structure given a `loglik`
refuses all thirteen.

It still does not solve item 7. `check_frame` receives the assembled
frame, and the frame carries no record of the taper: the taper is a
property of a transform that happened before the data frame existed.
The route is open for the frequency-grid check the freq lane wanted,
and that is now recorded correctly for whoever wants it.

## Item 8: the test

The assertion was a Hann-tapered exponent-3 periodogram fitted through
`whittle()`, at pinned seed 21, where the statistic sat 1.289 times the
trigger. Over 2000 fresh seeds the lag-1 rule refused 0.25% and the
smallest margin was 0.907: one seed in 400 failed the suite.

It is now two assertions.

1. **A rate over 1000 seeds**, in a new test, read off the shipped
   `valid_y()` and compared against the untapered response the
   threshold was calibrated on. The slack is three refusals, which is
   binomial rather than a tolerance: at the calibrated ceiling of one
   in 10000 the probability of exceeding it is below 1e-10, while the
   lag-1 rule expects 19. Checked over 30 different seed sets: the
   assertion holds 30 out of 30 on the shipped rule (0 refusals per
   1000 in every set) and fails 30 out of 30 on the lag-1 rule (10 to
   26 per 1000 against a control of 0).
2. **A deterministic margin** on the pinned seed, which does not
   rewrite the statistic in the test. The statistic is a variance of
   log differences, so `y^a` scales it by `a^2` exactly; asking the
   shipped refusal to accept `y^(1/sqrt(1.5))` asserts that the
   response clears the trigger by half again. It clears it by 111%.

Both were run against the unfixed rule and both fail there: the
deterministic one errors with the refusal, the rate one reports 16 > 3.

The rate test does not need a gate. Timed four times against the
baseline file on the same machine, the file costs about 2 seconds more
(2.3, 2.7, 3.5 and 4.2 s before against 4.2, 4.5, 4.9 and 6.5 s after,
the spread being other lanes on the machine), and it grows from 109 to
113 passing checks. It stays cheap because the rate is read from
`valid_y()` on a constructed periodogram rather than from a fitted
model: 3000 periodograms of 512 points cost about 0.5 s.

The third assertion in that test measures the window property the lag
rests on: the lag-1 autocorrelation of a Hann-tapered log periodogram
against the largest of its lag-2 and lag-3 values, asserted as a ratio
so that nothing absolute is pinned. Over 60 seeds the ratio never came
within a factor of 1.9 of its bound.

## What I did not do

- **I did not raise the 0.5 cap on the threshold fraction.** The null
  tail allows 0.65 at 384 ordinates and up, and raising it would lift
  the k = 2 detection gap and the exponent-2 leakage gap. It would also
  lift the false refusals on legitimate near-unit-root spectra, which
  are at 5 to 6% and are what the cap is holding down. That is a
  separate decision with its own table, and it is not item 7.
- **I did not make the lag adaptive.** Choosing the lag from the
  response's own autocorrelation would let a leakage floor, whose
  correlation is long-range, push the lag up until the check could no
  longer see it.
- **I did not touch the k = 2 limit or the near-unit-root limit.** Both
  are documented and both are unchanged by this lane, except that the
  AR(2) resonance case improved by a factor of three as a side effect.

## Verification

| what | result |
| --- | --- |
| `test-spectral.R`, own process | **113 pass**, 0 fail, 0 error, 0 skip |
| the same file at the base commit | 109 pass, 0 fail |
| the two new assertions, unfixed rule | 1 failure, 1 error, named above |
| `roxygen2::roxygenise()` | writes `whittle.Rd`, then nothing |
| `R CMD check --as-cran`, with manual | **2 NOTEs, 0 WARNINGs, 0 ERRORs** |
| that check's stages, manual PDF | 60 stages, 576 KB |

The two NOTEs:

1. `checking CRAN incoming feasibility ... [16s] NOTE`, "New
   submission". The package is not on CRAN. Environmental.
2. `checking HTML version of manual ... [26s] NOTE`, "package 'V8'
   unavailable". The expected environmental one.

Everything this lane could have broken is OK, including
`checking R code for possible problems ... [72s] OK`,
`checking Rd files ... OK`, `checking Rd cross-references ... OK`,
`checking examples ... [40s] OK`,
`checking examples with --run-donttest ... [62s] OK`,
`checking tests ... [290s] OK` (the whole suite in one process on the
tarball) and `checking re-building of vignette outputs ... [320s] OK`,
which is the stage that proves the vignette's two refusal blocks still
fire under the new rule. `checking PDF version of manual ... [10s] OK`.

The check was run twice. The first run, on a tree that was still being
edited, reported a third NOTE at `checking examples`, naming
`dharma_residuals` at 6.36 s elapsed. It is gone in the second run of
the same examples (40 s against 63 s for the stage), so it was machine
contention from three other lanes rather than content, which is what
the 0.54.0 round concluded about the same NOTE. No spectral topic
appeared in it.

I did NOT rebuild `docs/`. The pkgdown site still carries the 0.54.0
text of `?whittle` and the spectral vignette, and needs a rebuild
whenever this round's lanes are consolidated.

## Defects found and not fixed

- **`extract_y()` strips every attribute from the response**
  (`R/frame.R:247`). That is deliberate coercion and probably right,
  but it means no family can ever be handed a fact about its response
  by the code that built it. If a family ever needs one, this is the
  line, and the decision belongs to whoever owns `R/frame.R`.
- ~~**The 0.54.0 AR(2) resonance figure is not reproducible.**~~
  **WITHDRAWN 2026-09-08 by the review of this lane. The released
  figure is correct and needs no correction.** This lane read the
  record's AR(2) entry as a claim about a short record and measured it
  at 127 ordinates, where the lag-1 rate really is about 16%. The
  0.54.0 measurement was taken at 1023 ordinates, which
  `dev/reviews/2026-09-08-freq.md` states in its own table, and it
  reproduces there: on an unmodified 0.54.0 build, poles at radius
  0.995 and normalized frequency 0.1, 1000 replicates, the rate is
  16.2% at 127 ordinates, 5.9% at 255, 1.1% at 511, 0.1% at 1023 and 0
  in 1000 at 2047. It falls with length because a longer record
  resolves the pole and the periodogram stops being leakage-dominated,
  which is also why the burn-in sweep RISES: a longer burn-in gives the
  sharper stationary resonance. What the record was missing is the
  length beside the number, and `dev/freq-findings.md` now carries it.
