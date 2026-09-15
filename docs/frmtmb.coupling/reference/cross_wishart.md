# The complex Wishart family for a pair of signals

The likelihood of a two-channel cross-spectral matrix. The response is
the first auto-spectrum and the rest of the Hermitian matrix rides in
`vreal()`, with the degrees of freedom in `vint()`. Its four
distributional parameters are the two channel powers, the magnitude
squared coherence and the phase, each with its own linear predictor, so
a random effect on coherence or a smooth in coherence over frequency is
an ordinary
[`frmtmb::frm()`](https://aforren1.github.io/frmtmb/reference/frm.html)
formula.

## Usage

``` r
cross_wishart()
```

## Value

A `frmtmb_family` object.

## The formula this family needs

    frm(bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
           pow2 ~ 1, coh ~ 1, phase ~ 1),
        family = cross_wishart(), data = xs)

[`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md)
returns exactly those columns under exactly those names. Every one of
the four dpars needs a formula: a dpar left out is a constant, which is
right for `pow2` in some designs and almost never right for `coh`.

## The density

For `W` the summed Hermitian matrix at one frequency and `n` its degrees
of freedom,

    log f(W) = (n - 2) log det W - tr(S^-1 W) - n log det S - log CGamma_2(n)
    log CGamma_2(n) = log(pi) + lgamma(n) + lgamma(n - 1)

In this package's coordinates both hard parts collapse. Writing `S11`,
`S22` for the powers, `C` for the coherence and `phi` for the phase,
`det S = S11 S22 (1 - C)`, so the log determinant is a sum of three logs
and no determinant is taken; and `tr(S^-1 W)` is

    (S22 w11 + S11 w22 - 2 sqrt(C S11 S22) (cos(phi) w12r + sin(phi) w12i))
      / (S11 S22 (1 - C))

which is real arithmetic throughout. **At two channels the complex
Wishart needs no complex arithmetic on the tape.** RTMB 1.9's
`adcomplex` was used to derive and check this form at general `p`, and
it works, including under a Laplace approximation; but at `p = 2` it was
measured 5.8 times slower than the closed form for the same answer to
six decimals, so the closed form is what ships.

## The trap

a power formula too simple for the data: The likelihood couples the
powers and the coherence. A power model too simple for the spectrum does
not merely leave `mu` and `pow2` wrong: it reweights the rows each
coherence parameter pools, and moves the coherence with them. It
converges without complaint.

**Check it by comparing AIC with and without a frequency term on the
powers.** Nothing here does it for you. A draft of this package warned
automatically from a `fit_check` hook that measured the spread of the
power residuals, and that statistic is not monotone in the coherence
displacement it was meant to catch: it went silent at a 2.7-fold power
spread where the reported coherence was already nearly double the truth,
and it fired at a ratio of 20.5 on a correctly specified model whose
coherences were right. It was withdrawn rather than retuned, because a
check that punishes correct models teaches users to ignore it.

An AIC comparison does work, and works where the hook did not. Two bands
over 400 frequencies at `n = 8`, true coherences 0.35 and 0.05,
`coh ~ band`, varying the power spread across the record:

|              |                    |                       |           |
|--------------|--------------------|-----------------------|-----------|
| power spread | flat powers report | correct powers report | delta AIC |
| 1.6x         | 0.295, 0.083       | 0.330, 0.068          | 80        |
| 2.7x         | 0.268, 0.078       | 0.342, 0.041          | 487       |
| 4.4x         | 0.256, 0.136       | 0.344, 0.054          | 954       |
| 7.3x         | 0.240, 0.206       | 0.342, 0.047          | 1849      |
| 19.6x        | 0.229, 0.458       | 0.341, 0.052          | 3716      |
| 385x         | 0.222, 0.875       | 0.342, 0.059          | 11518     |

The flat fit reports the high band as nine times its true coherence at a
19.6-fold spread and as seventeen times at 385-fold, and AIC rejects it
by 80 or more even at the mildest spread tested. The correctly specified
fits recover 0.35 and 0.05 throughout.

So: fit the powers twice, once with a frequency term and once without,
and read the coherence off the one AIC prefers. **Model the powers
before reading any coherence.**

## What a condition contrast recovers, and what the interval costs

A contrast between two conditions measured in the same subjects needs a
random effect on the grouping the contrast varies WITHIN. `(1 | id)` is
not it: a subject effect shifts both of that subject's conditions
together and cancels out of their difference, so it widens the
intercept's interval and leaves the contrast's alone. Measured on three
seeds of a 16-subject design, `(1 | id)` multiplies the intercept's
standard error by 1.678, 1.937 and 2.096 and the contrast's by 1.005,
1.034 and 1.071. The term that carries the contrast's between-subject
spread is `(1 | id:cond)`. Omit it and the estimate hardly moves while
its interval collapses, which is the failure a single fit cannot show.

Measured over **148 replicates** of the design in
`dev/scale-findings.md`: 40 subjects, 2 conditions, 60 frequencies, 8
segments a row. The truth is on the logit scale, with intercept -0.6, a
condition contrast of **0.5**, a Gaussian bump of 0.8 over frequency,
`sd(id)` 0.35 and `sd(id:cond)` 0.20. Every model below gives the powers
and the phase an intercept, which is correct here. `width` is the mean
interval width relative to the last row's, and `sd(z)` is the spread of
`(estimate - truth) / se` over replicates, where 1.0 is an honest
interval:

|                                   |                      |       |       |
|-----------------------------------|----------------------|-------|-------|
| coherence formula                 | coverage             | width | sd(z) |
| condition only                    | 0.493 (0.414, 0.573) | 0.37  | 2.75  |
| and a smooth over frequency       | 0.486 (0.407, 0.566) | 0.38  | 2.75  |
| and a subject effect              | 0.527 (0.447, 0.606) | 0.37  | 2.77  |
| crossed effect, no subject effect | 1.000 (0.975, 1.000) | 1.90  | 0.56  |
| both, the correct model           | 0.953 (0.906, 0.977) | 1.00  | 1.06  |

The correct model covers at its nominal rate and recovers both
components, `sd(id)` 0.341 plus or minus 0.049 and `sd(id:cond)` 0.197
plus or minus 0.027. The model without `(1 | id:cond)` puts the contrast
within 0.0036 of the correct model's on the same data, and reports an
interval two and three quarter times too narrow: on 63 of the 148
replicates it missed where the correct model covered, and it never
covered where the correct model missed. **That is the rung a wider
interval would repair, and it is the only one.**

**The first two rows are a different failure.** A model with no random
effect at all is not estimating 0.5 badly. A logit contrast is not
collapsible, so a model that leaves variance out of the linear predictor
is consistent for the population-averaged contrast, which is smaller.
The factor `1 / sqrt(1 + 0.346 V)`, with `V` the omitted variance taken
from the truth and nothing fitted, predicts -0.0196 for the first row
against -0.0197 measured and -0.0135 for the second against -0.0133.
Their coverage is low for that reason and a wider interval does not
repair it: the remedy is to model the variance, or to say that the
marginal contrast is the quantity wanted. `dev/scale-findings.md`
recorded this attenuation, at 0.016 across the ladder, before the
replicate study ran.

**The direction depends on which term is missing.** Dropping `(1 | id)`
and keeping `(1 | id:cond)` goes the other way, 1.9 times too wide,
because the subject variance is then forced into the crossed term. So
the rule is not "more random effects" and not "fewer": a within-subject
contrast needs BOTH terms, the subject one for the intercept and the
crossed one for the contrast. The two one-term models fail in OPPOSITE
directions, which is what makes the pair worth naming: the crossed term
alone is conservative and covers 148 of 148, while the subject term
alone is anti-conservative and covers 78. An interval that is too wide
costs power; one that is too narrow reports a finding that is not there.

The width is what changes, not the arithmetic. On the same generator
with `sd(id:cond)` set to 0, over 36 replicates, the two models agree to
a paired width ratio of 1.000 at the median and 1.264 at its largest,
and the extra component is estimated at 0.015 plus or minus 0.018. One
number says it best: on data that HAS the component, the model without
it reports a mean width of 0.070, and 0.072 is what the correct model
reports on data that does not. Leaving the term out answers a different
experiment's question.

`dev/coh-findings.md` in the repository has the construction, the seeds
and the null arm in full.

## Why the links are the constraint, in floating point too

The spectral matrix must be positive definite at every frequency, in
every group, for every draw the optimizer tries. `det S > 0` needs only
`C < 1`, and `C` reaches the density through a logit link, so no value
of any linear predictor can produce an invalid matrix **in exact
arithmetic**. No eigenvalue floor and no projection step are applied,
and none is needed.

Exact arithmetic is not what runs. `plogis(eta)` is exactly 1 in double
precision from about `eta = 36.74` upward, so a density that formed
`1 - C` by subtraction would read `log(0)` and divide by zero there, and
would already be wrong in the third digit at `eta = 30`. That region is
reachable: a random effect on `coh` with one group near a true coherence
of 1 drives the estimate past `eta = 34`.

This family therefore never subtracts. Core keeps each dpar's linear
predictor beside it while the objective is taped, and
[`frmtmb::dpar_log1m()`](https://aforren1.github.io/frmtmb/reference/frmtmb-robust-dpars.html)
reads `log(1 - C)` off it exactly, so `1 / (1 - C)` is its exponential.
Both are exact to `eta = 709`, where the double range itself ends,
against `36.74` for the naive form. Measured against a high-precision
reference at `n = 16`: relative error 1.0e-03 at `eta = 30` and NaN at
`eta = 40` before the change, below 1e-15 at both after it.

One path keeps the plain round trip, and says so rather than being
floored: `residuals(type = "deviance")` runs off the tape, where core
does not store the linear predictor, so the accessor falls back to the
plain form. That is accurate until `C` rounds to exactly 1, which needs
`eta` past 36.74 and is further than a fit reaches in practice. Measured
on the degenerate case above, two signals differing by 1e-5 of noise:
the fit lands at `eta = 23.0`, `C` is 0.99999999989743915 rather than 1,
and all 127 response, Pearson and deviance residuals are finite. Past
36.74 they would be `NaN`, and no floor is applied to hide it.

## What it refuses

- More than two channels. This is a scope decision, and **core is not
  the obstacle**: `frm()` carries a matrix-valued response, and both a
  gaussian fit and a custom family indexing `y[, 1]` and `y[, 2]` work
  today. What is missing is here rather than there. A `p`-channel model
  needs `p^2` linear predictors written out by hand; the clean
  `power, power, coherence, phase` coordinates of this family do not
  survive past `p = 2`, so it would be written on raw Cholesky entries
  that mean nothing to a reader; and the log determinant and the trace
  stop being closed forms, which puts `adcomplex` back on the tape at
  about six times the cost. The density itself generalizes and was
  checked at `p = 4`. Fit pairs.

- `n < 2`. A single complex draw gives a rank-one matrix whose coherence
  is exactly 1 and whose log determinant is `-Inf`; the normalizing
  constant is `Inf` there as well, since `lgamma(n - 1)` reaches a
  non-positive argument. This is checked in `valid_y()`, before any tape
  is built.

- [`simulate()`](https://rdrr.io/r/stats/simulate.html). A draw from
  this family is a whole Hermitian matrix and the response slot carries
  one number of the four, so a simulated `w11` with the observed `w22`,
  `w12r` and `w12i` beside it is not a draw from anything.
  [`frm_cross_simulate()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_simulate.md)
  returns all four columns instead.

## See also

[`frm_cross_spectrum()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_cross_spectrum.md)
to build the response,
[`frm_coherence()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_coherence.md)
and
[`frm_phase()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/frm_coherence.md)
to read the fit

## Examples

``` r
set.seed(4)
src <- rnorm(4096)
xs <- frm_cross_spectrum(src + rnorm(4096), 0.9 * src + rnorm(4096),
                         segments = 16)
xs <- xs[xs$freq < 0.1, ]
fit <- frmtmb::frm(
  frmtmb::bf(w11 | vreal(w22, w12r, w12i) + vint(n) ~ 1,
             pow2 ~ 1, coh ~ 1, phase ~ 1),
  family = cross_wishart(), data = xs)
frm_coherence(fit)[1, ]
#>   .estimate       .se   .lower    .upper      .eta
#> 1 0.2031608 0.1568791 0.157874 0.2574666 -1.366655
```
