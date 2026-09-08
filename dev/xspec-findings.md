# Cross-spectra in frmtmb: feasibility probe

Date: 2026-09-07. Branch `wt-xspec`, frmtmb 0.53.0, RTMB 1.9, TMB 1.9.25,
R 4.6.1, Windows 11.

This lane owns item 5 of `dev/frequency-domain-todo.md`, the complex
cross-spectrum. Items 1 to 4, the univariate Whittle helpers, belong to a
sibling lane working on core. **No core file was read for writing and no
core file was changed.**

## Verdict

**Yes, and the package is built.** `extensions/frmtmb.coupling` 0.1.0
ships a complex Wishart family for a channel pair, a cross-spectrum
constructor, coherence and phase extractors with intervals, and a
hierarchical worked example. 337 expectations pass in 113 seconds.

The four phase 1 answers, in one line each:

* **(a)** Complex AD tapes and differentiates correctly, to 1e-10
  relative against `numDeriv`, at first and second order. It is not
  needed at two channels: the density there is real arithmetic in closed
  form, and the complex route is 5.8 times slower for the same answer.
* **(b)** The Wishart form needs `n >= p`; a raw cross-periodogram at
  one frequency has coherence exactly 1 and log determinant `-Inf`.
  Segments, sine tapers and frequency smoothing each buy their nominal
  count; **tapers and smoothing together do not** and are refused. The
  working minimum is `n >= 4` per unit with a grouping factor.
* **(c)** The complex Cholesky, which at `p = 2` is
  `log power, log power, logit coherence, phase`. Positive definiteness
  across the whole frequency range is free: the logit link IS the
  constraint.
* **(d)** Naive coherence at a true 0.2 with 4 segments averages 0.375
  over 4000 replicates; the model gives 0.209. Averaging per-subject
  coherences does not fix it, and concatenating subjects makes it worse
  (0.23 against a truth of 0.50, from phase cancellation). **But the
  bias goes as `(1 - C)^2 / n`, so at a true coherence of 0.9 the naive
  estimator is already unbiased (0.9007 at n = 4) and this package buys
  nothing.** The whole argument lives at low to moderate coherence. The
  original summary of this lane omitted that and the review was right to
  call it out.

The single most useful thing measured: **put a random effect on every
distributional parameter, not only on coherence.** The two channel
powers are incidental parameters, and left free they starve the
coherence variance component, which collapses to zero in 144 of 150
runs at 4 segments and drops interval coverage from 0.912 to 0.622. The
fit converges without a warning either way.

## 1. Does complex AD tape

RTMB 1.9 carries an S4 `adcomplex` class. What it actually provides,
read off the namespace rather than the news file:

| operation | present | how it is implemented |
| --- | --- | --- |
| `+ - * /`, `[`, `[<-`, `dim`, `t`, `rep` | yes | pairs of `advector` |
| `Re`, `Im`, `Mod`, `Arg`, `Conj` | yes | pairs of `advector` |
| `exp`, `log`, `sqrt` | yes | pairs of `advector` |
| `%*%` | yes | four real matrix products |
| `solve(a)`, `solve(a, b)` | yes | `solve_complex_atomic`, a real TMB atomic |
| `fft` | yes | `fft_complex` atomic |
| `det`, `determinant` | **no** | absent for `adcomplex` and for base complex |
| `chol`, `Cholesky` | **no** | absent for `adcomplex` |

`adcomplex(real, imag)` demands `advector` in **both** slots; passing a
plain numeric imaginary part fails `validObject` with a message about
slot `imag`. `solve_complex(x)` is not a two-argument solve, it is
`solve(x + 0i)`; the two-argument form is the S4 `solve(a, b)` method,
which inverts `a` and then multiplies.

### Gradient and Hessian against numDeriv

Every probe below is written **once** against a pluggable complex
constructor and then run twice, through base R `complex` and through
`adcomplex`, so the tape and the reference execute identical source.

| probe | expression | value agreement | grad max abs err | Hessian max abs err |
| --- | --- | --- | --- | --- |
| 1a scalar | `Mod(q)^2 + Re(q) - 2*Im(q) + Arg(z*Conj(w))`, `q = z*w + z/w` | exact | 2.38e-10 (3.7e-11 rel) | 4.56e-11 (5.3e-12 rel) |
| 1b 2x2 Hermitian | `2*sum(logdiag) + Re tr(S^-1 B) + Im((S^-1 B)[1,2])`, `S = L L^H` | exact | 3.36e-10 (3.2e-10 rel) | 8.39e-11 (1.4e-11 rel) |
| 1c 5x5 Hermitian | same shape, 25 parameters | 1.4e-14 | 2.63e-08 (3.6e-10 rel) | not taken |
| 1d `fft` | `sum(Mod(fft(z))^2)`, length 8 | 1.4e-14 | 5.71e-09 (1.6e-10 rel) | not taken |

Errors at this size are the finite-difference method's, not the tape's.
1d also reproduces Parseval exactly (121.287135939 both ways).

**Complex AD works, including second order, which is what Laplace needs.**

### The missing determinant is not a problem

There is no complex `det` on the tape, and base R refuses
`determinant()` on complex matrices as well, so a naive transcription of
the complex Wishart density does not run. The fix removes the need
rather than supplying the operation: parameterize the spectral matrix by
its **complex Cholesky factor** `S = L L^H`, with `L` lower triangular
and a real positive diagonal. Then

    log det S = 2 * sum(log(diag(L)))

which is a sum of parameters, exact and free. Checked against the
eigenvalues of the assembled matrix: -0.2 both ways, Hermitian residual
0.00e+00, minimum eigenvalue 0.393.

That single choice decides the parameterization for the rest of this
document, and it also disposes of the positive definiteness constraint
in section 3.

## 2. Is the complex Wishart log density right

The density used throughout, for `W` Hermitian positive definite and
`n` the segment count:

    log f(W) = (n - p) log det W - tr(S^-1 W) - n log det S - log CGamma_p(n)
    log CGamma_p(n) = p(p-1)/2 log(pi) + sum_{j=1..p} lgamma(n - j + 1)

`log det W` is data, so it is computed once outside the tape;
`log det S` is `2 * sum(theta[1:p])` by the Cholesky parameterization;
`tr(S^-1 W)` is the only term that needs complex arithmetic on the tape.

| check | result |
| --- | --- |
| 2a p = 1 must be a Gamma | `cw_lpdf` -3.7179638186 against `dgamma(3.7, shape = 5, scale = 2.3, log = TRUE)` -3.7179638186, difference 0 |
| 2b E[W] = n S, 20000 draws, n = 2 / 4 / 16 | max abs error 0.0112 / 0.0067 / 0.0011 |
| 2b E[log det W], same draws | -0.8125 / 1.4889 / 4.7300 against theory -0.8396 / 1.4937 / 4.7302 |
| 2c gradient, p = 2, n = 5, 4 parameters | max abs 1.42e-09, 1.9e-10 relative |
| 2c Hessian, same | 1.2e-10 relative |
| 2c gradient, p = 4, n = 7, 16 parameters | max abs 2.66e-09, 4.2e-10 relative |
| 2c Hessian, same | 2.2e-11 relative |
| 2e MLE from one W must be W/n, p = 3, n = 8 | max abs 1.68e-05, `nlminb` default tolerance |

**Yes.** The density tapes, and its first and second derivatives are
correct to finite-difference precision.

### The Wishart form and the segment form are the same model

Writing the likelihood on the `p x n` complex coefficient matrix `D`
instead of on `W = D D^H` gives

    log f(D) = -tr(S^-1 W) - n log det S - n p log(pi)

which differs from the Wishart form by a term free of the parameters.
Measured: the two scores agree to 9.1e-09 (p = 2) and 6.5e-09 (p = 4),
and the gap between the two log densities is 24.1052 at one `theta` and
24.1052 at another (p = 2), 54.2988 twice (p = 4). Same maximum, same
Hessian, same Wald intervals. Only the reported log likelihood differs.

This matters for section 3: below `n = p` the Wishart form does not
exist but the segment form still does.

## 3. What the model estimates, and its constraints

### Parameterization

`S(w) = L(w) L(w)^H`, `L` complex lower triangular with a real positive
diagonal held as its log. For `p x p` that is `p^2` free real numbers,
which is the dimension of the Hermitian positive definite cone, so the
map is a bijection with no redundancy and no constraint.

**Positive definiteness across the whole frequency range comes free.**
This is the point that makes a smooth-in-frequency spectral matrix
practical. Whatever a spline, a random effect, or a group contrast does
to `theta(w)`, `exp()` keeps the diagonal positive and `L L^H` is
positive definite at every `w`. No eigenvalue floor, no projection step,
no rejected proposal. A model written directly on the entries of `S`
would have to enforce positive definiteness at every frequency and every
level of every grouping factor, and would fail at the boundary; this one
cannot leave the cone.

### What comes out of it

For `p = 2`, write `L = [[a, 0], [c, b]]` with `a, b > 0` real and `c`
complex. Then

    S11 = a^2                     channel 1 power
    S22 = |c|^2 + b^2             channel 2 power
    S21 = a c                     the cross-spectrum
    coherence = |c|^2 / (|c|^2 + b^2) = plogis(2 * eta),  eta = log(|c| / b)
    phase     = Arg(S21) = Arg(c)

Checked over 500 random draws: max absolute error 5.55e-16 on all three
identities. Checked over 200 random Hermitian matrices: swapping the two
channels changes `L` but leaves the derived coherence unchanged to
4.44e-16, and it agrees with `|S12|^2 / (S11 S22)` to the same precision.

So the working coordinates are **power, power, coherence logit, phase**,
and every one of them is unbounded and modelable:

| quantity | coordinate | link back |
| --- | --- | --- |
| channel powers | `log a`, `log b` | `exp`, squared |
| coherence | `eta` | `plogis(2 * eta)`, so a Wald interval cannot leave (0, 1) |
| phase | `phi` | identity, wrapped |

This is what makes the extension worth writing. A random effect on `eta`
is a hierarchical model for coherence; a random effect on `phi` is a
hierarchical model for lag; a smooth in `eta` over frequency is a
coherence spectrum; and none of them can produce an invalid spectral
matrix.

## 4. Does the degrees-of-freedom parameter behave

### Below n = p the Wishart density does not exist

Measured at `p = 2`, one frequency, `W = D D^H`:

| n | eigenvalues of W | log det W | log CGamma_2(n) | naive coherence |
| --- | --- | --- | --- | --- |
| 1 | 4.175, -5.6e-16 | -Inf | Inf | **1.000000** |
| 2 | 4.935, 0.432 | 0.7572 | 1.14473 | 0.431683 |
| 3 | 8.553, 0.321 | 1.0100 | 1.83788 | 0.840693 |
| 4 | 3.564, 0.430 | 0.4257 | 3.62964 | 0.583161 |
| 5 | 6.038, 1.033 | 1.8300 | 6.11454 | 0.500187 |

A raw cross-periodogram at one frequency is rank one, so its coherence is
exactly 1 by arithmetic, its log determinant is `-Inf`, and
`log CGamma_p(n)` is `Inf` because `lgamma(n - p + 1)` reaches a
non-positive argument. **`n >= p` is a hard requirement of the Wishart
form.** The segment form of section 2 stays finite for any `n >= 1`, but
below `n = p` it cannot identify `S` either, because its maximizer is
singular.

### Smoothing, tapers or trials: all three work, and n_eff is the nominal n

At true coherence 0 the naive estimate has mean exactly `1/n`, so
`n_eff = 1 / mean(coherence)` reads the real degrees of freedom straight
off a simulation.

**Corrected after review.** The first version of this table read ONE
frequency bin per replicate and discarded the rest, which left it far
too noisy to support the conclusion drawn from it. Re-run with every
retained bin of 250 replicates pooled, series length 4096:

| route | nominal 2 | 4 | 8 | 16 |
| --- | --- | --- | --- | --- |
| disjoint segments, white | 1.998 | 4.002 | 7.998 | 15.829 |
| frequency smoothing, white | 2.002 | 3.992 | 8.011 | 15.965 |
| sine multitaper, white | 2.004 | 4.012 | 7.958 | 16.021 |
| disjoint segments, AR(1) phi = 0.9 | 1.997 | 3.989 | 8.000 | 16.070 |
| frequency smoothing, AR(1) phi = 0.9 | 1.994 | 3.993 | 8.008 | 16.069 |
| sine multitaper, AR(1) phi = 0.9 | 2.000 | 3.992 | 7.971 | 16.197 |

All three routes deliver what they promise, on a white spectrum and on a
steeply sloped one alike, to about 1 percent at nominal 16.

The earlier table reported frequency smoothing on AR(1) at nominal 16 as
15.13 and called it a 5 percent local-flatness cost, and that number
reached the shipped manual. **It was noise, and the conclusion drawn
from it was wrong.** The reviewer's independent measurement agrees with
the corrected numbers above. The cost is real but small and arrives
later: at `smooth = 32` on the same spectrum the measurement is 31.5,
about 1.6 percent. What smoothing genuinely cannot survive is a band
across which the COHERENCE varies, which is a modelling error rather
than a degrees-of-freedom one and which no count detects.

### How many

The Wishart form needs `n >= p`. That is the arithmetic floor and it is
nowhere near enough in practice. What sets the practical floor is not the
density but the variance component, measured in section 5: with random
effects on every parameter, `n = 4` per subject over 40 subjects reaches
a bias of +0.008 in coherence and 0.91 Wald coverage. At `n = 2` a
per-subject coherence has a standard deviation near 0.29 and the model
has almost nothing to shrink toward.

**Recommendation: refuse `n < p`, warn below `n = p + 2`, and document
`n >= 4` per unit with a grouping factor as the working minimum.**

## 5. The payoff: coherence bias, naive against the model

### Homogeneous subjects

20 subjects, identical true spectral matrix, 4000 replicates. "Naive
mean" is the mean over subjects of each subject's own coherence, which is
what a two-stage pipeline computes. "Model pooled" is the complete
pooling complex Wishart estimate, `sum_i W_i / (N n)`, verified against
the taped optimizer to 4.0e-06 in the matrix and 7 digits in the
coherence.

| true C | n | naive, one subject | naive mean of 20 | model pooled 20 | ratio |
| --- | --- | --- | --- | --- | --- |
| 0.2 | 2 | +0.3737 | +0.3714 | +0.0199 | 0.054 |
| 0.2 | 4 | +0.1762 | +0.1749 | +0.0085 | 0.049 |
| 0.2 | 8 | +0.0847 | +0.0837 | +0.0037 | 0.045 |
| 0.2 | 16 | +0.0394 | +0.0405 | +0.0016 | 0.041 |
| 0.2 | 32 | +0.0196 | +0.0200 | +0.0009 | 0.044 |
| 0.5 | 2 | +0.2034 | +0.1944 | +0.0081 | 0.042 |
| 0.5 | 4 | +0.0830 | +0.0804 | +0.0038 | 0.047 |
| 0.5 | 8 | +0.0379 | +0.0356 | +0.0018 | 0.051 |
| 0.5 | 16 | +0.0199 | +0.0161 | +0.0004 | 0.024 |
| 0.5 | 32 | +0.0089 | +0.0091 | +0.0013 | 0.147 |
| 0.8 | 2 | +0.0511 | +0.0506 | +0.0009 | 0.018 |
| 0.8 | 4 | +0.0184 | +0.0158 | +0.0005 | 0.033 |
| 0.8 | 8 | +0.0086 | +0.0068 | +0.0008 | 0.123 |
| 0.8 | 16 | +0.0036 | +0.0028 | +0.0001 | 0.043 |
| 0.8 | 32 | +0.0012 | +0.0013 | +0.0000 | 0.032 |

The naive bias goes as `1/n` and **does not improve with more subjects**,
because averaging 20 biased numbers returns the same bias. The model bias
goes as `1/(N n)`. At a true coherence of 0.2 with 4 segments, a routine
setting, the naive answer averages 0.375, nearly double the truth; the
model gives 0.209.

**Where this buys nothing.** Read the C = 0.8 rows: the naive bias at
n = 4 is +0.018 and at n = 32 it is +0.001. The reviewer measured C = 0.9
and got +0.0007 for one subject at n = 4, against +0.0003 pooled. The
bias is `(1 - C)^2 / n`, so it disappears as the coupling gets strong,
and a package sold on removing it should say plainly that at high
coherence there is nothing to remove. What survives at any coherence is
the rest: an interval with the right sampling distribution, partial
pooling, and a group comparison on a scale that cannot leave (0, 1).

An independent check of the same arithmetic, 40000 replicates at true
C = 0: the naive mean is 0.5021, 0.3301, 0.2481, 0.1252, 0.0625, 0.0315
at n = 2, 3, 4, 8, 16, 32, against the exact `1/n` of 0.5, 0.3333, 0.25,
0.125, 0.0625, 0.03125.

### Heterogeneous subjects, which is the real case

Pooling is legitimate only when every subject has the same spectral
matrix. Once coherence and phase vary across subjects, pooling is worse
than naive averaging. 40 subjects, coherence logit `eta ~ N(0, 0.4^2)`,
phase `phi ~ N(0.8, sd_phi^2)`, per-subject powers `~ N(0, 0.3^2)`,
150 replicates, true median coherence 0.500:

| sd_phi | n | naive mean | pooled | hierarchical | sd_eta collapsed to 0 |
| --- | --- | --- | --- | --- | --- |
| 0.0 | 4 | 0.5832 (+0.083) | 0.4139 (-0.086) | 0.5559 (+0.056) | 149/150 |
| 0.0 | 16 | 0.5176 (+0.018) | 0.4198 (-0.080) | 0.5058 (+0.006) | 8/150 |
| 0.8 | 4 | 0.5853 (+0.085) | **0.2326 (-0.267)** | 0.5666 (+0.067) | 144/150 |
| 0.8 | 16 | 0.5168 (+0.017) | **0.2324 (-0.268)** | 0.5066 (+0.007) | 12/150 |

Pooling across subjects with a spread of phases **destroys the
coherence**: 0.23 against a truth of 0.50, and the error does not shrink
with more segments, because it is not a variance problem but cancellation
of the cross terms. That is the argument for doing this hierarchically
instead of concatenating everyone's segments, which is the other thing
pipelines do.

### Why the hierarchical fit is still biased at n = 4, and the fix

The `+0.067` above is not the likelihood. It is the variance component
collapsing to zero in 144 of 150 runs, which turns the model back into
complete pooling on `eta`. Three variants, 40 subjects, sd_phi = 0.8,
150 replicates, true coherence 0.500:

| variant | n | coherence | bias | mean sd_eta | Wald 95% coverage |
| --- | --- | --- | --- | --- | --- |
| powers free per subject, sds free | 4 | 0.5683 | +0.0683 | 0.000 | **0.622** |
| sds fixed at the truth | 4 | 0.5042 | +0.0042 | 0.400 | 0.933 |
| powers also random effects, sds free | 4 | 0.5084 | +0.0084 | 0.291 | 0.912 |
| powers free per subject, sds free | 16 | 0.5125 | +0.0125 | 0.282 | 0.820 |
| sds fixed at the truth | 16 | 0.5036 | +0.0036 | 0.400 | 0.933 |
| powers also random effects, sds free | 16 | 0.5075 | +0.0075 | 0.366 | 0.893 |

**The per-subject powers are incidental parameters.** Two free numbers
per subject that gain no information as subjects are added starve the
variance component, and the coherence bias and the broken 0.622 coverage
follow from that, not from the complex Wishart. Give the powers a random
effect too and the bias falls by a factor of eight at n = 4 while
coverage rises from 0.622 to 0.912.

That is a design instruction for the package, and a natural one in this
grammar: put `(1 | id)` on **every** distributional parameter, not only
on the coherence.

## 6. Verdict on phase 1

**Yes, on all four questions, and the package is built.**

(a) Complex AD tapes, first and second order, to finite-difference
precision. It is not needed at two channels, which is the case that
matters; it is what let the general-`p` form be derived and checked, and
it is what a `p >= 3` version would use.

(b) The degrees of freedom behave. `n >= p` is a hard floor for the
Wishart form and `n >= 4` is the working one. Segments, sine tapers and
frequency smoothing each deliver their nominal count; tapers and
smoothing together do not, and that combination is refused.

(c) The parameterization is the complex Cholesky, which at two channels
is exactly `log power, log power, logit coherence, phase`. Positive
definiteness over the whole frequency range costs nothing, because the
logit link **is** the constraint.

(d) The payoff is large and it is not where it first appears. The naive
per-unit bias is real and goes as `1/n`, but the model does not remove
it for a single unit: with one unit the maximum likelihood estimate is
the naive estimate. What the model buys is that pooling is done **inside
a likelihood that knows the sampling distribution**, so partial pooling
across units and frequencies is available and the group answer is nearly
unbiased where the average of per-unit answers is not.

## 7. What was built

`extensions/frmtmb.coupling`, version 0.1.0, following
`extensions/frmtmb.spline` in layout, registration and test structure.
**Core lists it nowhere**: no `Suggests`, no `Imports`, no core file
mentions it. The three repository files touched outside the extension
are `README.md` (one table row and one install line), `dev/build-docs.R`
(one `PKGS` entry and a stale count in a comment), and `_pkgdown.yml`
(one Extensions menu entry), plus a new workflow
`.github/workflows/check-frmtmb-coupling.yaml`.

| export | what it is |
| --- | --- |
| `frm_cross_spectrum()` | signals to rows: drops bin zero and Nyquist, counts the degrees of freedom, accepts a matrix pair for many units |
| `cross_wishart()` | the family: four dpars, `log`, `log`, `logit`, `identity` |
| `frm_coherence()`, `frm_phase()` | extractors with intervals, built on `predict(se.fit = TRUE)` |
| `frm_cross_simulate()` | whole Hermitian draws, since `simulate()` refuses |

### At two channels the complex Wishart is real arithmetic

This is the finding that decided the implementation. In the coordinates
of section 3, with `C` the coherence and `phi` the phase,

    det S       = S11 S22 (1 - C)
    tr(S^-1 W)  = (S22 w11 + S11 w22
                   - 2 sqrt(C S11 S22) (cos(phi) w12r + sin(phi) w12i))
                  / (S11 S22 (1 - C))

No determinant, no complex number, no matrix inverse. Both identities
are checked against complex linear algebra over 50 random matrices in
`test-cross-wishart.R`, and the whole log likelihood is checked against
an independent evaluation that does use a complex inverse and an
eigenvalue log determinant.

Measured cost of the alternative: the same hierarchical model written
through `solve()` on `adcomplex` reaches the same objective
(110.596902 both ways) and the same standard error (0.0935 both ways)
in 1.63 s against 0.28 s, **5.8 times slower**. So RTMB's complex AD is
the reason this lane could be certain the general form is right, and is
not in the shipped code path.

### The response is one number of four

The frame carries a scalar response, so a Hermitian matrix reaches the
density as `w11` plus three `vreal()` columns and one `vint()`. That
works and is declared through `required_aterms`, without which an absent
column reaches the density as `NULL` and the fit returns a
log-likelihood of zero.

**This lane recorded a core seam here that does not exist, and the
review caught it.** The claim was that the frame carries no Hermitian
response. It carries a matrix response perfectly well:
`R/frame.R:241-247` preserves one explicitly, and the reviewer fitted
both a gaussian model with a 60 by 2 matrix response and a custom family
indexing `y[, 1]` and `y[, 2]`. **Core is not the obstacle and no lane
should be sent after this.**

The real reasons to stop at two channels are this package's own: `p`
channels need `p^2` linear predictors written out by hand; the
`power, power, coherence, phase` coordinates of section 3 do not survive
past `p = 2`, so the model would be written on raw Cholesky entries that
mean nothing to a reader; and the log determinant and the trace stop
being closed forms, which puts `adcomplex` back on the tape at about six
times the cost. The arithmetic is fine: probe 2 tapes `p = 4` with
gradients correct to 4.2e-10 relative.

**The seam that IS real is a different one.** Core stores each dpar's
linear predictor beside it in the `dpars` list as `.eta_<dpar>`, and the
link registry carries `logit_eta` and `log_eta` fields for exactly the
saturation problem section 9 describes. It works, core's own densities
use it through internal accessors, and this package now depends on it.
It is documented nowhere for extension authors, no sibling extension
uses it, and `robust_logit()` is `@noRd`. A promise would be worth
making.

### Validation

| claim | how it is checked | result |
| --- | --- | --- |
| the density is right | reduces to `dgamma` at p = 1 | difference 0 |
| the estimator is right | with intercepts only the estimate must be `sum(W)/(N n)` on all four coordinates | 1.4e-12 relative through `frm()`, at four different truths and `n` from 2 to 32 |
| the log likelihood is right | recomputed with a complex inverse and eigenvalue log det | agrees to 1e-6 |
| weights | unit weights change nothing, doubled weights double the log likelihood | exact |
| the simulator | first moment over 4000 draws against the fitted matrix | within 5 percent |
| the bias claim | naive against pooled at true coherence 0, 600 replicates | naive 1/n, pooled below naive/10 |
| interval coverage | see section 5 | 0.912 at n = 4 with all dpars random |

Suite: **337 expectations, 0 failures, 0 errors, 113 s** in one process
with `NOT_CRAN=true`.

### Two things found while building, both worth more than the code

**1. `tapers = 1` must mean no taper.** The first version applied the
first sine taper even at `tapers = 1`. A tapered transform has
correlated neighbouring bins, so `segments = 2, smooth = 4` bought 6.1
independent draws where it claimed 8. The bug was invisible except
through the effective degrees of freedom measurement, and every standard
error downstream would have been about 15 percent too small. Fixed by
making `tapers = 1` a boxcar; every accepted configuration now measures
its nominal count:

| configuration, nominal n = 8 | measured |
| --- | --- |
| `segments = 8` | 7.92 |
| `segments = 4, smooth = 2` | 8.17 |
| `segments = 2, smooth = 4` | 7.99 |
| `segments = 1, smooth = 8` | 7.98 |
| `segments = 4, tapers = 2` | 8.03 |
| `segments = 2, tapers = 4` | 7.76 |
| `segments = 1, tapers = 8` | 7.98 |

The same measurement is what condemned tapers and smoothing **together**:
4 tapers over 4 bins buys 5.80 against a nominal 16, and 2 over 2 buys
2.70 against 4. That combination is now refused rather than counted
wrong.

**2. A power formula too simple for the data moves the coherence.** This
is the sharpest edge in the package and it converges without complaint.
The likelihood couples the powers and the coherence, so misspecifying
`mu` and `pow2` reweights the rows each coherence parameter pools.
Measured on a pair driven by an AR(1) source whose power varies 236-fold
across frequency, with `coh ~ band`:

| power formula | low band | high band | AIC |
| --- | --- | --- | --- |
| `mu ~ 1, pow2 ~ 1` | **0.485** | **0.468** | 8636 |
| `mu ~ band, pow2 ~ band` | 0.641 | 0.050 | 5599 |
| `mu ~ s(freq), pow2 ~ s(freq)` | 0.350 | 0.050 | 4097 |

**The automatic check this lane shipped for it was wrong, and the review
caught that too.** A `post$fit_check` hook compared the spread of
`log(W11 / (n S11))` against `trigamma(n)` and warned above a ratio of
2. That measures residual spread in the POWER model, which is not
monotone in the coherence displacement the user is being warned about.
The reviewer's measurement, two bands at `n = 8`, true coherences 0.35
and 0.05:

* at a **2.7-fold** power spread the constant-power model reported a
  high-band coherence of 0.092 against a truth of 0.05, nearly double,
  and the hook stayed **silent**;
* on the **correct** band-varying model, which recovered 0.368 and
  0.042 accurately, the hook **fired** at a ratio of 20.5.

So it was silent where it mattered and loud where it did not. The
threshold of 2 had been tuned to the single 236-fold example in its own
docstring. **The hook is withdrawn rather than retuned**, because a
check that punishes correct models teaches users to ignore it.

What replaces it is a comparison that was measured before being
recommended: fit the powers with and without a frequency term and
compare AIC. Two bands over 400 frequencies at `n = 8`, true coherences
0.35 and 0.05, varying the power spread:

| power spread | flat powers report | correct powers report | delta AIC |
| --- | --- | --- | --- |
| 1.6x | 0.295, 0.083 | 0.330, 0.068 | 80 |
| 2.7x | 0.268, 0.078 | 0.342, 0.041 | 487 |
| 4.4x | 0.256, 0.136 | 0.344, 0.054 | 954 |
| 7.3x | 0.240, 0.206 | 0.342, 0.047 | 1849 |
| 19.6x | 0.229, 0.458 | 0.341, 0.052 | 3716 |
| 385x | 0.222, 0.875 | 0.342, 0.059 | 11518 |

AIC rejects the flat model by 80 or more at **every** spread, including
the two where the hook was silent, and by construction it never punishes
the correct model. That is what the manual and the vignette now say to
do.

### The hierarchical model, through the shipped package

40 units, coherence logit `~ N(0, 0.4^2)`, phase `~ N(0.8, 0.8^2)`,
per-unit powers `~ N(0, 0.3^2)`, true coherence 0.500, random effects on
all four dpars, 60 replicates each:

| n per unit | clean fits | warned | failed | estimate | bias | variance collapsed |
| --- | --- | --- | --- | --- | --- | --- |
| 4 | 58 | 2 | 0 | 0.5102 | +0.0102 | 0/60 |
| 8 | 60 | 0 | 0 | 0.5057 | +0.0057 | 0/60 |
| 16 | 60 | 0 | 0 | 0.5099 | +0.0099 | 0/60 |

Against the naive mean of per-unit coherences at the same settings,
+0.085 at n = 4 and +0.017 at n = 16, and against pooling, -0.267 at
both. **The variance component never collapsed once every dpar carried a
random effect**, where it collapsed in 144 of 150 runs with the powers
left free.

## 8. What was refused, and why

* **More than two channels.** The arithmetic is checked at `p = 4`.
  Core is NOT the obstacle: a matrix response fits, and so does a
  custom family indexing its columns. This is a scope decision; see section 7.
* **`simulate()`.** A draw is a matrix and the response slot holds one
  of its four numbers, so a simulated `w11` beside the observed rest is
  not a draw from anything and is usually not even positive definite.
  Refused through `sim_refusal`, with `frm_cross_simulate()` supplied in
  its place rather than leaving a hole.
* **`cens()` and `trunc()`.** No `lcdf`, no `lccdf`, and
  `accepts_aterms` does not list them, so they are refused by name at
  frame assembly. A censored cross-spectral matrix is not a thing anyone
  has.
* **Overlapping segments.** Not implemented at all rather than
  implemented with a wrong `n`.
* **Tapers and smoothing together.** Refused on the measurement above.
* **`n < 2`, at two places.** `frm_cross_spectrum()` refuses to build
  such rows and `valid_y()` refuses to fit them.
* **Touching core.** No core file was modified. The univariate Whittle
  helpers, items 1 to 4 of the specification, were left to the sibling
  lane; `frm_cross_spectrum()` duplicates none of them, because it
  returns a Hermitian matrix per row rather than a periodogram ordinate.

## 9. The numerical failure the review found, and the fix

The mathematics of section 3 is exact: `det S > 0` needs only `C < 1`,
and a logit link cannot produce `C = 1`. The arithmetic was not.
`plogis(eta)` is **exactly 1** in double precision from
`eta = 36.736800569677101` upward, so a density that forms `1 - C` by
subtraction reads `log(0)` and divides by zero there. This lane shipped
that subtraction and a code comment blessing it, and three documentation
claims asserting the parameterization could not fail.

Measured by the reviewer, one legal row at `n = 16`, shipped expression
against the stable one:

| eta | shipped | stable | relative error |
| --- | --- | --- | --- |
| 10 | -4299.993882 | -4299.993882 | 9.9e-13 |
| 25 | -14400919385 | -14400979522 | 4.2e-06 |
| 30 | -2139477257233 | -2137294915880 | **1.0e-03** |
| 36.74 | **NaN** | -1.807e+15 | undefined |
| 40 | **NaN** | -4.708e+16 | undefined |

And it is reachable. A random effect on `coh` with one group at a true
coherence of `1 - 1e-15` drove the estimate to `eta = 34.031` and
returned **NaN standard errors**, with `valid_y()` passing every row
because each row's own determinant was positive. Without any random
effect, two signals differing by 1e-5 of noise gave `frm_coherence()` a
**95 percent interval of width zero**, lower equal to estimate equal to
upper.

### The fix, and why it needed no floor

Core had already solved this and this lane had not read far enough.
`build_objective()` stores each dpar's linear predictor beside it as
`.eta_<dpar>`, and `R/links.R` carries `logit_eta` and `log_eta` fields
on exactly the links whose round trip saturates, with a comment
explaining the same failure in the same words. The logit's `logit_eta`
is the identity, so the log-odds are available exactly on the tape:

    log(1 - C) = -logspace_add(0, eta)
    1 / (1 - C) = exp(logspace_add(0, eta))

Both are exact to `eta = 709`, where the double range itself ends,
against 36.74 for the naive form. **No floor and no projection step was
added**, so the "no eigenvalue floor" property the package sells
survives intact; the complement is simply never formed by subtraction.

Measured after the fix, the reviewer's own reproductions:

| case | before | after |
| --- | --- | --- |
| `cw_lpdf` at eta = 40 | NaN | -3.77e+15, finite |
| `cw_lpdf` at eta = 100 | NaN | -4.30e+41, finite |
| `cw_lpdf` at eta = 700 | NaN | -1.62e+302, finite |
| nine groups, top group at eta = 16 | NaN standard errors | logLik -167.8234, interval (0.982216, 0.995904) |
| two signals differing by 1e-5 | interval of width 0 | refused by name |
| relative error against a stable reference, eta 10 to 700 | up to 4.5e-02 | below 1e-12 |

The zero-width interval is closed separately, at the point the number is
consumed: `cp_link_se()` now refuses a standard error that is zero or
not finite instead of pushing it through `plogis()` and returning a
point. That refusal is what the two-signal case now hits.

One path keeps the plain round trip and says so. `residuals(type =
"deviance")` runs off the tape, where core does not store the linear
predictor by design, so it uses `log1p(-C)`. That holds until `C` rounds
to exactly 1, which needs `eta` past 36.74 and is further than a fit
reaches: the deliberately degenerate two-signal fit lands at
`eta = 23.0` with `C = 0.99999999989743915` and all 127 response,
Pearson and deviance residuals finite.

## 10. Sharp edges a user will meet

1. **The power model, above.** The one that silently lies, and the one
   this lane's automatic check got wrong. Compare AIC with and without a
   frequency term on the powers; nothing warns for you.
2. **A coherence at 1 to machine precision costs the interval, not the
   fit.** The density is exact there now, but the Hessian is not
   informative, so `frm_coherence()` refuses by name rather than
   returning a point. More segments per unit, or a prior on `coh`, is
   the remedy.
3. **Phase wraps and nothing detects it.** `frm_phase()` returns radians
   and does not wrap. A group whose phase sits near the wrap point gets
   a wide interval for a parameterization reason rather than a data
   reason. Re-center with an offset.
4. **A small residual gradient on smooth fits.** Some `s(freq)` fits
   return with a maximum absolute gradient near 1e-3 and core's
   convergence warning. It is seed dependent, it appears on
   `Gamma(link = "log")` and `exponential()` fits of the same response
   as well, and it is core's absolute threshold against a large
   objective rather than anything in this family.
5. **`n = 4` is the floor, not a target.** Two of 60 hierarchical fits
   warned at `n = 4` and none did at 8 or 16.
6. **The pooled estimate is power weighted.** `sum(W) / (N n)` gives the
   coherence of the summed matrix, which is dominated by the
   highest-power rows. It is a different functional from the average of
   per-row coherences, and the difference is large on a sloped spectrum
   (0.641 against 0.390 on the AR(1) example). Neither is wrong; they
   answer different questions, and a model with a frequency-resolved
   power formula answers the one usually meant.

## 11. Verification

All measured on R 4.6.1, Windows 11, with the private library this lane
built. `R CMD check` was run on the tarball, `--as-cran`, with the
manual, pandoc and TinyTeX both on `PATH`. **These are the post-review
numbers**; the counts this lane reported before the review were
understated, which the reviewer caught.

| file | expectations, `NOT_CRAN=true` | gate on |
| --- | --- | --- |
| `test-coherence.R` | 58 | 51 pass, 3 skip |
| `test-cross-spectrum.R` | 61 | 56 pass, 1 skip |
| `test-cross-wishart.R` | 203 | 200 pass, 1 skip |
| `test-message-uniqueness.R` | 4 | 4 |
| `test-surface.R` | 58 | 58 |
| **total** | **384 in 63 s** | **369 pass, 5 skip, 16 s** |

46 `test_that` blocks.

| step | result |
| --- | --- |
| suite, one process, `NOT_CRAN=true` | 384 expectations, 0 failures, 0 errors |
| suite, gated tier off | 369 expectations, 5 skips |
| tests inside `R CMD check` | `FAIL 0 \| WARN 4 \| SKIP 1 \| PASS 380`, 78 s |
| roxygen | idempotent: a second run writes nothing |
| examples | OK |
| vignette rebuild | OK |
| PDF manual | OK |
| HTML manual | OK |
| **`R CMD check --as-cran`** | **1 WARNING, 0 NOTEs** |

The one in-check skip is `test-message-uniqueness.R`, which skips itself
when the sources are not there to scan. The 4 warnings are the two sharp
edges section 10 records, both inside tests that tolerate them: three
from the `n = 4` hierarchical sweep, where a minority of replicates fail
to converge, and one small residual gradient on a smooth fit.

### The one WARNING, and its three measured causes

It is the CRAN incoming feasibility check. All three parts are
structural rather than defects in this package, and the reviewer
confirmed the third independently.

1. **"New submission".** The first check of any package not already on
   CRAN.
2. **"Strong dependencies not in the CRAN or BioC software
   repositories: frmtmb".** `frmtmb` is in no repository. Every
   extension in this monorepo raises this, and it cannot clear until
   core is released. This lane's first report named only two causes and
   left this one out.
3. **A 404 on `https://aforren1.github.io/frmtmb/frmtmb.coupling`.**
   The subsite is unpublished, not mis-linked. The reviewer fetched the
   siblings live: `frmtmb.learn` and `frmtmb.spline` return 200 with and
   without a trailing slash, and this one returns 404 both ways. The
   URL form matches five of the six siblings and clears the first time
   `dev/build-docs.R` runs with the `PKGS` entry this lane added.

### Two NOTEs cleared earlier, and one that came back and went again

* **"Files 'README.md' or 'NEWS.md' cannot be checked without
  'pandoc'".** Environment, not package: `RSTUDIO_PANDOC` is read by
  rmarkdown and not by `R CMD check`, which wants pandoc on `PATH`.
* **"Namespace in Imports field not imported from: 'RTMB'".** Real when
  the density was pure base arithmetic, and `RTMB` was dropped. It is
  back in `Imports` now, and genuinely used: `logspace_add()` is what
  makes the coherence complement exact (section 9). The NOTE does not
  return, because the import is real.

## 12. Counts

* 1 new package, 5 R files, 5 test files, 1 vignette, 5 help topics.
* 3 repository files modified outside the extension: `README.md`,
  `_pkgdown.yml`, `dev/build-docs.R`. 1 new workflow. 1 findings
  document.
* **0 core files changed.** Core lists the extension in no `Imports` and
  no `Suggests`.
* 23 compatibility rows registered, of which 6 are `untested` and say so.
* Probes run and discarded: complex AD taping, the complex Wishart
  density and its gradients, the degrees-of-freedom sweep, the coherence
  bias Monte Carlo, the Laplace test through the complex atomic, and the
  hierarchical convergence sweep. None left in the worktree; every
  number they produced is in this document.

### What the punch round changed

| item | severity | outcome |
| --- | --- | --- |
| S1, `1 - ch` catastrophic cancellation | blocker | fixed on the log scale through core's `.eta_<dpar>` and `logspace_add`; 3 false doc claims corrected; `cp_link_se()` now refuses a degenerate standard error; 6 new tests above `eta = 3` |
| S2, the p > 2 refusal blamed core falsely | blocker | refusal restated in this package's own terms; the false core seam withdrawn from section 7 and replaced by the real one, `.eta_<dpar>` being undocumented |
| S3, `fit_check` fired on the wrong statistic | blocker | hook **withdrawn**; an AIC comparison measured across power spreads from 1.6x to 385x replaces it in the manual, the vignette and a test |
| S4, an unreproducible smoothing number | minor | table re-run with every bin pooled; the 5 percent AR(1) loss was noise and is gone from section 4 and the manual |
| S5, vignette prose covered one of two warnings | minor | prose rewritten; the surviving warning named and explained |
| S6, stale workflow counts | trivial | 384 / 369 / 380, audited per file |
| S7, the WARNING has three causes | trivial | all three now named in section 11 |
| framing, high coherence | not a blocker | section 5 and the vignette now say the bias vanishes as coherence rises and the package buys nothing at 0.9 |

Tests went from 341 expectations to 384, and the 43 new ones are almost
all about the region the review found: the density and its slope at
`eta` of 10, 20, 30, 36.74 and 40; the complement against an independent
evaluation to `eta = 700`; the interval invariant that a 95 percent
interval is never returned with zero width; a group at a coherence of 1
not poisoning a fit; and the AIC comparison that replaced the hook.

### What the review confirmed, and what that is worth

The likelihood was verified harder than this lane had verified it: 1.6e-12
against an independent complex-matrix implementation over 200 rows, and
an importance-sampling mass test that varies `n` between proposal and
target, which tests the normalizing constant's degrees-of-freedom
dependence rather than letting it cancel. The bias tables, the
degrees-of-freedom tables, the `tapers = 1` fix and its 15 percent, and
the refusals all reproduced. **The mathematics held; what failed was
arithmetic, one false claim about someone else's code, and a diagnostic
tuned to a single example.** Those are the three things a lane cannot
check by rederiving its own work, which is the argument for the review.
