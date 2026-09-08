# frmtmb.coupling: Coherence and Phase Between Signals for 'frmtmb' Models

Hierarchical models for the coupling between two recorded signals. A
spectral model of ONE series needs no new likelihood, because a
periodogram ordinate is exponential about the spectral density; two
series is different, since the cross-spectrum is complex and the pair's
periodogram at each frequency is complex Wishart. This package supplies
that likelihood, parameterized so that coherence and phase are
distributional parameters with their own linear predictors, which is
what makes a random effect on coherence or a smooth in coherence over
frequency an ordinary 'frmtmb' formula. The parameterization cannot
leave the positive definite cone, because a logit link on coherence IS
the constraint, and the complement it needs is computed on the log scale
from the linear predictor so that the guarantee survives in floating
point rather than only on paper. The reason to want it: coherence
estimated from few segments is biased upward by about the reciprocal of
the segment count, the bias does not shrink when subject coherences are
averaged, and concatenating subjects instead cancels the cross terms
whenever phase varies between them. Both failures are measured in the
vignette against a known truth.

## What this package is for

One recorded signal needs no new likelihood. A periodogram ordinate of a
stationary series is exponential about the spectral density, so a
spectral model of one channel is an
[`frm()`](https://aforren1.github.io/frmtmb/reference/frm.html) fit with
an
[`exponential()`](https://aforren1.github.io/frmtmb/reference/frmtmb-families.html)
family and a log link, and `dev/frequency-domain-todo.md` in the
repository measures that it already works.

Two signals is different. The cross-spectrum is complex and the pair's
periodogram at each frequency is complex Wishart, so coherence and phase
need a likelihood that does not exist elsewhere in the package. That
likelihood is
[`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md),
and everything else here exists to feed it or to read it.

## What it buys, measured

Coherence estimated from `n` segments is biased upward by roughly
`(1 - C)^2 / n`, and the two things analysis pipelines do about it both
fail:

- **Averaging per-subject coherences does not help.** The mean of 20
  biased numbers carries the same bias. At a true coherence of 0.2 with
  4 segments the naive per-subject estimate averages 0.375 over 4000
  replicates, against a model estimate of 0.209.

- **Concatenating everyone's segments is worse.** Pooling is valid only
  when every subject has the same spectral matrix. With a phase spread
  of 0.8 radians between subjects the pooled estimate reads 0.23 against
  a truth of 0.50, and the error does not shrink with more segments,
  because it is cancellation of the cross terms rather than variance.

A hierarchical fit is what is left, and it reaches a bias of +0.008 and
0.91 interval coverage at 4 segments over 40 subjects. The full tables
are in `dev/xspec-findings.md`.

## The one thing to get right

Put a random effect on **every** distributional parameter, not only on
`coh`. The two channel powers are incidental parameters: two free
numbers per subject that gain no information as subjects are added. Left
free, they starve the coherence variance component, which collapses to
zero in 144 of 150 replicates at 4 segments, and interval coverage falls
from 0.912 to 0.622. This is the failure that looks most like success,
because the fit converges without a warning.

## What core supplies and what it does not

Everything read off a fitted object here is a documented seam:
[`stats::predict()`](https://rdrr.io/r/stats/predict.html) with
`se.fit = TRUE` and `dpar =`, which is what both extractors are built
out of, and [`stats::family()`](https://rdrr.io/r/stats/family.html) to
refuse a fit of the wrong family.

One thing this package uses is real but not promised: core stores each
dpar's linear predictor beside it in the `dpars` list as `.eta_<dpar>`,
which is what makes the coherence complement exact (see
[`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md),
"Why the links are the constraint"). It works, it is what core's own
binomial and beta densities rely on through internal accessors, and no
sibling extension uses it. It is undocumented for extension authors, and
a promise would be worth making.

What is NOT missing, and this package said otherwise in 0.1.0: a
matrix-valued response. `R/frame.R` preserves one explicitly, and a
custom family indexing `y[, 1]` and `y[, 2]` fits today. The reason this
package stops at two channels is its own, and it is set out in
[`cross_wishart()`](https://aforren1.github.io/frmtmb/frmtmb.coupling/reference/cross_wishart.md)
under "What it refuses".

## See also

Useful links:

- <https://aforren1.github.io/frmtmb/frmtmb.coupling>

- <https://github.com/aforren1/frmtmb>

- Report bugs at <https://github.com/aforren1/frmtmb/issues>

## Author

**Maintainer**: Alex Forrence <alex.forrence@gmail.com>
([ORCID](https://orcid.org/0000-0002-9728-6337))

Authors:

- Alex Forrence <alex.forrence@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-9728-6337))
