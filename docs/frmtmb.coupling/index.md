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
