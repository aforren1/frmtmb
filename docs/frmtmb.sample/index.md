Samples the posterior of a model assembled by 'frmtmb', by handing its
'RTMB' objective to 'tmbstan' and so to Stan's NUTS sampler. Applies the
'brms' default priors, samples the qualifying random-effect blocks in
their non-centered form, and returns draws under 'frmtmb' parameter
names. Provides the posterior method surface those draws support:
expectations and predictions, pointwise log-likelihood and approximate
leave-one-out cross-validation, hypothesis tests, conditional effects,
posterior predictive checks, and conversion to the 'posterior' and
'coda' draws formats.
