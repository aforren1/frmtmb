# frmtmb: Formula-Based Regression Models via 'RTMB'

Fits regression models specified with a 'brms'-style formula grammar by
maximum likelihood, with the Laplace approximation for random effects.
Model objectives are generated as R closures and differentiated with
'RTMB', so no compilation happens at run time. Supports 'lme4'-style
random effects with structured covariances, and (in later versions)
distributional regression, smooths, nonlinear formulas, and multivariate
models.

## Author

**Maintainer**: Alex Forrence <alex.forrence@gmail.com>

Authors:

- Alex Forrence <alex.forrence@gmail.com>
