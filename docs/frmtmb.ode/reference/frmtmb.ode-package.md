# frmtmb.ode: Ordinary Differential Equation Dynamics for 'frmtmb' Models

Adds 'frm_ode()' to the nonlinear formula grammar of 'frmtmb'. A model
writes its dynamics as a plain R function, and 'frm_ode()' solves one
system per group and returns the states as a predictor. The solve is
differentiated exactly, so the fit stays a maximum-likelihood fit with
the Laplace approximation for random effects. Repeated doses,
steady-state dosing and time-varying parameters are supported.
Compartment models in pharmacokinetics are the main application.

## See also

Useful links:

- <https://aforren1.github.io/frmtmb/>

- <https://github.com/aforren1/frmtmb>

- Report bugs at <https://github.com/aforren1/frmtmb/issues>

## Author

**Maintainer**: Alex Forrence <alex.forrence@gmail.com>
([ORCID](https://orcid.org/0000-0002-9728-6337))

Authors:

- Alex Forrence <alex.forrence@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-9728-6337))
