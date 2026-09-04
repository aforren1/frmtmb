Adds 'frm_ode()' to the nonlinear formula grammar of 'frmtmb'. A model
writes its dynamics as a plain R function, and 'frm_ode()' solves one
system per group and returns the states as a predictor. The solve is
differentiated exactly, so the fit stays a maximum-likelihood fit with
the Laplace approximation for random effects. Repeated doses,
steady-state dosing and time-varying parameters are supported.
Compartment models in pharmacokinetics are the main application.
