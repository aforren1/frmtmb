
## Bug fixes

* The drawn route for an exact `gp()` needs long chains at low noise, and the
  numbers are stated. At true sigma 1e-3, 15 rows, 4 chains of 12000 read sigma
  5.9e-4, 1.05e-3, 2.3e-3 (Rhat 1.025, sdgp Rhat 1.07, 58 minutes), and brms's
  read 5.2e-4, 1.03e-3, 2.4e-3 (1.031, 1.15, 53 minutes, 20000 iterations).
  Neither reaches Rhat 1.01 at those lengths. `frm_sample(laplace = TRUE)` serves
  the outer parameters only: it refuses `posterior_epred()`,
  `posterior_linpred()`, `posterior_predict()`, `fitted()`, `predict()`,
  `log_lik()`, `loo()` and `conditional_effects()`.
