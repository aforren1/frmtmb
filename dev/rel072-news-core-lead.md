The main changes:

- brms's `categorical()` and `multinomial()` signatures, `refcat`,
  `fitted()` and `posterior_linpred()` of every location predictor,
  and `logistic_normal()` in `mvbf()`.
- Per-coefficient `sd` priors, brms's labels and order of nested and
  multi-membership groups, and one Student-t latent scale for the
  terms of one grouping factor.
- `uniform()` and `lognormal()` priors, and support guards on the
  tape.
- Every addition term on `newdata`, coded as the fit coded it.
- Accurate standard errors for Cook's distance, D1 and D2, REML in the
  design's coordinates, and bounded Newton steps.
- The exact Hessian for every fit without random effects.
- `drop1()`, `step()`, `update()` and `influence()` refit on the data
  the fit holds.
