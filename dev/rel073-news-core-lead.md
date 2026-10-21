The main changes:

- `cens()` and `trunc()` for 23 more families, every one brms 2.23.0
  censors, with brms's reading of a discrete bound.
- `brmsfamily("negbinomial2")` and `dirichlet_multinomial()`.
- The generics of brms's `loo_*()` functions, `kfold()`, `reloo()`,
  `loo_moment_match()` and `loo_subsample()`, with a refusal on a fit.
- `quadrature = TRUE` refuses a row that reads more than one
  random-effect block, and `frm()` refuses a weighted likelihood with
  no maximum.
