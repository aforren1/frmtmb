Needs frmtmb 0.73.0. It registers its draws methods on the `loo_*()`,
`kfold()`, `reloo()`, `loo_moment_match()` and `loo_subsample()`
generics that frmtmb now defines, and re-exports them. It calls the
sampling-API export `check_weights_unbounded()`.
