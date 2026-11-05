The main changes:

- BREAKING: the exact `gp()` is non-centered, brms's `zgp`, at a nugget of
  `1e-12 * sdgp^2`, and starts at the maximum of the model's marginal.
- BREAKING: an autocorrelation term on a family with no residual `y - mu`
  is brms's latent residual, as `ar(cov = TRUE)` and more.
- BREAKING: a smooth's standard deviation has brms's class `sds`.
- New: a nonlinear body under `cumulative()`, `sratio()`, `cratio()` and
  `acat()`, and on a `Gamma()` or `inverse.gaussian()` default link.
- `fitted()` on a new row of `me()` or of a latent residual is the mean
  over that draw, and `Est.Error` holds the draw.
- `frm(importance =)` reads the Laplace fit's own scale, in fewer rounds.
- Update the extensions together with frmtmb. frmtmb.sample 0.26.0 and
  frmtmb.learn 0.7.4 need frmtmb 0.79.0.
