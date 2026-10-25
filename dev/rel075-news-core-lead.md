The main changes:

- BREAKING: coefficient names and orders follow brms 2.23.0. A term
  that calls `exp()`, or another function named like a covariance
  structure, has brms's name (`expx`, not `.frm_expx`).
  `default_prior()` lists each coefficient by brms's name (`IxE2`,
  `polyx21`), and the prior tables sort their rows as brms sorts them.
  `fixef()` orders nonlinear and multivariate rows as brms does.
- BREAKING: brms's `bhaz()` addition term is in frmtmb.survival. Without
  that package, `bhaz()` stops with a message that names it.
- BREAKING: `set_prior()` accepts a `dirichlet()` row for every class.
  `validate_prior()` and `frm()` refuse it when the class is not a
  simplex of the model.
- New: brms's `fcor()` residual term is in frmtmb, and brms's `sar()`
  term is in frmtmb.spatial.
- New seams for extension packages: `frmtmb_register_aterm(option =)`
  for an addition term that sets an argument of the family, and
  `frmtmb_register_autocor()` for a residual correlation term. With
  them come `frmtmb_autocor_pars()`, `frmtmb_covstruct_info()` and
  `frmtmb_register_covstruct(full_rank =)`.
- A saved fit that needs extension packages this session cannot load
  names every one of them in one refusal.
- Update the extensions together with frmtmb. frmtmb.survival 0.2.0,
  frmtmb.spatial 0.2.0, frmtmb.sample 0.22.0 and frmtmb.latent 0.7.1
  need frmtmb 0.75.0.
