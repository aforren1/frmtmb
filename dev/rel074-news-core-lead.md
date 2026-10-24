The main changes:

- Four features moved to extension packages. Install the package and
  load it to use them:
  - `cox()` and `cox_baseline()`: `library(frmtmb.survival)`.
  - `royston_parmar()` and `rp_floored()`, from the retired
    frmtmb.spline: `library(frmtmb.survival)`.
  - `mixture_mvn()`: `library(frmtmb.latent)`.
  - `car()` and `spde()`: `library(frmtmb.spatial)`.

  Without the package, `family = "cox"`, `brmsfamily("cox")`, `car()`
  and `spde()` stop with a message that names the package. The README
  shows how to install every extension at once.
- frmtmb.spline is retired. Curve inference, `frm_curve()`,
  `frm_curve_deriv()` and `frm_curve_feature()`, is in frmtmb now.
  Remove an installed copy with `remove.packages("frmtmb.spline")`:
  it still loads with frmtmb 0.74.0, masks these three functions and
  replaces their print methods.
- New seams for extension packages: `frmtmb_register_covstruct()`,
  `frmtmb_family(simplex_pars =, aterm_data =)`, the family hook
  `post$curve_check`, `us_sigma()` and `us_theta()`, and
  `frmtmb_local_registries()`.
- Every brms family, addition term and formula special resolves in
  core, resolves in a named extension, or has a refusal that says why.
- Update the extensions together with frmtmb. frmtmb.latent 0.7.0,
  frmtmb.survival 0.1.0, frmtmb.spatial 0.1.0 and frmtmb.sample 0.21.1
  need frmtmb 0.74.0, and frmtmb.latent 0.6.2 does not load with
  frmtmb 0.74.0.
