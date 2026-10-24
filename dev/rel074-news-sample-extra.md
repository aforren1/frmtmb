## Bug fixes

* The note of default priors from `frm_sample()` names every simplex
  that a family declares with `frmtmb_family(simplex_pars =)`, by its
  brms name. It named only the `sbhaz` simplex of `cox()`, so the
  `dirichlet(1)` default on the simplex of another family was not
  announced. For `cox()` the note now reads "sbhaz: dirichlet(1),
  brms's default, on the simplex the cox() family declares".

## Other changes

* frmtmb.survival is a suggested package. The tests of `cox()` draws
  read `cox()` from it, and skip without it.
