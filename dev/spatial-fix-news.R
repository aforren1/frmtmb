# Fix round 1: core NEWS for the A and B fixes, on the tree named as the
# argument (the text edited is the same at the AB stage and after C).
root <- commandArgs(trailingOnly = TRUE)[1]
p <- file.path(root, "NEWS.md")
s <- paste(readLines(p, warn = FALSE), collapse = "\n")
rep1 <- function(s, old, new) {
  if (!grepl(old, s, fixed = TRUE)) stop("no match: ", substr(old, 1, 60))
  sub(old, new, s, fixed = TRUE)
}
s <- rep1(s, '  are entries of the same shape. The structure joins `frm_compat()` as
  a feature of kind `"covstruct"`. `?frmtmb-extension-api` has the
  steps.', '  are entries of the same shape. The structure joins `frm_compat()` as
  a feature of kind `"covstruct"`, where quadrature and `importance =`
  refuse it, as `frm()` does. `?frmtmb-extension-api` has the steps.
  The call checks the entry it builds, and each block is checked
  against its entry when the frame is built: the length of `start`, the
  positions of `sd_idx`, the dimension and names of `vcov`, and for a
  registered structure that `lpdf` is a finite log density, larger at
  zero than far from it. A name the rest of frmtmb, brms or another
  extension owns is refused (`poly`, `ns`, `gaussian`, `dec`, ...), and
  a name another package registered is replaced only with
  `replace = TRUE`. A fit records the package each structure came from:
  read where that package is not installed, the fit is refused by the
  package\'s name, with class `frmtmb_extension_missing`, before a method
  reads the structure.

* `frmtmb_local_registries()` saves every registry a package fills from
  its `.onLoad()` and restores it when a test ends, for an extension\'s
  own test suite.

* `frm_compat()` resolves the rules that packages contribute in the
  order of the packages\' names, so its answer no longer depends on the
  order of the `library()` calls.')
s <- rep1(s, '  frmtmb spelling: `bhaz()` (the basis is set in `cox(df =, degree =,
  intercept =)`) and `brmsfamily("custom")` (use `custom_family()`).
  Every family, addition term and formula special of brms 2.23.0 now
  resolves in core, resolves in a named extension, or has such a
  refusal; before, three reached a generic "not supported" list.', '  frmtmb spelling: `bhaz()` (the basis is set in `cox(df =, degree =,
  intercept =)`, and `bhaz(gr = )` has no equivalent) and
  `brmsfamily("custom")` (use `custom_family()`). Every family, addition
  term and formula special of brms 2.23.0, in its plainest spelling and
  in the others brms documents, now resolves in core, resolves in a
  named extension, or has such a refusal; before, `custom`, `bhaz()`,
  `dec()` and bare `gr(g)` did not.')
s <- rep1(s, '## Bug fixes

* A call around a bar term that names no covariance structure is
  refused. Through 0.73.0 the call was dropped and the term fitted with
  the default `us` structure with no message: `foo(1 | g)` and
  `us_t(1 | g)` both gave the log-likelihood of `(1 | g)`.', '## Bug fixes

* A call around a bar term that names no covariance structure is
  refused, at any depth and in `re_formula =` too. Through 0.73.0 the
  call was dropped and the term fitted with the default `us` structure
  with no message: `foo(1 | g)`, `foo(bar(1 | g))` and `us_t(1 | g)` all
  gave the log-likelihood of `(1 | g)`, and `re_formula = ~ foo(1 | g)`
  was ignored.

* `(1 | gr(g))`, brms\'s plainest spelling of a grouping factor, is
  `(1 | g)`, as it is in brms. Through 0.73.0 it was refused.

* A gaussian model whose field has one row per level, such as `car()`
  over areas with one observation each, no longer warns that its
  variance is confounded with the residual sd. A fixed matrix
  correlates the levels, as for `gr(cov = A)`, which the check already
  exempted; through 0.73.0 it named three structures and missed the
  rest.')
writeLines(s, p)
cat("fix-news applied to", root, "\n")
