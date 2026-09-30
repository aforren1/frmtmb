# Lane ceplot punch 1: frmtmb.sample NEWS (m5, m6).
f <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/extensions/frmtmb.sample/NEWS.md"
s <- paste(readLines(f), collapse = "\n")
rep1 <- function(s, old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("expected one match, got ", n)
  sub(old, new, s, fixed = TRUE)
}
s <- rep1(s, "* `posterior_samples(pars = )` returns the coefficients in brms's
  order, every intercept first (`b_Intercept`, `b_sigma_Intercept`,
  `b_x`, `b_sigma_x`), as brms's `variables()` lists them.",
"* `posterior_samples(pars = )` returns the coefficients in brms's
  order, as brms's `variables()` lists them: the intercept of every
  distributional parameter first (`b_Intercept`, `b_sigma_Intercept`,
  `b_x`, `b_sigma_x`), a nonlinear parameter's intercept staying with
  its own coefficients (`b_a_Intercept`, `b_a_z`, `b_b_Intercept`).")
s <- rep1(s, "* `parnames()` is frmtmb's generic, re-exported here; the method for
  draws is unchanged.",
"* `parnames()` is frmtmb's generic, re-exported here. With brms loaded
  it warns once, as brms's does; brms's generic and the method each
  warned before.")
writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], f)
cat("edited\n")
