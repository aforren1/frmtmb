# Review of lane code8: the SE analysis of the poly tight-sd stop, run
# directly on each build (so the base's verdict at a stop it calls short
# is seen too), and on the lane with fd_hessian()'s retake as 0.80.0
# kept it (only when the matrix does not factor) and with the Ritz step
# off, to say which change loses the fixed effects.
#   RV_ARM=before|after Rscript dev/review-code8-tightsd4.R <variant>
variant <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
if (variant == "noritz") {
  assignInNamespace("se_ridge_ritz", function(...) NULL, "frmtmb")
}
env <- new.env(parent = asNamespace("frmtmb"))
for (ex in parse("tests/testthat/test-brms-names-specials.R")) {
  if (is.call(ex) && identical(ex[[1]], as.name("<-"))) eval(ex, env)
}
d <- env$bns_prior_data()
form <- bf(y ~ poly(x, 2) + (1 + poly(x, 2) | g))
p <- set_prior("normal(2, 0.001)", class = "sd", coef = "polyx22", group = "g")
f <- quiet(frm(form, data = d, prior = p))
f$obj$env$.frmtmb_outer_h <- NULL
for (nm in c("se_fit_an", "hessian_fixed", "se_analysis", "sdr")) {
  if (exists(nm, f$cache)) rm(list = nm, envir = f$cache)
}
an <- quiet(ns$se_fit_analysis0(f, TRUE))
cat(variant, "tier", an$tier, "lost:", paste(names(an$lost), an$lost), "\n")
h <- ns$fit_outer_hessian(f)
Hh <- h$H; D <- sqrt(abs(diag(Hh)))
cat("FD unit-diag ev:", signif(eigen(Hh / outer(D, D), TRUE, TRUE)$values, 3), "\n")
cat("DONE\n")
