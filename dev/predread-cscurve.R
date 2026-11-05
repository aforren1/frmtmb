# D4 (pair p013, seed 13): the cs() offsets on every reader of the
# ordinal linear predictor, frmtmb against brms 2.23.0 at the same
# parameter values (a Fixed_param brmsfit through brms_fixed_fit()).
#   Rscript dev/predread-cscurve.R <lib> > dev/predread-log/cscurve-<tag>.txt
# <lib> is the library put before rellib-r16 (rellib-r16 itself for the
# "before" arm).
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else "C:/Users/adf44/source/r/rellib-r16"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r16",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
spec <- E$eq_specs[["cs_pred() x frm_curve"]]()
fit <- E$eq_frm_fit(spec)
nd <- spec$margs$newdata
tr <- function(e) tryCatch(e, error = function(err) {
  paste("ERROR:", conditionMessage(err))
})
bfit <- E$eq_fixed(spec, brms::bf(y ~ x + cs(z)), brms::sratio(),
                   spec$data, fit)
bl <- brms::posterior_linpred(bfit, newdata = nd, re_formula = NA)
cat("brms posterior_linpred dim:", dim(bl), "\n")
cat("brms dimnames:\n"); str(dimnames(bl))
b1 <- bl[1, , ]
rel <- function(a, b) E$eq_rel(a, b)
cv <- tr(frm_curve(fit, newdata = nd, simultaneous = FALSE))
cat("frm_curve nrow:", tr(nrow(cv)), " grid rows:", nrow(nd), "\n")
cat("frm_curve .estimate vs brms (flattened):",
    tr(rel(cv$.estimate, b1)), "\n")
if (is.data.frame(cv)) print(utils::head(cv, 3))
lk <- tr(frm_linpred(fit, newdata = nd, type = "link"))
cat("frm_linpred(link) dim:", tr(paste(dim(lk), collapse = " x ")),
    " length:", length(lk), "\n")
cat("frm_linpred(link) vs brms:", tr(rel(lk, b1)), "\n")
ft <- fitted(fit, newdata = nd, scale = "linear")
cat("fitted(linear) Estimate vs brms:", rel(ft[, "Estimate", ], b1), "\n")
se <- tr(frm_linpred(fit, newdata = nd, type = "link", se.fit = TRUE))
if (is.list(se)) {
  cat("frm_linpred se.fit vs fitted(linear) Est.Error:",
      rel(se$se.fit, ft[, "Est.Error", ]), "\n")
}
lb <- tr(frm_lp_basis(fit, newdata = nd, re_formula = NA))
if (is.list(lb)) {
  cat("frm_lp_basis eta vs brms:", rel(lb$eta, b1), "; coef_names:",
      paste(lb$coef_names, collapse = " "), "\n")
}
if (is.data.frame(cv)) {
  cat("frm_curve .se vs fitted(linear) Est.Error:",
      tr(rel(cv$.se, ft[, "Est.Error", ])), "\n")
}
cvs <- tr(frm_curve(fit, newdata = nd, simultaneous = TRUE, nsim = 2000,
                    seed = 1))
if (is.data.frame(cvs)) cat("simultaneous crit:", cvs$.crit_sim[1], "\n")
dv <- tr(frm_curve_deriv(fit, var = "x", newdata = nd,
                         simultaneous = FALSE))
cat("frm_curve_deriv nrow:", tr(nrow(dv)), "\n")
if (is.data.frame(dv)) {
  # the cs() part does not move with x, so every threshold's slope is b_x
  cat("deriv .estimate range:", range(dv$.estimate), " b_x:",
      fixef(fit)["x", "Estimate"], "\n")
}
dz <- tr(frm_curve_deriv(fit, var = "z", newdata = nd,
                         simultaneous = FALSE))
if (is.data.frame(dz)) {
  # d/dz of eta + cs_k is the cs() coefficient of threshold k
  fz <- fixef(fit)
  cat("deriv in z, by row:\n")
  print(utils::head(dz[, intersect(c("x", "z", ".threshold", ".estimate",
                                     ".se"), names(dz))], 12))
  cat("cs coefficients:\n")
  print(fz[grepl("z", rownames(fz)), c("Estimate", "Est.Error")])
}
# brms's posterior_linpred in sample, same layout
bl0 <- brms::posterior_linpred(bfit, re_formula = NA)
lk0 <- tr(frm_linpred(fit, type = "link", re_formula = NA))
cat("in sample: brms dim", dim(bl0), "; frm_linpred vs brms:",
    tr(rel(lk0, bl0[1, , ])), "\n")
