# Reviewer of lane seflags: under OpenBLAS LAPACK the lane also names
# theta_1 and theta_2 of fixture 1 "flat", which the base does not (the
# lane's own ext-ob3/dh7diag-lane.txt; its findings say the verdicts are
# unchanged). Is the likelihood flat along that direction? Objective
# changes along each stored flat column scaled so its largest theta
# move is h, both ways; base SEs for comparison; which rule removed it.
#   Rscript dev/seflags-rev-fx1theta.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-common.R")
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages(library(testthat))
env <- new.env(parent = ns)
root <- "C:/Users/adf44/source/r/frmtmb-wt-seflags/tests/testthat"
for (h in list.files(root, "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = env)
}
if (exists("se_ridge_curvature", ns)) {
  trace("se_ridge_curvature", exit = quote(cat(sprintf(
    "  RIDGE k %d ev %.3g -> c0 %.3g\n", k, e$values[k], returnValue()))),
    print = FALSE, where = ns)
  trace("se_tier3", exit = quote(if (exists("rescued")) cat(sprintf(
    "  TIER3 tau %.3g flat_tol %s | removed ev %s | try_k %s\n", tau,
    if (exists("flat_tol")) format(flat_tol, digits = 3) else "-",
    paste(signif(ev[!keep], 3), collapse = " "),
    paste(signif(ev[try_k], 3), collapse = " ")))), print = FALSE, where = ns)
}
withr::local_seed(3L)
fit <- suppressWarnings(suppressMessages(env$brms_fixture(1)))
on <- ns$outer_par_names(fit)
sdr <- suppressWarnings(ns$sdr_of(fit))
cat("code", fit$opt$convergence, "| lost", lost_of(fit), "\n")
se <- sqrt(diag(as.matrix(sdr$cov.fixed)))
cat("theta estimates:", paste(sprintf("%s=%.4g", on[grepl("^theta", on)],
                                       fit$opt$par[grepl("^theta", on)]),
                               collapse = " "), "\n")
cat("SE theta_1..4:", format(se[grepl("^theta_", on)], digits = 4), "\n")
p <- fit$opt$par
f0 <- fit$obj$fn(p)
N <- sdr$se_null
th <- grepl("^theta_[12]$", on)
for (j in seq_len(NCOL(N))) {
  v <- N[, j]
  if (max(abs(v[th])) < 1e-3 * max(abs(v))) next
  v <- v / max(abs(v[th]))
  ch <- vapply(c(-1, -0.3, -0.1, 0.1, 0.3, 1), function(h) {
    fit$obj$fn(p + h * v) - f0
  }, 0)
  cat(sprintf("null col %d (theta_1 %.3f theta_2 %.3f): objective change at theta move -1 -0.3 -0.1 0.1 0.3 1: %s\n",
              j, v[on == "theta_1"], v[on == "theta_2"],
              paste(format(ch, digits = 3), collapse = " ")))
}
fit$obj$fn(p)
# each of theta_1, theta_2 alone
for (nm in c("theta_1", "theta_2")) {
  ch <- vapply(c(-1, -0.3, 0.3, 1), function(h) {
    q <- p; q[on == nm] <- q[on == nm] + h; fit$obj$fn(q) - f0
  }, 0)
  cat(nm, "alone: change at -1 -0.3 0.3 1:", format(ch, digits = 3), "\n")
}
fit$obj$fn(p)
