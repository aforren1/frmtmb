# Reviewer of lane seflags, claim 1, the cap: on ported fixture 1 under
# OpenBLAS LAPACK flat_tol is 0.5 * the largest projection (0.359). A
# hypothesis b_Age + t * bs_sAge_1 is estimable only at t = 0.190041
# (the exact alias: fx1 = -0.0464 + 0.1900 Age). Sweep t; print the
# reported Est.Error, the warnings, the cosine of the hypothesis with the
# EXACT alias and its projection onto the stored flat subspace, both in
# the unit-diagonal metric jc_nonest() uses.
#   Rscript dev/seflags-rev2-pred3.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev2-common.R")
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages(library(testthat))
env <- new.env(parent = ns)
root <- "C:/Users/adf44/source/r/frmtmb-wt-seflags/tests/testthat"
for (h in list.files(root, "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = env)
}
withr::local_seed(3L)
fit <- suppressWarnings(suppressMessages(env$brms_fixture(1)))
cat("fixture 1 | code", fit$opt$convergence, "| lost", lost_of(fit), "\n")
d <- fit$data
X <- fit$frame$linpreds[[1]]$X
cn <- colnames(X)
cf <- coef(stats::lm(as.numeric(X[, grep("^s[(]Age[)][.]?fx", cn)[1]]) ~
                       d$Age))
on <- ns$outer_par_names(fit)
nex <- numeric(length(on))
nex[on == "(Intercept)"] <- cf[1]
nex[on == "Age"] <- cf[2]
nex[on == "s(Age).fx1"] <- -1
sdr <- suppressWarnings(ns$sdr_of(fit))
fl <- sdr$se_null_flat
H <- ns$fit_outer_hessian(fit)$H
sc <- if (!is.null(fl)) fl$scale else sqrt(abs(diag(H)))
r <- which(is.finite(sc) & sc > 0)
if (!is.null(fl)) {
  B <- sdr$se_null[r, fl$cols, drop = FALSE] * sc[r]
  B <- qr.Q(qr(B))
  cat("flat tol", fl$tol, "| flat columns", length(fl$cols), "\n")
} else B <- NULL
n_s <- nex[r] * sc[r]
n_s <- n_s / sqrt(sum(n_s^2))
if (!is.null(B)) {
  cat("projection of the exact alias onto the stored subspace:",
      sqrt(sum(crossprod(B, n_s)^2)), "\n")
}
vn <- variables(fit)
cat("variables:", paste(vn[1:8], collapse = " "), "\n")
for (t in c(cf[2], 0.15, 0.1, 0.05, 0, -0.05, -0.1, -0.2, -0.5, -1)) {
  G <- numeric(length(on))
  G[on == "Age"] <- 1
  G[on == "s(Age).fx1"] <- t
  gs <- G[r] / sc[r]
  ctrue <- abs(sum(gs * n_s)) / sqrt(sum(gs^2))
  pB <- if (!is.null(B)) {
    sqrt(sum(crossprod(B, gs)^2)) / sqrt(sum(gs^2))
  } else NA
  h <- sprintf("b_Age + %.6f * bs_sAge_1 = 0", t)
  z <- cap(hypothesis(fit, h, class = NULL))
  est <- if (inherits(z$v, "rev_err")) z$v else z$v$hypothesis$Est.Error
  cat(sprintf("t %9.6f cos(exact alias) %.4f proj(stored) %.4f | %s | warnings %d\n",
              t, ctrue, pB, format(est, digits = 6), length(z$w)))
}
