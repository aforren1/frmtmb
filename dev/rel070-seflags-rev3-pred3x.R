# Reviewer of lane seflags, re-check of B1 on fixture 1 (OpenBLAS
# LAPACK): the hypothesis exactly on the alias (t at full precision),
# at relative offsets 1e-9 to 1e-3, b_Intercept alone, and
# fitted(newdata) with and without re_formula = NA.
#   Rscript dev/seflags-rev3-pred3x.R <lane|base|op2>
source("C:/Users/adf44/source/r/frmtmb-wt-release/dev/rel070-seflags-rev3-common.R")
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
X <- fit$frame$linpreds[[1]]$X
cn <- colnames(X)
cf <- coef(stats::lm(as.numeric(X[, grep("^s[(]Age[)][.]?fx", cn)[1]]) ~
                       fit$data$Age))
for (rel in c(0, 1e-9, 1e-7, 1e-5, 1e-3)) {
  t <- cf[[2]] * (1 + rel)
  h <- sprintf("b_Age + %.17g * bs_sAge_1 = 0", t)
  z <- cap(hypothesis(fit, h, class = NULL))
  cat(sprintf("rel offset %-6g Est.Error %-10s warnings %d\n", rel,
              format(z$v$hypothesis$Est.Error, digits = 6), length(z$w)))
}
z <- cap(hypothesis(fit, "b_Intercept = 0", class = NULL))
cat("b_Intercept Est.Error", z$v$hypothesis$Est.Error, "warnings",
    length(z$w), "\n")
nd <- fit$data[1:5, ]
for (rf in list(NULL, NA)) {
  p <- cap(fitted(fit, newdata = nd, re_formula = rf))
  cat("fitted(newdata, re_formula =", format(rf), ") se",
      paste(signif(p$v[, 2], 4), collapse = " "), "| warnings",
      length(p$w), "\n")
}
