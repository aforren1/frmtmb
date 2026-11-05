# Item 10: is.numeric() on an advector, the jitter on gpgap's tape, and gpgap's two test
# files against this build (a copy of this tree's tests/testthat holds them).
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
library(RTMB)
res <- NULL
tp <- RTMB::MakeTape(function(p) { res <<- list(num = is.numeric(p), cls = class(p), adv = inherits(p, "advector")); sum(p) }, c(1, 2))
str(res)
# gpgap's line on a tape, at a nugget where chol() fails at the build
gp_chol_num_g <- function(K) {
  R <- tryCatch(chol(K), error = function(e) NULL)
  jit <- 1e-12
  while (is.null(R) && jit <= 1e-6) { Kj <- K; diag(Kj) <- diag(Kj) + jit; R <- tryCatch(chol(Kj), error = function(e) NULL); jit <- jit * 10 }
  if (is.null(R)) chol(K) else R
}
set.seed(2); u <- sort(runif(30, 0, 6)); D2 <- outer(u, u, "-")^2
for (nug in c(1e-12, 1e-17)) {
  kf <- function(lr) exp(-D2 / (2 * exp(2 * lr))) + diag(nug, 30)
  for (lr0 in log(c(5, 500))) {
    tape1 <- tryCatch(RTMB::MakeADFun(function(p) { K <- kf(p$lr); L <- t(if (is.numeric(K)) gp_chol_num_g(K) else chol(K)); sum(log(RTMB::diag(L))) }, list(lr = lr0), silent = TRUE), error = function(e) e)
    tape2 <- tryCatch(RTMB::MakeADFun(function(p) { K <- kf(p$lr); L <- t(if (inherits(K, "advector")) chol(K) else gp_chol_num_g(K)); sum(log(RTMB::diag(L))) }, list(lr = lr0), silent = TRUE), error = function(e) e)
    a <- if (inherits(tape1, "error")) paste("BUILD ERROR", substr(conditionMessage(tape1), 1, 60)) else paste(vapply(log(c(5, 500, 5e4)), function(v) format(tape1$fn(v), digits = 7), ""), collapse = " ")
    b <- if (inherits(tape2, "error")) paste("BUILD ERROR", substr(conditionMessage(tape2), 1, 60)) else paste(vapply(log(c(5, 500, 5e4)), function(v) format(tape2$fn(v), digits = 7), ""), collapse = " ")
    cat(sprintf("RMG nug %g start rho %.0f | is.numeric line: %s | inherits line (rho 500): %s\n", nug, exp(lr0), a, b))
  }
}
td <- "C:/Users/adf44/source/r/rev-ncgp-src/tests"
unlink(td, recursive = TRUE); dir.create(td, recursive = TRUE)
for (f in list.files(file.path(REV_WT, "tests/testthat"), pattern = "^helper-", full.names = TRUE)) file.copy(f, td)
gg <- "C:/Users/adf44/source/r/frmtmb-wt-gpgap/tests/testthat"
for (f in c("test-gp-nugget-floor.R", "test-brms-lik-gp-exact.R")) if (file.exists(file.path(gg, f))) file.copy(file.path(gg, f), td)
Sys.setenv(FRMTMB_STAN_CACHE = file.path(REV_WT, "dev/review-ncgp-stan-cache"))
suppressPackageStartupMessages(library(testthat))
for (f in list.files(td, pattern = "^test-", full.names = TRUE)) {
  r <- as.data.frame(testthat::test_file(f, package = "frmtmb", env = testthat::test_env("frmtmb"), reporter = "silent", stop_on_failure = FALSE))
  cat(sprintf("RMT %s tests %d pass %d fail %d err %d skip %d\n", basename(f), nrow(r), sum(r$nb - r$failed - r$warning - r$skipped), sum(r$failed), sum(r$error), sum(r$skipped)))
  print(r[r$failed > 0 | r$error, c("test", "failed", "error")], row.names = FALSE)
}
cat("DONE\n")
