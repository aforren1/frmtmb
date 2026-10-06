# Punch round 1: test-case-studies.R's fmv_id (identity cov, genetic and
# residual variances one sum). Does the SE check defer, and what does it
# report at first SE use?
.libPaths(c("C:/Users/adf44/source/r/wt-nanse-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
src <- readLines("tests/testthat/test-case-studies.R")
eval(parse(text = src[8:41]))
ped <- build_pedigree(12, 2, 6)
n <- nrow(ped)
A <- build_A(ped)
set.seed(4)
G <- matrix(c(1.0, 0.6, 0.6, 0.8), 2, 2)
U <- t(chol(A)) %*% matrix(rnorm(n * 2), n, 2) %*% chol(G)
long <- data.frame(
  id = factor(rep(ped$id, times = 2), levels = ped$id),
  trait = factor(rep(c("y1", "y2"), each = n)),
  value = c(3 + U[, 1] + rnorm(n, 0, 0.7),
            1 + U[, 2] + rnorm(n, 0, 0.9)))
I <- diag(n)
dimnames(I) <- dimnames(A)
t0 <- proc.time()[[3]]
f <- withCallingHandlers(
  frm(bf(value ~ 0 + trait + (0 + trait | gr(id, cov = I)),
         sigma ~ 0 + trait) + gaussian(), data = long, data2 = list(I = I)),
  warning = function(w) {
    cat("FIT WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })
cat("fit time", proc.time()[[3]] - t0, "\n")
cat("deferred:", !is.null(f$cache$se_deferred), " hess cached:",
    !is.null(f$cache$hessian_fixed), "\n")
t0 <- proc.time()[[3]]
H <- ns$fit_outer_hessian(f)
cat("FD Hessian time", proc.time()[[3]] - t0, "n par", length(f$opt$par), "\n")
withCallingHandlers(print(sqrt(diag(vcov(f)))),
  warning = function(w) {
    cat("USE WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })
