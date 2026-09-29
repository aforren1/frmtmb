# REVIEW script 13: rescor_row_loglik() itself across nu, on BOTH
# builds, on one fit and one parameter vector per nu. This is the
# reference build's side of script 12, which its draws_chain_id() bug
# stops before the density is reached.
#
#   Rscript dev/arcovsample-rev-13-rrl.R <lane|ref>

a <- commandArgs(trailingOnly = TRUE)
arm <- a[1L]
LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
suppressMessages(library(frmtmb))
newx <- "arma_cond_resp" %in% getNamespaceExports("frmtmb")
cat("ARM ", arm, " newexports=", newx, "\n", sep = "")
stopifnot(identical(newx, identical(arm, "lane")))

set.seed(777L)
N <- 50L
x <- rnorm(N)
E <- matrix(rnorm(2 * N), N) %*% chol(matrix(c(1, 0.4, 0.4, 1), 2)) /
  sqrt(rchisq(N, 5) / 5)
dd <- data.frame(x = x, y1 = 0.4 + 0.9 * x + E[, 1],
                 y2 = -0.1 + 0.5 * x + E[, 2])
fit <- frm(bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE) + student(),
           data = dd)
bn <- names(fit$frame$par_template$betad)
jnu <- match("nu_(Intercept)", bn)
stopifnot(!is.na(jnu))
nus <- c(3, 10, 100, 1e3, 1e5, 1e7, 1e9, 1e11, 1e15, 1e20, 1e50, 1e100,
         1e200, 1e300, 1e306)
cat("\n", sprintf("%10s %18s %18s %14s", "nu", "sum rescor_row_loglik",
                  "-obj$fn", "abs resid"), "\n", sep = "")
nbad <- 0L
for (nu in nus) {
  p <- fit$opt$par
  p[which(names(p) == "betad")[jnu]] <- log(nu - 1)
  objv <- -as.numeric(fit$obj$fn(p))
  f2 <- fit
  f2$estimates <- fit$obj$env$parList(p)
  f2$cache <- new.env(parent = emptyenv())
  rl <- sum(frmtmb::rescor_row_loglik(f2, frmtmb::eval_dpars(f2)))
  bad <- !is.finite(rl) || abs(rl - objv) > 1e-6 * abs(objv)
  if (bad) nbad <- nbad + 1L
  cat(sprintf("%10.1e %18.8f %18.8f %14.6g%s\n", nu, rl, objv,
              abs(rl - objv), if (bad) "  BAD" else ""))
}
cat("\nnu values where the pointwise density disagrees with the ",
    "objective by more than 1e-6 relative, or is not finite: ", nbad,
    " of ", length(nus), "\n", sep = "")
cat("DONE\n")
