# Reviewer: cost of sar() at a realistic areal size, a 50 x 50 rook
# lattice (n = 2500), symmetric-normalized W. Seed 8.
# Rscript dev/review-sarfcor-scale.R
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb.spatial); library(Matrix)})
cat("lib:", find.package("frmtmb"), "\n")
k <- 50
n <- k * k
id <- matrix(1:n, k)
ii <- c(id[-k, ], id[, -k]); jj <- c(id[-1, ], id[, -1])
A <- sparseMatrix(i = c(ii, jj), j = c(jj, ii), x = 1, dims = c(n, n))
dg <- 1 / sqrt(rowSums(A))
W <- Diagonal(x = dg) %*% A %*% Diagonal(x = dg)
set.seed(8)
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(Diagonal(n) - 0.5 * W, 1 + 0.5 * d$x + rnorm(n)))
tm <- function(lab, expr) {
  t0 <- proc.time()[["elapsed"]]
  r <- tryCatch(force(expr), error = function(e) conditionMessage(e))
  cat(sprintf("%-28s %8.2f s\n", lab, proc.time()[["elapsed"]] - t0))
  invisible(r)
}
fit <- tm("frm() lag, n = 2500", frm(y ~ x + sar(W), data = d,
                                      data2 = list(W = W)))
cat("lagsar", format(confint_varcorr(fit)$estimate, digits = 4), "\n")
tm("fitted()", fitted(fit))
tm("residuals(pearson)", residuals(fit, type = "pearson"))
tm("autocor_matrix()", autocor_matrix(fit))
tm("simulate(nsim = 20)", simulate(fit, nsim = 20))
tm("summary()", summary(fit))
