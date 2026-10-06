# The coefficient variance of the A grid of dev/cifix-scan2.R's model,
# as frm_lp_basis() returns it, on each fit named, one fit per process
# call. Usage: Rscript dev/cifix-qrange.R <lib or "base"> <n> <seed> ...
args <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r6"
user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(args[1], "base")) c(base, user) else
  c(args[1], base, user))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
cases <- matrix(as.integer(args[-1]), ncol = 2, byrow = TRUE)
for (k in seq_len(nrow(cases))) {
  n <- cases[k, 1]; s <- cases[k, 2]
  set.seed(s)
  d <- data.frame(x = sort(stats::runif(n, 0, 6)),
                  fac = factor(rep(c("A", "B"), length.out = n)))
  d$y <- sin(d$x) + ifelse(d$fac == "B", 0.3 * d$x, 0) +
    stats::rnorm(n, 0, 0.3)
  fit <- suppressWarnings(
    frm(bf(y ~ fac + s(x, by = fac, k = 6) + gp(x)), data = d))
  gx <- d$x[-1] - diff(d$x) / 2
  A <- data.frame(x = gx, fac = factor("A", levels = levels(d$fac)))
  lb <- frm_lp_basis(fit, newdata = A)
  C <- as.matrix(lb$A)
  q <- rowSums((C %*% lb$V) * C)
  se <- frm_linpred(fit, newdata = A, se.fit = TRUE)$se.fit
  cat(sprintf("n = %d seed = %d theta1 = %.4f: q in [%.4g, %.4g],",
              n, s, fit$estimates$theta[1], min(q), max(q)),
      sprintf("se.fit in [%.4g, %.4g]\n", min(se), max(se)))
}
