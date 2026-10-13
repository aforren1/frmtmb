# Reviewer of lane optrap, item 1: is the fit of y * 2^k the fit of y
# scaled, bit for bit, on one build in one process? ri20 design
# (dev/optrap-scale.R), seeds 1 to 10, k in -10, 10, 20.
#   Rscript dev/optrap-rev-pow2.R <lib or none> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
lab <- a[2]
cat(lab, "lib", dirname(find.package("frmtmb")), "\n")
X <- list()
for (s in 1:10) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.5)[d$g] + rnorm(100)
  f1 <- suppressWarnings(frm(y ~ x + (1 | g), data = d))
  for (k in c(-10, 10, 20)) {
    dk <- d
    dk$y <- d$y * 2^k
    fk <- suppressWarnings(frm(y ~ x + (1 | g), data = dk))
    b <- fk$estimates$beta / 2^k
    X[[length(X) + 1L]] <- data.frame(
      seed = s, k = k,
      beta_identical = identical(unname(b), unname(f1$estimates$beta)),
      beta_maxrel = max(abs(b / f1$estimates$beta - 1)),
      theta_diff = fk$estimates$theta[1] - k * log(2) -
        f1$estimates$theta[1],
      ll_diff = as.numeric(logLik(fk)) + 100 * k * log(2) -
        as.numeric(logLik(f1)),
      evals_same = identical(fk$opt$evaluations, f1$opt$evaluations))
  }
}
X <- do.call(rbind, X)
print(X)
cat("beta bitwise identical:", sum(X$beta_identical), "of", nrow(X),
    "| same evaluation counts:", sum(X$evals_same), "| max beta rel diff",
    signif(max(X$beta_maxrel), 3), "| max |ll diff|",
    signif(max(abs(X$ll_diff)), 3), "\n")
