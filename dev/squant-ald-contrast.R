# Lane squant: the 5449 design (80 rows, groups of 4, seed 20261114, the
# test's data) under asym_laplace() on rellib-r22 and quantile_elf() on
# the lane build. Run: Rscript dev/squant-ald-contrast.R <lib>
.libPaths(c(commandArgs(TRUE)[1],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(20261114)
n <- 80
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / 4), each = 4)))
d$y <- 1 + 0.5 * d$x + rnorm(n / 4, 0, 0.7)[as.integer(d$g)] + rnorm(n)
r <- tryCatch(frm(bf(y ~ x + (1 | g), quantile = 0.5),
                  family = asym_laplace(), data = d),
              error = function(e) conditionMessage(e))
cat("asym_laplace:", if (is.character(r)) substr(r, 1, 200) else
  paste("fit, logLik", logLik(r)), "\n")
if (exists("quantile_elf")) {
  f <- frm(bf(y ~ x + (1 | g), quantile = 0.5),
           family = quantile_elf(smoothing = 0.5), data = d)
  cat("quantile_elf: conv", f$opt$convergence, "logLik",
      as.numeric(logLik(f)), "\n")
}
