# Lane optrap: test-open-issues.R's zero weights against the subset
# (lme4#880), both fits' units, starts and estimates.
a <- commandArgs(TRUE)
.libPaths(c(if (length(a) && nzchar(a[1])) a[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(1)
gm <- rnorm(10, sd = 5) + c(rep(0, 8), -50, 50)
dw <- data.frame(group = factor(rep(LETTERS[1:10], each = 3)), x = rnorm(30))
dw$y <- gm[as.integer(dw$group)] + 10 * dw$x + rnorm(30)
dw$w <- as.integer(dw$group %in% LETTERS[1:8])
sub <- frm(bf(y ~ x + (1 | group)) + gaussian(), data = subset(dw, w == 1))
wt <- frm(bf(y | weights(w) ~ x + (1 | group)) + gaussian(), data = dw)
for (f in list(sub, wt)) {
  cat("units", f$par_units, "| code", f$opt$convergence, f$opt$message,
      "| obj", format(f$opt$objective, digits = 15), "\n")
  cat("  par", format(f$opt$par, digits = 12), "\n")
}
cat("varcorr rel diff", varcorr_matrices(wt)[[1]][1, 1] /
      varcorr_matrices(sub)[[1]][1, 1] - 1, "\n")
