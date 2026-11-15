# Short check of parity2's fix round: does a fit's boot band at an
# OBSERVED level carry the level's conditional variance (what
# "old_levels" would lend)? Compare with the Wald band and fitted() at
# that level. Data as F1 (seed 100).
#   Rscript dev/reviews/parity2-short-obslevel.R [base] > <out> 2>&1
args <- commandArgs(TRUE)
base <- length(args) && identical(args[1], "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(100)
d <- data.frame(x = rnorm(60), g = factor(rep(1:30, each = 2)))
d$y <- 1 + 0.5 * d$x + rnorm(30)[d$g] + rnorm(60, sd = 2)
fb <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = d)
ic <- list(x = c(-1, 0, 1))
wid <- function(ce) ce$x$upper__ - ce$x$lower__
for (l in c("10", "3")) {
  cw <- conditional_effects(fb, "x", re_formula = NULL, int_conditions = ic,
                            conditions = data.frame(g = l))
  cb <- suppressMessages(conditional_effects(
    fb, "x", re_formula = NULL, int_conditions = ic, band = "boot",
    boot = 400, seed = 2, conditions = data.frame(g = l)))
  cat("level", l, "wald width", format(wid(cw), digits = 4),
      "boot 400 width", format(wid(cb), digits = 4), "\n")
}
cat("done\n")
