# Review H: G3's seed 5 (group log sd -3.04, draw sd 6.84, flat side
# 0.003, not held): predict() at a NEW level, against the plug-in.
# Usage: Rscript ppdraw-review-h.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (length(args) && args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("C:/Users/adf44/source/r/frmtmb-wt-ppdraw/dev/reviews/ppdraw-review-common.R")
set.seed(5)
d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, 0.1)[d$g] + stats::rnorm(120)
fit <- q(frm(bf(y ~ x + (1 | g)), data = d))
nd <- data.frame(x = 0, g = factor("new"))
for (pe in c(TRUE, FALSE)) {
  set.seed(1)
  p <- q(predict(fit, newdata = nd, ndraws = 4000, summary = FALSE,
                 allow_new_levels = TRUE, propagate_error = pe))
  cat(sprintf(paste("propagate_error %-5s sd %.4g  q(0.01, 0.5, 0.99)",
                    "%s  share |y - 1| > 10: %.4f\n"), pe, stats::sd(p),
              paste(format(stats::quantile(p, c(0.01, 0.5, 0.99)),
                           digits = 4), collapse = " "),
              mean(abs(p - 1) > 10)))
}
