# Lane postfix, fix round 2: the default-start error of the review's
# sinusoid (P 3, seed 1001, dev/postfix-rev2-sinus.R), printed whole.
.libPaths(c("C:/Users/adf44/source/r/wt-postfix-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(1001)
xs <- seq(0, 15, length.out = 120)
ds <- data.frame(x = xs, y = 2 * sin(6.283185307179586 * xs / 3) +
                   rnorm(120, 0, 0.6))
f1 <- bf(y ~ A * sin(6.283185307179586 * x / P), A + P ~ 1, nl = TRUE)
cat(tryCatch(frm(f1, data = ds), error = conditionMessage), "\n")
