# Does brms 2.23.0 fit se() beside ar(cov = TRUE)? And frmtmb?
.libPaths(c("C:/Users/adf44/source/r/wt-sigmainn-lib", "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
d <- data.frame(t = rep(1:6, 10), g = rep(1:10, each = 6), x = rnorm(60), s = runif(60, 0.2, 0.5))
d$y <- rnorm(60)
code <- brms::stancode(brms::bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)), data = d)
cat(grep("target [+]=.*time", strsplit(code, "\n")[[1]], value = TRUE), sep = "\n")
r <- tryCatch(frm(bf(y | se(s, sigma = TRUE) ~ x + ar(t, g, cov = TRUE)), data = d), error = function(e) conditionMessage(e))
cat("frmtmb:", if (is.character(r)) r else "fits", "\n")
