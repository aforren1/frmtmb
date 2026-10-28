# mointerop x sigmainn: the saved-fit refusal through marginaleffects,
# in full. Rscript dev/rel077-ms-probe.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r15", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(60), h = factor(rep(1:5, each = 12)), t = rep(1:12, 5))
d$y <- 1 + 0.5 * d$x + rnorm(60)
fit0 <- suppressWarnings(frm(y ~ x + ar(time = t, gr = h, cov = TRUE), data = d))
old <- fit0
old$frame$autocor[[1L]]$sigma_is <- NULL
nd <- data.frame(x = c(-1, 0, 1), t = 1, h = factor(1, levels = 1:5))
e <- tryCatch(suppressWarnings(marginaleffects::predictions(old, newdata = nd, dpar = "sigma")), error = function(e) e)
cat(conditionMessage(e), "\n")
cat("frmtmb's own words in it:", grepl("0.76.0 or earlier", conditionMessage(e)), "\n")
e2 <- tryCatch(marginaleffects::get_predict(old, newdata = nd, dpar = "sigma"), error = function(e) e)
cat("get_predict():", conditionMessage(e2), "\n")
