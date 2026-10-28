# Short check: would `.get_sigma.frmtmb_fit` calling autocor_check_fit()
# then NextMethod() keep insight's default for fresh fits and refuse a
# saved 0.76.0 fit? Registered in this session only.
.libPaths(c("C:/Users/adf44/source/r/wt-sigmainn-lib", "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
registerS3method(".get_sigma", "frmtmb_fit", function(x, verbose = TRUE, ...) {
  frmtmb:::autocor_check_fit(x)
  NextMethod()
}, envir = asNamespace("insight"))
set.seed(3)
d <- data.frame(x = rnorm(120), x2 = rnorm(120))
d$y <- 1 + d$x + rnorm(120, sd = exp(0.2 * d$x2))
f <- frm(bf(y ~ x, sigma ~ x2), data = d)
cat("sigma ~ x2:", suppressWarnings(as.numeric(insight::get_sigma(f))), "\n")
cat("gaussian:", as.numeric(insight::get_sigma(frm(y ~ x, data = d))), "\n")
old <- readRDS("dev/review-sigmainn-out/saved-r14-fits.rds")$covT
cat("saved covT:", tryCatch(insight::get_sigma(old), error = function(e) substr(conditionMessage(e), 1, 60)), "\n")
