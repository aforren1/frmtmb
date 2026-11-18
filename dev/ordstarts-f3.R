# Lane ordstarts, build round: the F3 test case of test-ordinal-mixture.R
# (sratio + acat with cs(x), the review sweep's cs_sr_acat seed 1) under
# arm lane or base; prints the warnings, the starts record and the
# threshold standard errors.
#   Rscript dev/ordstarts-f3.R <lane|base> [starts]
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- stats::rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
d$y <- 1L + rowSums(outer(lat, stats::quantile(lat, 1:3 / 4), ">"))
ctl <- if (length(a) > 1) frmtmb_control(starts = as.integer(a[2])) else
  frmtmb_control()
w <- character()
m <- character()
fit <- withCallingHandlers(
  frm(bf(y ~ cs(x) + z), family = mixture(sratio(), acat()), data = d,
      control = ctl),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, substr(conditionMessage(x), 1, 200))
    invokeRestart("muffleMessage")
  })
cat("logLik", as.numeric(logLik(fit)), "code", fit$opt$convergence, "\n")
print(fit$opt$mixture_starts)
cat("WARNINGS:\n", paste(substr(w, 1, 300), collapse = "\n"), "\n")
cat("MESSAGES:\n", paste(m, collapse = "\n"), "\n")
fx <- suppressWarnings(fixef(fit))
print(fx[grep("^mu2_Intercept", rownames(fx)), ])
