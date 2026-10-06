# Lane fixes, punch round 2: the Gamma("identity") update of test-nl.R
# without the pre-fit check, and the ported :955 update.
.libPaths(c("C:/Users/adf44/source/r/wt-fixes-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
set.seed(955)
n <- 60
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- exp((3 + 0.5 * d$x + 0.3 * d$z + rnorm(n, 0, 0.4)) / 3)
f0 <- frm(bf(y ~ 1 / (1 + exp(-a)) * exp(b * z), a ~ 1 + x, b ~ 1 + x,
             nl = TRUE), data = d, family = Gamma("identity"),
          start = list(beta = c(2, 0, 0.5, 0)))
r <- withCallingHandlers(tryCatch(
  suppressMessages(update(f0, formula. = bf(y ~ a + b, nl = TRUE))),
  error = function(e) conditionMessage(e)),
  warning = function(w) {cat("WARN:", conditionMessage(w), "\n");
    invokeRestart("muffleWarning")})
if (is.character(r)) cat("ERROR:", r, "\n") else {
  cat("conv", r$opt$convergence, "\n"); print(fixef(r)[, 1:2]) }
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
fit2 <- brms_fixture(2)
r <- withCallingHandlers(tryCatch(
  suppressMessages(update(fit2, formula. = bf(count ~ a + b, nl = TRUE))),
  error = function(e) conditionMessage(e)),
  warning = function(w) {cat("WARN:", substr(conditionMessage(w), 1, 300), "\n");
    invokeRestart("muffleWarning")})
cat(":955:", if (is.character(r)) paste("ERROR", substr(r, 1, 300)) else
  paste("fit, conv", r$opt$convergence), "\n")
