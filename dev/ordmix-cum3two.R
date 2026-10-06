# Three cumulative components on two-class data (seed 20261016, the
# lpcheck case cum3 before it got three classes): on the first lane
# build the optimizer stopped at "NA/NaN gradient evaluation" with
# component 2's increment at exp(-37.66); the gap floor of
# ord_log_interior() is meant to let it finish.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261005 + 11)
n <- 400
x <- rnorm(n); z <- rnorm(n)
g <- factor(sample(c("a", "b"), n, TRUE))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
d <- data.frame(y, x)
f <- withCallingHandlers(
  tryCatch(frm(bf(y ~ x), family = mixture(cumulative(), cumulative(),
                                           cumulative()), data = d),
           error = function(e) e),
  warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })
if (inherits(f, "error")) cat("ERROR:", conditionMessage(f), "\n") else {
  cat("logLik", as.numeric(logLik(f)), "\n")
  print(round(f$opt$par, 3))
  f2 <- frm(bf(y ~ x), family = mixture(cumulative(), cumulative()), data = d)
  cat("two components: logLik", as.numeric(logLik(f2)), "\n")
}
