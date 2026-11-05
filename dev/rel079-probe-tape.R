# The MakeADFun calls that carry no gp_tape_scope(): can each tape an exact gp()?
# R/fit.R:1539 is the Laplace tape behind integrate = and importance =,
# R/fit.R:2018 quadrature, R/importance.R:1812 the saved proposal. A call that
# stops before taping has nothing to scope.
.libPaths(c("C:/Users/adf44/source/r/rellib-r17",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
set.seed(1)
d <- data.frame(x = seq(0, 6, length.out = 25))
d$y <- sin(d$x) + rnorm(25, 0, 0.3)
d$k <- rpois(25, exp(1 + sin(d$x)))
d$g <- factor(rep(1:5, 5))
tt <- function(lab, e) {
  r <- tryCatch(withCallingHandlers(e, warning = function(w) invokeRestart("muffleWarning")),
                error = function(e) e)
  cat(sprintf("%-28s %s\n", lab, if (inherits(r, "error")) {
    paste("STOPS:", substr(conditionMessage(r), 1, 130))
  } else {
    paste("FITS, logLik", format(as.numeric(logLik(r)), digits = 8))
  }))
}
tt("plain", frm(y ~ gp(x), data = d))
tt("importance", frm(k ~ gp(x), family = poisson(), data = d, importance = 100))
tt("integrate", frm(k ~ gp(x), family = poisson(), data = d, integrate = 5))
tt("quadrature", frm(k ~ gp(x), family = poisson(), data = d, quadrature = TRUE))
tt("mi retape", {
  d2 <- d
  d2$z <- d2$y + rnorm(25)
  d2$y[3:5] <- NA
  frm(bf(y | mi() ~ gp(x)) + bf(z ~ 1 + x) + gaussian(), data = d2)
})
