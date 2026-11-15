# Lane famend, fix round F6: every warning (class and text) the fits of
# test-family-fitend.R raise, under the BLAS of the R that runs this.
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("R", R.home(), "| LAPACK", La_library(), "\n")
show <- function(lab, expr) {
  w <- list()
  withCallingHandlers(expr, warning = function(x) {
    w[[length(w) + 1L]] <<- x; invokeRestart("muffleWarning")
  }, message = function(x) invokeRestart("muffleMessage"))
  cat(lab, ":", length(w), "warnings\n")
  for (x in w) cat("   [", class(x)[1], "]", substr(conditionMessage(x), 1, 90), "\n")
}
fd <- function(seed, n = 300, nu0 = 1.1) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(n)^(-1 / nu0)
  d
}
mk <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
show("frechet 9031", frm(y ~ x, family = frechet(), data = fd(9031)))
show("frechet 9001", frm(y ~ x, family = frechet(), data = fd(9001)))
show("cs 3", frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), sratio()),
                 data = mk(3)))
show("cs 4", frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), sratio()),
                 data = mk(4)))
set.seed(4)
d <- data.frame(x = rnorm(300))
d$y <- rpois(300, exp(1 + 0.3 * d$x))
for (k in c(0, 9)) {
  dk <- d
  dk$x <- dk$x * (1 + 2e-15 * k)
  show(paste("nb k", k), frm(y ~ x, family = negbinomial(), data = dk))
  show(paste("nb2 k", k), frm(y ~ x, family = brmsfamily("negbinomial2"),
                              data = dk))
}
