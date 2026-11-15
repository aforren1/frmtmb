# Lane famend: d logLik / d nu at the stop of seed 9031 (3060), and along
# the ridge log(mu) + log(nu - 1) = const.
.libPaths(c("C:/Users/adf44/source/r/wt-famend-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
nu0 <- 1.1
set.seed(9000 + 31)
d <- data.frame(x = rnorm(300))
d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(300)^(-1 / nu0)
f <- suppressWarnings(frm(y ~ x, family = frechet(), data = d))
lp <- f$frame[["linpreds"]][[ns$linpred_key(f$spec$responses[[1]]$resp_name,
                                             "nu")]]
cat("opt message:", f$opt$message, " iterations", f$opt$iterations, "\n")
cat("par:", format(f$opt$par, digits = 10), "\n")
cat("gradient:", format(f$obj$gr(f$opt$par), digits = 4), "\n")
cat("par_units:", format(f$par_units), "\n")
cat("dll:", ns$gev_xi_slope(f, lp), "\n")
p <- f$opt$par
for (de in c(-2, -1, 0, 1, 2, 4)) {
  q <- p
  q[3] <- p[3] + de
  q[1] <- p[1] - de
  cat(sprintf("ridge step %+d: objective %.10f\n", de, f$obj$fn(q)))
}
# the profile in log(nu - 1), the other two maximized
for (e in c(-25, -19.138, -15, -10, -6, -4, -3)) {
  o <- function(b) f$obj$fn(c(b, e))
  m <- optim(p[1:2] + c(-(e - p[3]), 0), o, method = "BFGS")
  cat(sprintf("profile log(nu - 1) = %.3f: -logLik %.8f at mu-icpt %.4f\n",
              e, m$value, m$par[1]))
}
