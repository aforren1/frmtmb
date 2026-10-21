# Reviewer, fix round 1: the censored Gamma at shape 3000 now reports
# false convergence at the maximum. AD vs FD at the stop, and the fit's
# gradient there. ARM=lane or base-rev (the pre-fix lane build is gone,
# so only the lane).
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
for (shp in c(300, 3000, 3e4)) {
  set.seed(3)
  d <- data.frame(x = rnorm(200))
  d$y <- rgamma(200, shape = shp, rate = shp / exp(1 + 0.1 * d$x))
  qq <- quantile(d$y, 0.8)
  d$cc <- as.integer(d$y > qq)
  d$y[d$cc == 1] <- qq
  w <- character()
  fit <- withCallingHandlers(frm(y | cens(cc) ~ x, family = Gamma(link = "log"), data = d),
                             warning = function(x) { w <<- c(w, substr(conditionMessage(x), 1, 60)); invokeRestart("muffleWarning") })
  p <- fit$opt$par
  g <- as.numeric(fit$obj$gr(p))
  fd <- numDeriv::grad(fit$obj$fn, p)
  cat(sprintf("shape %g: conv %d msg '%s' | AD %s | FD %s\n   warnings: %s\n", shp, fit$opt$convergence,
              fit$opt$message, paste(signif(g, 4), collapse = ","), paste(signif(fd, 4), collapse = ","),
              paste(w, collapse = " || ")))
}
