# Punch round 1, m2: the residuals two tests of test-ordinal-mixture.R
# bound with absolute tolerances, measured on the lane build, so the
# bounds can be ratios to measured quantities. (1) The equidistant
# component's second difference of thresholds, against the size of the
# thresholds. (2) The two slopes of the shared-threshold mixture fitted
# from equal starts: identical, or how far apart.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
omx_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) +
    stats::rlogis(n) / exp(0.3 * d$z * cls)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d
}
d <- omx_data(20261008)
f3 <- frm(bf(y ~ x), family = mixture(cumulative(threshold = "equidistant"),
                                      acat(threshold = "equidistant")),
          data = d)
for (k in 1:2) {
  tk <- fixef(f3)[paste0("mu", k, "_Intercept[", 1:3, "]"), "Estimate"]
  r <- abs(diff(diff(tk)))
  cat(sprintf("equidistant mu%d: |second difference| = %.3g = %.3g eps * max|tau|\n",
              k, r, r / (.Machine$double.eps * max(abs(tk)))))
}
set.seed(20261014)
d <- data.frame(x = rnorm(400))
cls <- stats::rbinom(400, 1, 0.4)
lat <- ifelse(cls == 1, 2 * d$x, -1.5 * d$x) + stats::rlogis(400)
d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
fam <- mixture(cumulative(), cumulative(), order = "mu")
stuck <- suppressWarnings(frm(bf(y ~ x), family = fam, data = d,
                              start = list(beta = c(0, 0))))
bs <- stuck$estimates$beta
cat(sprintf("stuck slopes %s and %s, identical %s, theta %s\n",
            format(bs[[1]], digits = 17), format(bs[[2]], digits = 17),
            identical(bs[[1]], bs[[2]]),
            paste(format(unlist(stuck$estimates[grep("theta|mix",
                  names(stuck$estimates))]), digits = 17), collapse = " ")))
