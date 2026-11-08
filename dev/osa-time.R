# Lane osa: the cost of osa_method = "cdf" against the route it replaced,
# in one process of the lane build.
#
#   Rscript dev/osa-time.R
#
# "old" is the call 0.79.1 made for a count, RTMB::oneStepPredict(method =
# "oneStepGeneric", discrete = TRUE, range = c(0, Inf)) on the fit's own
# objective; "cdf" is residuals(type = "osa") now. "control" runs the cdf
# route a second time and must read 1.0. The arms are interleaved, each
# round's time is the elapsed time of one call, and the reported time is
# the minimum over 5 rounds. Designs: seed 11, 20 groups of 10 (as in
# dev/osa-re.R), poisson and bernoulli, with and without (1 | g).
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(11)
n <- 200
d <- data.frame(g = factor(rep(1:20, each = 10)), x = rnorm(n))
e <- rnorm(20, 0, 0.7)[d$g]
d$yp <- rpois(n, exp(0.2 + 0.5 * d$x + e))
d$yb <- rbinom(n, 1, plogis(0.2 + 0.5 * d$x + e))
old <- function(fit) {
  RTMB::oneStepPredict(fit$obj, observation.name = ".frm_obs",
                       method = "oneStepGeneric", discrete = TRUE,
                       range = c(0, Inf), trace = FALSE)$residual
}
new <- function(fit) residuals(fit, type = "osa")[, "Estimate"]
tm <- function(f) {
  t0 <- proc.time()[["elapsed"]]
  suppressWarnings(f())
  proc.time()[["elapsed"]] - t0
}
cases <- list(
  poisson = list(yp ~ x, poisson()),
  poisson_re = list(yp ~ x + (1 | g), poisson()),
  bernoulli = list(yb ~ x, bernoulli()),
  bernoulli_re = list(yb ~ x + (1 | g), bernoulli()))
cat("BEGIN osa-time\n")
for (nm in names(cases)) {
  fit <- suppressWarnings(frm(cases[[nm]][[1]], family = cases[[nm]][[2]],
                              data = d))
  t <- matrix(NA, 5, 3, dimnames = list(NULL, c("old", "cdf", "control")))
  for (k in 1:5) {
    t[k, "old"] <- tm(function() old(fit))
    t[k, "cdf"] <- tm(function() new(fit))
    t[k, "control"] <- tm(function() new(fit))
  }
  m <- apply(t, 2, min)
  cat(sprintf("%-13s old %6.2fs  cdf %6.2fs  control %6.2fs  cdf/old %5.2f  control/cdf %5.2f\n",
              nm, m["old"], m["cdf"], m["control"], m["cdf"] / m["old"],
              m["control"] / m["cdf"]))
}
cat("END osa-time\n")
