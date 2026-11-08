# Review lane osa: where the "cdf" route's cost with (1 | g) comes from.
# 500 rows, 25 groups of 20, one call per arm, bernoulli and negbinomial
# beside dev/review-osa-time.R's poisson: a family whose CDF is cheap on
# the tape against one whose CDF is a sum.
#   Rscript dev/review-osa-time2.R
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
old <- function(fit) {
  RTMB::oneStepPredict(fit$obj, observation.name = ".frm_obs",
                       method = "oneStepGeneric", discrete = TRUE,
                       range = c(0, Inf), trace = FALSE)$residual
}
tm <- function(f) {
  t0 <- proc.time()[["elapsed"]]
  suppressWarnings(f())
  proc.time()[["elapsed"]] - t0
}
cat("BEGIN review-osa-time2\n")
G <- 25L
set.seed(601)
n <- G * 20L
d <- data.frame(g = factor(rep(seq_len(G), each = 20)), x = rnorm(n))
e <- rnorm(G, 0, 0.7)[d$g]
d$yp <- rpois(n, exp(0.5 + 0.4 * d$x + e))
d$yb <- rbinom(n, 1, plogis(0.5 + 0.4 * d$x + e))
d$yn <- rnbinom(n, size = 3, mu = exp(0.5 + 0.4 * d$x + e))
for (cs in list(list("bernoulli", yb ~ x + (1 | g), bernoulli()),
                list("negbinomial", yn ~ x + (1 | g), negbinomial()))) {
  f1 <- suppressWarnings(frm(cs[[2]], family = cs[[3]], data = d))
  t_cdf <- tm(function() residuals(f1, type = "osa"))
  t_old <- tm(function() old(f1))
  cat(sprintf("%-12s rows %d (1 | g): cdf %.1fs old %.1fs cdf/old %.2f\n",
              cs[[1]], n, t_cdf, t_old, t_cdf / t_old))
}
cat("END review-osa-time2\n")
