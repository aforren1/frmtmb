# Review lane osa: cost of the "cdf" route on a larger design, poisson with
# (1 | g), 2000 rows in 100 groups of 20 (and 500 rows in 25 groups, for
# the scaling). One call per arm (the large calls are minutes long);
# "old" is the 0.79.1 call made on the same objective, as in
# dev/osa-time.R. Also a design with no group-level term.
#
#   Rscript dev/review-osa-time.R
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
new <- function(fit) residuals(fit, type = "osa")[, "Estimate"]
tm <- function(f) {
  t0 <- proc.time()[["elapsed"]]
  suppressWarnings(f())
  proc.time()[["elapsed"]] - t0
}
cat("BEGIN review-osa-time\n")
for (G in c(25L, 100L)) {
  set.seed(601)
  n <- G * 20L
  d <- data.frame(g = factor(rep(seq_len(G), each = 20)), x = rnorm(n))
  d$y <- rpois(n, exp(0.5 + 0.4 * d$x + rnorm(G, 0, 0.7)[d$g]))
  f0 <- suppressWarnings(frm(y ~ x, family = poisson(), data = d))
  fit_t <- tm(function() frm(y ~ x + (1 | g), family = poisson(), data = d))
  f1 <- suppressWarnings(frm(y ~ x + (1 | g), family = poisson(), data = d))
  t_cdf0 <- tm(function() new(f0))
  t_cdf <- tm(function() new(f1))
  t_old <- tm(function() old(f1))
  cat(sprintf(paste("rows %4d groups %3d | fit %.1fs | no group term: cdf",
                    "%.2fs | (1 | g): cdf %.1fs old %.1fs cdf/old %.2f\n"),
              n, G, fit_t, t_cdf0, t_cdf, t_old, t_cdf / t_old))
}
cat("END review-osa-time\n")
