# Review probe, arm A: fits made on 0.83.0 (rellib-r22) and saved, for
# arm B to load on the lane's core.
#   Rscript dev/reviews/cdfname-review-rds-a.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
set.seed(3088)
n <- 120
x <- rnorm(n)
d <- data.frame(x = x, y = rweibull(n, 2, exp(0.5 + 0.3 * x)))
d$cc <- ifelse(d$y > 2, "right", "none")
d$y <- pmin(d$y, 2)
dp <- data.frame(x = x, yc = rpois(n, exp(1 + 0.2 * x)))
dp$cc <- ifelse(dp$yc > 4, "right", "none")
dp$yc <- pmin(dp$yc, 4)
pois_lcdf <- frmtmb_family(
  "pois_old", dpars = "mu", links = list(mu = "log"),
  type = "discrete", accepts_aterms = c("weights", "cens", "trunc"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dpois(y, dpars[["mu"]], log = TRUE)
  },
  lcdf = function(q, dpars, aterms) RTMB::ppois(q, dpars[["mu"]]),
  init_dpars = list(mu = function(y, aterms) mean(y)))
fits <- list(
  weibull = frm(y | cens(cc) ~ x, family = weibull(), data = d),
  poisson = frm(yc | cens(cc) ~ x, family = poisson(), data = dp),
  custom = frm(yc | cens(cc) ~ x, family = pois_lcdf, data = dp),
  weibull_plain = frm(y ~ x, family = weibull(), data = d))
for (nm in names(fits)) {
  cat(sprintf("A %-14s logLik=%.10f has_lcdf_elt=%s\n", nm,
              as.numeric(logLik(fits[[nm]])),
              "lcdf" %in% names(fits[[nm]]$family %||% list())))
}
saveRDS(list(fits = fits, d = d, dp = dp),
        "dev/reviews/cdfname-review-old-fits.rds")
cat("saved\n")
