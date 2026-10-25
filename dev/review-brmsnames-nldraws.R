## fixef() row order: fit vs frmtmb.sample draws, nonlinear and
## distributional models.  Rscript dev/review-brmsnames-nldraws.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
cat("arm", arm, dirname(system.file(package = "frmtmb")), "\n")
set.seed(12)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- 1 + 0.5 * exp(d$z / 2) + (0.3 + 0.2 * d$z) * d$x +
  rnorm(n, 0, exp(0.1 * d$z))
for (f in list(bf(y ~ a + b * x, a ~ exp(z / 2), b ~ 1 + z, sigma ~ z,
                  nl = TRUE),
               bf(y ~ x, sigma ~ z))) {
  fit <- q(frm(f, data = d))
  ds <- q(frm_sample(fit, chains = 1, iter = 200, refresh = 0, seed = 2))
  cat(deparse1(f$formula), "\n  fit  :", rownames(fixef(fit)),
      "\n  draws:", rownames(fixef(ds)), "\n")
}
