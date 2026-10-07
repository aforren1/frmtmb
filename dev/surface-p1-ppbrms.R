# Lane surface, punch round 1, B2: pp_mixture() on a fit, Est.Error as
# the SD of the logit-normal law of its interval, against brms 2.23.0's
# pp_mixture() on the same data (the reviewer's fit,
# dev/surface-rev-out/brms-ppmix.rds: 4 chains x 1500 post-warmup
# draws, seed 1, brms's default priors) and against a Monte Carlo of
# the parameters' normal law. All 300 rows, unfiltered.
#
#   Rscript dev/surface-p1-ppbrms.R lane|base > dev/surface-out/p1-ppbrms-<arm>.txt
arm <- commandArgs(TRUE)[1]
source("dev/surface-env.R")
surface_env(arm)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "| frmtmb", find.package("frmtmb"), "\n")
b <- readRDS("dev/surface-rev-out/brms-ppmix.rds")
set.seed(4)
dm <- data.frame(y = c(rnorm(150, -1.5), rnorm(150, 1.5)))
stopifnot(identical(dm$y, b$data$y))
fmx <- frm(y ~ 1, family = mixture(gaussian(), gaussian()), data = dm)
pf <- pp_mixture(fmx)
pb <- b$pm
k <- if (cor(pb[, "Estimate", 1], pf[, "Estimate", 1]) < 0) 2 else 1
p <- pf[, "Estimate", 1]
bins <- cut(pmin(p, 1 - p), c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, 0.2, 0.5),
            right = FALSE)
rat <- pf[, "Est.Error", 1] / pb[, "Est.Error", k]
# the Monte Carlo of the parameters' normal law (the lane's original
# construction, all rows): 4000 draws, seed 1
ds <- frmtmb:::fit_draw_space(fmx)
v0 <- frmtmb:::fit_outer_vector(fmx, ds$map)
set.seed(1)
Z <- MASS::mvrnorm(4000, v0, ds$V)
Pm <- vapply(seq_len(4000), function(r) {
  mixture_probs(frmtmb:::fit_set_outer(fmx, Z[r, ], ds$map))[, 1]
}, numeric(nrow(dm)))
mc <- apply(Pm, 1, sd)
cat("\nEst.Error / brms posterior SD, and / Monte Carlo SD, by min(p, 1 - p):\n")
print(data.frame(n = tapply(rat, bins, length),
                 to_brms = signif(tapply(rat, bins, median), 3),
                 to_mc = signif(tapply(pf[, "Est.Error", 1] / mc, bins,
                                       median), 3),
                 brms_sd = signif(tapply(pb[, "Est.Error", k], bins,
                                         median), 3)))
cat("all 300 rows: median Est.Error / brms SD", signif(median(rat), 3),
    "| / MC SD", signif(median(pf[, "Est.Error", 1] / mc), 3), "\n")
mid <- pmin(p, 1 - p) >= 1e-3
cat("rows with p in [0.001, 0.999]:", sum(mid), "| Est.Error / brms SD",
    "range of the bin medians",
    paste(signif(range(tapply(rat[mid], droplevels(bins[mid]), median)), 3),
          collapse = " to "), "\n")
wf <- pf[, "Q97.5", 1] - pf[, "Q2.5", 1]
wb <- pb[, "Q97.5", k] - pb[, "Q2.5", k]
cat("interval width / brms's: median", signif(median(wf / wb), 3), "\n")
inside <- p >= pb[, "Q2.5", k] & p <= pb[, "Q97.5", k]
cat("estimates inside brms's 95% interval:", sum(inside), "of", length(p),
    "\n")
cat("Est.Error where brms's SD is above 0.002 / brms SD: median",
    signif(median(rat[pb[, "Est.Error", k] > 0.002]), 3), "\n")

# the integral itself, against stats::integrate() on the same law
lse <- (frmtmb:::fit_fd_se(fmx, mixture_probs)[, 1]) / (p * (1 - p))
check <- c(which.min(abs(p - 0.5)), which.min(p), which.max(p),
           which.min(abs(p - 0.01)), which.min(abs(p - 0.999)))
for (i in check) {
  s <- if (p[i] > 0.5) -1 else 1  # the complement, against cancellation
  f1 <- function(z) stats::plogis(s * z) * stats::dnorm(z, qlogis(p[i]), lse[i])
  f2 <- function(z) stats::plogis(s * z)^2 * stats::dnorm(z, qlogis(p[i]), lse[i])
  lo <- qlogis(p[i]) - 12 * lse[i]
  hi <- qlogis(p[i]) + 12 * lse[i]
  m1 <- integrate(f1, lo, hi, rel.tol = 1e-12)$value
  m2 <- integrate(f2, lo, hi, rel.tol = 1e-12)$value
  ref <- sqrt(max(m2 - m1^2, 0))
  cat(sprintf("p %.6g lse %.4g: Gauss-Hermite %.6g, integrate() %.6g, ratio %.8f\n",
              p[i], lse[i], pf[i, "Est.Error", 1], ref,
              pf[i, "Est.Error", 1] / ref))
}
