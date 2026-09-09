# ndt lane: is the recovered between-subject spread in `ndt` the random
# effect, or is it the FLOORS?
#
# ndt_i = fraction_i * floor_i, and floor_i tracks ndt_i because a
# slower non-decision time pushes the whole response-time distribution
# right. So a fit could "recover sd(ndt)" while the random effect did
# nothing. Two designs settle it.
#
#   A. ndt varies between subjects (the scale design's shape).
#   B. ndt is CONSTANT and the boundary separation varies instead, so
#      the floors vary and the truth does not. A fit that reports a
#      between-subject spread in `ndt` here is reporting the floors.
#
# 15 subjects x 250 trials, seeds 771 + 0:7 for each design.
# Run: Rscript --vanilla dev/ndt-scripts/ndt-floor-artifact.R
.libPaths(c("C:/Users/adf44/source/r/ndt-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

ns <- 15L
nt <- 250L

one <- function(seed, sd_lndt, sd_lbs) {
  set.seed(seed)
  u_nd <- rnorm(ns, 0, sd_lndt)
  u_bs <- rnorm(ns, 0, sd_lbs)
  s <- rep(seq_len(ns), each = nt)
  t0 <- 0.25 * exp(u_nd)
  d <- ddm_simulate(ns * nt, mu = 1.1, bs = 1.4 * exp(u_bs),
                    ndt = t0[s], bias = 0.5)
  d$s <- factor(s)
  fit <- frm(bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1 + (1 | s),
                ndt ~ 1 + (1 | s), bias = 0.5),
             family = wiener(), data = d)
  key <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]
  fl <- as.numeric(tapply(d$rt, d$s, min))
  hat <- as.numeric(ndt_time(fit, newdata = key))
  frac_pop <- as.numeric(predict(fit, newdata = key[1L, ], dpar = "ndt",
                                 type = "response", re.form = NA))
  c(conv = fit$opt$convergence,
    sd_true = sd(t0), sd_hat = sd(hat), sd_floor_only = sd(frac_pop * fl),
    sd_floor = sd(fl),
    cor_hat_true = stats::cor(hat, t0),
    cor_floor_true = stats::cor(fl, t0),
    bias_ms = 1000 * mean(hat - t0),
    rmse_ms = 1000 * sqrt(mean((hat - t0)^2)),
    rmse_floor_only_ms = 1000 * sqrt(mean((frac_pop * fl - t0)^2)))
}

for (design in c("A: ndt varies", "B: ndt constant, bs varies")) {
  sd_lndt <- if (startsWith(design, "A")) 0.12 else 0
  sd_lbs <- if (startsWith(design, "A")) 0 else 0.25
  m <- t(vapply(0:7, function(k) one(771L + k, sd_lndt, sd_lbs),
                numeric(10)))
  cat("\n==", design, " (", ns, "x", nt, ", 8 seeds ) ==\n")
  cat("converged:", sum(m[, "conv"] == 0), "of", nrow(m), "\n")
  for (nm in colnames(m)[-1]) {
    cat(sprintf("  %-20s mean %10.5f  sd %10.5f  range %10.5f %10.5f\n",
                nm, mean(m[, nm]), sd(m[, nm]), min(m[, nm]),
                max(m[, nm])))
  }
}
