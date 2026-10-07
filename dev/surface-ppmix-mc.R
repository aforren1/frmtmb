# Lane surface: pp_mixture() on a fit, its delta-method Est.Error
# against a Monte Carlo of the same asymptotic law (parameters drawn
# from N(estimate, vcov(full = TRUE)), the probabilities recomputed at
# each draw). Same data as dev/surface-check.R (seed 4).
#
#   Rscript dev/surface-ppmix-mc.R > dev/surface-out/ppmix-mc.txt
source("dev/surface-env.R")
surface_env("lane")
suppressPackageStartupMessages(library(frmtmb))
set.seed(4)
dm <- data.frame(y = c(rnorm(150, -1.5), rnorm(150, 1.5)))
fmx <- frm(y ~ 1, family = mixture(gaussian(), gaussian()), data = dm)
pf <- pp_mixture(fmx)
ds <- frmtmb:::fit_draw_space(fmx)
v0 <- frmtmb:::fit_outer_vector(fmx, ds$map)
set.seed(1)
R <- 4000
Z <- MASS::mvrnorm(R, v0, ds$V)
P <- vapply(seq_len(R), function(r) {
  mixture_probs(frmtmb:::fit_set_outer(fmx, Z[r, ], ds$map))[, 1]
}, numeric(nrow(dm)))
mc_sd <- apply(P, 1, sd)
ratio <- pf[, "Est.Error", 1] / mc_sd
cat("outer parameters:", names(v0), "\n")
cat("delta Est.Error / Monte Carlo SD (R = 4000, seed 1), rows with MC SD",
    "> 0.01:\n")
print(summary(ratio[mc_sd > 0.01]))
cat("the same, rows with MC SD > 0.05:\n")
print(summary(ratio[mc_sd > 0.05]))
q <- apply(P, 1, quantile, c(0.025, 0.975))
cover <- mean(pf[, "Q2.5", 1] <= q[2, ] & pf[, "Q97.5", 1] >= q[1, ])
cat("rows whose logit-Wald interval overlaps the MC 95% range:", cover, "\n")
w_wald <- pf[, "Q97.5", 1] - pf[, "Q2.5", 1]
w_mc <- q[2, ] - q[1, ]
cat("interval width, Wald / MC, rows with MC width > 0.02:\n")
print(summary((w_wald / w_mc)[w_mc > 0.02]))
