# PUNCH ROUND 2, N13. The truth of the correlated block the grouped
# rlddm row actually estimates.
#
#   Rscript dev/rlddm-scripts/rlddm-cortrue.R <lib>
#
# Seed 20260908, the tier's own design and truths, 100 learners by 200
# trials. NO FIT: this is arithmetic on the drawn deviations and the
# drawn data's own floors, which is the point. If the truth of the block
# can be computed without fitting, then recording it as zero is a
# statement the design does not support.
#
# The design draws four independent deviations. Three of them are the
# ones the model estimates: `alpha` on the logit, `drift` on the
# identity, `bs` on the log. The fourth is not. Under a per-group bound
# the model estimates qlogis(ndt_i / floor_i), and floor_i is that
# learner's own fastest response, a draw from its WHOLE parameter
# vector.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})
cat("lib:", lib, "\n\n")

seed <- 20260908L
ns <- 100L
nt <- 200L
tr <- list(alpha = 0.35, drift = 2.5, bs = 1.5, ndt = 0.25,
           sd_alpha = 0.5, sd_drift = 1.0, sd_bs = 0.2, sd_ndt = 0.15)

d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                     seed = seed)
set.seed(seed + 2L)
i <- as.integer(d$id)
ua <- stats::rnorm(ns, 0, tr$sd_alpha)
ud <- stats::rnorm(ns, 0, tr$sd_drift)
ub <- stats::rnorm(ns, 0, tr$sd_bs)
un <- stats::rnorm(ns, 0, tr$sd_ndt)
s <- frm_task_simulate(
  rlddm(subject = id, trial = trial), d,
  pars = list(alpha = stats::plogis(stats::qlogis(tr$alpha) + ua[i]),
              drift = tr$drift + ud[i],
              bs = tr$bs * exp(ub[i]),
              ndt = tr$ndt * exp(un[i]),
              bias = 0.5),
  seed = seed)[[1L]]

ndt_true <- tr$ndt * exp(un)
floor <- as.numeric(tapply(s$rt, s$id, min))
cat("rows:", nrow(s), " learners:", ns, "\n")
cat("mean(ndt truth):", format(mean(ndt_true), digits = 9), "\n")
cat("mean(own floor):", format(mean(floor), digits = 9), "\n\n")

drawn <- cbind(alpha = ua, drift = ud, bs = ub, ndt = un)
fitted_par <- cbind(alpha = ua, drift = ud, bs = ub,
                    ndt = stats::qlogis(ndt_true / floor))

mx <- function(m) {
  cr <- stats::cor(m)
  max(abs(cr[lower.tri(cr)]))
}
show <- function(m, lbl) {
  cr <- stats::cor(m)
  cat(lbl, "\n")
  for (a in seq_len(ncol(m) - 1L)) {
    for (b in (a + 1L):ncol(m)) {
      cat(sprintf("  %-6s %-6s %+.4f\n", colnames(m)[a], colnames(m)[b],
                  cr[b, a]))
    }
  }
  cat(sprintf("  max abs off-diagonal: %.4f\n\n", mx(m)))
}
show(drawn, "-- the four deviations AS DRAWN --")
show(fitted_par, "-- the same truths in the parameterization now fitted --")

cat("-- the mechanism --\n")
cat(sprintf("  cor(log floor, bs deviation) : %+.4f\n",
            stats::cor(log(floor), ub)))
cat(sprintf("  cor(log floor, ndt deviation): %+.4f\n",
            stats::cor(log(floor), un)))
cat("\ncor_true for the learn-rlddm row:",
    sprintf("%.6f", mx(fitted_par)), "\n")
cat("what the row recorded before punch round 2: 0\n")
