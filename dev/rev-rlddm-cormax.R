# REVIEW round 2, item 1.0b: WHY `cor_max_abs` moves from 0.396 to
# 0.916 against a stated truth of 0.
#
#   Rscript dev/rev-rlddm-cormax.R <mode>
#
# mode = "truth"  arithmetic only, no fit: the drawn truths expressed in
#                 the parameterization the model now estimates
# mode = "fit"    the scale row refitted, saving the whole correlation
#                 matrix and the per-learner random effects
#
# Seed 20260908, the tier's own, 100 learners by 200 trials, the data
# function copied from
# extensions/frmtmb.learn/tests/testthat/test-scale.R at
# scale_small() FALSE.
#
# THE HYPOTHESIS UNDER TEST. `ndt` is now a FRACTION of each learner's
# own fastest response, and that floor is itself a draw from the
# learner's whole parameter vector: a slow learner (low drift, high
# boundary) has a later fastest response. So the deviation the model
# estimates on `ndt` is not the deviation the design drew; it is
# `qlogis(ndt_i / floor_i)`, and if `floor_i` moves with `drift_i` or
# `bs_i` then that deviation is correlated with them BY CONSTRUCTION,
# with a diagonal truth. If the truths reparameterized this way already
# carry a correlation near 0.92, the fitted 0.916 is recovery and not a
# defect.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1L]] else "truth"

.libPaths(c("C:/Users/adf44/source/r/rev-rlddm-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})
cat("mode :", mode, "\n")
cat("learn:", format(packageVersion("frmtmb.learn")), "at",
    dirname(system.file(package = "frmtmb.learn")), "\n")
cat("eam  :", format(packageVersion("frmtmb.eam")), "at",
    dirname(system.file(package = "frmtmb.eam")), "\n")

truth <- list(alpha = 0.35, drift = 2.5, bs = 1.5, ndt = 0.25,
              sd_alpha = 0.5, sd_drift = 1.0, sd_bs = 0.2,
              sd_ndt = 0.15)
seed <- 20260908L
ns <- 100L
nt <- 200L

d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                     seed = seed)
set.seed(seed + 2L)
i <- as.integer(d$id)
ua <- stats::rnorm(ns, 0, truth$sd_alpha)
ud <- stats::rnorm(ns, 0, truth$sd_drift)
ub <- stats::rnorm(ns, 0, truth$sd_bs)
un <- stats::rnorm(ns, 0, truth$sd_ndt)
s <- frm_task_simulate(
  rlddm(subject = id, trial = trial), d,
  pars = list(alpha = stats::plogis(stats::qlogis(truth$alpha) + ua[i]),
              drift = truth$drift + ud[i],
              bs = truth$bs * exp(ub[i]),
              ndt = truth$ndt * exp(un[i]),
              bias = 0.5),
  seed = seed)[[1L]]
lv <- levels(s$id)
own <- as.numeric(tapply(s$rt, s$id, min))[match(lv, lv)]
tru <- truth$ndt * exp(un)

# the deviations the DESIGN drew, on the links the row states truths on
dev_drawn <- cbind(alpha = ua, drift = ud, bs = ub, ndt_log = un)
# the deviation the model NOW estimates: ndt as a fraction of the
# learner's own floor, on a plain logit
frac <- tru / own
eta <- stats::qlogis(frac)
dev_fit <- cbind(alpha = ua, drift = ud, bs = ub,
                 ndt_frac = eta - mean(eta))

show <- function(lab, m) {
  cr <- stats::cor(m)
  cat("\n--", lab, "--\n")
  print(round(cr, 4))
  cat("max abs off-diagonal:",
      format(max(abs(cr[lower.tri(cr)])), digits = 6), "\n")
  invisible(cr)
}

cat("\nrows", nrow(s), " learners", ns, "\n")
show("the deviations the design drew", dev_drawn)
show("the same truths, in the parameterization the row now fits",
     dev_fit)

cat("\n-- the mechanism: what a learner's own floor is made of --\n")
for (nm in c("alpha", "drift", "bs", "ndt")) {
  v <- switch(nm, alpha = ua, drift = ud, bs = ub, ndt = un)
  cat(sprintf("  cor(log floor, %-6s deviation) = %+.6f\n", nm,
              stats::cor(log(own), v)))
}
cat(sprintf("  cor(log floor, log true ndt)     = %+.6f\n",
            stats::cor(log(own), log(tru))))
cat(sprintf("  cor(qlogis(frac), drift dev)     = %+.6f\n",
            stats::cor(eta, ud)))
cat(sprintf("  sd of the drawn ndt deviation, log scale  = %.6f\n",
            stats::sd(un)))
cat(sprintf("  sd of the fitted-parameterization deviation = %.6f\n",
            stats::sd(eta)))

# the OLD parameterization, for the same truths: a fraction of the
# GLOBAL fastest response
gmin <- min(s$rt)
old_frac <- tru / gmin
cat(sprintf("\n  learners whose truth is ABOVE the global floor %.6f: %d",
            gmin, sum(old_frac >= 1)))
cat(" of ", ns, "\n", sep = "")
cat("  (those have no representable value under one bound, which is\n")
cat("   why the 0.396 was measured on a fit that did not converge)\n")

if (!identical(mode, "fit")) quit(save = "no")

form <- bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(id) ~
             1 + (1 | p | id),
           drift ~ 1 + (1 | p | id), bs ~ 1 + (1 | p | id),
           ndt ~ 1 + (1 | p | id), bias = 0.5)
t0 <- Sys.time()
fit <- frm(form, family = rlddm(subject = id, trial = trial), data = s,
           se = TRUE)
cat("\nfit wall s:",
    format(as.numeric(difftime(Sys.time(), t0, units = "secs")),
           digits = 6), "\n")
vc <- VarCorr(fit)[[1L]]
cr <- stats::cov2cor(vc)
cat("\n-- the FITTED correlation matrix of the random-effect block --\n")
print(round(cr, 4))
cat("max abs off-diagonal:",
    format(max(abs(cr[lower.tri(cr)])), digits = 6), "\n")
j <- which(abs(cr) == max(abs(cr[lower.tri(cr)])), arr.ind = TRUE)
cat("the pair that carries it:", paste(rownames(cr)[j[1L, 1L]],
                                       colnames(cr)[j[1L, 2L]]), "\n")
re <- ranef(fit)
saveRDS(list(cor_fit = cr, vc = vc, ranef = re, dev_fit = dev_fit,
             dev_drawn = dev_drawn, own = own, tru = tru,
             logLik = as.numeric(stats::logLik(fit))),
        "dev/rev-rlddm-cormax.rds")
cat("wrote dev/rev-rlddm-cormax.rds\n")
