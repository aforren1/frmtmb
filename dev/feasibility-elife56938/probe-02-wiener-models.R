# Probe 2: the three analytically-tractable models from the paper, through frm().
#
# M1 (Eq. 14)  3-parameter DDM      drift proportional to coherence
# M2 (Eq. 16)  8-parameter DDM      free drift per coherence level
# M4 (Eq. 17)  18-parameter DDM     free drift, bound and non-decision time per level
#
# Part A is a recovery study over repeated simulations at Roitman dimensions; a single
# replicate cannot distinguish a biased estimator from an unlucky draw.
# Part B fits the real Roitman & Shadlen monkey 1 data that ships with PyDDM.
#
# PyDDM parameterizes bounds as +/-B with an unbiased start. wiener() uses brms
# conventions: boundary SEPARATION bs = 2B, relative start bias = 0.5. Drift is the same
# quantity in both. fixef() returns a list keyed by dpar; ndt comes back on a scaled
# logit whose ceiling is min(y), so back-transform with plogis(eta) * min(y).

sp <- Sys.getenv("EL_SCRATCH")
.libPaths(c(file.path(sp, "el-lib"), .libPaths()))
suppressPackageStartupMessages({
  library(frmtmb); library(frmtmb.ddm); library(RWiener)
})
cat("frmtmb", as.character(packageVersion("frmtmb")),
    "/ frmtmb.ddm", as.character(packageVersion("frmtmb.ddm")), "\n\n")

COH <- c(0, 0.032, 0.064, 0.128, 0.256, 0.512)

sim_ddm <- function(n_per, mu0, bs, ndt, coh = COH) {
  do.call(rbind, lapply(coh, function(cc) {
    d <- RWiener::rwiener(n = n_per, alpha = bs, tau = ndt, beta = 0.5, delta = mu0 * cc)
    data.frame(rt = d$q, upper = as.integer(d$resp == "upper"),
               coh = cc, cohf = factor(cc, levels = coh))
  }))
}

natural <- function(fit, dat) {
  fe <- fixef(fit)
  c(mu0 = unname(fe$mu[1]),
    bs = unname(exp(fe$bs[1])),
    ndt = unname(plogis(fe$ndt[1]) * min(dat$rt)))
}

# ---- Part A: recovery study --------------------------------------------------------
truth <- c(mu0 = 10.5, bs = 1.6, ndt = 0.28)
NREP <- 25
cat("=== M1 recovery,", NREP, "replicates of", 440 * length(COH), "trials ===\n")
t0 <- proc.time()[3]
rec <- t(vapply(seq_len(NREP), function(r) {
  set.seed(1000 + r)
  d <- sim_ddm(440, truth[["mu0"]], truth[["bs"]], truth[["ndt"]])
  f <- frm(bf(rt | vint(upper) ~ 0 + coh, bias = 0.5), family = wiener(), data = d)
  natural(f, d)
}, numeric(3)))
cat("elapsed", round(proc.time()[3] - t0, 1), "s (", round((proc.time()[3] - t0) / NREP, 2),
    "s per fit )\n\n")
summ <- data.frame(
  truth = truth,
  mean_est = colMeans(rec),
  bias = colMeans(rec) - truth,
  rel_bias_pct = 100 * (colMeans(rec) - truth) / truth,
  mc_se = apply(rec, 2, sd) / sqrt(NREP),
  sd_across_reps = apply(rec, 2, sd))
print(round(summ, 4))
cat("\nbias in units of its own Monte Carlo standard error:\n")
print(round(summ$bias / summ$mc_se, 2))

# ---- one replicate, all three models -----------------------------------------------
set.seed(1001)
dat <- sim_ddm(440, truth[["mu0"]], truth[["bs"]], truth[["ndt"]])
cat("\n=== all three models on one replicate (", nrow(dat), "trials ) ===\n")

m1 <- frm(bf(rt | vint(upper) ~ 0 + coh, bias = 0.5), family = wiener(), data = dat)
m2 <- frm(bf(rt | vint(upper) ~ 0 + cohf, bias = 0.5), family = wiener(), data = dat)
m4 <- frm(bf(rt | vint(upper) ~ 0 + cohf, bs ~ 0 + cohf, ndt ~ 0 + cohf, bias = 0.5),
          family = wiener(), data = dat)

cat("\nM1 natural scale:\n"); print(round(natural(m1, dat), 4))
cat("\nM2 drift per coherence vs truth (mu0 * C):\n")
print(round(data.frame(coh = COH, est = unname(fixef(m2)$mu),
                       truth = truth[["mu0"]] * COH), 3))
cat("\nM4 free parameters:", length(coef(m4)), "\n")

cat("\nmodel comparison (data generated under M1, so BIC should prefer M1):\n")
print(data.frame(model = c("M1 (3 par)", "M2 (8 par)", "M4 (18 par)"),
                 npar = c(length(coef(m1)), length(coef(m2)), length(coef(m4))),
                 logLik = round(c(logLik(m1), logLik(m2), logLik(m4)), 2),
                 BIC = round(c(BIC(m1), BIC(m2), BIC(m4)), 2)))

# ---- Part B: the real Roitman & Shadlen monkey 1 data ------------------------------
csv <- file.path(sp, "el-pyddm", "roitman_rts.csv")
if (file.exists(csv)) {
  cat("\n\n================ real data: Roitman & Shadlen monkey 1 ================\n")
  rd <- read.csv(csv)
  rd <- subset(rd, monkey == 1 & rt > 0.1 & rt < 1.65)   # the paper's own trimming
  rd$upper <- as.integer(rd$correct == 1)
  rd$cohf <- factor(rd$coh, levels = COH)
  cat(nrow(rd), "trials; min RT", min(rd$rt), "s; accuracy", round(mean(rd$upper), 3), "\n")

  r1 <- frm(bf(rt | vint(upper) ~ 0 + coh, bias = 0.5), family = wiener(), data = rd)
  r2 <- frm(bf(rt | vint(upper) ~ 0 + cohf, bias = 0.5), family = wiener(), data = rd)
  r4 <- frm(bf(rt | vint(upper) ~ 0 + cohf, bs ~ 0 + cohf, ndt ~ 0 + cohf, bias = 0.5),
            family = wiener(), data = rd)

  cat("\nM1 on real data, natural scale:\n"); print(round(natural(r1, rd), 4))
  cat("\nM2 drift per coherence:\n")
  print(round(data.frame(coh = COH, drift = unname(fixef(r2)$mu)), 3))

  cat("\nmodel comparison, real monkey 1 data:\n")
  print(data.frame(model = c("M1 (3 par)", "M2 (8 par)", "M4 (18 par)"),
                   npar = c(length(coef(r1)), length(coef(r2)), length(coef(r4))),
                   logLik = round(c(logLik(r1), logLik(r2), logLik(r4)), 2),
                   BIC = round(c(BIC(r1), BIC(r2), BIC(r4)), 2)))

  # The ndt link ceiling is a real constraint, not a formality: PyDDM's tutorial fit puts
  # non-decision time ABOVE the fastest observed RT, which only the lapse mixture makes
  # admissible. Without a contaminant component wiener() cannot go there by construction.
  cat("\nndt ceiling check: min(RT) =", min(rd$rt),
      "; PyDDM tutorial nondectime = 0.211\n")
  cat("  fitted ndt =", round(natural(r1, rd)[["ndt"]], 4),
      "-- capped below min(RT) by the scaled_logit link.\n")
} else {
  cat("\n[real data not found at", csv, "- Part B skipped]\n")
}
