# ndt lane, measurement 1: what a per-subject floor is, on the plan's
# own eam design, and how far it sits from the truth it bounds.
#
# Run: Rscript --vanilla dev/ndt-scripts/ndt-floors.R
# Seeds: 20260908 (the tier's own), then 20260908 + 1:19 for the bias
# study.

.libPaths(c("C:/Users/adf44/source/r/ndt-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

cat("RWiener available:", requireNamespace("RWiener", quietly = TRUE),
    "\n")

tr <- list(mu0 = 0.4, mu_cond = 0.9, bs = 1.4, ndt = 0.25,
           sd_mu = 0.35, sd_lbs = 0.20, sd_lndt = 0.12, sv = 0.4)

make <- function(seed, sv = 0, ns = 30L, nt = 400L) {
  set.seed(seed)
  u_mu <- stats::rnorm(ns, 0, tr$sd_mu)
  u_bs <- stats::rnorm(ns, 0, tr$sd_lbs)
  u_nd <- stats::rnorm(ns, 0, tr$sd_lndt)
  s <- rep(seq_len(ns), each = nt)
  cond <- rep(rep(0:1, each = nt / 2L), times = ns)
  d <- ddm_simulate(ns * nt,
                    mu = tr$mu0 + tr$mu_cond * cond + u_mu[s],
                    bs = tr$bs * exp(u_bs[s]),
                    ndt = tr$ndt * exp(u_nd[s]),
                    bias = 0.5, sv = sv)
  d$s <- factor(s)
  attr(d, "ndt_subject") <- tr$ndt * exp(u_nd)
  d
}

# ------------------------------------------------- the tier's own draw
d <- make(20260908L)
truth <- attr(d, "ndt_subject")
floors <- tapply(d$rt, d$s, min)
gl <- min(d$rt)

cat("\n== the tier's design, seed 20260908, 30 x 400 ==\n")
cat("global floor min(rt)     :", format(gl, digits = 9), "\n")
cat("true ndt range           :", format(range(truth), digits = 4), "\n")
cat("subjects with true ndt > global floor:",
    sum(truth > gl), "of", length(truth), "\n")
cat("subjects with true ndt > own floor   :",
    sum(truth > floors), "of", length(truth), "\n")
cat("per-subject floor range  :", format(range(floors), digits = 4),
    "\n")
cat("floor - true ndt (ms)    : min", round(1000 * min(floors - truth), 2),
    " median", round(1000 * median(floors - truth), 2),
    " max", round(1000 * max(floors - truth), 2), "\n")
cat("truth / own floor        : min", round(min(truth / floors), 4),
    " median", round(median(truth / floors), 4),
    " max", round(max(truth / floors), 4), "\n")
cat("sd(true ndt)             :", format(sd(truth), digits = 4), "\n")
cat("sd(own floor)            :", format(sd(floors), digits = 4), "\n")
cat("sd(truth / own floor)    :", format(sd(truth / floors), digits = 4),
    "\n")

# The plan's claim: at the bound-lifted optimum every subject's ndt sits
# below its own fastest response with 22 ms to spare. The quantity that
# claim is about is floor - ndt, so the smallest of these over subjects
# is what "22 ms of spare" has to be measured against.
cat("\nsmallest margin floor - true ndt, ms:",
    round(1000 * min(floors - truth), 2), "\n")

# ------------------------------------------- is the floor a biased ceiling
#
# floor_i - ndt_i is the fastest DECISION time of that subject's 400
# trials, so it is positive by construction and shrinks with the trial
# count. If the fit is asked for ndt as a fraction of the floor, the
# quantity it can never exceed is floor_i, and the question is how far
# above ndt_i that is, and how much that overshoot VARIES between
# subjects: a constant overshoot is absorbed by the intercept, a
# varying one is what the random effect has to carry.
cat("\n== the floor as a ceiling on ndt, 20 seeds x 30 subjects ==\n")
rows <- list()
for (k in 0:19) {
  dk <- make(20260908L + k)
  tk <- attr(dk, "ndt_subject")
  fk <- as.numeric(tapply(dk$rt, dk$s, min))
  rows[[length(rows) + 1L]] <- data.frame(
    seed = 20260908L + k, over_ms = 1000 * (fk - tk),
    ratio = tk / fk, above_global = tk > min(dk$rt))
}
a <- do.call(rbind, rows)
cat("replicates                 :", nrow(a), "subject-draws over 20 seeds\n")
cat("floor - ndt, ms            : min", round(min(a$over_ms), 2),
    " q25", round(quantile(a$over_ms, 0.25), 2),
    " median", round(median(a$over_ms), 2),
    " q75", round(quantile(a$over_ms, 0.75), 2),
    " max", round(max(a$over_ms), 2), "\n")
cat("mean overshoot, ms         :", round(mean(a$over_ms), 2),
    " sd", round(sd(a$over_ms), 2), "\n")
cat("ndt / own floor            : mean", round(mean(a$ratio), 4),
    " sd", round(sd(a$ratio), 4), "\n")
cat("subjects above GLOBAL floor: ", sum(a$above_global), "of", nrow(a),
    "\n")

# How much of the per-subject spread in the ratio is the overshoot and
# how much is the truth? If ndt were constant across subjects the ratio
# would still vary, purely from the floor. Decompose on one seed set.
cat("\nsd of log(ndt) over subjects, truth   :",
    round(sd(log(tr$ndt * exp(rnorm(0)))), 4), "(0 by construction)\n")
b <- do.call(rbind, lapply(rows, function(r) {
  data.frame(sd_ratio = sd(r$ratio), sd_over = sd(r$over_ms))
}))
cat("per-seed sd(ndt/floor)     : mean", round(mean(b$sd_ratio), 4),
    " range", round(min(b$sd_ratio), 4), round(max(b$sd_ratio), 4), "\n")
cat("per-seed sd(floor - ndt) ms: mean", round(mean(b$sd_over), 2),
    " range", round(min(b$sd_over), 2), round(max(b$sd_over), 2), "\n")

# 400 trials against 100: how fast does the overshoot fall?
cat("\n== overshoot against trials per subject, seed 20260908 ==\n")
for (nt in c(50L, 100L, 200L, 400L, 800L)) {
  dn <- make(20260908L, ns = 30L, nt = nt)
  tn <- attr(dn, "ndt_subject")
  fn <- as.numeric(tapply(dn$rt, dn$s, min))
  cat(sprintf("nt=%4d  mean overshoot %6.2f ms  sd %5.2f  max %6.2f\n",
              nt, mean(1000 * (fn - tn)), sd(1000 * (fn - tn)),
              max(1000 * (fn - tn))))
}
