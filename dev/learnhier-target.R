# Lane `learnhier`: is the TARGET the biased quantity?
#
#   Rscript dev/learnhier-target.R
#
# NO FITS for the decomposition. Everything comes from each learner's
# drawn `ndt` and observed floor, which `dev/learnhier-run2.R` saves.
#
# THE QUESTION. Every other component of the block has a population
# truth on the scale the model fits it on. `ndt` under `ndt_group(id)`
# does not: the scale tier's own record says `ndt=no link truth`,
# because the link's ceiling is each learner's own fastest response. Its
# only target is the realized `sd(qlogis(ndt_i / floor_i))`, and
# `floor_i` is an OBSERVED MINIMUM. So `ndt` is the one component whose
# target is estimated rather than known, and it is the one component
# whose bias is six times the others. That coincidence is worth a test.
#
# THE DECOMPOSITION IS AN EXACT IDENTITY, not an approximation, and it
# says what the fitted deviation actually is. A learner's floor is its
# own non-decision time plus the fastest DECISION it happened to make,
#
#   floor_i = ndt_i + m_i,
#
# where `m_i` is that learner's smallest first-passage time over its
# trials. Then
#
#   qlogis(ndt_i / floor_i) = log(ndt_i / (floor_i - ndt_i))
#                           = log(ndt_i) - log(m_i).
#
# So what `ndt ~ (1 | p | id)` estimates under `ndt_group(id)` is NOT
# the learner's non-decision time on some link. It is the log
# non-decision time MINUS the log of an observed minimum. The first term
# is the random effect the design drew. The second is a nuisance: it
# depends on how fast that learner happened to decide, which is a
# function of their drift and their boundary and of 200 draws of luck,
# and has nothing to do with their non-decision time.
#
#   Var(eta) = Var(log ndt) + Var(log m) - 2 Cov(log ndt, log m)
#
# A first attempt held one factor at its mean instead. That is wrong:
# `ndt_i / mean(floor)` leaves the unit interval, `qlogis` saturates,
# and the pieces came back ten times the whole, 4.44 and 5.63 against
# 0.43. The identity above needs no such surgery, and the run below
# checks it numerically rather than asserting it.
source("dev/learnhier-env.R")

fs <- list.files("dev/learnhier-rec", pattern = "^rlddm-[0-9]+[.]rds$",
                 full.names = TRUE)
R <- Filter(function(r) !is.null(r$per_learner), lapply(fs, readRDS))
cat("replicates carrying per-learner vectors: ", length(R), "\n",
    sep = "")
if (!length(R)) quit(save = "no")

row <- function(r) {
  p <- r$per_learner
  nd <- p$ndt_true
  m <- p$margin
  eta <- log(nd) - log(m)
  data.frame(seed = r$seed,
             full = stats::sd(eta),
             from_ndt = stats::sd(log(nd)),
             from_min = stats::sd(log(m)),
             cor_nm = stats::cor(log(nd), log(m)),
             est = r$vc[["sd_ndt.est"]],
             est_lwr = r$vc[["sd_ndt.lwr"]],
             est_upr = r$vc[["sd_ndt.upr"]],
             stringsAsFactors = FALSE)
}
d <- do.call(rbind, lapply(R, row))

cat("\n== the realized target, decomposed exactly ==\n")
q <- function(x) c(mean = mean(x), sd = stats::sd(x))
print(round(t(vapply(d[c("full", "from_ndt", "from_min", "cor_nm",
                         "est")], q, numeric(2))), 4))

err <- vapply(R, function(r) {
  p <- r$per_learner
  a <- stats::var(log(p$ndt_true) - log(p$margin))
  b <- stats::var(log(p$ndt_true)) + stats::var(log(p$margin)) -
    2 * stats::cov(log(p$ndt_true), log(p$margin))
  (a - b) / a
}, numeric(1))
cat("\nidentity check over ", nrow(d),
    " replicates, worst relative error ",
    format(max(abs(err)), digits = 3), "\n", sep = "")

cat("\n== how much of the target is the NUISANCE term? ==\n")
cat("sd(log ndt), the random effect the design drew:          ",
    round(mean(d$from_ndt), 4), "\n", sep = "")
cat("sd(log m), the observed-minimum term it did not:         ",
    round(mean(d$from_min), 4), "\n", sep = "")
cat("their correlation across learners:                       ",
    round(mean(d$cor_nm), 4), "\n", sep = "")
cat("sd(eta), the realized target:                            ",
    round(mean(d$full), 4), "\n", sep = "")
cat("variance share of the target from the observed minimum:  ",
    round(100 * mean(d$from_min^2) / mean(d$full^2), 1),
    " percent\n", sep = "")

cat("\n== where does the estimate sit? ==\n")
cat("estimate ", round(mean(d$est), 4), " against the target ",
    round(mean(d$full), 4), ": ",
    round(100 * mean(d$est - d$full) / mean(d$full), 1),
    " percent, below on ", sum(d$est < d$full), " of ", nrow(d), "\n",
    sep = "")
cat("estimate against sd(log ndt) alone, ", round(mean(d$from_ndt), 4),
    ": ", round(100 * mean(d$est - d$from_ndt) / mean(d$from_ndt), 1),
    " percent\n", sep = "")
cat("the fit's interval covers the target on ",
    sum(d$est_lwr <= d$full & d$est_upr >= d$full), " of ", nrow(d),
    ", and sd(log ndt) on ",
    sum(d$est_lwr <= d$from_ndt & d$est_upr >= d$from_ndt), " of ",
    nrow(d), "\n", sep = "")

cat("\n== how noisy is the target itself? ==\n")
cat("spread across replicates, target ", round(stats::sd(d$full), 4),
    ", estimate ", round(stats::sd(d$est), 4),
    "; the sampling spread of a standard deviation of 100 draws would",
    " be about ", round(mean(d$full) / sqrt(200), 4), "\n", sep = "")
saveRDS(d, "dev/learnhier-target.rds")
cat("\nwrote dev/learnhier-target.rds\n")
