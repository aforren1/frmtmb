# Where the Phase 0 design's time goes, at the Phase 0 design's size.
#
# The whole-fit numbers are in lincmt-phase0.R. This one decomposes
# them: tape build, and one gradient at the starting values, for both
# paths in ONE process, interleaved, with a control arm built from the
# same code that must report 1.0. Tape build is
# `frm(dry_run = "objective")` less `frm(dry_run = "frame")`, which is
# `dev/scale-findings.md`'s own construction.
.libPaths(c("C:/Users/adf44/source/r/lincmt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(RTMB); library(frmtmb); library(frmtmb.ode)
})
args <- commandArgs(trailingOnly = TRUE)
NS <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
ROUNDS <- if (length(args) >= 2L) as.integer(args[[2]]) else 3L

TR <- list(ka = 1.0, ke = 0.15, V = 20, sigma = 0.3, sd_lka = 0.3,
           sd_lke = 0.25, amt = 100, ii = 12)
ode_dyn <- function(t, y, p) {
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}
ode_doses <- data.frame(time = c(0, 12), state = "depot",
                        value = TR$amt, ii = c(12, 12),
                        addl = c(0L, 12L), ss = c(TRUE, FALSE))
conc_ss <- function(t, ka, ke, V, amt, ii) {
  u <- t %% ii
  amt * ka / (V * (ka - ke)) *
    (exp(-ke * u) / (1 - exp(-ke * ii)) -
       exp(-ka * u) / (1 - exp(-ka * ii)))
}
set.seed(20260908L)
lka <- log(TR$ka) + stats::rnorm(NS, 0, TR$sd_lka)
lke <- log(TR$ke) + stats::rnorm(NS, 0, TR$sd_lke)
tt <- 144 + c(0.5, 1, 2, 4, 6, 8, 10, 12)
d <- data.frame(id = factor(rep(seq_len(NS), each = length(tt))),
                time = rep(tt, times = NS))
i <- as.integer(d$id)
d$conc <- conc_ss(d$time, exp(lka[i]), exp(lke[i]), TR$V, TR$amt,
                  TR$ii) + stats::rnorm(nrow(d), 0, TR$sigma)
ST <- list(beta = c(log(0.8), log(0.2), log(15)))

f_ode <- function() {
  doses <- ode_doses
  bf(conc ~ frm_ode(ode_dyn, init = list(0, 0), times = time,
                    parms = list(exp(lka), exp(lke), exp(lV)),
                    group = id, states = c("depot", "central"),
                    output = "central", events = doses),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
f_l20 <- function() {
  doses <- ode_doses
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, events = doses, n_ss = 20L),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
f_inf <- function() {
  doses <- ode_doses
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, events = doses),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
el <- function(e) {
  t0 <- proc.time()[["elapsed"]]
  force(e)
  proc.time()[["elapsed"]] - t0
}
arms <- list(ode = f_ode(), lin20 = f_l20(), linInf = f_inf(),
             control = f_l20())

cat("subjects:", NS, " rows:", nrow(d), " rounds:", ROUNDS, "\n\n")

# one call through frm() first, so that nobody pays a process's
# one-time costs inside a timed block
invisible(frm(f_inf() + gaussian(), data = d, start = ST,
              dry_run = "frame"))

build <- matrix(NA_real_, ROUNDS, length(arms),
                dimnames = list(NULL, names(arms)))
for (r in seq_len(ROUNDS)) {
  for (a in names(arms)) {
    tf <- el(frm(arms[[a]] + gaussian(), data = d, start = ST,
                 dry_run = "frame"))
    to <- el(frm(arms[[a]] + gaussian(), data = d, start = ST,
                 dry_run = "objective"))
    build[r, a] <- to - tf
  }
}
cat("tape build (s), minimum over", ROUNDS, "rounds:\n")
print(signif(apply(build, 2, min), 4))
cat("control / lin20:",
    format(min(build[, "control"]) / min(build[, "lin20"]),
           digits = 4), "\n")

# The control is a SECOND BLOCK OF THE SAME OBJECT, not a fourth
# objective built from the same formula. The first version of this
# script did the latter and the control reported 4.409: a process
# holding four Laplace objectives does not time the fourth the way it
# times the second, so a duplicate object measures the process and not
# the code. That is the instrument failing its own check, which is
# what the check is for.
objs <- lapply(arms[c("ode", "lin20", "linInf")], function(f)
  frm(f + gaussian(), data = d, start = ST, dry_run = "objective")$obj)
objs$control <- objs$lin20
p0 <- objs$ode$par
cat("\nobjective at the starting values, ode",
    format(objs$ode$fn(p0), digits = 14), " lincmt(20)",
    format(objs$lin20$fn(p0), digits = 14), " rel",
    format(abs(objs$lin20$fn(p0) - objs$ode$fn(p0)) /
             abs(objs$ode$fn(p0))), "\n")

# Every call is at a DIFFERENT parameter. Repeating `gr()` at one point
# measures the second and later gradient there, and TMB's Laplace
# objective caches the inner Newton solution on the parameter, so that
# is not the operation an optimizer performs. This script's third
# version repeated one point and reported the arm with FEWER tape nodes
# as twice as slow, which is the reverse of what its whole-fit time
# says.
set.seed(1L)
PJIT <- matrix(rnorm(64 * length(p0), 0, 0.02), 64,
               byrow = TRUE) + rep(p0, each = 64)
blk <- function(o, n) {
  t0 <- proc.time()[["elapsed"]]
  for (k in seq_len(n)) o$gr(PJIT[(k - 1L) %% 64L + 1L, ])
  proc.time()[["elapsed"]] - t0
}
grow <- function(o) {
  # WARM IT FIRST. The first gradient after a tape is built pays a
  # cold inner Newton solve, and a block-growing loop that starts at
  # n = 1 mistakes that for the per-call cost and stops at n = 1. That
  # is how an arm with FEWER tape nodes came out slower than one with
  # more, in this script's own second version.
  o$gr(p0)
  n <- 1L
  repeat {
    if (blk(o, n) > 1.2 || n >= 2048L) return(n)
    n <- n * 2L
  }
}
nb <- vapply(objs, grow, 1L)
cat("\nblock sizes:", paste(names(nb), nb, collapse = "  "), "\n")
gr <- matrix(NA_real_, ROUNDS, length(objs),
             dimnames = list(NULL, names(objs)))
for (r in seq_len(ROUNDS)) {
  for (a in names(objs)) gr[r, a] <- blk(objs[[a]], nb[[a]]) / nb[[a]]
}
best <- apply(gr, 2, min)
cat("\nseconds per gradient at the starting values, minimum over",
    ROUNDS, "rounds:\n")
print(signif(best, 4))
cat("control / lin20 (must be 1.0):",
    format(best[["control"]] / best[["lin20"]], digits = 4), "\n")
cat("speedup ode / lin20 :", format(best[["ode"]] / best[["lin20"]],
                                    digits = 5), "\n")
cat("speedup ode / linInf:", format(best[["ode"]] / best[["linInf"]],
                                    digits = 5), "\n")
cat("round spread (max/min):\n")
print(signif(apply(gr, 2, max) / apply(gr, 2, min), 4))
