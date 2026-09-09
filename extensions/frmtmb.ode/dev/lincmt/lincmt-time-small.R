# Development-scale timing: the Phase 0 schedule at a smaller subject
# count, both arms in ONE process, interleaved, with a control arm built
# from the same code that must report 1.0.
#
# Usage: Rscript lincmt-time-small.R <n_subject> <rounds>
args <- commandArgs(trailingOnly = TRUE)
NS <- if (length(args) >= 1L) as.integer(args[[1]]) else 10L
ROUNDS <- if (length(args) >= 2L) as.integer(args[[2]]) else 3L

.libPaths(c("C:/Users/adf44/source/r/lincmt-lib",
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(RTMB); library(frmtmb); library(frmtmb.ode)
})
cat("frmtmb.ode from:", dirname(getNamespaceInfo("frmtmb.ode",
                                                 "path")), "\n")

TR <- list(ka = 1.0, ke = 0.15, V = 20, sigma = 0.3, sd_lka = 0.3,
           sd_lke = 0.25, amt = 100, ii = 12)
ode_dyn <- function(t, y, p) {
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}
ode_doses <- data.frame(time = c(0, 12), state = "depot",
                        value = TR$amt, ii = c(12, 12),
                        addl = c(0L, 12L), ss = c(TRUE, FALSE))
ode_times <- 144 + c(0.5, 1, 2, 4, 6, 8, 10, 12)
conc_ss <- function(t, ka, ke, V, amt, ii) {
  u <- t %% ii
  amt * ka / (V * (ka - ke)) *
    (exp(-ke * u) / (1 - exp(-ke * ii)) -
       exp(-ka * u) / (1 - exp(-ka * ii)))
}
mk_data <- function(ns, seed = 20260908L) {
  set.seed(seed)
  lka <- log(TR$ka) + stats::rnorm(ns, 0, TR$sd_lka)
  lke <- log(TR$ke) + stats::rnorm(ns, 0, TR$sd_lke)
  tt <- ode_times
  d <- data.frame(id = factor(rep(seq_len(ns), each = length(tt))),
                  time = rep(tt, times = ns))
  i <- as.integer(d$id)
  d$conc <- conc_ss(d$time, exp(lka[i]), exp(lke[i]), TR$V, TR$amt,
                    TR$ii) + stats::rnorm(nrow(d), 0, TR$sigma)
  d
}

form_ode <- function() {
  doses <- ode_doses
  bf(conc ~ frm_ode(ode_dyn, init = list(0, 0), times = time,
                    parms = list(exp(lka), exp(lke), exp(lV)),
                    group = id, states = c("depot", "central"),
                    output = "central", events = doses),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
# `n_ss` has to be a literal in a nonlinear body: a bare name there is
# a request for a column of the data, so a scalar held in a variable is
# looked for in the frame and not found.
form_lin20 <- function() {
  doses <- ode_doses
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, events = doses, n_ss = 20L),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
form_linInf <- function() {
  doses <- ode_doses
  bf(conc ~ frm_lincmt(parms = list(ka = exp(lka), ke = exp(lke),
                                    V = exp(lV)),
                       times = time, group = id, ncmt = 1,
                       depot = TRUE, events = doses),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE)
}
ST <- list(beta = c(log(0.8), log(0.2), log(15)))

d <- mk_data(NS)
cat("subjects:", NS, " rows:", nrow(d), "\n")

# ---- values agree, at the same parameters -------------------------
o_ode <- frm(form_ode() + gaussian(), data = d, start = ST,
             dry_run = "objective")$obj
o_lin <- frm(form_lin20() + gaussian(), data = d, start = ST,
             dry_run = "objective")$obj
o_inf <- frm(form_linInf() + gaussian(), data = d, start = ST,
             dry_run = "objective")$obj
p0 <- o_ode$par
cat("objective ode", format(o_ode$fn(p0), digits = 14),
    " lincmt(n_ss=20)", format(o_lin$fn(p0), digits = 14),
    " rel", format(abs(o_lin$fn(p0) - o_ode$fn(p0)) /
                     abs(o_ode$fn(p0))), "\n")
cat("objective lincmt(n_ss=Inf)", format(o_inf$fn(p0), digits = 14),
    " rel to ode",
    format(abs(o_inf$fn(p0) - o_ode$fn(p0)) / abs(o_ode$fn(p0))), "\n")
g_ode <- o_ode$gr(p0); g_lin <- o_lin$gr(p0)
cat("gradient max rel diff",
    format(max(abs(g_lin - g_ode)) / max(abs(g_ode))), "\n")

# ---- tape sizes, of the PREDICTOR rather than of the marginal
# objective: the Laplace objective is itself an atomic node, so taping
# it counts nothing about either path
d1 <- d[as.integer(d$id) == 1L, ]
ndf <- function(f) nrow(MakeTape(f, numeric(3))$data.frame())
nd_lin <- ndf(function(p) sum(frm_lincmt(
  parms = list(ka = exp(p[1]), ke = exp(p[2]), V = exp(p[3])),
  times = d1$time, ncmt = 1, depot = TRUE, events = ode_doses,
  n_ss = 20L)))
nd_inf <- ndf(function(p) sum(frm_lincmt(
  parms = list(ka = exp(p[1]), ke = exp(p[2]), V = exp(p[3])),
  times = d1$time, ncmt = 1, depot = TRUE, events = ode_doses)))
nd_ode <- ndf(function(p) sum(frm_ode(
  ode_dyn, init = list(0, 0), times = d1$time,
  parms = list(exp(p[1]), exp(p[2]), exp(p[3])),
  states = c("depot", "central"), output = "central",
  events = ode_doses)))
cat("
nodes for ONE subject's predictor: ode", nd_ode,
    " lincmt(20)", nd_lin, " lincmt(Inf)", nd_inf, "
")
cat("(the ode number is not comparable on its own: RTMBode::ode() is",
    "ONE ADjoint atomic node per solve, whose reverse pass is a whole",
    "numerical solve)
")

# ---- interleaved gradient timing ----------------------------------
# Every call is at a DIFFERENT parameter, and every object is warmed
# before its block size is chosen. Repeating one point measures the
# second and later gradient there, which a Laplace objective serves
# from a cached inner solve; see lincmt-instrument.R for the size of
# that error.
set.seed(1L)
PJIT <- matrix(rnorm(64 * 3, 0, 0.02), 64, byrow = TRUE)
block <- function(obj, p, n) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(n)) obj$gr(p + PJIT[(i - 1L) %% 64L + 1L, ])
  proc.time()[["elapsed"]] - t0
}
grow <- function(obj, p) {
  obj$gr(p)
  n <- 1L
  repeat {
    s <- block(obj, p, n)
    if (s > 1.2 || n >= 4096L) return(list(n = n, s = s))
    n <- n * 2L
  }
}
n_ode <- grow(o_ode, p0)$n
n_lin <- grow(o_lin, p0)$n
n_inf <- grow(o_inf, p0)$n
# the control is a second block of the SAME object, not a duplicate
# objective: see lincmt-time-scale.R for why that matters
cat("\nblock sizes: ode", n_ode, " lincmt(20)", n_lin, " lincmt(Inf)",
    n_inf, "\n")
res <- matrix(NA_real_, ROUNDS, 4,
              dimnames = list(NULL, c("ode", "lin20", "linInf",
                                      "control")))
for (r in seq_len(ROUNDS)) {
  res[r, "ode"] <- block(o_ode, p0, n_ode) / n_ode
  res[r, "lin20"] <- block(o_lin, p0, n_lin) / n_lin
  res[r, "linInf"] <- block(o_inf, p0, n_inf) / n_inf
  res[r, "control"] <- block(o_lin, p0, n_lin) / n_lin
}
best <- apply(res, 2, min)
cat("\nseconds per gradient (minimum over", ROUNDS, "rounds):\n")
print(signif(best, 4))
cat("control / lin20 (must be 1.0):",
    format(best[["control"]] / best[["lin20"]], digits = 4), "\n")
cat("speedup ode / lincmt(20):",
    format(best[["ode"]] / best[["lin20"]], digits = 5), "\n")
cat("speedup ode / lincmt(Inf):",
    format(best[["ode"]] / best[["linInf"]], digits = 5), "\n")
cat("round spread (max/min) per arm:\n")
print(signif(apply(res, 2, max) / apply(res, 2, min), 4))
