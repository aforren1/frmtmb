# Probe A1: hand-rolled RTMB population PK with per-subject RTMBode::ode()
# solves, WITHOUT frmtmb. Establishes that the per-group solve pattern
# tapes, differentiates and optimizes, so that any later frmtmb failure is
# attributable to the nl machinery rather than to the pattern itself.
library(RTMB)
library(RTMBode)
source("dev/ode/pk-common.R")

d <- sim_pk()
cat("rows:", nrow(d), " subjects:", nlevels(d$id), "\n")

## --- 1. numeric-mode check: ODE solution vs closed form ---------------
ka <- exp(PK_TRUTH$lka + d$b_ka); ke <- exp(PK_TRUTH$lke + d$b_ke)
num <- pk_ode(ka, ke, rep(exp(PK_TRUTH$lV), nrow(d)), d$time, d$id, d$dose)
cat("max |ode - analytic| (numeric mode):",
    format(max(abs(num - d$mu_true)), digits = 3), "\n")

## --- 2. taped: MakeADFun over the same helper -------------------------
gi <- as.integer(d$id); ng <- nlevels(d$id)
nll <- function(p) {
  getAll(p)
  n <- length(conc)
  ka <- exp(lka + u_ka[gi])
  ke <- exp(lke + u_ke[gi])
  V  <- exp(rep(lV, n))
  mu <- pk_ode(ka, ke, V, time, id, dose)
  nll <- -sum(dnorm(u_ka, 0, exp(lsd_ka), log = TRUE)) -
    sum(dnorm(u_ke, 0, exp(lsd_ke), log = TRUE))
  nll - sum(dnorm(conc, mu, exp(lsigma), log = TRUE))
}
conc <- d$conc; time <- d$time; id <- d$id; dose <- d$dose

p0 <- list(lka = 0, lke = log(0.25), lV = log(8), lsigma = 0,
           lsd_ka = log(0.3), lsd_ke = log(0.3),
           u_ka = numeric(ng), u_ke = numeric(ng))

t_tape <- system.time(
  obj <- MakeADFun(nll, p0, random = c("u_ka", "u_ke"), silent = TRUE))
cat("tape build:", t_tape[["elapsed"]], "s\n")
cat("fn at start:", obj$fn(obj$par), "\n")
cat("gr at start:", format(obj$gr(obj$par), digits = 4), "\n")

t_opt <- system.time(
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 500, eval.max = 500)))
cat("optimize:", t_opt[["elapsed"]], "s  conv:", opt$convergence,
    " nll:", opt$objective, "\n")
sdr <- sdreport(obj)
est <- as.list(sdr, "Est"); sd_ <- as.list(sdr, "Std")
res <- data.frame(
  par = c("lka", "lke", "lV", "log sigma", "log sd_lka", "log sd_lke"),
  truth = c(PK_TRUTH$lka, PK_TRUTH$lke, PK_TRUTH$lV, log(PK_TRUTH$sigma),
            log(PK_TRUTH$sd_lka), log(PK_TRUTH$sd_lke)),
  est = c(est$lka, est$lke, est$lV, est$lsigma, est$lsd_ka, est$lsd_ke),
  se = c(sd_$lka, sd_$lke, sd_$lV, sd_$lsigma, sd_$lsd_ka, sd_$lsd_ke))
print(res, digits = 4)
cat("PROBEA1 OK\n")
