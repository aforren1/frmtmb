# Probe 4: can the paper's headline GDDM (Eq. 13) be fit by gradient-based ML on an
# RTMB tape? This is the decisive question, because everything else in the paper either
# reduces to the analytic Wiener family or is a quadrature away from it.
#
# Three obstacles in PyDDM's own solver are value-dependent branches, which a tape cannot
# record. Each needs a replacement, not a translation:
#
#   1. MOVING DOMAIN. The collapsing bound shrinks the spatial grid, and PyDDM sandwiches
#      the bound between two integer grids with weights weight_inner / weight_outer
#      (model.py, solve_numerical). Those weights are continuous but kinked in B0 and tau,
#      and the grid indices themselves are parameter-dependent. Replacement: substitute
#      y = x / B(t). The walls sit at +/-1 for all t, the grid is fixed, and no index
#      depends on a parameter:
#
#        dx = (mu0 (C/Cmax)^alpha - l x) dt + dW,   B(t) = B0 exp(-t/tau)
#        dy = [mu0 (C/Cmax)^alpha / B(t) - l y + y/tau] dt + (1/B(t)) dW
#
#      Every coefficient is now a smooth function of every parameter.
#
#   2. EARLY EXIT. PyDDM breaks out of the time loop when the surviving mass drops below
#      1e-4. Replacement: none needed, just run all steps. It costs time, not correctness.
#
#   3. NON-DECISION TIME LOOKUP. Indexing the density at round((rt - t_nd)/dt) makes an
#      integer index depend on a parameter. Replacement: convolve the density with a
#      kernel of FIXED length spanning the whole admissible t_nd range, with
#      parameter-dependent weights. Indices fixed, weights smooth. See gddm-core.R for
#      why the kernel ended up as a truncated-power cubic B-spline.
#
# PyDDM uses backward Euler whenever bounds are time-varying (Crank-Nicolson is listed in
# the paper's Table 1 as requiring time-independent bounds). The rescaling removes that
# restriction, so Crank-Nicolson is available here and is used below.
#
# READ probe-05 FOR THE HEADLINE NUMBERS, not this file. Two caveats apply here:
#   - the f_ms / grad_ms columns below are unreliable; they report 0 for a tape that
#     plainly takes longer than that. probe-05 times obj$fn and obj$gr, the functions the
#     optimizer actually calls, and gets 0.4 ms and 28 ms.
#   - this fit does NOT renormalize the defective density, and it is fit to
#     Euler-Maruyama data, so its leak and B0 estimates are biased. probe-05 isolates
#     both causes and recovers all six parameters.
# What this file does establish, and it is the point: the GDDM likelihood tapes, its AD
# gradient is correct to 1e-9, and nlminb converges on it.

sp <- Sys.getenv("EL_SCRATCH")
.libPaths(c(file.path(sp, "el-lib"), .libPaths()))
suppressPackageStartupMessages(library(RTMB))

source(file.path("dev", "feasibility-elife56938", "gddm-core.R"))

truth <- list(mu0 = 8.0, alpha = 0.8, leak = 1.2, B0 = 1.5, tau = 1.0, t_nd = 0.25)
start <- list(lmu0 = log(6), lalpha = log(1), leak_raw = 0.5, lB0 = log(1.2),
              ltau = log(1.5), qt_nd = qlogis(0.2 / NDT_MAX))

dat <- sim_gddm(600, truth)
cat("simulated GDDM:", nrow(dat), "trials over", length(COH), "coherences\n")
cat("accuracy by coherence:\n")
print(round(tapply(dat$resp, dat$coh, mean), 3))

# ---- cost at increasing grid resolution --------------------------------------------
cat("\n=== tape build and evaluation cost ===\n")
cat("PyDDM's tutorial GDDM uses dx=0.01, dt=0.01, T_dur=2, i.e. about 360 spatial by\n",
    "200 temporal nodes, solved once per coherence.\n\n", sep = "")
grids <- list(c(nt = 50, ny = 51), c(nt = 100, ny = 101), c(nt = 200, ny = 201),
              c(nt = 200, ny = 361))
bench <- do.call(rbind, lapply(grids, function(g) {
  ob <- make_obj(dat, g[["nt"]], g[["ny"]])
  t0 <- proc.time()[3]
  tp <- MakeTape(ob$f, start)
  t_tape <- proc.time()[3] - t0
  pv <- unlist(start)
  R <- 200L
  t0 <- proc.time()[3]; for (i in 1:R) v <- tp(pv); t_f <- (proc.time()[3] - t0) / R
  tg <- tp$jacfun()
  t0 <- proc.time()[3]; for (i in 1:R) g2 <- tg(pv); t_g <- (proc.time()[3] - t0) / R
  data.frame(nt = g[["nt"]], ny = g[["ny"]],
             tape_s = round(t_tape, 2), f_ms = round(1000 * t_f, 3),
             grad_ms = round(1000 * t_g, 3), nll = round(as.numeric(v), 3))
}))
print(bench)

# ---- gradient correctness ------------------------------------------------------------
cat("\n=== AD gradient vs finite differences (nt=100, ny=101) ===\n")
ob <- make_obj(dat, 100, 101)
tp <- MakeTape(ob$f, start)
pv <- unlist(start)
ad <- as.numeric(tp$jacfun()(pv))
fd <- numeric(length(pv))
for (i in seq_along(pv)) {
  e <- 1e-5; p1 <- pv; p2 <- pv; p1[i] <- p1[i] + e; p2[i] <- p2[i] - e
  fd[i] <- (tp(p1) - tp(p2)) / (2 * e)
}
print(data.frame(par = names(pv), ad = signif(ad, 7), fd = signif(fd, 7),
                 rel_err = signif((ad - fd) / pmax(abs(fd), 1e-8), 3)))

# ---- fit and recovery ----------------------------------------------------------------
cat("\n=== ML fit of the 6-parameter GDDM (nt=200, ny=201) ===\n")
ob <- make_obj(dat, 200, 201)
cat("aggregated to", ob$ncell, "counted cells from", nrow(dat), "trials\n")
t0 <- proc.time()[3]
obj <- MakeADFun(ob$f, start, silent = TRUE)
cat("tape build", round(proc.time()[3] - t0, 1), "s\n")
t0 <- proc.time()[3]
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(iter.max = 300, eval.max = 400, trace = 0))
cat("optimization", round(proc.time()[3] - t0, 1), "s;",
    opt$iterations, "iterations; convergence code", opt$convergence, "\n")
cat("message:", opt$message, "\n")

est <- opt$par
fit <- c(mu0 = exp(est[["lmu0"]]), alpha = exp(est[["lalpha"]]), leak = est[["leak_raw"]],
         B0 = exp(est[["lB0"]]), tau = exp(est[["ltau"]]),
         t_nd = NDT_MAX * plogis(est[["qt_nd"]]))
tru <- unlist(truth)[names(fit)]
hs <- try(solve(obj$he(est)), silent = TRUE)
se <- if (inherits(hs, "try-error")) rep(NA_real_, 6) else sqrt(diag(hs))

cat("\nGDDM parameter recovery:\n")
print(data.frame(truth = round(tru, 4), estimate = round(fit, 4),
                 se_link = round(se, 4)))
