# Probe 5: is the GDDM estimator of probe 4 actually right, and what does it cost?
#
# Probe 4 fit the 6-parameter GDDM by gradient-based ML and recovered mu0, alpha, tau and
# t_nd well but returned leak = 2.16 against a simulating value of 1.2 and B0 = 1.27
# against 1.5. Two candidate explanations, and they call for different verdicts:
#
#   (a) the estimator is biased, e.g. by the grid or by the B-spline shift, or
#   (b) the SIMULATOR is biased. Euler-Maruyama with discrete monitoring misses boundary
#       excursions that occur between steps, so simulated first passages are late. That
#       is a property of the data generator, not of the fit.
#
# The clean discriminator is to simulate from the solver's own density. If recovery is
# then tight, the estimator is sound and probe 4's gap was (b). This also measures how
# much the answer moves with grid resolution, which is the number that decides whether a
# real implementation can afford PyDDM's default grid.

sp <- Sys.getenv("EL_SCRATCH")
.libPaths(c(file.path(sp, "el-lib"), .libPaths()))
suppressPackageStartupMessages(library(RTMB))
source(file.path("dev", "feasibility-elife56938", "gddm-core.R"))

truth <- list(mu0 = 8.0, alpha = 0.8, leak = 1.2, B0 = 1.5, tau = 1.0, t_nd = 0.25)
start <- list(lmu0 = log(6), lalpha = log(1), leak_raw = 0.5, lB0 = log(1.2),
              ltau = log(1.5), qt_nd = qlogis(0.2 / NDT_MAX))

to_nat <- function(p) c(mu0 = exp(p[["lmu0"]]), alpha = exp(p[["lalpha"]]),
                        leak = p[["leak_raw"]], B0 = exp(p[["lB0"]]),
                        tau = exp(p[["ltau"]]), t_nd = NDT_MAX * plogis(p[["qt_nd"]]))

# ---- 1. solver convergence in the grid ---------------------------------------------
# Plain R, no AD, so this is cheap and says how fine a grid the likelihood actually needs.
cat("=== solver convergence: choice probability and mean RT vs grid ===\n")
ref <- NULL
tab <- do.call(rbind, lapply(list(c(50, 51), c(100, 101), c(200, 201), c(200, 361),
                                  c(400, 401), c(800, 801)), function(g) {
  nt <- g[1]; ny <- g[2]; dt <- T_DUR / nt
  s <- gddm_solve(truth$mu0 * (0.128 / CMAX)^truth$alpha, truth$leak, truth$B0,
                  truth$tau, nt, ny, dt)
  pu <- sum(s$up) * dt; pl <- sum(s$lo) * dt
  tg <- seq(0, T_DUR, by = dt)
  data.frame(nt = nt, ny = ny, p_upper = pu, p_lower = pl, mass = pu + pl,
             mean_rt = sum(tg * (s$up + s$lo)) * dt / (pu + pl))
}))
tab$d_p_upper <- c(NA, diff(tab$p_upper))
print(round(tab, 6))
cat("\nThe change from ny=201 to the finest grid is the discretization error you accept\n",
    "by fitting on the coarse grid.\n", sep = "")

# ---- 2. recovery against the solver's own density ------------------------------------
# Draw counts from the model the estimator assumes. Any remaining gap is estimator error.
cat("\n=== recovery from the solver's own density (n = 24000) ===\n")
NT <- 200; NY <- 201; dt <- T_DUR / NT
sim_from_solver <- function(p, n_per, seed = 11) {
  set.seed(seed)
  do.call(rbind, lapply(seq_along(COH), function(j) {
    drift <- if (COH[j] == 0) 0 else p$mu0 * (COH[j] / CMAX)^p$alpha
    s <- gddm_solve(drift, p$leak, p$B0, p$tau, NT, NY, dt)
    sh <- function(v) {
      out <- numeric(length(v))
      for (k in 0:(as.integer(ceiling(NDT_MAX / dt)) + 2L))
        out[(k + 1):length(v)] <- out[(k + 1):length(v)] +
          b3_kernel(p$t_nd / dt - k) * v[1:(length(v) - k)]
      out
    }
    pu <- sh(s$up) * dt; pl <- sh(s$lo) * dt
    m <- sum(pu) + sum(pl)
    pu <- (1 - LAPSE) * pu / m + LAPSE * 0.5 / NT
    pl <- (1 - LAPSE) * pl / m + LAPSE * 0.5 / NT
    prob <- c(pu, pl); prob <- prob / sum(prob)
    cnt <- as.vector(rmultinom(1, n_per, prob))
    nb <- length(pu)
    data.frame(coh = COH[j],
               bin = c(seq_len(nb), seq_len(nb)) - 1L,
               resp = rep(c(1L, 0L), each = nb),
               n = cnt)
  }))
}
cells <- sim_from_solver(truth, 4000)
cells <- cells[cells$n > 0 & cells$bin >= 1L & cells$bin <= NT, ]
cat("nonzero cells:", nrow(cells), " total trials:", sum(cells$n), "\n")

# make_obj expects trials; feed it pre-aggregated cells instead
# renorm divides each condition's defective density by its own total mass. The
# discretized solve loses a little probability (see the convergence table above), and that
# loss DEPENDS ON THE PARAMETERS: configurations that absorb faster lose less. Maximizing
# an un-renormalized sum(n log p) therefore pays a bonus for fast absorption, which biases
# leak up and B0 down. PyDDM exposes the same choice.
make_obj_cells <- function(cells, nt, ny, renorm = TRUE) {
  dt <- T_DUR / nt
  wmax <- as.integer(ceiling(NDT_MAX / dt)) + 2L
  ci <- match(cells$coh, COH)
  f <- function(par) {
    getAll(par)
    mu0 <- exp(lmu0); alpha <- exp(lalpha); B0 <- exp(lB0)
    tau <- exp(ltau); t_nd <- NDT_MAX * plogis(qt_nd); leak <- leak_raw
    nll <- 0
    for (j in seq_along(COH)) {
      drift <- if (COH[j] == 0) 0 * mu0 else mu0 * (COH[j] / CMAX)^alpha
      s <- gddm_solve(drift, leak, B0, tau, nt, ny, dt)
      pu <- shift_density(s$up, t_nd, dt, wmax)
      pl <- shift_density(s$lo, t_nd, dt, wmax)
      if (renorm) {
        m <- (sum(pu) + sum(pl)) * dt
        pu <- pu / m; pl <- pl / m
      }
      pu <- (1 - LAPSE) * pu + LAPSE * 0.5 / T_DUR
      pl <- (1 - LAPSE) * pl + LAPSE * 0.5 / T_DUR
      k <- which(ci == j)
      ku <- k[cells$resp[k] == 1L]; kl <- k[cells$resp[k] == 0L]
      if (length(ku)) nll <- nll - sum(cells$n[ku] * log(pu[cells$bin[ku] + 1L]))
      if (length(kl)) nll <- nll - sum(cells$n[kl] * log(pl[cells$bin[kl] + 1L]))
    }
    nll
  }
  f
}

t0 <- proc.time()[3]
obj <- MakeADFun(make_obj_cells(cells, NT, NY), start, silent = TRUE)
tape_s <- proc.time()[3] - t0
cat("tape build", round(tape_s, 1), "s\n")

# ---- honest cost measurement, using the very functions the optimizer calls -----------
R <- 50L
pv <- obj$par
t0 <- proc.time()[3]; for (i in 1:R) v <- obj$fn(pv); t_f <- (proc.time()[3] - t0) / R
t0 <- proc.time()[3]; for (i in 1:R) g <- obj$gr(pv); t_g <- (proc.time()[3] - t0) / R
cat(sprintf("objective %.1f ms per call; gradient %.1f ms per call (6 conditions, %dx%d grid)\n",
            1000 * t_f, 1000 * t_g, NT, NY))

t0 <- proc.time()[3]
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = 400, eval.max = 500))
opt_s <- proc.time()[3] - t0
cat("optimization", round(opt_s, 1), "s;", opt$iterations, "iterations; code",
    opt$convergence, "-", opt$message, "\n")
cat("max |gradient| at optimum:", format(max(abs(obj$gr(opt$par))), digits = 3), "\n")

est <- to_nat(opt$par)
hs <- try(solve(obj$he(opt$par)), silent = TRUE)
se_link <- if (inherits(hs, "try-error")) rep(NA_real_, 6) else sqrt(diag(hs))
# delta method onto the natural scale
jac <- c(est[["mu0"]], est[["alpha"]], 1, est[["B0"]], est[["tau"]],
         NDT_MAX * dlogis(opt$par[["qt_nd"]]))
se_nat <- se_link * abs(jac)
cat("\nrecovery from the model's own density:\n")
print(data.frame(truth = unlist(truth), estimate = round(est, 4),
                 se = round(se_nat, 4),
                 z = round((est - unlist(truth)) / se_nat, 2)))

# ---- 2b. the same fit WITHOUT renormalization, to isolate the mass-loss bias ---------
cat("
=== same data, same grid, renorm = FALSE ===
")
objn <- MakeADFun(make_obj_cells(cells, NT, NY, renorm = FALSE), start, silent = TRUE)
optn <- nlminb(objn$par, objn$fn, objn$gr, control = list(iter.max = 400, eval.max = 500))
print(data.frame(truth = unlist(truth),
                 renorm_TRUE = round(est, 4),
                 renorm_FALSE = round(to_nat(optn$par), 4)))

# ---- 3. the same fit at a finer grid, to separate grid error from estimator error ----
cat("\n=== refit at nt=200, ny=361 ===\n")
obj2 <- MakeADFun(make_obj_cells(cells, 200, 361, renorm = TRUE), start, silent = TRUE)
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, control = list(iter.max = 400, eval.max = 500))
cat("code", opt2$convergence, "-", opt2$message, "\n")
print(data.frame(truth = unlist(truth),
                 ny201 = round(est, 4), ny361 = round(to_nat(opt2$par), 4)))
