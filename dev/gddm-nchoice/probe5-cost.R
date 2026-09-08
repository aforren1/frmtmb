# Probe 5: what the grid route costs as the number of alternatives grows,
# and what one solve of the SHIPPED 1-D family costs, for scale.
#
# Run from the worktree root:
#   Rscript dev/gddm-nchoice/probe5-cost.R
#
# n alternatives means an (n-1)-dimensional domain. This measures the tape
# build and one gradient for the same Douglas ADI scheme in one, two and
# three dimensions, at matched resolution per axis and matched step count,
# so that the exponent in the dimension is measured rather than assumed.
#
# The three-dimensional case here uses a CUBE, not the regular tetrahedron
# that four alternatives actually needs. That is deliberate: the cost of
# the scheme is set by the node count and the number of directional sweeps,
# both of which the cube reproduces, and building a correct tetrahedral
# domain would add code without changing the number this probe is for. The
# cube is not a fiction either: it is the domain of the n-accumulator
# diffusion with ABSOLUTE thresholds and correlated noise, which is the
# race with its independence assumption removed.

source("dev/gddm-nchoice/common.R")
source("dev/gddm-nchoice/pde2d.R")
suppressPackageStartupMessages(library(frmtmb.eam))

# --- for scale: one solve and one tape of the shipped 1-D family ----------
fam <- gddm()
comp <- fam[["gddm"]]$comp
dpn <- fam[["gddm"]]$dpars
ctl1 <- function(dt, ny) {
  nt <- as.integer(round(2 / dt))
  list(dt = dt, ny = ny, t_max = nt * dt, nt = nt, renormalize = TRUE,
       max_ndt = 0.3, wmax = as.integer(ceiling(0.3 / dt)) + 2L,
       dpars = dpn, tridiagonal = "recorded")
}
cat("-- shipped 1-D gddm, one solve in plain doubles --\n")
for (cfg in list(c(0.01, 201), c(0.005, 201), c(0.01, 401))) {
  ct <- ctl1(cfg[1], as.integer(cfg[2]))
  tt <- gd_time(frmtmb.eam:::gd_solve(
    list(mu = 2, bs = 2, bias = 0.5, ndt = 0), 0, comp, ct), reps = 3L)
  cat(sprintf("dt = %.3f, ny = %3d, nt = %3d, %6d node-steps: %.4f s\n",
              cfg[1], cfg[2], ct$nt, cfg[2] * ct$nt, tt$sec))
}

# --- the same ADI in d dimensions on a cube ------------------------------
#
# One explicit stage carrying every term, then d implicit corrections, each
# a set of tridiagonal systems along one axis. All coefficients constant,
# so the recurrence vectorizes across lines exactly as pde2d.R does.
gd_cube_grid <- function(n, d) {
  np <- n + 2L                       # one wall node each side, per axis
  gr <- do.call(expand.grid, rep(list(seq_len(n)), d))
  flat <- function(m) as.integer(as.matrix(m) %*% np^(seq_len(d) - 1L)) + 1L
  ctr <- flat(gr)
  sh <- function(k, s) { m <- gr; m[[k]] <- m[[k]] + s; flat(m) }
  lines <- lapply(seq_len(d), function(k)
    lapply(seq_len(n), function(q) which(gr[[k]] == q)))
  list(n = n, d = d, np = np, N = nrow(gr), ctr = ctr,
       ip = lapply(seq_len(d), sh, s = 1L),
       im = lapply(seq_len(d), sh, s = -1L), lines = lines)
}

gd_cube_solve <- function(b, h, g, nt, dt, theta = 0.5) {
  "[<-" <- RTMB::ADoverload("[<-")
  d <- g$d; n <- g$n
  P <- numeric(g$np^d) + b[1L] * 0
  p <- numeric(g$N) + b[1L] * 0
  p[seq_len(g$N)] <- 1 / g$N
  one <- rep(1, g$N)
  lo <- lapply(seq_len(d), function(k) one * (-theta * dt *
                 (b[k] / (2 * h) + 1 / (2 * h * h))))
  up <- lapply(seq_len(d), function(k) one * (-theta * dt *
                 (-b[k] / (2 * h) + 1 / (2 * h * h))))
  di <- lapply(seq_len(d), function(k) one * (1 + theta * dt / (h * h)))
  for (k in seq_len(d)) {
    lo[[k]][g$lines[[k]][[1L]]] <- lo[[k]][g$lines[[k]][[1L]]] * 0
    up[[k]][g$lines[[k]][[n]]] <- up[[k]][g$lines[[k]][[n]]] * 0
  }
  out <- numeric(nt + 1L) + b[1L] * 0
  for (st in seq_len(nt)) {
    P[g$ctr] <- p
    a <- lapply(seq_len(d), function(k)
      -b[k] * (P[g$ip[[k]]] - P[g$im[[k]]]) / (2 * h) +
        (P[g$ip[[k]]] - 2 * P[g$ctr] + P[g$im[[k]]]) / (2 * h * h))
    y <- p + dt * Reduce(`+`, a)
    for (k in seq_len(d))
      y <- gd_pde_thomas(lo[[k]], di[[k]], up[[k]], y - theta * dt * a[[k]],
                         g$lines[[k]])
    p <- y
    out[st + 1L] <- sum(p[g$lines[[1L]][[n]]])
  }
  out
}

cat("\n-- the same scheme in d dimensions, cube of n nodes per axis --\n")
cat(sprintf("%3s %5s %5s %9s %14s %12s %12s %14s\n",
            "d", "n", "nt", "nodes", "build s", "fn s", "grad s",
            "us/node-step"))
cfgs <- list(c(1, 51, 100), c(1, 51, 200), c(2, 25, 100), c(2, 37, 100),
             c(2, 37, 200), c(3, 13, 100), c(3, 17, 100), c(3, 21, 100))
for (cf in cfgs) {
  d <- as.integer(cf[1]); n <- as.integer(cf[2]); nt <- as.integer(cf[3])
  g <- gd_cube_grid(n, d)
  mk <- function() {
    f <- function(pp) {
      v <- gd_cube_solve(rep(pp$b, d), 2 / (n + 1), g, nt, 2 / nt)
      -sum(log(frmtmb.eam:::ddm_floor(v[-1L], 1e-300)))
    }
    RTMB::MakeADFun(f, list(b = 0.3), silent = TRUE)
  }
  tb <- gd_time(mk())
  ob <- tb$value
  tf <- gd_time(ob$fn(ob$par), reps = 3L)
  tgr <- gd_time(ob$gr(ob$par), reps = 3L)
  cat(sprintf("%3d %5d %5d %9d %14.2f %12.4f %12.4f %14.2f\n",
              d, n, nt, g$N, tb$sec, tf$sec, tgr$sec,
              1e6 * tb$sec / (g$N * nt)))
}
