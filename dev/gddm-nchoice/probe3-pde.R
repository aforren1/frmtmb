# Probe 3: the two-dimensional Fokker-Planck route, against the closed form
# and against the shipped one-dimensional solver's cost.
#
# Run from the worktree root:
#   Rscript dev/gddm-nchoice/probe3-pde.R
#
# probe0 numbers to compare against: the shipped 1-D gddm at dt = 0.01,
# ny = 201, t_max = 2 builds a tape in 8.56 s per condition, evaluates a
# gradient over 2000 trials in 6.7 ms, and matches the Wiener closed form to
# 0.0018 in the log density past 0.2 s.

source("dev/gddm-nchoice/common.R")
source("dev/gddm-nchoice/images.R")
source("dev/gddm-nchoice/pde2d.R")
suppressPackageStartupMessages({library(frmtmb.eam); library(numDeriv)})

tl <- gd_tiles(250L)
U <- gd_basis(3L)
mu <- c(0.9, 0.1, -1.0)
cc <- 0.7071
b <- mu - mean(mu)                    # drift in s coordinates
b <- b[1:2]
a <- as.numeric(t(U) %*% mu)          # drift in z coordinates
r <- sqrt(1.5) * cc

cat(sprintf("mu = (%s), c = %.4f\n", paste(mu, collapse = ", "), cc))

# --- convergence of the PDE to the closed form ----------------------------
t_max <- 2
cat("\n-- PDE against the image series, max |log f_k| error over t in [0.15, 2] --\n")
cat(sprintf("%5s %8s %8s %10s %10s %10s %10s %10s\n",
            "N", "nt", "nodes", "solve s", "mass", "err k=1", "err k=2", "err k=3"))
res <- list()
for (gh in c(FALSE, TRUE)) {
cat(sprintf("ghost correction: %s
", gh))
for (N in c(24L, 36L, 54L, 81L)) {
  g <- gd_pde_grid(N)
  for (nt in c(100L, 400L)) {
    dt <- t_max / nt
    t0 <- proc.time()[["elapsed"]]
    fl <- gd_pde_solve(b, cc, c(0, 0), g, nt, dt, ghost = gh)
    el <- proc.time()[["elapsed"]] - t0
    tg <- seq(0, t_max, by = dt)
    sel <- tg >= 0.15
    ex <- lapply(1:3, function(k) gd_img_density(tg[sel], a, c(0, 0), r, tl, k))
    er <- vapply(1:3, function(k)
      max(abs(log(pmax(fl[[k]][sel], 1e-300)) - log(ex[[k]]))), numeric(1))
    m <- sum(vapply(fl, sum, numeric(1))) * dt
    cat(sprintf("%5d %8d %8d %10.3f %10.6f %10.5f %10.5f %10.5f\n",
                N, nt, g$nin, el, m, er[1], er[2], er[3]))
    res[[paste(gh, N, nt)]] <- list(N = N, nt = nt, fl = fl, err = er)
  }
}
}

# --- is the three-fold symmetry preserved? --------------------------------
cat("\n-- symmetry: zero drift, centred start, the three fluxes should agree --\n")
for (N in c(36L, 54L, 81L)) {
  g <- gd_pde_grid(N)
  fl <- gd_pde_solve(c(0, 0), cc, c(0, 0), g, 200L, t_max / 200)
  s <- vapply(fl, sum, numeric(1)) * (t_max / 200)
  cat(sprintf("N = %3d   P = %.6f %.6f %.6f   spread %.2e   total %.6f\n",
              N, s[1], s[2], s[3], max(s) - min(s), sum(s)))
}

# --- the tape ------------------------------------------------------------
mk_obj <- function(ntrial, N, nt, seed = 31L) {
  set.seed(seed)
  y <- runif(ntrial, 0.35, 1.6)
  ch <- sample.int(3L, ntrial, TRUE)
  dt <- t_max / nt
  g <- gd_pde_grid(N)
  nb <- nt + 1L
  s <- (y - 0.25) / dt
  k0 <- as.integer(floor(s)); wq <- s - k0
  base <- (ch - 1L) * nb
  i1 <- base + k0 + 1L; i2 <- base + k0 + 2L
  f <- function(p) {
    fl <- gd_pde_solve(c(p$b1, p$b2), exp(p$logc), c(p$s01, p$s02), g, nt, dt)
    v <- do.call(c, fl)
    -sum(log(frmtmb.eam:::ddm_floor((1 - wq) * v[i1] + wq * v[i2], 1e-300)))
  }
  RTMB::MakeADFun(f, list(b1 = b[1], b2 = b[2], logc = log(cc),
                          s01 = 0, s02 = 0), silent = TRUE)
}

cat("\n-- tape build and one gradient, 2000 trials, one condition --\n")
cat(sprintf("%5s %6s %8s %12s %10s %10s\n",
            "N", "nt", "nodes", "build s", "fn s", "grad s"))
for (cfg in list(c(24, 100), c(36, 100), c(36, 200), c(54, 200))) {
  N <- as.integer(cfg[1]); nt <- as.integer(cfg[2])
  tb <- gd_time(mk_obj(2000L, N, nt))
  ob <- tb$value
  tf <- gd_time(ob$fn(ob$par), reps = 3L)
  tg2 <- gd_time(ob$gr(ob$par), reps = 3L)
  cat(sprintf("%5d %6d %8d %12.2f %10.4f %10.4f\n",
              N, nt, gd_pde_grid(N)$nin, tb$sec, tf$sec, tg2$sec))
}

cat("\n-- PDE gradient against numDeriv, N = 24, nt = 100, 200 trials --\n")
ob <- mk_obj(200L, 24L, 100L)
x <- ob$par + c(0.05, -0.03, 0.02, 0.01, -0.01)
g1 <- ob$gr(x)
g2 <- numDeriv::grad(function(v) ob$fn(v), x)
cat(sprintf("%8s %16s %16s %12s\n", "par", "RTMB", "numDeriv", "rel"))
for (i in seq_along(x))
  cat(sprintf("%8s %16.6f %16.6f %12.2e\n", names(ob$par)[i], g1[i], g2[i],
              abs(g1[i] - g2[i]) / max(abs(g2[i]), 1e-8)))
