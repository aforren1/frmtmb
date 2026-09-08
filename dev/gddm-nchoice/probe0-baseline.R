# Probe 0: what the shipped two-choice gddm() costs, so every n-choice
# number later has something to be a multiple of.
#
# Run from the worktree root: Rscript dev/gddm-nchoice/probe0-baseline.R
#
# Tape build and one gradient for the 1-D Fokker-Planck likelihood at the
# shipped grid, swept over the number of conditions, at a realistic trial
# count; plus the 1-D solver's accuracy against the closed-form Wiener
# density, which is the bar any n-choice route has to match.

source("dev/gddm-nchoice/common.R")
suppressPackageStartupMessages(library(frmtmb.eam))
gdens <- frmtmb.eam:::gd_densities
gsolve <- frmtmb.eam:::gd_solve

fam <- gddm()
comp <- fam[["gddm"]]$comp
dpn <- fam[["gddm"]]$dpars
cat("gddm dpars:", paste(dpn, collapse = ", "), "\n")

mk_ctl <- function(dt = 0.01, ny = 201L, t_max = 2, ndt_ub = 0.3) {
  nt <- as.integer(round(t_max / dt))
  list(dt = dt, ny = ny, t_max = nt * dt, nt = nt, renormalize = TRUE,
       max_ndt = ndt_ub, wmax = as.integer(ceiling(ndt_ub / dt)) + 2L,
       dpars = dpn, tridiagonal = "recorded")
}

# ---- accuracy of the shipped 1-D solve against the Wiener closed form ----
ctl <- mk_ctl()
tg <- seq(0, ctl$t_max, by = ctl$dt)
pj <- list(mu = 2, bs = 2, bias = 0.5, ndt = 0)
s <- gsolve(pj, 0, comp, ctl)
w_up <- frmtmb.eam:::ddm_lpdf_both(tg[-1L], 2, 2, 0.5, rep(1, ctl$nt))
w_lo <- frmtmb.eam:::ddm_lpdf_both(tg[-1L], 2, 2, 0.5, rep(0, ctl$nt))
sel <- tg[-1L] >= 0.2
err <- max(abs(log(pmax(s$up[-1L], 1e-300))[sel] - w_up[sel]),
           abs(log(pmax(s$lo[-1L], 1e-300))[sel] - w_lo[sel]))
cat(sprintf("1-D solver vs Wiener, worst |log| error past 0.2 s: %.5f\n", err))
mass <- sum(s$up + s$lo) * ctl$dt
cat(sprintf("1-D unrenormalized mass over [0, %.1f]: %.6f\n", ctl$t_max, mass))

# ---- tape build and gradient, swept over conditions ----------------------
bench_1d <- function(ncond, ntrial = 2000L, ctl = mk_ctl()) {
  set.seed(1)
  gindex <- sort(rep_len(seq_len(ncond), ntrial))
  first <- match(seq_len(ncond), gindex)
  y <- runif(ntrial, 0.35, 1.6)
  up <- rbinom(ntrial, 1L, 0.7)
  nb <- ctl$nt + 1L
  ss <- y / ctl$dt; k0 <- as.integer(floor(ss)); w <- ss - k0
  base <- (gindex - 1L) * 2L * nb + (1L - up) * nb
  d <- list(ncond = ncond, first = first,
            cov = matrix(0, ncond, 1L), gindex = gindex, w = w,
            i1 = base + k0 + 1L, i2 = base + k0 + 2L)
  pars <- list(mu = rep(1.5, ncond), bs = rep(2, ncond),
               bias = 0.5, ndt = 0.25)
  f <- function(p) {
    dp <- list(mu = p$mu[d$gindex], bs = p$bs[d$gindex],
               bias = p$bias, ndt = p$ndt)
    pv <- gdens(dp, comp, ctl, d)
    -sum(log(frmtmb.eam:::ddm_floor((1 - d$w) * pv[d$i1] + d$w * pv[d$i2], 1e-300)))
  }
  tb <- gd_time(RTMB::MakeADFun(f, pars, silent = TRUE))
  obj <- tb$value
  x <- obj$par
  tf <- gd_time(obj$fn(x), reps = 3L)
  tgr <- gd_time(obj$gr(x), reps = 3L)
  list(ncond = ncond, build = tb$sec, fn = tf$sec, gr = tgr$sec,
       nll = tf$value, ngrad = length(x))
}

cat("\n-- shipped 1-D gddm, dt = 0.01, ny = 201, t_max = 2, 2000 trials --\n")
cat(sprintf("%6s %10s %10s %10s %14s\n",
            "ncond", "build s", "fn s", "grad s", "nll"))
for (k in c(1L, 2L, 4L, 8L)) {
  r <- bench_1d(k)
  cat(sprintf("%6d %10.3f %10.4f %10.4f %14.4f\n",
              r$ncond, r$build, r$fn, r$gr, r$nll))
}
