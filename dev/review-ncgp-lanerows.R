# Item 1: the lane's three short rows and noise-free row, read by my closed form at 256 bits, by
# the local polish, and with the marginal warm start. Seeds as in dev/ncgp-truth.R.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
source(file.path(REV_WT, "dev/review-ncgp-margin.R"))
rows <- list()
for (seed in c(5701, 5702)) {
  set.seed(seed); n <- 100
  x <- sort(runif(n, 0, 10)); e0 <- rnorm(n)
  rows[[length(rows) + 1L]] <- list(lab = sprintf("exp(0.8x) seed %d", seed), x = x, y = exp(0.8 * x) + 0.3 * e0)
}
set.seed(131); u <- sort(runif(20, 0, 6)); x <- rep(u, each = 30); e0 <- rnorm(length(x))
rows[[length(rows) + 1L]] <- list(lab = "rep 20x30 s=1e-2 seed 131", x = x, y = sin(2 * x) + 1e-2 * e0)
set.seed(20261102); x <- round(runif(60, 0, 6), 3)
rows[[length(rows) + 1L]] <- list(lab = "noise-free sin(2x)", x = x, y = sin(2 * x))
set.seed(1080); x <- sort(round(runif(80, 0, 10), 4)); e0 <- rnorm(80)
rows[[length(rows) + 1L]] <- list(lab = "logx s=1e-3 n80", x = x, y = log1p(x) + 1e-3 * e0)
for (r0 in rows) {
  y <- r0$y; x <- r0$x
  t0 <- proc.time()[[3]]; r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x)); t1 <- proc.time()[[3]] - t0
  g <- rev_gp_read(r$fit)
  s <- cf_setup(y, x)
  m <- cf_max(s, extra = list(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd)))), ngrid = 30L)
  ef <- cf_ll_mpfr(s, g$beta, g$sigma, g$sd, g$rho); em <- cf_ll_mpfr(s, m$beta, m$sigma, m$sd, m$rho)
  fp <- function(p) { v <- cf_prof(s, p[1], p[2])$ll; if (is.finite(v)) -v else 1e300 }
  lo <- nlminb(c(log(g$rho), 2 * (log(g$sigma) - log(g$sd))), fp, control = list(rel.tol = 1e-13, iter.max = 400))
  t0 <- proc.time()[[3]]; mo <- marg_opt(r$fit$frame); tm <- proc.time()[[3]] - t0
  r2 <- rev_fit(y ~ gp(x), data.frame(y = y, x = x), start = mo$start)
  g2 <- rev_gp_read(r2$fit)
  e2 <- cf_ll_mpfr(s, g2$beta, g2$sigma, g2$sd, g2$rho)
  cat(sprintf(paste0("RLR %-26s frm ll %.5f sd %.4g rho %.4g sigma %.4g conv %d w %d %.1fs | 256 bits at frm's %.5f at cf's max %.5f (sd %.4g rho %.4g sigma %.4g) gap %.4f | local polish of the marginal from frm's point %.5f | ",
                     "marginal start %.1fs: ll %.5f sd %.4g conv %d w %d, 256 bits %.5f gap %.2e\n"),
              r0$lab, g$ll, g$sd, g$rho, g$sigma, g$conv, length(r$warn), t1, ef, em, m$sd, m$rho, m$sigma, ef - em, -lo$objective, tm, g2$ll, g2$sd, g2$conv, length(r2$warn), e2, e2 - em))
  flush.console()
}
cat("DONE\n")
