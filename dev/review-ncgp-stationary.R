# Item 1: is the short stop a second stationary point of the marginal, or a stop
# of nlminb on the tape? For each design: frm()'s point P, the marginal's gradient
# at P by central differences, a local polish of the marginal from P, the profile of
# the marginal along sd, and frm() restarted at the marginal's maximum.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ll3 <- function(s, p) {
  # p = log sigma, log sd, log rho; mu profiled
  r <- exp(2 * (p[1] - p[2]))
  a <- cf_prof(s, p[3], log(r))
  if (!is.finite(a$ll)) return(-Inf)
  cf_ll_at(s, a$beta, exp(p[1]), exp(p[2]), exp(p[3]))
}
num_grad <- function(s, p, h = 1e-4) {
  vapply(seq_along(p), function(i) {
    e <- replace(numeric(3), i, h)
    (ll3(s, p + e) - ll3(s, p - e)) / (2 * h)
  }, 0)
}
designs <- list()
for (seed in 1:3) {
  designs[[length(designs) + 1L]] <- local({
    set.seed(1080 + seed - 1L + if (seed == 1) 0L else 500L)
    x <- sort(round(runif(80, 0, 10), 4))
    list(lab = sprintf("logx s=1e-3 n80 seed %d", seed), x = x, y = log1p(x) + 1e-3 * rnorm(80))
  })
}
for (cs in list(c(20, 30, 0.01), c(20, 30, 0.02), c(30, 10, 0.01), c(15, 40, 0.01))) for (seed in 1:3) {
  designs[[length(designs) + 1L]] <- local({
    set.seed(7000 + seed)
    u <- sort(round(runif(cs[1], 0, 6), 3))
    x <- rep(u, each = cs[2])
    list(lab = sprintf("rep %dx%d s=%g seed %d", cs[1], cs[2], cs[3], seed), x = x, y = sin(2 * x) + cs[3] * rnorm(length(x)))
  })
}
for (d in designs) {
  x <- d$x; y <- d$y
  r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x))
  f <- r$fit
  if (inherits(f, "error")) { cat("RS", d$lab, "ERROR\n"); next }
  g <- rev_gp_read(f)
  s <- cf_setup(y, x)
  P <- c(log(g$sigma), log(g$sd), log(g$rho))
  m <- cf_max(s, extra = list(c(P[3], 2 * (P[1] - P[2]))), ngrid = 28L)
  gr <- num_grad(s, P)
  tg <- f$obj$gr(f$opt$par)
  # local polish from P on the marginal in full coordinates
  o <- nlminb(P, function(p) { v <- ll3(s, p); if (is.finite(v)) -v else 1e300 }, control = list(rel.tol = 1e-14, iter.max = 500, eval.max = 800))
  # restart on the z tape at the marginal's max
  r2 <- rev_fit(y ~ gp(x), data.frame(y = y, x = x),
                start = list(theta = c(log(m$sd), log(m$rho)), betad = log(m$sigma), beta = m$beta[1]))
  g2 <- if (!inherits(r2$fit, "error")) rev_gp_read(r2$fit)
  cat(sprintf(paste0("RS %-24s frm ll %.5f sd %.4g rho %.4g sigma %.4g conv %d w %d msg %s | cf max %.5f sd %.4g rho %.4g | gap %.3e | ",
                     "cf grad at P (sigma,sd,rho) %.2e %.2e %.2e | tape grad max %.2e | local cf polish from P -> %.5f sd %.4g rho %.4g (gain %.3e) | ",
                     "tape from cf max: ll %.5f sd %.4g conv %d w %d\n"),
              d$lab, g$ll, g$sd, g$rho, g$sigma, g$conv, length(r$warn), gsub(" ", "_", f$opt$message), m$ll, m$sd, m$rho, g$ll - m$ll,
              gr[1], gr[2], gr[3], max(abs(tg)), -o$objective, exp(o$par[2]), exp(o$par[3]), -o$objective - g$ll,
              if (is.null(g2)) NA else g2$ll, if (is.null(g2)) NA else g2$sd, if (is.null(g2)) NA else g2$conv, length(r2$warn)))
  if (abs(g$ll - m$ll) > 1e-3) {
    # the profile of the marginal along sd, sigma and rho maximized at each
    sds <- exp(seq(log(min(g$sd, m$sd)) - 0.3, log(max(g$sd, m$sd)) + 0.3, length.out = 12))
    prof <- vapply(sds, function(sdv) {
      fp <- function(q) { v <- ll3(s, c(q[1], log(sdv), q[2])); if (is.finite(v)) -v else 1e300 }
      -nlminb(c(log(m$sigma), log(m$rho)), fp, control = list(rel.tol = 1e-12))$objective
    }, 0)
    cat("RSP sd  ", paste(signif(sds, 4), collapse = " "), "\n")
    cat("RSP ll-m", paste(sprintf("%.4f", prof - m$ll), collapse = " "), "\n")
  }
  flush.console()
}
cat("DONE\n")
