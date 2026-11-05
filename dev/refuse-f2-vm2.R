# N1: dev/review-refuse-f1-vm2.R and dev/review-refuse-f1-vm.R's seed 11400 on
# the build with the window back. Same draws. The line arm is gone: frm() is
# the window. Reports shifts against an exact maximum on (-pi, pi), the range of
# the values, the warning, and L1 through the survey's instrument.
source("dev/review-refuse-lib.R")
rvm <- function(mu, k) {
  vapply(mu, function(m) { repeat { w <- runif(1, -pi, pi); if (runif(1) < exp(k * (cos(w) - 1))) break }; ((w + m + pi) %% (2 * pi)) - pi }, 0)
}
exact <- function(d, start) {
  rows <- which(is.na(d$y)); ob <- setdiff(seq_len(nrow(d)), rows)
  vals <- seq(-pi, pi, length.out = 1001); dg <- diff(vals)[1]
  ldv <- function(y, eta, lk) exp(lk) * cos(y - 2 * atan(eta)) - log(2 * pi) - log(besselI(exp(lk), 0))
  nll <- function(p) {
    ey <- p[1] + p[2] * d$x
    s <- sum(ldv(d$y[ob], ey[ob], p[3])) + sum(dnorm(d$z[ob], p[4] + p[5] * d$y[ob] + p[6] * d$x[ob], exp(p[7]), log = TRUE))
    for (i in rows) {
      li <- ldv(vals, ey[i], p[3]) + dnorm(d$z[i], p[4] + p[5] * vals + p[6] * d$x[i], exp(p[7]), log = TRUE)
      mx <- max(li); s <- s + mx + log(sum(exp(li - mx)) * dg)
    }
    if (!is.finite(s)) 1e10 else -s
  }
  o <- nlminb(start, nll, control = list(iter.max = 2000, eval.max = 4000, rel.tol = 1e-10))
  list(par = o$par, ll = -o$objective, se = sqrt(diag(solve(stats::optimHess(o$par, nll)))), fn = nll)
}
pars <- function(a) {
  e <- a$estimates; c(e$beta[grep("^y_", names(e$beta))], e$betad[grep("^y_", names(e$betad))],
                      e$beta[grep("z_.Intercept", names(e$beta))], e$beta[grep("miy", names(e$beta))], e$beta[grep("^z_x", names(e$beta))], e$betad[grep("^z_", names(e$betad))])
}
set.seed(11200 + 7)
n <- 200
for (r in 0:12) {
  x <- rnorm(n); eta <- 0.3 + 0.5 * x
  y <- rvm(2 * atan(eta), 3)
  d <- data.frame(x = x, y = y); d$z <- 0.5 + 0.7 * y + 0.3 * x + rnorm(n, 0, 0.7)
  rows <- sample(n, 60); d$y[rows] <- NA
  a <- tr(frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y) + x) + gaussian(), data = d))
  if (is_err(a)) { pr("r=%d %s", r, msg(a, 90)); next }
  pa <- pars(a)
  ex <- tr(exact(d, pa))
  if (is_err(ex)) { pr("r=%d exact %s", r, msg(ex, 60)); next }
  v <- frmtmb:::mi_values(a, "y")[sort(rows)]
  sh <- (pa - ex$par) / ex$se
  pr("VM2 r=%2d chart %s; values %.2f %.2f; max|shift| %.2f (zy %+.2f); logLik %.2f exact %.2f; zy %.3f exact %.3f; total %.2f warned %s",
     r, a$frame$mi_map$y$chart, min(v), max(v), max(abs(sh)), sh[5], as.numeric(logLik(a)), ex$ll, pa[5], ex$par[5],
     a$mi_laplace$total, any(grepl("Laplace approximation of the", attr(a, "w"))))
}
# seed 11400: three shapes, L1 at each
set.seed(11400)
n <- 150; x <- rnorm(n)
for (cfg in list(c(0, 3), c(3.05, 3), c(3.14, 20))) {
  th <- cfg[1] + 0.05 * x
  y0 <- rvm(th, cfg[2])
  d <- data.frame(x = x, y = y0)
  d$z <- 0.5 + 0.5 * y0 + rnorm(n, 0, 0.5)
  rows <- seq(5, n, by = 10); d$y[rows] <- NA
  fit <- tr(frm(bf(y | mi() ~ x) + von_mises() + bf(z ~ mi(y)) + gaussian(), data = d))
  if (is_err(fit)) { pr("mean %.2f kappa %g: %s", cfg[1], cfg[2], msg(fit, 150)); next }
  v <- frmtmb:::mi_values(fit, "y")[rows]
  pr("VM mean %.2f kappa %g: chart %s; latent %.3f %.3f (true %.3f %.3f); total %.2f, warned %s",
     cfg[1], cfg[2], fit$frame$mi_map$y$chart, min(v), max(v), min(y0[rows]), max(y0[rows]),
     fit$mi_laplace$total, any(grepl("Laplace approximation of the", attr(fit, "w"))))
  spec <- list(bf = "mvbf(bf(y | mi() ~ x, family = von_mises()), bf(z ~ mi(y)), rescor = FALSE)", family = NA, data = d)
  r <- tr(meas(spec))
  if (is_err(r)) pr("   L1: %s", msg(r, 160)) else pr("   L1: verdict %s const %s pts [%s] grad %s err %s", r$verdict, format(r$l1_const %||% NA, digits = 3),
    paste(format(unlist(r$l1_extra_const) %||% NA, digits = 3), collapse = ","), format(r$l1_grad %||% NA, digits = 3), substr(r$l1_error %||% "", 1, 110))
}
