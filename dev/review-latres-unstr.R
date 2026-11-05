source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
set.seed(4242)
d <- expand.grid(t = 1:4, g = factor(1:8))
d$x <- rnorm(nrow(d))
e <- rnorm(nrow(d), 0, 0.5)
d$y <- rpois(nrow(d), exp(0.6 + 0.4 * d$x + e))
fit <- frm(y ~ x + unstr(t, g), data = d, family = poisson())
print(fit$estimates$theta)
print(names(fit$obj$par))
print(fit$frame$re_blocks[[1]]$theta_idx)
par <- fit$obj$env$last.par.best[fit$obj$env$lfixed()]
print(par)
f0 <- fit$obj$fn(par)
cat("at est", f0, "\n")
for (k in which(names(par) == "theta")) {
  p <- par
  p[k] <- p[k] + 0.5
  cat(k, fit$obj$fn(p), "\n")
}
C <- 0.3^abs(outer(1:4, 1:4, "-"))
ucc <- fi("us_chol_cor")
f <- function(th) sum((ucc(th, 4)[lower.tri(C)] - C[lower.tri(C)])^2)
o <- nlminb(rep(0, 6), f, control = list(rel.tol = 1e-20, x.tol = 1e-14))
print(o$par); print(o$objective); print(ucc(o$par, 4))
