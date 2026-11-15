# Lane gpfield fix round: gp_field_newton()'s field at a replicate's
# parameters against the tape's inner mode there (obj$fn()), on
# dev/gpfield-cost.R's rows (set.seed(9300 + n)), 20 parameter draws from
# vcov(full = TRUE) (set.seed(5)). Args: fam n
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
fam <- a[1]; n <- as.integer(a[2])
ns <- asNamespace("frmtmb")
set.seed(9300 + n)
d <- data.frame(x = sort(round(stats::runif(n, 0, 50), 4)))
eta <- sin(d$x / 3)
d$y <- switch(fam, gaussian = 2 * eta + stats::rnorm(n, 0, 0.3),
              poisson = stats::rpois(n, exp(1 + eta)),
              bernoulli = stats::rbinom(n, 1, plogis(2 * eta)))
ff <- switch(fam, gaussian = gaussian(), poisson = poisson(),
             bernoulli = bernoulli())
fit <- gf_fit(bf(y ~ gp(x)), data = d, family = ff)$fit
map <- ns$outer_par_map(fit)
v0 <- ns$fit_outer_vector(fit, map)
V <- vcov(fit, full = TRUE)
R <- chol(V)
gov <- ns$gp_drawn_b(fit)
plan <- ns$gp_field_draw_plan(fit, gov)
cat("newton plan", !is.null(plan$newton), "\n")
bk <- fit$frame$re_blocks[[1]]
set.seed(5)
for (s in 1:20) {
  v <- v0 + as.vector(t(R) %*% rnorm(length(v0)))
  fs <- ns$fit_set_outer(fit, v, map)
  zn <- ns$gp_field_newton(fs, plan$newton)$z
  par <- ns$outer_from_template(fs$estimates, fit$obj, fit$frame)
  fit$obj$fn(par)
  ze <- fit$obj$env$last.par[which(names(fit$obj$env$last.par) == "b")]
  th <- fs$estimates$theta
  fn <- ns$gp_field(zn, th, bk)
  fe <- ns$gp_field(ze, th, bk)
  sdpost <- sqrt(diag(solve(as.matrix(fit$obj$env$spHess(
    fit$obj$env$last.par, random = TRUE)))))
  cat(sprintf(paste0("NEWTON %s s %2d | max |f newton - f exact| %.3e | ",
                     "over the field's sd %.3e | z %.3e\n"), fam, s,
              max(abs(fn - fe)), max(abs(fn - fe)) / exp(th[1]),
              max(abs(zn - ze) / sdpost)))
}
