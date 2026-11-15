# Lane gpfield fix round: the parts of gp_field_newton() on the poisson
# n 300 fit of dev/gpfield-cost.R, 200 calls each, and how many Newton
# steps the 100 replicates of dev/gpfield-newton-prof.R take.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
ns <- asNamespace("frmtmb")
n <- 300
set.seed(9300 + n)
d <- data.frame(x = sort(round(stats::runif(n, 0, 50), 4)))
d$y <- stats::rpois(n, exp(1 + sin(d$x / 3)))
fit <- gf_fit(bf(y ~ gp(x)), data = d, family = poisson())$fit
bk <- fit$frame$re_blocks[[1]]
th <- fit$estimates$theta
tm <- function(lab, expr) {
  e <- substitute(expr)
  t0 <- proc.time()[[3]]
  for (i in 1:200) eval(e)
  cat(sprintf("BENCH %-22s %.2f ms\n", lab, 1000 * (proc.time()[[3]] - t0) / 200))
}
K <- ns$gp_corr(th, bk)
tm("gp_corr", ns$gp_corr(th, bk))
tm("chol 300", chol(K + diag(0.1, n)))
tm("eval_dpars b NULL", ns$eval_dpars(fit, b = NULL))
tm("row_lpdf", ns$row_lpdf(fit$spec$responses[[1]]$family, d$y, d$y,
                           list(mu = exp(rep(1, n))), list(), NULL))
tm("expand_b", ns$expand_b(fit$frame, fit$estimates$b, th))
# steps per replicate
map <- ns$outer_par_map(fit)
v0 <- ns$fit_outer_vector(fit, map)
R <- chol(vcov(fit, full = TRUE))
plan <- ns$gp_field_draw_plan(fit, ns$gp_drawn_b(fit))
cnt <- 0L
trace("chol", quote(cnt <<- cnt + 1L), print = FALSE, where = ns)
set.seed(5)
for (s in 1:100) {
  fs <- ns$fit_set_outer(fit, v0 + as.vector(t(R) %*% rnorm(length(v0))), map)
  ns$gp_field_newton(fs, plan$newton)
}
untrace("chol", where = ns)
cat("STEPS chol calls per replicate (steps + 1 for L)", cnt / 100, "\n")
