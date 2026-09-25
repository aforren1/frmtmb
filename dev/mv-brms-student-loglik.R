# Lane mv, 2026-09-25: student() rescor against brms 2.23.0 itself.
#
# brms is run with algorithm = "fixed_param" from an init equal to a
# frmtmb parameter point, so its one draw IS that point, and
# brms::log_lik() (brms's own R implementation of the multivariate t,
# log_lik_student_mv) is compared with frmtmb's taped objective there.
# brms centers the design, so its Intercept_<resp> is b0 + b1 * mean(x).
# Needs Stan: /opt/rlib/stan first on .libPaths(), /opt/r/bin on PATH.
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-mv")
.libPaths(c("/opt/rlib/stan", lib, "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
set.seed(20260925)
n <- 150
x <- rnorm(n)
E <- (matrix(rnorm(2 * n), n) %*% chol(matrix(c(1, 0.4, 0.4, 1), 2))) /
  sqrt(rchisq(n, 5) / 5)
dd <- data.frame(x = x, y1 = 1 + 0.5 * x + E[, 1], y2 = -1 + 0.8 * E[, 2])

fit <- frm(bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE) + student(), data = dd)
set.seed(3)
p <- fit$opt$par + rnorm(length(fit$opt$par), 0, 0.05)
pl <- fit$obj$env$parList(p)
bn <- names(fit$frame$par_template$betad)
bd <- function(nm) pl$betad[match(nm, bn)]
C <- us_chol_cor(pl$thetar, 2L)
mx <- mean(x)
init <- list(
  b_y1 = array(pl$beta[2]), Intercept_y1 = pl$beta[1] + pl$beta[2] * mx,
  sigma_y1 = exp(bd("y1_sigma_(Intercept)")),
  b_y2 = array(pl$beta[4]), Intercept_y2 = pl$beta[3] + pl$beta[4] * mx,
  sigma_y2 = exp(bd("y2_sigma_(Intercept)")),
  Lrescor = t(chol(C)), nu = 1 + exp(bd("nu_(Intercept)")))
suppressMessages(library(brms))
bfit <- brm(bf(mvbind(y1, y2) ~ x) + set_rescor(TRUE), data = dd,
            family = student(), chains = 1, iter = 1, warmup = 0,
            algorithm = "fixed_param", init = list(init), refresh = 0,
            seed = 1)
ll_brms <- sum(log_lik(bfit))
ll_frm <- -fit$obj$fn(p)
cat(sprintf("brms log_lik summed   %.12f\n", ll_brms))
cat(sprintf("frmtmb objective      %.12f\n", ll_frm))
cat(sprintf("relative difference   %.3e\n",
            abs(ll_brms - ll_frm) / abs(ll_brms)))
cat("brms draw nu:", as_draws_df(bfit)$nu, " frmtmb nu:", init$nu, "\n")
