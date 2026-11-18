# Lane gpselimit: why gpselimit-exact.R's rep 3 seed 2 (exact gp(), length
# scale 0.18 of the smallest gap) loses theta_2 as "concave" with no
# condition. Prints the fit-time analysis, the cache's explanations and
# the conditions of summary().
# Usage: Rscript gpselimit-exact2.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-gpselimit/dev/gpselimit-lib.R")
a <- commandArgs(TRUE)
gs_libs(a[1])
rep <- 3; seed <- 2
set.seed(7700 + 10 * rep + seed)
np <- 40
xs <- sort(runif(np, 0, 10))
u <- rnorm(np)
x <- rep(xs, each = rep)
y <- rep(u, each = rep) + 0.3 * rnorm(np * rep)
fit <- withCallingHandlers(
  frm(bf(y ~ gp(x)), data = data.frame(x = x, y = y), family = gaussian()),
  warning = function(w) {
    cat("WARNING", class(w)[1], ":", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  },
  message = function(m) {
    cat("MESSAGE", class(m)[1], ":", conditionMessage(m), "\n")
    invokeRestart("muffleMessage")
  })
cat("se_fit_an lost:\n"); print(fit$cache$se_fit_an$lost)
cat("se_fit_an tier:", fit$cache$se_fit_an$tier, "\n")
cat("se_explained:", format(fit$cache$se_explained), "\n")
cat("se_deferred:", format(fit$cache$se_deferred), "\n")
cat("check_at_fit:", frmtmb:::se_check_at_fit(fit), "evals",
    format(fit$opt$evals), "\n")
cat("grad:", format(fit$obj$gr(fit$opt$par)), "\n")
withCallingHandlers(print(summary(fit)),
  warning = function(w) {
    cat("WARNING", class(w)[1], ":", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  })
print(fit$cache$sdr$se_lost)
H <- fit$cache$hessian_fixed$H
print(H)
