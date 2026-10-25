# Upstream repro (dev/upstream-bugs.md): brms 2.23.0's
# posterior_epred() and fitted() fail on every fit with
# sar(type = "error"). prepare_predictions_ac() reads both lagsar and
# errorsar with prepare_draws(); with no lagsar column that gives a
# matrix with no columns, not NULL, so posterior_epred_gaussian() takes
# the lag branch and solve() stops. posterior_predict() and log_lik()
# work. Brms only, no frmtmb. Rscript dev/sarfcor-upstream-epred.R
# (compiles one Stan program)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
k <- 5
g <- expand.grid(r = 1:k, c = 1:k)
n <- nrow(g)
A <- outer(1:n, 1:n, function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
# the symmetric normalization, whose eigenvalues are real, so that the
# fit itself runs (see the entry on complex eigenvalues)
dg <- 1 / sqrt(rowSums(A))
Ws <- A * outer(dg, dg)
set.seed(5)
d <- data.frame(x = rnorm(n))
d$y <- 1 + 0.5 * d$x + as.vector(solve(diag(n) - 0.4 * Ws, rnorm(n)))
fit <- brm(y ~ x + sar(Ws, type = "error"), data = d,
           data2 = list(Ws = Ws), chains = 1, iter = 400, refresh = 0,
           seed = 1)
pd <- prepare_predictions(fit)
cat("prep$ac$lagsar: NULL", is.null(pd$ac$lagsar), "| dim",
    paste(dim(pd$ac$lagsar), collapse = " x "), "\n")
for (call in c("posterior_epred(fit)", "fitted(fit)",
               "posterior_predict(fit)", "log_lik(fit)")) {
  r <- tryCatch({eval(parse(text = call)); "ok"},
                error = function(e) conditionMessage(e))
  cat(sprintf("%-24s %s\n", call, r))
}
