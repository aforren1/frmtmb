# Reviewer: frmtmb.sample on a sar() fit: what frm_sample() says about
# the default prior of lagsar, and whether a registered term whose
# natural parameter is not named lagsar/errorsar reaches the draws.
# Rscript dev/review-sarfcor-sample.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb.spatial); library(frmtmb.sample); library(Matrix)})
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.sample"), "\n")
frmtmb_local_registries(environment())
set.seed(3)
g <- expand.grid(r = 1:6, c = 1:6)
n <- nrow(g)
A <- outer(1:n, 1:n, function(i, j)
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
dg <- 1 / sqrt(rowSums(A)); W <- A * outer(dg, dg)
d <- data.frame(x = rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.5 * W, 1 + 0.5 * d$x + rnorm(n)))
fit <- frm(y ~ x + sar(W), data = d, data2 = list(W = W))
msgs <- character(0)
s <- withCallingHandlers(
  frm_sample(fit, chains = 1, iter = 400, warmup = 200, seed = 3, refresh = 0),
  message = function(m) {msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")},
  warning = function(w) invokeRestart("muffleWarning"))
cat("frm_sample messages:\n"); cat(paste0("  > ", msgs), sep = "")
cat("mentions lagsar/rho/thetaac in messages:",
    any(grepl("lagsar|thetaac|rho", msgs)), "\n")
cat("draws columns:", colnames(as.matrix(s)), "\n")
# prior class lagsar
r <- tryCatch({frm_sample(fit, prior = set_prior("uniform(-1, 1)", class = "lagsar"),
                          chains = 1, iter = 50, warmup = 25, refresh = 0); "accepted"},
              error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 200)))
cat("prior class lagsar:", r, "\n")

# a registered term whose natural parameter is 'toyphi'
ar1_build <- function(spec, rows, n_all, lookup) {
  n <- length(rows)
  list(S = sparseMatrix(i = 2:n, j = 1:(n - 1), x = 1, dims = c(n, n)),
       first = c(1, numeric(n - 1)))
}
ar1_whiten <- function(y, mu, sigma, theta, aux) {
  phi <- tanh(theta[1]); e <- (y - mu) / sigma
  u <- e - phi * as.vector(aux$S %*% e) + (sqrt(1 - phi^2) - 1) * e * aux$first
  lsig <- if (length(sigma) == 1L) length(y) * log(sigma) else sum(log(sigma))
  list(u = u, logjac = 0.5 * log(1 - phi^2) - lsig)
}
ar1_unwhiten <- function(u, mu, sigma, theta, aux) {
  phi <- tanh(theta[1]); u <- as.numeric(u); u[1] <- u[1] / sqrt(1 - phi^2)
  mu + sigma * as.numeric(stats::filter(u, phi, method = "recursive"))
}
frmtmb_register_autocor("toyar", args = "M", build = ar1_build,
                        whiten = ar1_whiten, unwhiten = ar1_unwhiten, npar = 1L,
                        natural = function(theta, aux) c(toyphi = tanh(theta[1])),
                        sigma_varies = FALSE)
ft <- frm(y ~ x + toyar(1), data = d)
cat("toy variables(fit):", variables(ft), "\n")
st <- suppressWarnings(suppressMessages(
  frm_sample(ft, chains = 1, iter = 300, warmup = 150, seed = 4, refresh = 0)))
cat("toy draws columns:", colnames(as.matrix(st)), "\n")
