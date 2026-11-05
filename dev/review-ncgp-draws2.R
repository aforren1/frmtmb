# Item 4, 5: draws read by hand. From the labelled columns alone (b_Intercept, sigma, sdgp,
# lscale, b[i]), with my own Cholesky over gp_order, the field, log_lik(), fitted(), and the
# kriging at new positions. A test that compares labels against the kernel, not against the
# same labels. Then a bound, an autocor chart, the formula route, loo and check_laplace.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
suppressPackageStartupMessages(library(frmtmb.sample))
ns <- asNamespace("frmtmb")
set.seed(3003)
n <- 30L
d <- data.frame(x = sort(runif(n, 0, 6)), z = rnorm(n)); d$y <- 1.5 + sin(2 * d$x) + 0.3 * d$z + rnorm(n, sd = 0.1)
fit <- suppressWarnings(frm(y ~ z + gp(x), d))
dr <- suppressWarnings(suppressMessages(frm_sample(fit, chains = 1, iter = 400, warmup = 200, seed = 5, refresh = 0)))
m <- dr$draws
bk <- Filter(function(k) identical(k$covstruct, "gp"), fit$frame$re_blocks)[[1]]
print(colnames(m)[!grepl("^b\\[", colnames(m))])
print(fit$frame$linpreds[[1]]$param_colnames)
xs <- sort(unique(d$x)); p <- bk$gp_order
fld <- function(i, X0 = NULL) {
  z <- m[i, paste0("b[", bk$b_idx, "]")]
  sd <- m[i, "sdgp_gpx"]
  # lscale is reported on scaled inputs: brms's scale; kernel's rho is lscale * div
  div <- bk$gp_lscale_div %||% 1
  rho <- m[i, "lscale_gpx"] * div
  K <- exp(-outer(xs, xs, "-")^2 / (2 * rho^2)); diag(K) <- diag(K) + 1e-12
  L <- t(chol(K[p, p]))
  f <- sd * as.vector(L %*% z)
  fs <- numeric(length(xs)); fs[p] <- f
  list(f = fs, rho = rho, sd = sd, K = K)
}
ll <- suppressWarnings(log_lik(dr))
err <- numeric(); errf <- numeric(); errk <- numeric(); ep <- suppressWarnings(fitted(dr, summary = FALSE))
x0 <- c(0.33, 2.2, 5.1)
ek <- suppressWarnings(fitted(dr, newdata = data.frame(x = x0, z = 0), summary = FALSE))
for (i in seq(1, nrow(m), length.out = 12)) {
  fi <- fld(i)
  eta <- m[i, "b_Intercept"] + m[i, "b_z"] * d$z + fi$f[match(d$x, xs)]
  lh <- dnorm(d$y, eta, m[i, "sigma"], log = TRUE)
  err <- c(err, max(abs(lh - ll[i, ])))
  errf <- c(errf, max(abs(eta - ep[i, ])))
  Ks <- exp(-outer(x0, xs, "-")^2 / (2 * fi$rho^2))
  mk <- m[i, "b_Intercept"] + as.vector(Ks %*% solve(fi$K, fi$f))
  errk <- c(errk, max(abs(mk - ek[i, ])))
}
cat(sprintf("RH labelled columns by hand: log_lik max diff %.2e, fitted %.2e, kriging mean (a draw has its residual, so read the spread) %.2e\n", max(err), max(errf), max(errk)))
cat(sprintf("RH intercept column mean %.3f (fixef %.3f), z %.3f (%.3f)\n", mean(m[, "b_Intercept"]), fixef(fit, flatten = TRUE)[1], mean(m[, "b_z"]), fixef(fit, flatten = TRUE)[2]))
rg <- suppressWarnings(tryCatch(ranef(dr), error = function(e) e))
cat("RH ranef(draws) class:", class(rg)[1], if (inherits(rg, "error")) conditionMessage(rg), "\n")
tr <- function(lab, e) cat(sprintf("RH %-28s %s\n", lab, tryCatch({ v <- suppressWarnings(suppressMessages(e)); paste(class(v)[1], paste(dim(v), collapse = "x"), paste(utils::head(signif(suppressWarnings(as.numeric(unlist(v))), 4), 4), collapse = " ")) }, error = function(err) paste("ERR", substr(conditionMessage(err), 1, 140)))))
tr("loo", loo(dr)$estimates)
tr("posterior_predict new", posterior_predict(dr, newdata = data.frame(x = x0, z = 0)))
tr("conditional_effects", conditional_effects(dr, effects = "x"))
tr("bayes_R2", bayes_R2(dr))
tr("VarCorr", VarCorr(dr))
tr("hypothesis", hypothesis(dr, "z > 0")$hypothesis$Estimate)
tr("prior_draws", { d0 <- frm_sample(fit, chains = 1, iter = 200, warmup = 100, seed = 5, refresh = 0, sample_prior = "yes"); colnames(d0$draws)[grepl("^prior", colnames(d0$draws))] })
tr("check_laplace", frmtmb.sample::check_laplace(fit, chains = 1, iter = 200, warmup = 100, seed = 5, refresh = 0))
cat("== a bound on sigma, tape holds b first\n")
pb <- set_prior("uniform(0, 0.2)", class = "sigma", ub = 0.2)
dd <- suppressWarnings(suppressMessages(tryCatch(frm_sample(fit, prior = pb, chains = 1, iter = 300, warmup = 150, seed = 6, refresh = 0), error = function(e) e)))
if (inherits(dd, "error")) cat("RH bound ERR", conditionMessage(dd), "\n") else cat(sprintf("RH bound: sigma range %.4f %.4f; expected at most 0.2\n", min(dd$draws[, "sigma"]), max(dd$draws[, "sigma"])))
cat("== autocor chart beside the field\n")
da <- data.frame(t = 1:40, x = sort(runif(40, 0, 6))); da$y <- sin(da$x) + as.numeric(arima.sim(list(ar = 0.5), 40, sd = 0.2))
fa <- suppressWarnings(frm(y ~ gp(x) + ar(time = t), da))
dra <- suppressWarnings(suppressMessages(tryCatch(frm_sample(fa, chains = 1, iter = 300, warmup = 150, seed = 7, refresh = 0), error = function(e) e)))
if (inherits(dra, "error")) cat("RH ar ERR", conditionMessage(dra), "\n") else {
  mc <- dra$draws
  print(round(apply(mc[, !grepl("^b\\[|^lp__", colnames(mc)), drop = FALSE], 2, median), 3))
  print(round(unlist(fa$estimates[c("beta", "betad", "thetaac")]), 3))
}
cat("== the formula route\n")
df1 <- suppressWarnings(suppressMessages(tryCatch(frm_sample(y ~ z + gp(x), data = d, chains = 1, iter = 300, warmup = 150, seed = 8, refresh = 0), error = function(e) e)))
if (inherits(df1, "error")) cat("RH formula ERR", conditionMessage(df1), "\n") else {
  print(round(apply(df1$draws[, c("b_Intercept", "b_z", "sigma", "sdgp_gpx", "lscale_gpx")], 2, median), 3))
}
cat("DONE\n")
