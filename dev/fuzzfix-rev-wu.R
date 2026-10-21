# Reviewer: does check_weights_unbounded() fire on a model whose weighted
# likelihood IS bounded, and does it miss ones that are not?
#   Rscript dev/fuzzfix-rev-wu.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")

run <- function(name, expr) {
  ws <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) e)
  if (inherits(r, "error")) {
    m <- conditionMessage(r)
    tag <- if (grepl("has no maximum", m)) "REFUSED" else "ERROR"
    cat(sprintf("ARM %-28s %s | %s\n", name, tag, substr(m, 1, 150)))
  } else {
    bd <- r$estimates$betad
    cat(sprintf("ARM %-28s ok logLik %.6g | betad %s | warnings %d%s\n",
                name, as.numeric(logLik(r)),
                paste(sprintf("%.4g", bd), collapse = " "), length(ws),
                if (length(ws)) paste0(" (", substr(ws[1], 1, 70), ")")
                else ""))
  }
  invisible(r)
}

set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), xs = rnorm(n),
                w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
cat(sprintf("sum(w) - n = %.4f, xs range %.3f %.3f\n", sum(d$w) - n,
            min(d$xs), max(d$xs)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, exp(0.3 * d$xs))
d$yg <- rgamma(n, shape = 3, rate = 3 / exp(eta))
d$k <- rpois(n, exp(eta))
d$yb <- rbinom(n, 1, plogis(eta))
d$nt <- 10L
d$kb <- rbinom(n, 10, plogis(eta))
d$lsig <- log(0.7) + 0 * d$x
d$ws <- d$w / mean(d$w)

G <- gaussian()
cat("--- controls: the refusal must fire\n")
run("gauss_obs", frm(bf(y | weights(w) ~ x + (1 | obs)), family = G,
                     data = d))
run("gauss_obs_sigx", frm(bf(y | weights(w) ~ x + (1 | obs), sigma ~ xs),
                          family = G, data = d))
run("gauss_obs_trunc", frm(bf(y | weights(w) + trunc(lb = -20) ~ x +
                                (1 | obs)), family = G, data = d))

cat("--- the spike direction is not in the dpar's span\n")
run("FP1_sigma_0+xs", frm(bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + xs),
                          family = G, data = d))
run("FP1_sigma_0+xs_unw", frm(bf(y ~ x + (1 | obs), sigma ~ 0 + xs),
                              family = G, data = d))
run("FP2_sigma_offset", frm(bf(y | weights(w) ~ x + (1 | obs),
                               sigma ~ 0 + offset(lsig)),
                            family = G, data = d))
run("FP2b_sigma_const", frm(bf(y | weights(w) ~ x + (1 | obs), sigma = 0.7),
                            family = G, data = d))

cat("--- bounded families and non-owned effects\n")
run("bern_obs", frm(bf(yb | weights(w) ~ x + (1 | obs)), family = bernoulli(),
                    data = d))
run("binom_obs", frm(bf(kb | weights(w) + trials(nt) ~ x + (1 | obs)),
                     family = binomial(), data = d))
run("pois_obs", frm(bf(k | weights(w) ~ x + (1 | obs)), family = poisson(),
                    data = d))
run("negbin_obs", frm(bf(k | weights(w) ~ x + (1 | obs)),
                      family = negbinomial(), data = d))
run("gauss_group", frm(bf(y | weights(w) ~ x + (1 | g)), family = G,
                       data = d))
run("gamma_group", frm(bf(yg | weights(w) ~ x + (1 | g)), family = Gamma("log"),
                       data = d))

cat("--- frequency weights on aggregated rows\n")
# 20 groups x 4 distinct rows, integer counts; the group effect is shared
dagg <- d; dagg$fw <- sample(1:3, n, replace = TRUE)
run("agg_group_fw", frm(bf(y | weights(fw) ~ x + (1 | g)), family = G,
                        data = dagg))
# each aggregated row its own effect: the same as replicated identical y
# sharing one effect, whose within-cell variance is 0
run("agg_obs_fw", frm(bf(y | weights(fw) ~ x + (1 | obs)), family = G,
                      data = dagg))
dis <- dagg[rep(seq_len(n), dagg$fw), ]
run("disagg_cell_unitw", frm(bf(y ~ x + (1 | obs)), family = G, data = dis))

cat("--- priors\n")
run("prior_sigma_t", frm(bf(y | weights(w) ~ x + (1 | obs)), family = G,
                         data = d,
                         prior = set_prior("student_t(3, 0, 2.5)",
                                           class = "sigma")))
run("prior_sigma_gamma2", frm(bf(y | weights(w) ~ x + (1 | obs)), family = G,
                              data = d,
                              prior = set_prior("gamma(2, 1)",
                                                class = "sigma")))
run("prior_b_only", frm(bf(y | weights(w) ~ x + (1 | obs)), family = G,
                        data = d, prior = set_prior("normal(0, 5)",
                                                    class = "b")))
run("prior_shape_gamma", frm(bf(yg | weights(w) ~ x + (1 | obs)),
                             family = Gamma("log"), data = d,
                             prior = set_prior("gamma(0.01, 0.01)",
                                               class = "shape")))

cat("--- sum(w) = n exactly\n")
run("gauss_obs_scaled", frm(bf(y | weights(ws) ~ x + (1 | obs)), family = G,
                            data = d))
run("gamma_obs_scaled", frm(bf(yg | weights(ws) ~ x + (1 | obs)),
                            family = Gamma("log"), data = d))
run("gauss_obs_scaleTRUE", frm(bf(y | weights(w, scale = TRUE) ~ x +
                                    (1 | obs)), family = G, data = d))

cat("--- families outside the table\n")
run("exgauss_obs", frm(bf(y | weights(w) ~ x + (1 | obs)),
                       family = exgaussian(), data = d))
d2 <- d; d2$y2 <- d$y + rnorm(n)
run("mv_gauss_obs", frm(mvbf(bf(y | weights(w) ~ x + (1 | obs)),
                             bf(y2 | weights(w) ~ x + (1 | obs))),
                        family = G, data = d2))

# The gaussian marginal with one effect per row has a closed form, so
# boundedness of FP1 does not rest on the optimizer: per row,
# (1 - w)/2 log(2 pi s^2) - log(w)/2 + log N(r; 0, tau^2 + s^2 / w).
cat("--- closed-form gaussian marginal, sigma = exp(gamma * xs)\n")
cf <- function(p, X = cbind(1, d$x), S = cbind(d$xs)) {
  r <- d$y - X %*% p[1:2]
  ls <- S %*% p[3:(2 + ncol(S))]
  tau2 <- exp(2 * p[length(p)])
  s2 <- exp(2 * ls)
  sum((1 - d$w) / 2 * log(2 * pi * s2) - log(d$w) / 2 +
        dnorm(r, 0, sqrt(tau2 + s2 / d$w), log = TRUE))
}
best <- -Inf
for (st in list(c(0, 0, 0, 0), c(0.5, 0.4, 2, -1), c(0.5, 0.4, -2, -1),
                c(0, 0, 5, 0), c(0, 0, -5, 0))) {
  o <- optim(st, function(p) -cf(p), method = "BFGS",
             control = list(maxit = 2000, reltol = 1e-14))
  cat(sprintf("  start %s -> logLik %.8f gamma %.5f conv %d\n",
              paste(st, collapse = ","), -o$value, o$par[3], o$convergence))
  best <- max(best, -o$value)
}
# the profile over gamma, tau and beta maximized at each gamma
for (gm in c(-50, -10, -3, 3, 10, 50)) {
  o <- optim(c(0.5, 0.4, -1), function(p) -cf(c(p[1:2], gm, p[3])),
             method = "BFGS", control = list(maxit = 2000))
  cat(sprintf("  gamma fixed %6.1f profile logLik %.6g\n", gm, -o$value))
}
cat(sprintf("  intercept direction check: with sigma ~ 1 + xs the same form at
  log sigma = -20: %.6g\n", cf(c(0.5, 0.4, -20, -1), S = cbind(1)) ))
cat("DONE\n")
