# Release 0.70.0: is frmtmb or lmer() short on the probe's (1 | g)
# polynomials at x on [1, 2]? Both are compared with an exact profile
# maximum: for each (sigma, tau) the GLS coefficients come from a
# Householder QR of the whitened design and the marginal log-likelihood
# from the Cholesky factor of V, maximized over log sigma and log tau.
# frmtmb's and lmer's own parameters are also evaluated in that exact
# arithmetic, so a gap in a package's reported logLik shows up separately
# from a gap in its estimates.
#   Rscript dev/rel070-optrap-diag-llref.R <lib> <label>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", dirname(find.package("frmtmb")), "frmtmb",
    as.character(packageVersion("frmtmb")), "label", a[2], "\n")
ll_at <- function(y, X, Z, sig, tau, beta = NULL) {
  V <- sig^2 * diag(length(y)) + tau^2 * tcrossprod(Z)
  R <- chol(V)
  yw <- backsolve(R, y, transpose = TRUE)
  Xw <- backsolve(R, X, transpose = TRUE)
  # the residual is the trailing part of Q'y, which a LAPACK QR gives
  r <- if (is.null(beta)) {
    qr.qty(qr(Xw, LAPACK = TRUE), yw)[-seq_len(ncol(X))]
  } else yw - Xw %*% beta
  -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(r^2))
}
for (deg in 5:8) {
  set.seed(deg)
  d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
  d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                         "+ (1 | g)"))
  X <- cbind(1, as.matrix(d[paste0("x", 1:deg)]))
  Z <- model.matrix(~ 0 + g, d)
  f <- suppressWarnings(suppressMessages(frm(fo, data = d)))
  m <- suppressWarnings(suppressMessages(lmer(fo, data = d, REML = FALSE)))
  sf <- exp(f$estimates$betad[[1]]); tf <- exp(f$estimates$theta[[1]])
  vl <- as.data.frame(VarCorr(m))
  sl <- sigma(m); tl <- vl$sdcor[vl$grp == "g"]
  prof <- function(p) -ll_at(d$y, X, Z, exp(p[1]), exp(p[2]))
  o <- optim(log(c(sl, tl)), prof, method = "Nelder-Mead",
             control = list(reltol = 1e-14, maxit = 5000))
  o <- optim(o$par, prof, method = "BFGS", control = list(reltol = 1e-15))
  best <- -o$value
  fb <- suppressWarnings(fixef(f))[, "Estimate"]
  cat(sprintf(paste0("deg %d | exact max %.6f at sd %.6g/%.6g | frmtmb ",
                     "reported %.3g, own pars exact %.3g, own sds profiled ",
                     "%.3g | lmer reported %.3g, own pars exact %.3g, own ",
                     "sds profiled %.3g (all relative to the exact max)\n"),
              deg, best, exp(o$par[1]), exp(o$par[2]),
              as.numeric(logLik(f)) - best,
              ll_at(d$y, X, Z, sf, tf, fb) - best,
              ll_at(d$y, X, Z, sf, tf) - best,
              as.numeric(logLik(m)) - best,
              ll_at(d$y, X, Z, sl, tl, fixef(m)) - best,
              ll_at(d$y, X, Z, sl, tl) - best))
}
