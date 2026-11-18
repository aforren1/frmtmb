# Review of lane squant: follow-up of squant-review-density.txt 2d/2b.
# (i) where the tape's d lF / d mu of quantile_elf()'s cdf_log goes to
# zero (x = 1 / (1 + exp(-t)) on the tape: its derivative
# e^-t / (1 + e^-t)^2 overflows the square once |t| > 354.9);
# (ii) elf_quantile() at tiny shapes against uniroot on the CDF;
# (iii) a censored fit with rows in that band: is the stop stationary?
#
# Run (worktree root):
#   Rscript dev/reviews/squant-review-tape.R > dev/reviews/squant-review-tape.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
dp <- function(mu, sigma, p) list(mu = mu, sigma = sigma, quantile = p)

## (i) the band ------------------------------------------------------------
cat("\n== (i) tape d lF/d mu against -f/F, lambda 0.25, sigma 0.7, p 0.5\n")
lam <- 0.25
fam <- quantile_elf(smoothing = lam)
for (t in c(-50, -300, -354, -355, -356, -400, -600, -689, -691, -800,
            355, 400, 689)) {
  q <- 0.3 + t * lam * 0.7
  which <- if (t < 0) "lF" else "lS"
  tp <- RTMB::MakeTape(function(m) {
    fam$cdf_log(q, dp(m, 0.7, 0.5), list())[[which]]
  }, 0.3)
  g <- as.numeric(tp$jacobian(0.3))
  v <- fam$cdf_log(q, dp(0.3, 0.7, 0.5), list())[[which]]
  an <- exp(fam$lpdf(q, dp(0.3, 0.7, 0.5), list()) - v) *
    (if (which == "lF") -1 else 1)
  cat(sprintf("  t %5g %s: value %.6f | tape d/dmu %.8g | analytic %.8g\n",
              t, which, v, g, an))
}
# the bare expression, on the tape
tp <- RTMB::MakeTape(function(t) 1 / (1 + exp(-t)), -400)
cat(sprintf("  d/dt of 1 / (1 + exp(-t)) at t = -400 on the tape: %.6g (true %.6g)\n",
            as.numeric(tp$jacobian(-400)), stats::dlogis(-400)))

## (ii) elf_quantile at tiny shapes ----------------------------------------
cat("\n== (ii) elf_quantile() against uniroot on the log CDF\n")
for (p in c(0.01, 0.5, 0.99)) {
  for (lam in c(0.001, 0.02, 0.05, 0.25)) {
    f <- quantile_elf(smoothing = lam)
    for (pr in c(0.1, 0.5, 0.9)) {
      got <- ns$elf_quantile(pr, 0, 1, p, lam)
      root <- stats::uniroot(function(z) {
        exp(f$cdf_log(z, dp(0, 1, p), list())$lF) - pr
      }, c(-2000, 2000), tol = 1e-13)$root
      if (abs(got - root) > 1e-6 * max(1, abs(root))) {
        cat(sprintf(paste0("  p %.2f lam %.3f pr %.1f: elf_quantile %.8g,",
                           " uniroot %.8g, qbeta(pr, a, b) %.6g\n"),
                    p, lam, pr, got, root,
                    stats::qbeta(pr, lam * (1 - p), lam * p)))
      }
    }
  }
}

## (iii) a censored fit with rows in the band ------------------------------
cat("\n== (iii) left-censored fit, smoothing 0.02: is the stop stationary?\n")
fd_grad <- function(fn, par) {
  vapply(seq_along(par), function(k) {
    h <- 1e-5 * max(1, abs(par[k]))
    u <- par
    d <- par
    u[k] <- u[k] + h
    d[k] <- d[k] - h
    (fn(u) - fn(d)) / (2 * h)
  }, 0)
}
for (s in 1:5) {
  for (lam in c(0.02, 0.25)) {
    set.seed(20261117 + s)
    n <- 300
    d <- data.frame(x = stats::runif(n))
    d$y <- 1 + 2 * d$x + 1.5 * stats::rt(n, 2)
    L <- -1
    d$c <- ifelse(d$y < L, -1, 0)
    d$yc <- pmax(d$y, L)
    ws <- character()
    fit <- withCallingHandlers(
      frm(bf(yc | cens(c) ~ x, quantile = 0.5),
          family = quantile_elf(smoothing = lam), data = d),
      warning = function(w) {
        ws <<- c(ws, conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    par <- fit$opt$par
    gt <- as.numeric(fit$obj$gr(par))
    gf <- fd_grad(fit$obj$fn, par)
    dpv <- lapply(ns$eval_dpars(fit)[[1L]], as.numeric)
    tt <- (L - dpv$mu) / (lam * dpv$sigma)
    band <- sum(d$c == -1 & tt < -354.9 & tt > -690)
    # optimize the same objective from the stop with the finite-difference
    # gradient only (the value fn is right)
    o <- stats::optim(par, fit$obj$fn, function(p) fd_grad(fit$obj$fn, p),
                      method = "BFGS", control = list(reltol = 1e-14,
                                                      maxit = 500))
    cat(sprintf(paste0("  seed %d lam %.2f: conv %d warn %d | censored %d, in",
                       " the band %d | max |tape grad| %.3g, max |fd grad|",
                       " %.3g | logLik %.6f, BFGS on fn from the stop %.6f |",
                       " coef %s -> %s\n"),
                s, lam, fit$opt$convergence, length(ws), sum(d$c == -1), band,
                max(abs(gt)), max(abs(gf)), -fit$opt$objective, -o$value,
                paste(sprintf("%.4f", par[1:2]), collapse = " "),
                paste(sprintf("%.4f", o$par[1:2]), collapse = " ")))
  }
}
cat("DONE\n")
