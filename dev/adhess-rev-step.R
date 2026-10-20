# Review of lane adhess, item 3: the RE step rule on columns RMS does
# not describe. Poisson GLMM with (1 | g), n = 400, 20 groups; per
# column type the fit's covariance (fit_outer_hessian()) against a
# Richardson Jacobian of the Laplace gradient (two base steps, whose
# agreement is "ref self"), plus an alternative rule, 1e-3 / max|x|.
#   REV_ARM=lane|base Rscript dev/adhess-rev-step.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
mk <- function(seed, xfun, b = 0.3, fam = "poisson") {
  set.seed(seed)
  n <- 400
  g <- factor(rep(1:20, length.out = n))
  x <- xfun(n)
  xc <- (x - mean(x)) / stats::sd(x)
  u <- stats::rnorm(20, 0, 0.4)
  eta <- 0.3 + b * xc + u[g]
  y <- if (fam == "poisson") stats::rpois(n, exp(eta)) else
    stats::rbinom(n, 1, stats::plogis(eta))
  data.frame(x = x, g = g, y = y, z = factor(sample(c("a", "b"), n, TRUE)))
}
heavy <- function(rms) function(n) {
  x <- exp(stats::rnorm(n, 0, 1.6))
  x * rms / sqrt(mean(x^2))
}
cases <- list(
  tiny_1e2 = list(x = function(n) stats::rnorm(n, 0, 1e-2)),
  tiny_1e4 = list(x = function(n) stats::rnorm(n, 0, 1e-4)),
  shift_1e3 = list(x = function(n) 1000 + stats::rnorm(n)),
  shift_5 = list(x = function(n) 5 + stats::rnorm(n, 0, 0.1)),
  heavy_rms0.9 = list(x = heavy(0.9)),
  heavy_rms3 = list(x = heavy(3)),
  heavy_rms50 = list(x = heavy(50)),
  x500 = list(x = function(n) stats::rnorm(n, 500, 100)),
  inter_x500 = list(x = function(n) stats::rnorm(n, 500, 100),
                    f = y ~ x * z + (1 | g)),
  dummy = list(x = function(n) stats::rnorm(n), f = y ~ z + x + (1 | g)),
  bern_x500 = list(x = function(n) stats::rnorm(n, 500, 100),
                   fam = "bernoulli")
)
only <- commandArgs(TRUE)
if (length(only)) cases <- cases[only]
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fam <- cs$fam %||% "poisson"
  for (seed in 1:3) {
    d <- mk(seed, cs$x, fam = fam)
    f <- cs[["f"]] %||% (y ~ x + (1 | g))
    r <- cap(frm(f, data = d, family = if (fam == "poisson") poisson() else
      bernoulli()))
    fit <- r$v
    if (inherits(fit, "rev_err")) {
      cat(sprintf("%-13s s%d ERROR %s\n", nm, seed, fit))
      next
    }
    h <- ns$fit_outer_hessian(fit)
    at <- h$at
    V <- suppressWarnings(vcov(fit, full = TRUE))
    Hr1 <- rich_hess(fit$obj, at, d = 1e-2)
    Hr2 <- rich_hess(fit$obj, at, d = 1e-3)
    Vr <- solve(Hr2)
    self <- cov_err(solve(Hr1), Vr)
    eH <- cov_err(solve(h$H), Vr)
    eV <- cov_err(V, Vr)
    alt <- NA
    stp <- NA
    if (arm != "base") {
      pn <- names(fit$opt$par)
      stp <- ns$outer_fd_steps(fit$frame, pn, fit$par_units %||% 1)
      X <- fit$frame$linpreds[[1]]$X
      mx <- apply(abs(as.matrix(X)), 2, max)
      st <- ns$obj_state_save(fit$obj)
      u <- fit$par_units
      if (is.null(u) || all(u == 1)) {
        ea <- rep(1e-3, length(at))
        ea[seq_along(mx)] <- 1e-3 / pmax(1, mx)
        Ha <- ns$fd_hessian(at, fit$obj$gr, eps = ea)$H
        alt <- cov_err(solve(Ha), Vr)[["se"]]
      }
      ns$obj_state_restore(fit$obj, st)
    }
    xm <- as.matrix(fit$frame$linpreds[[1]]$X)
    cat(sprintf(paste0("%-13s s%d units %-5s rms %s max %s | step %s | ",
                       "ref self %.1e | H se %.2e cor %.2e | vcov se %.2e ",
                       "| alt(max) se %.2e | lost %s warn %d\n"),
                nm, seed, !is.null(fit$par_units),
                paste(signif(sqrt(colMeans(xm^2)), 3), collapse = "/"),
                paste(signif(apply(abs(xm), 2, max), 3), collapse = "/"),
                paste(signif(stp, 3), collapse = "/"),
                self[["se"]], eH[["se"]], eH[["cor"]], eV[["se"]], alt,
                paste(names(sdr_of(fit)$se_lost %||% character(0)),
                      collapse = ",") %||% "-", length(r$w)))
  }
}
