# Lane arstep, fix round: how a constant sigma is held in the frame, and
# whether the base's gaussian marginal reads it (review open question 3).
# bf(y ~ gp(x), sigma = 0.3), n 50, seed 12. Args: arm
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
rev_libs(args[1])
ns <- asNamespace("frmtmb")
set.seed(12); n <- 50
x <- sort(runif(n, 0, 10)); y <- sin(x) + 0.3 * rnorm(n)
d <- data.frame(x = x, y = y)
for (sg in c(0.3, 0.05)) {
  r <- rev_fit(bf(y ~ gp(x), sigma = sg), d)
  fit <- r$fit
  fr <- fit$frame
  cat("SC sigma", sg, "betad_fixed_idx", fr$betad_fixed_idx, "betad template", fr$par_template$betad,
      "names(par)", names(fit$obj$par), "\n")
  sl <- Filter(function(l) identical(l$dpar, "sigma"), fr$linpreds)
  cat("SC sigma linpreds:", length(sl), if (length(sl)) paste(names(sl[[1]]), collapse = ","), "\n")
  if (length(sl)) cat("SC constant:", format(sl[[1]]$constant), "par", sl[[1]]$par, "\n")
  s <- ns$gp_marginal_setup(fr)
  cat("SC marginal setup null:", is.null(s), "has_sigma", s$has_sigma, "se2[1]", s$se2[1], "\n")
  m <- fit$obj$env$.frmtmb_gp_marginal
  cat(sprintf("SC proof stored %s, fit objective %.6f, ll %.6f, conv %d, warnings %d\n",
              if (is.null(m)) "none" else sprintf("%.6f", m), fit$opt$objective, as.numeric(logLik(fit)),
              fit$opt$convergence, length(r$warn)))
  if (!is.null(s)) {
    nll <- ns$gp_marginal_nll(s)
    pl <- fit$obj$env$parList(fit$opt$par)
    p <- list(beta = pl$beta, theta = pl$theta)
    if (s$has_sigma) p$betad <- pl$betad
    cat(sprintf("SC at the fit's optimum: marginal %.6f objective %.6f\n", nll(p), fit$obj$fn(fit$opt$par)))
  }
  # the true maximum from a few starts
  best <- as.numeric(logLik(fit))
  for (rho in c(0.3, 1, 3)) for (sd0 in c(0.5, 2)) {
    r2 <- rev_fit(bf(y ~ gp(x), sigma = sg), d, start = list(theta = c(log(sd0), log(rho))))
    if (!inherits(r2$fit, "error")) best <- max(best, as.numeric(logLik(r2$fit)))
  }
  cat(sprintf("SC sigma %g default ll %.6f best of 6 starts %.6f\n", sg, as.numeric(logLik(fit)), best))
}
cat("DONE\n")
