# Reviewer of lane drawdir, check 2: on the random-slope fits the
# per-direction reading fires on (design S0.1 / S0.2 of
# dev/predraw2-cases.R, 15 groups of 8), is the normal law wrong against
# a posterior, or only against the oracle (known-parameter) law?
# For each id: the new level's (x = 2) predictive under
#   norm   the normal law (the verdict's fire set FALSE: 0.83.0's draws)
#   lane   the lane's reweighted draws
#   mlik   MCMC (tmbstan, laplace = TRUE) on the fit's own objective, the
#          density the reweighting targets (flat in the outer parameters),
#          its draws put through predict() in place of the law's
#   brms   frm_sample() with brms's default priors (what the warning names)
# and the oracle law's coverage of each arm's central 90% and 98%
# intervals, with the widths. 4000 draws each, seed 1.
# Usage: Rscript dev/reviews/drawdir-review-mcmc.R <id> [...]
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("lib:", dirname(find.package("frmtmb")), "\n")
source("dev/drawdir-cases.R")
ns <- asNamespace("frmtmb")
real <- ns$draw_reweight
nd <- data.frame(x = 2, g = factor("new"))
summ <- function(lab, s, tr) {
  out <- character(0)
  for (a in c(0.1, 0.02)) {
    q <- stats::quantile(s, c(a / 2, 1 - a / 2), names = FALSE)
    out <- c(out, sprintf("cov%02d %.3f w%02d %6.3f", round(100 * (1 - a)),
                          stats::pnorm(q[2], tr$mean, tr$sd) -
                            stats::pnorm(q[1], tr$mean, tr$sd),
                          round(100 * (1 - a)), q[2] - q[1]))
  }
  cat(sprintf("  %-5s %s | mean %.3f sd %.3f\n", lab,
              paste(out, collapse = " "), mean(s), stats::sd(s)))
}
pp <- function(fit) {
  set.seed(1)
  suppressMessages(suppressWarnings(
    predict(fit, newdata = nd, allow_new_levels = TRUE, ndraws = 4000,
            summary = FALSE)))[, 1]
}
for (id in args) {
  case <- case_fit(id)
  fit <- case$fit
  tr <- case$truth
  cat("==", id, "truth mean", tr$mean, "sd", round(tr$sd, 4), "\n")
  print(round(fit$opt$par, 4))
  print(round(sqrt(diag(vcov(fit, full = TRUE))), 4))
  s_lane <- pp(fit)
  chk <- fit$cache[["draw_check"]]
  cat(sprintf("  verdict fire %s excess %.3f dir %.3f (%s)\n", chk$fire,
              chk$excess, chk$dir, chk$dir_par))
  chk$fire <- FALSE
  assign("draw_check", chk, envir = fit$cache)
  s_norm <- pp(fit)
  summ("norm", s_norm, tr)
  summ("lane", s_lane, tr)
  # the fit's own objective, sampled
  t0 <- proc.time()[[3]]
  sm <- suppressWarnings(tmbstan::tmbstan(
    fit$obj, chains = 4, iter = 3000, warmup = 2000, seed = 1,
    laplace = TRUE, init = list(fit$opt$par, fit$opt$par, fit$opt$par,
                                fit$opt$par), refresh = 0, cores = 1))
  M <- as.matrix(sm)
  M <- M[, setdiff(colnames(M), "lp__"), drop = FALSE]
  cat(sprintf("  tmbstan %.0f s, %d draws, divergent %d, max rhat %.3f\n",
              proc.time()[[3]] - t0, nrow(M),
              sum(rstan::get_divergent_iterations(sm)),
              max(rstan::summary(sm)$summary[, "Rhat"], na.rm = TRUE)))
  qq <- apply(M, 2, stats::quantile, c(0.005, 0.5, 0.995))
  print(round(qq, 3))
  Dm <- t(M)
  assignInNamespace("draw_reweight", function(object, dw) {
    stopifnot(nrow(dw$D) == nrow(Dm), ncol(dw$D) == ncol(Dm))
    Dm
  }, ns = "frmtmb")
  s_mlik <- pp(fit)
  assignInNamespace("draw_reweight", real, ns = "frmtmb")
  summ("mlik", s_mlik, tr)
  # brms's default priors
  t0 <- proc.time()[[3]]
  ds <- suppressMessages(suppressWarnings(
    frm_sample(fit, chains = 4, iter = 3000, warmup = 2000, seed = 1,
               refresh = 0, cores = 1)))
  # posterior_predict() on draws refuses a new level, so the new level's
  # predictive is drawn by hand from each posterior draw: y = b0 + 2 b1 +
  # N(0, sd0^2 + 4 sd1^2 + 4 cor sd0 sd1 + sigma^2)
  A <- as.matrix(ds)
  cat("  draws columns:", paste(head(colnames(A), 8), collapse = " "), "\n")
  col <- function(pat) A[, grep(pat, colnames(A))[1]]
  b0 <- col("^b_Intercept$")
  b1 <- col("^b_x$")
  s0 <- col("^sd_g__Intercept$")
  s1 <- col("^sd_g__x$")
  rho <- col("^cor_g__Intercept__x$")
  sg <- col("^sigma$")
  cat(sprintf("  frm_sample quantiles: sd_g__x %s | cor %s\n",
              paste(round(stats::quantile(s1, c(0.01, 0.5, 0.99)), 3),
                    collapse = " "),
              paste(round(stats::quantile(rho, c(0.01, 0.5, 0.99)), 3),
                    collapse = " ")))
  set.seed(1)
  s_brms <- b0 + 2 * b1 + stats::rnorm(length(b0)) *
    sqrt(s0^2 + 4 * s1^2 + 4 * rho * s0 * s1 + sg^2)
  cat(sprintf("  frm_sample %.0f s, %d draws\n", proc.time()[[3]] - t0,
              length(s_brms)))
  summ("brms", s_brms, tr)
}
cat("DONE\n")
