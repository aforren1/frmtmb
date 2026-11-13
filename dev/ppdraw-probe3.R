# Why predict()'s mean leaves fitted()'s on some sweep fits
# (dev/ppdraw-sweep.R): the log-likelihood at the parameter draws, and
# predict()'s mean against the mean over the same draws of the expected
# response at each draw (fitted_point() at the drawn parameters). The
# second difference is the simulator's; the rest is the draw law's.
# Usage: Rscript dev/ppdraw-probe3.R <lib or "base">
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
n <- 300L
source("dev/ppdraw-sweep-gen.R")
cases <- list(c("hurdle_poisson", "interior"),
              c("hurdle_negbinomial", "interior"),
              c("mixture_gaussian", "boundary"),
              c("mixture_poisson", "boundary"),
              c("mixture_cumulative", "interior"),
              c("zero_inflated_poisson", "interior"))
B <- 400L
for (cs in cases) {
  fam <- cs[1]
  case <- cs[2]
  d <- gen(fam, case)
  s <- spec(fam)
  fit <- suppressMessages(suppressWarnings(frm(s[[1]], family = s[[2]],
                                               data = d)))
  cat("\n==", fam, case, "\n")
  V <- vcov(fit, full = TRUE)
  est <- fit$opt$par
  print(round(cbind(est = est, sd = sqrt(diag(V))), 3))
  set.seed(5)
  dr <- suppressMessages(frmtmb:::predict_par_drawer(fit, TRUE, B))
  f0 <- fit$obj$fn(est)
  ll <- vapply(seq_len(B), function(s) {
    v <- frmtmb:::fit_outer_vector(dr(s))
    tryCatch(fit$obj$fn(v) - f0, error = function(e) NA_real_)
  }, 0)
  p <- length(est)
  cat("drop in log-likelihood at the draws: quantiles 0.5, 0.9, 0.99:",
      format(quantile(ll, c(0.5, 0.9, 0.99), na.rm = TRUE), digits = 3),
      "; chi2_p / 2 says", format(qchisq(c(0.5, 0.9, 0.99), p) / 2,
                                    digits = 3), "\n")
  em <- vapply(seq_len(B), function(s) {
    e <- frmtmb:::fitted_point(dr(s), NULL, NULL, "response", NULL, NULL,
                               FALSE, NULL)
    if (is.matrix(e)) NA_real_ else mean(e)
  }, 0)
  set.seed(5)
  pp <- suppressMessages(suppressWarnings(
    predict(fit, ndraws = B, summary = FALSE)))
  fe <- suppressWarnings(fitted(fit))
  if (length(dim(fe)) == 3L) {
    # category shares: the draw-averaged probability of each category
    ep <- Reduce(`+`, lapply(seq_len(B), function(s) {
      frmtmb:::fitted_point(dr(s), NULL, NULL, "response", NULL, NULL,
                            FALSE, NULL)
    })) / B
    K <- ncol(ep)
    sh <- t(apply(pp, 2L, function(col) tabulate(col, K) / B))
    P <- fe[, "Estimate", ]
    cat("max |z| shares vs fitted:",
        round(max(abs((sh - P) / sqrt(P * (1 - P) / B))), 2),
        "; vs the draw-averaged probabilities:",
        round(max(abs((sh - ep) / sqrt(ep * (1 - ep) / B))), 2),
        "; mean |P - draw-averaged|:", signif(mean(abs(P - ep)), 3), "\n")
  } else {
    m <- rowMeans(pp)
    se <- stats::sd(m) / sqrt(B)
    cat("predict mean", signif(mean(m), 6), " fitted", signif(mean(fe[, 1]), 6),
        " draw-averaged epred", signif(mean(em), 6), "\n")
    cat("z vs fitted", round((mean(m) - mean(fe[, 1])) / se, 2),
        "; z vs draw-averaged epred",
        round((mean(m) - mean(em)) / sqrt(se^2 + var(em) / B), 2), "\n")
  }
}
