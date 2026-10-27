# Reviewer of lane loglik: every loo_*() of frmtmb.sample on a
# covariance block, at hand-set draws. Seed 96.
# Output: dev/review-loglik-log/loofuns.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
set.seed(96)
d <- data.frame(g = rep(1:20, each = 6), t = rep(1:6, 20))
d$x <- rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + rnorm(nrow(d))
K <- 200L
for (f in list(y ~ x + ar(t, g, cov = TRUE), y ~ x + cosy(t, g))) {
  fit <- frm(f, data = d, dry_run = "objective")
  lab <- c(brms_par_labels(fit), "lp__")
  m <- matrix(0, K, length(lab), dimnames = list(NULL, lab))
  m[, "b_Intercept"] <- rnorm(K, 1, 0.1)
  m[, "b_x"] <- rnorm(K, 0.5, 0.1)
  m[, "sigma"] <- exp(rnorm(K, 0, 0.05))
  m[, "thetaac_1"] <- rnorm(K, 0.3, 0.05)
  ds <- structure(list(stanfit = NULL, draws = m, fit = fit),
                  class = "frmtmb_draws")
  for (fn in c("loo", "waic", "loo_epred", "loo_linpred", "loo_predict",
               "loo_predictive_interval", "loo_R2", "loo_pit",
               "loo_moment_match", "reloo", "kfold")) {
    if (!exists(fn)) {
      cat(sprintf("%-28s %-24s not available\n", deparse1(f), fn))
      next
    }
    r <- tryCatch({
      v <- suppressWarnings(get(fn)(ds))
      paste("OK, finite:", all(is.finite(unlist(
        if (inherits(v, "loo")) v$estimates else v))))
    }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 120)))
    cat(sprintf("%-28s %-24s %s\n", deparse1(f), fn, r))
  }
}
