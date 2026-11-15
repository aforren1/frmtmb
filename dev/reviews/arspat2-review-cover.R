# Review of lane arspat2, attack 4: coverage of sar()'s rho interval at
# true rho = 0.9 (lag, row-standardized rook W), 200 seeds per lattice
# (5 x 5 as the lane, and 10 x 10). From each lane fit both intervals:
# the new one (confint_varcorr(), logit of (lo, hi)) and the 0.82.0 one,
# rebuilt as rho +/- z * Est.Error (the delta method on rho, which the
# lane says is unchanged). On the base build (REVLIB="") the script reads
# confint_varcorr() as the old interval, to check the rebuild.
source("dev/reviews/arspat2-review-common.R")
nseed <- as.integer(Sys.getenv("REV_NSEED", "200"))
rho <- as.numeric(Sys.getenv("REV_RHO", "0.9"))
ks <- as.integer(strsplit(Sys.getenv("REV_K", "5,10"), ",")[[1]])
z <- qnorm(0.975)
for (k in ks) {
  W <- rook(k)
  W <- W / rowSums(W)
  n <- nrow(W)
  res <- data.frame()
  for (seed in seq_len(nseed)) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n))
    d$y <- as.vector(solve(diag(n) - rho * W, 1 + 0.5 * d$x + rnorm(n)))
    msgs <- character(0)
    fit <- withCallingHandlers(
      frm(y ~ x + sar(W), data = d, data2 = list(W = W)),
      warning = function(w) {
        msgs <<- c(msgs, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    ac <- fit$frame$autocor[[1L]]
    cv <- quiet(confint_varcorr(fit))
    s <- summary(fit)$cor_pars
    res <- rbind(res, data.frame(
      seed = seed, est = cv$estimate, lwr = cv$lwr, upr = cv$upr,
      type = cv$type, ee = s[1, "Est.Error"],
      old_lwr = cv$estimate - z * s[1, "Est.Error"],
      old_upr = cv$estimate + z * s[1, "Est.Error"],
      lo = ac$aux$lo, hi = ac$aux$hi, conv = fit$opt$convergence,
      nmsg = length(msgs),
      msg = if (length(msgs)) substr(paste(msgs, collapse = " | "), 1, 80) else ""))
  }
  cov_new <- mean(res$lwr <= rho & rho <= res$upr, na.rm = TRUE)
  cov_old <- mean(res$old_lwr <= rho & rho <= res$old_upr, na.rm = TRUE)
  cat(sprintf("\n%d x %d lattice, %d seeds, true rho %.2f, hi %.4f\n", k, k,
              nseed, rho, res$hi[1]))
  cat(sprintf("  type: %s\n", paste(unique(res$type), collapse = ",")))
  cat(sprintf("  coverage: confint_varcorr %.3f (NA %d), raw-scale delta %.3f\n",
              cov_new, sum(is.na(res$lwr)), cov_old))
  cat(sprintf("  past hi: confint_varcorr %d, raw-scale delta %d\n",
              sum(res$upr >= res$hi, na.rm = TRUE),
              sum(res$old_upr >= res$hi, na.rm = TRUE)))
  cat(sprintf("  miss low (upr < 0.9): new %d old %d; miss high (lwr > 0.9): new %d old %d\n",
              sum(res$upr < rho, na.rm = TRUE), sum(res$old_upr < rho, na.rm = TRUE),
              sum(res$lwr > rho, na.rm = TRUE), sum(res$old_lwr > rho, na.rm = TRUE)))
  cat(sprintf("  mean width new %.4f old %.4f; est max %.6f; nonzero conv %d; fits with messages %d\n",
              mean(res$upr - res$lwr, na.rm = TRUE),
              mean(res$old_upr - res$old_lwr, na.rm = TRUE), max(res$est),
              sum(res$conv != 0), sum(res$nmsg > 0)))
  if (any(res$nmsg > 0)) print(utils::head(res[res$nmsg > 0, c("seed", "est", "msg")], 5))
  saveRDS(res, sprintf("dev/reviews/arspat2-review-out/cover-%d-%s-%s.rds", k, rho,
                       if (nzchar(lib)) "lane" else "base"))
}
