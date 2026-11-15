# Lane wcheck2, fix round F1: the check's time where the pair rule runs
# in full and refuses nothing. Weights uniform in (0.9, 1.1), so some
# uncensored weight is above 1 and the early exit does not apply, and
# normal(0, 1) on the sigma Intercept, so the one-cell rule passes the
# model on to the pair rule; every unit gains at most 0.1 per row and
# loses at least 0.9. The design of dev/reviews/wcheck2-review-time2.R
# (seed 3, sigma ~ 1 + fac). Each time is the minimum of 3 calls.
#   Rscript dev/wcheck2-time-pair.R <base|lane|lib>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
ns <- asNamespace("frmtmb")
pic <- set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma")
for (cfg in list(c(2000, 100), c(4000, 200), c(6000, 300),
                 c(10000, 1000))) {
  nn <- cfg[1]; L <- cfg[2]
  set.seed(3)
  dd <- data.frame(obs = factor(seq_len(nn)), x = rnorm(nn),
                   fac = factor(sample(sprintf("l%04d", seq_len(L)), nn,
                                       replace = TRUE)))
  dd$y <- rnorm(nn)
  dd$w <- runif(nn, 0.9, 1.1)
  fl <- bf(y | weights(w) ~ x + (1 | obs), sigma ~ 1 + fac)
  fr <- q(frm(fl, data = dd, dry_run = "frame"))
  tt <- vapply(1:3, function(i) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(ns$check_weights_unbounded(fr$spec, fr, pic),
                  error = function(e) "REFUSED")
    proc.time()[["elapsed"]] - t0
  }, 0)
  r <- tryCatch(ns$check_weights_unbounded(fr$spec, fr, pic),
                error = function(e) "REFUSED")
  cat(sprintf("weights (0.9, 1.1), n %d, %d levels: check %.2fs, %s\n", nn,
              L, min(tt), if (is.null(r)) "pass" else r))
}
cat("DONE\n")
