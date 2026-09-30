# Re-check: frmtmb.sample posterior_epred()/posterior_predict() fill of a
# missing newdata response under ar(cov = FALSE) (arma_cond_fill_epred),
# against brms's .predictor_arma() on the SAME draws. By value: brms run
# one draw at a time (S = 1) in draw order after the same set.seed(),
# which is frmtmb.sample's RNG order. By law: brms on all draws at once.
# Data seed 31, sampler seed 3, 1 chain, 4000 iter; RNG seed 2.
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-formrobust/dev/formrobust-rev-sample-arma.R")
eval(parse(text = src[1:(grep("^set.seed[(]2[)]", src)[1] - 1)]))
brms1 <- function(k, fill_only = TRUE) {
  pr <- prep; pr$ndraws <- 1; pr$dpars$sigma <- sig[k]
  sh <- pred_arma(eta[k, , drop = FALSE], ar = matrix(ar1[k], 1, 1),
                  Y = nd$y, J_lag = c(rep(1, nrow(nd) - 1), 0), fprep = pr)
  if (fill_only) return(sh[1, ])
  vapply(seq_len(ncol(sh)), function(i) rnorm(1, sh[1, i], sig[k]), 0)
}
set.seed(2); ef <- posterior_epred(ds, newdata = nd)
set.seed(2); eb <- t(vapply(seq_len(S), brms1, numeric(nrow(nd))))
cat("epred by value: identical", identical(unname(ef), unname(eb)),
    " max ulp", max(abs(ef - eb) / abs(eb)) / .Machine$double.eps, "\n")
set.seed(2); pf <- posterior_predict(ds, newdata = nd)
set.seed(2); pb <- t(vapply(seq_len(S), brms1, numeric(nrow(nd)),
                            fill_only = FALSE))
cat("predict by value: identical", identical(unname(pf), unname(pb)),
    " max abs diff", max(abs(pf - pb)), "\n")
set.seed(5)
shifted <- pred_arma(eta, ar = matrix(ar1, S, 1), Y = nd$y,
                     J_lag = c(rep(1, nrow(nd) - 1), 0), fprep = prep)
set.seed(6); ef2 <- posterior_epred(ds, newdata = nd)
rows <- 4:8
cat("epred sd brms (all draws):", format(apply(shifted, 2, sd)[rows],
                                         digits = 4), "\n")
cat("epred sd frmtmb.sample   :", format(apply(ef2, 2, sd)[rows],
                                         digits = 4), "\n")
cat("KS p rows 4..8:", format(vapply(rows, function(j) if (sd(shifted[, j]) == 0)
  NA_real_ else suppressWarnings(ks.test(ef2[, j], shifted[, j])$p.value), 0),
  digits = 3), "\n")
cat("rows 1..4 identical to full-response epred:",
    identical(ef2[, 1:4], posterior_epred(ds, newdata = d[d$g == "1", ])[, 1:4]),
    "\n")
