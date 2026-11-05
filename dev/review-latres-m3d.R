# Short check, item 5: the block at a limit with another coefficient that lost
# its SE for a reason of its own: arma on a random-intercept DGP. Is `ma[1]` named?
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
for (s in 1:25) {
  set.seed(s)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  d$x <- rnorm(nrow(d))
  d$y <- rpois(nrow(d), exp(.5 + .3 * d$x + rep(rnorm(10, 0, .5), each = 6)))
  ms <- character(0); ws <- character(0)
  f <- withCallingHandlers(frm(y ~ x + arma(t, g, cov = TRUE), family = poisson(), data = d),
    message = function(m) { ms <<- c(ms, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  bkx <- fi("frame_latres")(f$frame)[[1]]; a <- bkx$aux
  th <- f$estimates$theta
  lim <- fi("latres_block_limit")(f, bkx)$at
  if (!any(lim)) next
  sm <- summary(f)$cor_pars
  v <- sqrt(diag(as.matrix(vcov(f, full = TRUE))))[paste0("theta_", 1:3)]
  txt <- paste(c(ms, ws), collapse = " | ")
  cat(sprintf("seed %2d lim %s  ar %.4f ma %.3f sderr %.2g  SE chart %s  summary ma Est.Error %s  names ma in messages: %s\n", s, paste(as.integer(lim), collapse = ""),
              sm["ar[1]", 1], sm["ma[1]", 1], sm["sderr", 1], paste(signif(v, 2), collapse = ","), format(sm["ma[1]", 2], digits = 3), grepl("ma[1]", txt, fixed = TRUE)))
}
