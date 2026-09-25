source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# M1: benign scale, engaged vs FALSE. Seeds 11..13.
mk <- function(seed, sx, ng = 20, per = 15) {
  set.seed(seed)
  g <- factor(rep(seq_len(ng), each = per))
  x <- rnorm(ng * per) * sx
  u0 <- rnorm(ng, 0, 0.8); u1 <- rnorm(ng, 0, 0.5)
  xs <- x / sx
  y <- 1 + 0.7 * xs + u0[g] + u1[g] * xs + rnorm(ng * per)
  data.frame(y, x, g)
}
ctlF <- frmtmb_control(autoscale = FALSE)
ctlT <- frmtmb_control(autoscale = TRUE)
for (sx in c(1, 0.04)) for (seed in 11:13) {
  d <- mk(seed, sx)
  f <- y ~ x + (1 + x | g)
  a <- fitw(f, data = d)
  b <- fitw(f, data = d, control = ctlF)
  t <- fitw(f, data = d, control = ctlT)
  summ(sprintf("sx=%g seed=%d default vs F", sx, seed), a, b)
  summ(sprintf("sx=%g seed=%d TRUE vs F", sx, seed), t, b)
  cat("   template move (default):", tpl_move(a), "\n")
  cat("   template move (TRUE):   ", tpl_move(t), "\n")
  if (!is.null(t$tpl)) {
    # the template scored on the ORIGINAL objective vs the final optimum
    fb <- b$fit
    cat("   theta tpl:", format(t$tpl$theta, digits = 10), "\n")
    cat("   theta F  :", format(fb$estimates$theta, digits = 10), "\n")
  }
}
str(names(b$fit$estimates))
