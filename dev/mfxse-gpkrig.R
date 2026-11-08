# How much of the average slope's variance on an exact gp() is the field's
# conditional (kriging) variance at the unseen positions x +- h/2, which no
# parameter holds: w' C w with C = gp_krig_cov() over those rows. Data seed
# 8150 as dev/mfxse-gpslope.R; lane library.
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
ns <- asNamespace("frmtmb")
set.seed(8150)
m <- 90
d <- data.frame(x = runif(m, 0, 6), z = runif(m, 0, 6))
d$y <- sin(d$x) + 0.3 * d$z + rnorm(m, sd = 0.3)
for (h in c(1e-2, 1e-3, 1e-4)) for (sp in c("y ~ z + gp(x)",
                                           "y ~ gp(x, gr = FALSE)")) {
  fit <- quiet(frm(as.formula(sp), data = d))
  nd <- rbind(transform(d, x = x + h / 2), transform(d, x = x - h / 2))
  lp <- fit$frame$linpreds[[1]]
  ed <- ns$lp_eta_design(fit, lp, nd, TRUE, FALSE)
  sp1 <- ed$sm_parts[[1]]
  C <- matrix(0, 2 * m, 2 * m)
  if (!is.null(sp1$krig)) {
    Ck <- ns$gp_krig_cov(sp1$krig)
    C[sp1$krig$rows, sp1$krig$rows] <- as.matrix(Ck)
  }
  w <- c(rep(1, m), rep(-1, m)) / (m * h)
  sl <- quiet(avg_slopes(fit, variables = "x"))
  vk <- sum(w * (C %*% w))
  cat(sprintf("GPKRIG h %g %-24s coef se %.5f  kriging sd %.5f  total %.5f\n",
              h, sp, sl$std.error, sqrt(max(vk, 0)),
              sqrt(sl$std.error^2 + max(vk, 0))))
}
