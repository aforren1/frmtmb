# avg_slopes() on a gaussian exact gp() field (the filed MINOR), one arm,
# against the delta method by hand through frm_lp_basis() and against 150
# draws of the joint covariance. Data seed 8150 (the ncgp review's design),
# draw seed 8151.
#   MFX_ARM=before|after Rscript dev/mfxse-gpslope.R
source("C:/Users/adf44/source/r/frmtmb-wt-mfxse/dev/mfxse-lib.R")
set.seed(8150)
m <- 90
d <- data.frame(x = runif(m, 0, 6), z = runif(m, 0, 6))
d$y <- sin(d$x) + 0.3 * d$z + rnorm(m, sd = 0.3)
for (sp in c("y ~ z + gp(x)", "y ~ gp(x, gr = FALSE)", "y ~ gp(x, k = 15)")) {
  fit <- quiet(frm(as.formula(sp), data = d))
  sl <- quiet(avg_slopes(fit, variables = "x"))
  h <- 1e-4
  bh <- frm_lp_basis(fit, newdata = transform(d, x = x + h / 2))
  bl <- frm_lp_basis(fit, newdata = transform(d, x = x - h / 2))
  J <- (colMeans(as.matrix(bh$A)) - colMeans(as.matrix(bl$A))) / h
  ref <- sqrt(sum(J * (bh$V %*% J)))
  jc <- frmtmb:::get_joint_cov(fit)
  V <- as.matrix(jc$V)
  ee <- eigen((V + t(V)) / 2, symmetric = TRUE)
  Rt <- ee$vectors %*% diag(sqrt(pmax(ee$values, 0)), length(ee$values))
  idx <- split(seq_along(jc$names), factor(jc$names, levels = unique(jc$names)))
  center <- unlist(fit$estimates[names(idx)], use.names = FALSE)
  set.seed(8151)
  dr <- vapply(1:150, function(i) {
    p <- center + as.numeric(Rt %*% rnorm(length(center)))
    f <- fit
    for (k in names(idx)) f$estimates[[k]][] <- p[idx[[k]]]
    f$cache <- new.env(parent = emptyenv())
    mean((frm_linpred(f, newdata = transform(d, x = x + h / 2)) -
            frm_linpred(f, newdata = transform(d, x = x - h / 2))) / h)
  }, 0)
  cat(sprintf("GPSLOPE %s %-24s estimate %.5f  se %.3g  basis %.5f  draws %.5f\n",
              mfx_arm, sp, sl$estimate, sl$std.error, ref, sd(dr)))
}
