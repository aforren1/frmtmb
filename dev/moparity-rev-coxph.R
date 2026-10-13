# Reviewer: cox() ML slope on the lane build against survival::coxph's
# partial-likelihood slope, data of dev/moparity-cox-sample.R, seeds 1-10.
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
for (seed in 1:10) {
  set.seed(seed)
  n <- 300
  x <- rnorm(n)
  tt <- stats::rexp(n, exp(-0.5 + 0.7 * x))
  ct <- stats::rexp(n, 0.2)
  dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
  fit <- suppressWarnings(frm(bf(time | cens(cens) ~ x), family = cox(),
                              data = dd))
  cp <- survival::coxph(survival::Surv(time, 1 - cens) ~ x, data = dd)
  b <- fixef(fit)["x", ]
  cat(sprintf("seed %2d frm b_x %.4f (se %.4f) | coxph %.4f (se %.4f) | diff/se %.2f | code %d\n",
              seed, b[1], b[2], coef(cp), sqrt(vcov(cp)), (b[1] - coef(cp)) / sqrt(vcov(cp)),
              fit$opt$convergence))
}
