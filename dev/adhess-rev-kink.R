# Review of lane adhess: asym_laplace, whose log density is piecewise
# linear in mu, so its AD Hessian has no curvature along mu.
#   REV_ARM=lane|base Rscript dev/adhess-rev-kink.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-rev-common.R")
for (q in c("fixed", "free")) for (seed in 1:5) {
  set.seed(seed)
  n <- 400
  x <- rnorm(n)
  y <- 1 + 0.5 * x + rnorm(n)
  d <- data.frame(x, y)
  f <- if (q == "fixed") bf(y ~ x, quantile = 0.5) else bf(y ~ x)
  r <- cap(frm(f, family = asym_laplace(), data = d))
  fit <- r$v
  if (inherits(fit, "rev_err")) {
    cat(q, seed, "ERROR", fit, "\n")
    next
  }
  V <- suppressWarnings(vcov(fit, full = TRUE))
  at <- the_point(fit)
  st <- ns$obj_state_save(fit$obj)
  Hx <- fit$obj$he(at)
  ns$obj_state_restore(fit$obj, st)
  lost <- sdr_of(fit)$se_lost
  cat(sprintf(paste0("%-5s s%d code %d | grad %.2g | diag(AD H) %s | SE %s ",
                     "| lost %s | warn: %s\n"),
              q, seed, fit$opt$convergence,
              max(abs(fit$obj$gr(fit$opt$par))),
              paste(signif(diag(Hx), 3), collapse = "/"),
              paste(signif(sqrt(pmax(diag(V), 0)), 4), collapse = "/"),
              if (length(lost)) paste(names(lost), lost, sep = "=",
                                      collapse = ",") else "-",
              paste(substr(r$w, 1, 70), collapse = " || ")))
}
