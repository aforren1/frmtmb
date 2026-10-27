# Lane loglik: frmtmb's autocor_row_loglik() at the ML estimate against
# the leave-one-out conditional of a dense multivariate normal (or t,
# group by group) built directly from frmtmb's fitted covariance, for
# every residual covariance term, gaussian and student, with ragged
# groups, a missing time point and sigma ~ x. Seed 3.
# Output: dev/loglik-check.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.spatial)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spatial"), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-loglik/dev/loglik-dense.R")
set.seed(3)
sizes <- c(7, 5, 6, 4, 7, 3, 6, 5)
d <- do.call(rbind, lapply(seq_along(sizes), function(g) {
  data.frame(g = factor(g), t = seq_len(sizes[g]))
}))
d <- d[!(d$g == 2 & d$t == 3), ]   # one interior gap
d$x <- rnorm(nrow(d))
d$y <- 1 + 0.5 * d$x + rnorm(nrow(d)) + rnorm(8)[d$g]
d <- d[sample(nrow(d)), ]
n <- nrow(d)
V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)
P <- matrix(0, n, n)
P[cbind(1:(n - 1), 2:n)] <- 1
P <- P + t(P)
W <- P / rowSums(P)
forms <- list(
  "ar(cov)" = y ~ x + ar(t, g, cov = TRUE),
  "ar(p=2,cov)" = y ~ x + ar(t, g, p = 2, cov = TRUE),
  "ma(cov)" = y ~ x + ma(t, g, cov = TRUE),
  "arma(cov)" = y ~ x + arma(t, g, cov = TRUE),
  "cosy" = y ~ x + cosy(t, g),
  "unstr" = y ~ x + unstr(t, g),
  "ar(cov)+(1|g)" = y ~ x + ar(t, g, cov = TRUE) + (1 | g),
  "fcor" = y ~ x + fcor(V),
  "sar lag" = y ~ x + sar(W),
  "sar error" = y ~ x + sar(W, type = "error"))
for (nm in names(forms)) {
  for (fam in c("gaussian", "student")) {
    f <- frm(forms[[nm]], family = get(fam)(), data = d,
             data2 = list(V = V, W = W))
    fit <- f
    ll <- frmtmb::autocor_row_loglik(fit, eval_dpars(fit), "y")
    ref <- dense_loo(fit)
    nu <- eval_dpars(fit)[[1]]$nu
    cat(sprintf(paste0("%-16s %-8s max|frmtmb - dense| / max|dense| = ",
                       "%.3e  row sum %.4f  joint %.4f  nu %s\n"),
                nm, fam, max(abs(ll - ref)) / max(abs(ref)), sum(ll),
                dense_joint(fit),
                if (is.null(nu)) "-" else format(nu[1], digits = 4)))
  }
}
fh <- frm(bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x), data = d)
ll <- frmtmb::autocor_row_loglik(fh, eval_dpars(fh), "y")
cat(sprintf("%-16s %-8s max|frmtmb - dense| / max|dense| = %.3e\n",
            "ar sigma~x", "gaussian", max(abs(ll - dense_loo(fh))) /
              max(abs(dense_loo(fh)))))
fh <- frm(bf(y ~ x + fcor(V), sigma ~ x), data = d, data2 = list(V = V))
ll <- frmtmb::autocor_row_loglik(fh, eval_dpars(fh), "y")
cat(sprintf("%-16s %-8s max|frmtmb - dense| / max|dense| = %.3e\n",
            "fcor sigma~x", "gaussian", max(abs(ll - dense_loo(fh))) /
              max(abs(dense_loo(fh)))))
