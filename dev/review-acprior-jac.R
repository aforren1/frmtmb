# Reviewer: is the sampled default density flat on the NATURAL
# coefficients, measured independently of core's ac_trans_logjac()?
# The lane's test adds core's own log-Jacobian to core's prior, which
# is circular. Here: -log prior(theta) from the resolved default, plus
# log |det d phi / d theta| by central differences of the map
# autocor_natural(); flat on phi means the sum is constant in theta.
# Usage: Rscript dev/review-acprior-jac.R
# Log: dev/review-acprior-log/jac.txt
.libPaths(c("C:/Users/adf44/source/r/wt-acprior-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(4)
d <- data.frame(t = rep(1:6, 12), g = rep(1:12, each = 6), x = rnorm(72))
d$y <- 1 + 0.5 * d$x + rnorm(72)
numjac <- function(f, th, h = 1e-5) {
  J <- vapply(seq_along(th), function(j) {
    e <- replace(numeric(length(th)), j, h)
    (f(th + e) - f(th - e)) / (2 * h)
  }, numeric(length(f(th))))
  as.matrix(J)
}
cases <- list(y ~ x + ar(t, g, cov = TRUE), y ~ x + ma(t, g, cov = TRUE),
              y ~ x + ar(t, g, p = 2, cov = TRUE),
              y ~ x + ar(t, g, p = 3, cov = TRUE),
              y ~ x + ma(t, g, q = 2, cov = TRUE),
              y ~ x + arma(t, g, cov = TRUE),
              y ~ x + arma(t, g, p = 2, q = 2, cov = TRUE),
              y ~ x + cosy(t, g))
for (fo in cases) {
  r <- tryCatch({
    fit <- frm(fo, data = d, dry_run = "objective")
    ac <- fit$frame$autocor[[1L]]
    ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, NULL))$ri
    nlp <- frmtmb:::neg_log_prior_fn(ri$entries)
    nat <- function(th) unlist(frmtmb:::autocor_natural(th, ac))
    set.seed(1)
    v <- vapply(1:6, function(i) {
      th <- stats::rnorm(ac$npar, 0, 1.2)
      p <- fit$frame$par_template
      p$thetaac[] <- th
      c(nlp(p), log(abs(det(numjac(nat, th)))))
    }, c(0, 0))
    s <- colSums(v)
    sprintf("spread %.2e of scale %.2f; natural at last point: %s",
            diff(range(s)), max(abs(v)),
            paste(signif(nat(stats::rnorm(ac$npar)), 4), collapse = " "))
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(deparse(fo), ":", r, "\n")
}
# before: the same quantity on rellib-r13 would be the spread of the
# log-Jacobian itself; printed by the lane's test, not repeated here.
# Is the box of an order-k region large enough? Sample internal points
# widely and report the largest |phi_j| against choose(k, j).
for (k in 2:4) {
  fit <- frm(as.formula(paste0("y ~ x + ar(t, g, p = ", k, ", cov = TRUE)")),
             data = d, dry_run = "objective")
  ac <- fit$frame$autocor[[1L]]
  set.seed(2)
  m <- t(vapply(1:20000, function(i) {
    unlist(frmtmb:::autocor_natural(stats::rnorm(k, 0, 4), ac))
  }, numeric(k)))
  cat(sprintf("order %d: max |phi_j| %s; box %s\n", k,
              paste(signif(apply(abs(m), 2, max), 5), collapse = " "),
              choose(k, k %/% 2)))
}
