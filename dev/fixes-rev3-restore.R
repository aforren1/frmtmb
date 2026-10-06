# Reviewer of lane fixes, final check: does nl_flat_message() leave the
# fit's state as it found it? The same identified nonlinear fits, with
# the check and with it replaced by a no-op, then every downstream
# output compared bitwise. Also run on rellib-r5, which has no check.
#   Rscript dev/fixes-rev3-restore.R <lib> <arm: check|noop> <out.rds>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
if (identical(a[2], "noop") && exists("nl_flat_message", ns)) {
  assignInNamespace("nl_flat_message", function(obj, opt, frame) NULL,
                    ns = "frmtmb")
}
cat("LIB", find.package("frmtmb"), "arm", a[2], "\n")
set.seed(77)
n <- 240
d <- data.frame(x = runif(n, 0, 3), z = rnorm(n),
                g = factor(rep(1:12, 20)))
ag <- rnorm(12, 0, 0.3)
d$y <- (2 + ag[d$g]) * exp((0.4 + 0.15 * d$z) * d$x) + rnorm(n, 0, 0.4)
d$c <- rpois(n, exp(0.5 + 0.6 * plogis(3 * (d$x - 1.5)) + ag[d$g]))
nd <- data.frame(x = c(0.5, 2), z = c(0, 1), g = factor(c(1, 3),
                                                       levels = 1:12))
models <- list(
  re_nl = list(bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1 + z, nl = TRUE),
               gaussian(), list(beta = c(2, 0.4, 0))),
  pois_nl = list(bf(c ~ b0 + h * plogis(s * (x - m)), b0 ~ 1 + (1 | g),
                    h ~ 1, s ~ 1, m ~ 1, nl = TRUE),
                 poisson(), list(beta = c(0.5, 0.6, 3, 1.5))))
out <- list()
for (k in names(models)) {
  m <- models[[k]]
  fit <- suppressMessages(frm(m[[1]], data = d, family = m[[2]],
                              start = m[[3]]))
  fit_se <- suppressMessages(frm(m[[1]], data = d, family = m[[2]],
                                 start = m[[3]], se = TRUE))
  e <- fit$obj$env
  out[[k]] <- list(
    par = fit$opt$par, obj = fit$opt$objective,
    lpb = e$last.par.best, vb = e$value.best,
    vfull = vcov(fit, full = TRUE), vfull_se = vcov(fit_se, full = TRUE),
    fe = fixef(fit), re = ranef(fit),
    fitted = fitted(fit), pred = {set.seed(1); predict(fit, newdata = nd,
                                                      ndraws = 50)},
    sim = simulate(fit, nsim = 2, seed = 3),
    prof = suppressWarnings(confint(fit, parm = if (k == "re_nl") "b_Intercept" else "m",
                                    method = "profile")),
    boot = tryCatch({
      set.seed(9)
      b <- suppressWarnings(frm_bootstrap(fit, nsim = 4))
      b$t
    }, error = function(e) conditionMessage(e)))
}
saveRDS(out, a[3])
