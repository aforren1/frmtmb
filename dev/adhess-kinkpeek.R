# Lane adhess, fix round 1: the two sweep cases that keep finite SEs
# with no SE warning (dev/adhess-kinksweep.R), their estimates and SEs.
#   ADHESS_ARM=lane|base Rscript dev/adhess-kinkpeek.R
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
set.seed(1)
n <- 400
d <- data.frame(x = stats::rnorm(n))
d$y <- 1 + 0.5 * d$x + stats::rnorm(n)
d$c <- 0L
d$cen <- "right"
for (f in list(quote(frm(c ~ x, family = poisson(), data = d)),
               quote(frm(y | cens(cen) ~ x, data = d)))) {
  r <- cap(eval(f))
  fit <- r$v
  cat(deparse(f), "\n  est", format(fit$opt$par, digits = 5), "\n  SE ",
      format(sqrt(diag(suppressWarnings(vcov(fit, full = TRUE)))),
             digits = 5), "\n  grad", format(max(abs(fit$obj$gr(fit$opt$par))),
                                             digits = 3),
      "| logLik", format(as.numeric(logLik(fit)), digits = 8), "| warn",
      wtag(r$w), "\n")
}
