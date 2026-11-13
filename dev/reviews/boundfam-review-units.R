# Review of lane boundfam: se_edge_limit() reads natural units from
# fit$par_units whenever that is non-NULL, and from se_limit_units() only
# when it is NULL. A smooth (smooth_fx_units()) makes par_units non-NULL
# under autoscale = FALSE, leaving every other coefficient in raw units.
# Here the fit is started AT the maximum the autoscaled fit found, so the
# stop is the real maximum.
#
#   Rscript dev/reviews/boundfam-review-units.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-boundfam-lib",
            "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/reviews/boundfam-review-fp.R")
eval(parse(text = src[grep("^conds <- ", src):(grep("^ok <- ", src) - 1L)]))
for (sc in c(1e-5, 1e-6)) for (fam in c("poisson", "bernoulli")) {
  set.seed(9)
  n <- 300
  d <- data.frame(x = rnorm(n), z = runif(n))
  eta <- 0.5 * (fam == "poisson") + 0.3 * d$x + sin(2 * pi * d$z)
  d$y <- if (fam == "poisson") rpois(n, exp(eta)) else rbinom(n, 1, plogis(eta))
  d$xs <- d$x * sc
  famf <- if (fam == "poisson") poisson() else bernoulli()
  f1 <- suppressWarnings(frm(y ~ xs + s(z), family = famf, data = d))
  cat(sprintf("\n## %s sc=%g autoscaled fit: obj=%.8f xs=%.6g\n", fam, sc,
              f1$opt$objective, f1$estimates$beta[2]))
  st <- list(beta = f1$estimates$beta, theta = f1$estimates$theta)
  r <- conds(frm(y ~ xs + s(z), family = famf, data = d, start = st,
                 control = frmtmb_control(autoscale = FALSE)))
  show(sprintf("%s sc=%g autoscale=FALSE started at the max, y ~ xs + s(z)", fam, sc), r)
  cat("  par_units:", format(r$value$par_units, digits = 3), "\n")
  r2 <- conds(frm(y ~ xs, family = famf, data = d,
                  start = list(beta = f1$estimates$beta[1:2]),
                  control = frmtmb_control(autoscale = FALSE)))
  show(sprintf("%s sc=%g autoscale=FALSE started near the max, y ~ xs (no smooth)", fam, sc), r2)
  g <- glm(y ~ xs, family = if (fam == "poisson") poisson else binomial, data = d)
  cat(sprintf("  glm y ~ xs: xs=%.6g se=%.6g\n", coef(g)[2], sqrt(vcov(g)[2, 2])))
}
