# Reviewer of lane seflags: dev/seflags-rev-cost.R's glmm_1e5_year lost
# theta_1 as "flat" on both builds although the group sd is 0.7 over 500
# groups. Estimate, gradient, logLik against glmer(), and the same model
# with year centered.
#   Rscript dev/seflags-rev-year.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-common.R")
set.seed(1)
n <- 1e5
d <- data.frame(year = 2000 + sample(0:20, n, TRUE), z = rnorm(n),
                g = factor(sample(500, n, TRUE)))
d$y <- rbinom(n, 1, plogis(-0.5 + 0.05 * (d$year - 2010) + 0.3 * d$z +
                             rnorm(500, 0, 0.7)[d$g]))
d$yc <- d$year - 2010
for (fo in list(y ~ year + z + (1 | g), y ~ yc + z + (1 | g))) {
  r <- cap(frm(fo, family = bernoulli(), data = d))
  f <- r$v
  p <- f$opt$par
  cat(deparse(fo), "| code", f$opt$convergence, "| logLik",
      format(-f$opt$objective, digits = 10), "| par",
      paste(signif(p, 5), collapse = " "), "| max|grad|",
      signif(max(abs(f$obj$gr(p))), 3), "| lost", lost_of(f), "| warn",
      length(r$w), "\n")
  for (w in r$w) cat("  W:", substr(w, 1, 200), "\n")
}
m <- lme4::glmer(y ~ yc + z + (1 | g), family = binomial(), data = d)
cat("glmer (centered) logLik", format(as.numeric(logLik(m)), digits = 10),
    "sd", attr(lme4::VarCorr(m)$g, "stddev"), "\n")
