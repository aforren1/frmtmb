# Reviewer, fix round 1, item 2: the sigma default now reads brms's
# innovation sd, a user sigma prior frmtmb's marginal sd. Does the
# sampling table, passed back as prior =, still give the default
# density? And does writing the default's own text by hand?
# Compares -log prior at three internal points (theta, log sigma).
# Usage: Rscript dev/review-acprior-r1-trip.R
# Log: dev/review-acprior-log/r1-trip.txt
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
fo <- y ~ x + ar(t, g, cov = TRUE)
fit <- frm(fo, data = d, dry_run = "objective")
tb <- default_prior(fo, data = d, route = "sample")
print(as.data.frame(tb)[tb$class %in% c("sigma", "ar"),
                        c("prior", "class", "lb", "ub")])
nlp <- function(pr) {
  ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pr))$ri
  f <- frmtmb:::neg_log_prior_fn(ri$entries)
  vapply(list(c(0.3, -0.5), c(2, 0.2), c(-1, 1)), function(z) {
    p <- fit$frame$par_template
    p$thetaac[] <- z[1]
    p$sigma[] <- z[2]
    f(p)
  }, 0)
}
cat("sigma internal name(s):", names(unlist(fit$frame$par_template)), "\n")
r0 <- nlp(NULL)
r1 <- nlp(tb)
r2 <- nlp(set_prior("student_t(3, 0, 2.5)", class = "sigma"))
cat("defaults          :", format(r0, digits = 8), "\n")
cat("table as prior =  :", format(r1, digits = 8), "\n")
cat("default text typed:", format(r2, digits = 8), "\n")
cat("max |table - defaults|", max(abs(r1 - r0)), "\n")
