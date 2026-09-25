# Punch round 1, minor 4: item 5's REML/profile standard errors against
# glm() per family. frmtmb uses observed information; glm() expected.
#   PREDFIX_ARM=lane Rscript dev/predfix-p1-glmse.R > dev/predfix-log/p1-glmse.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
set.seed(5)
n <- 200
x <- rnorm(n)
eta <- 0.3 + 0.4 * x
d <- data.frame(x = x, pois = rpois(n, exp(eta)), bern = rbinom(n, 1, plogis(eta)),
                geom = rnbinom(n, mu = exp(eta), size = 1), expo = rexp(n, 1 / exp(eta)))
ctl <- stats::glm.control(epsilon = 1e-12, maxit = 100)
cases <- list(
  poisson = list("pois", poisson(), function() glm(pois ~ x, stats::poisson(), d, control = ctl)),
  bernoulli = list("bern", bernoulli(), function() glm(bern ~ x, stats::binomial(), d, control = ctl)),
  geometric = list("geom", geometric(), function() glm(geom ~ x, MASS::negative.binomial(1), d, control = ctl)),
  exponential = list("expo", exponential(), function() glm(expo ~ x, stats::Gamma("log"), d, control = ctl)))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  g <- cs[[3]]()
  seg <- sqrt(diag(summary(g, dispersion = 1)$cov.scaled))
  for (m in c("REML", "profile")) {
    f <- frm(bf(stats::as.formula(paste(cs[[1]], "~ x"))), family = cs[[2]], data = d,
             REML = m == "REML", control = frmtmb_control(profile = m == "profile"))
    se <- sqrt(diag(vcov(f)))
    cat(sprintf("%-12s %-8s max relative SE difference from glm %.2e\n", nm, m,
                max(abs(se / seg - 1))))
  }
}
