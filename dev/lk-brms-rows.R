# Candidate brms identity rows for the new links, run outside the tier
# file so a failure costs one row rather than the whole file.
suppressMessages({
  library(testthat); library(frmtmb); library(brms); library(rstan)
})
setwd("C:/Users/adf44/source/r/frmtmb-wt-links/tests/testthat")
for (h in c("helper-brms.R")) source(h)
options(frmtmb.brms_lp_report = TRUE)

run <- function(label, bform, family, data, fit, ...) {
  cat("\n#### ", label, "\n", sep = "")
  t0 <- Sys.time()
  r <- tryCatch(brms_lp_check(bform, family, data, fit, ...),
                error = function(e) {
                  cat("ERROR: ", conditionMessage(e), "\n", sep = ""); NULL
                })
  cat(sprintf("   %.0f s\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(r)
}

set.seed(17)
n <- 250
d <- data.frame(x = rnorm(n))

# 22a bernoulli(probit)
d$y <- rbinom(n, 1, pnorm(0.2 + 0.8 * d$x))
run("22a bernoulli(probit)", brms::bf(y ~ x), brms::bernoulli("probit"), d,
    frm(frmtmb::bf(y ~ x) + bernoulli("probit"), data = d))

# 22b bernoulli(probit_approx)
run("22b bernoulli(probit_approx)", brms::bf(y ~ x),
    brms::bernoulli("probit_approx"), d,
    frm(frmtmb::bf(y ~ x) + bernoulli("probit_approx"), data = d))

# 22c bernoulli(softit)
run("22c bernoulli(softit)", brms::bf(y ~ x), brms::bernoulli("softit"), d,
    frm(frmtmb::bf(y ~ x) + bernoulli("softit"), data = d))

# 22d binomial(cauchit) with trials
db <- data.frame(x = rnorm(n), nt = sample(3:12, n, TRUE))
db$y <- rbinom(n, db$nt, pcauchy(0.3 + 0.7 * db$x))
run("22d binomial(cauchit)", brms::bf(y | trials(nt) ~ x),
    brms::brmsfamily("binomial", link = "cauchit"), db,
    frm(frmtmb::bf(y | trials(nt) ~ x) + stats::binomial(link = "cauchit"),
        data = db))

# 22e beta(probit)
dbe <- data.frame(x = rnorm(n))
mu <- pnorm(0.2 + 0.5 * dbe$x)
dbe$y <- rbeta(n, mu * 8, (1 - mu) * 8)
run("22e Beta(probit)", brms::bf(y ~ x), brms::Beta("probit"), dbe,
    frm(frmtmb::bf(y ~ x) + Beta("probit"), data = dbe))

# 22f poisson(sqrt)
dp <- data.frame(x = rnorm(n))
dp$y <- rpois(n, (2 + 0.4 * pmax(pmin(dp$x, 3), -3))^2)
run("22f poisson(sqrt)", brms::bf(y ~ x), stats::poisson(link = "sqrt"), dp,
    frm(frmtmb::bf(y ~ x) + stats::poisson(link = "sqrt"), data = dp))

# 22g negbinomial(softplus)
dn <- data.frame(x = rnorm(n))
dn$y <- rnbinom(n, size = 3, mu = log1p(exp(1.5 + 0.5 * dn$x)))
run("22g negbinomial(softplus)", brms::bf(y ~ x),
    brms::negbinomial("softplus"), dn,
    frm(frmtmb::bf(y ~ x) + negbinomial("softplus"), data = dn))

# 22h negbinomial(squareplus)
run("22h negbinomial(squareplus)", brms::bf(y ~ x),
    brms::negbinomial("squareplus"), dn,
    frm(frmtmb::bf(y ~ x) + negbinomial("squareplus"), data = dn))

# 22i inverse.gaussian(1/mu^2)
di <- data.frame(x = rnorm(n))
di$eta <- 2 + 0.3 * pmax(pmin(di$x, 3), -3)
di$y <- rgamma(n, 8, 8 * sqrt(di$eta))
run("22i inverse.gaussian(1/mu^2)", brms::bf(y ~ x),
    stats::inverse.gaussian(), di,
    frm(frmtmb::bf(y ~ x) + stats::inverse.gaussian(), data = di))

cat("\nDONE\n")
