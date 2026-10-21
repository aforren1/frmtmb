# Reviewer, fix round 1: brms 2.23.0's discrete conventions. The SAME
# unshifted data frame goes to frmtmb and to brms's Stan program (flat
# priors, evaluated at frmtmb's estimates through brms_lp_check() of
# helper-brms.R): lp - logLik should be a constant (0 here) and brms's
# gradient zero. Then the edge cases, in both packages.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-cdfs"
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev/stan-cache"),
           NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages({library(testthat); library(frmtmb); library(brms)})
bf <- frmtmb::bf
options(frmtmb.notices = FALSE)
testthat::source_test_helpers(file.path(wt, "tests/testthat"), env = globalenv())
cat("lib:", find.package("frmtmb"), "\n")
set.seed(20261022)
n <- 160
x <- rnorm(n)
codes <- sample(c(-1, 0, 0, 0, 1, 2), n, TRUE)
mk <- function(y, top = Inf) {
  d <- data.frame(y = y, x = x, cc = codes)
  d$y2 <- pmin(d$y + 2, top)
  d$cc[d$cc == 2 & d$y2 <= d$y] <- 0
  # brms's P(Y > n) = 0 at a binomial's trials: not a censoring event
  if (is.finite(top)) d$cc[d$cc == 1 & d$y >= top] <- 0
  d
}
run <- function(lab, fb, fam_f, fam_b, d) {
  r <- tryCatch({
    fit <- frm(fb, family = fam_f, data = d)
    bb <- brms::bf(fb$formula)
    rr <- brms_lp_check(bb, fam_b, d, fit, tol_grad = Inf)
    sprintf("lp - logLik %.3e | brms max |grad| %.2e | logLik %.6f | codes %s", rr$lp - rr$ours,
            rr$max_grad, rr$ours, paste(table(d$cc), collapse = "/"))
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 150)))
  cat(sprintf("%-40s %s\n", lab, r))
}
mu <- exp(1 + 0.3 * x)
run("poisson, four codes", bf(y | cens(cc, y2) ~ x), poisson(), poisson(), mk(rpois(n, mu)))
run("negbinomial, four codes", bf(y | cens(cc, y2) ~ x), negbinomial(), brms::negbinomial(),
    mk(rnbinom(n, mu = mu, size = 2)))
run("geometric, four codes", bf(y | cens(cc, y2) ~ x), geometric(), brms::geometric(),
    mk(rnbinom(n, mu = mu, size = 1)))
db <- mk(rbinom(n, 12, plogis(0.3 * x)), top = 12); db$tr <- 12
run("binomial(12), four codes", bf(y | trials(tr) + cens(cc, y2) ~ x), binomial(), binomial(), db)
run("zero_inflated_poisson, four codes", bf(y | cens(cc, y2) ~ x), zero_inflated_poisson(),
    brms::zero_inflated_poisson(), mk(rpois(n, mu) * rbinom(n, 1, 0.8)))
run("hurdle_poisson, four codes", bf(y | cens(cc, y2) ~ x), hurdle_poisson(),
    brms::hurdle_poisson(), mk(ifelse(runif(n) < 0.2, 0, 1 + rpois(n, mu))))
# truncation composed with censoring: window [2, 9]
yt <- rpois(3 * n, rep(mu, 3)); xt <- rep(x, 3)
keep <- which(yt >= 2 & yt <= 9)[1:n]
dt <- data.frame(y = yt[keep], x = xt[keep], cc = codes)
dt$y2 <- pmin(dt$y + 2, 9)
dt$cc[dt$cc == 2 & dt$y2 <= dt$y] <- 0
dt$cc[dt$cc == 1 & dt$y >= 9] <- 0
dt$cc[dt$cc == -1 & dt$y <= 2] <- 0
run("poisson trunc(lb = 2, ub = 9) + cens", bf(y | cens(cc, y2) + trunc(lb = 2, ub = 9) ~ x),
    poisson(), poisson(), dt)
run("negbinomial trunc(lb = 2, ub = 9) + cens", bf(y | cens(cc, y2) + trunc(lb = 2, ub = 9) ~ x),
    negbinomial(), brms::negbinomial(), dt)
run("poisson trunc(lb = 2) alone", bf(y | trunc(lb = 2) ~ x), poisson(), poisson(), dt)
run("poisson trunc(lb = 0, ub = 9)", bf(y | trunc(lb = 0, ub = 9) ~ x), poisson(), poisson(), dt)

cat("\n== edge cases: frmtmb vs brms (standata)\n")
edge <- function(lab, f, fam_f, fam_b, d) {
  a <- tryCatch({ fit <- frm(f, family = fam_f, data = d)
                  sprintf("fits, logLik %.6f", as.numeric(logLik(fit))) },
                error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 110)))
  b <- tryCatch({ brms::standata(brms::bf(f), data = d, family = fam_b); "accepts" },
                error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 80)))
  cat(sprintf("%-44s frmtmb: %s\n%-44s brms:   %s\n", lab, a, "", b))
}
d0 <- data.frame(x = x, y = rpois(n, mu), cc = 0)
d0$y2 <- d0$y + 1
edge("poisson trunc(lb = 0), zeros present", y | trunc(lb = 0) ~ x, poisson(), poisson(), d0)
cat(sprintf("   untruncated logLik %.6f\n", as.numeric(logLik(frm(y ~ x, family = poisson(), data = d0)))))
edge("poisson trunc(lb = -2)", y | trunc(lb = -2) ~ x, poisson(), poisson(), d0)
d1 <- d0; d1$cc[1:5] <- 1; d1$y[1:5] <- d1$y[1:5] + 0.5
edge("poisson right-censored at 2.5 (fractional)", y | cens(cc) ~ x, poisson(), poisson(), d1)
edge("poisson trunc(lb = 0.5) (fractional)", y | trunc(lb = 0.5) ~ x, poisson(), poisson(), d0)
d2 <- d0; d2$cc[1:5] <- 2; d2$y2[1:5] <- d2$y[1:5]
edge("poisson interval with y2 == y", y | cens(cc, y2) ~ x, poisson(), poisson(), d2)
d3 <- d0; d3$cc[1:5] <- 2; d3$y2[1:5] <- d3$y[1:5] - 1
edge("poisson interval with y2 < y", y | cens(cc, y2) ~ x, poisson(), poisson(), d3)
d4 <- data.frame(x = x, y = rnorm(n), cc = 0); d4$y2 <- d4$y + 1
d4$cc[1:5] <- 2; d4$y2[1:5] <- d4$y[1:5] - 0.5
edge("gaussian interval with y2 < y", y | cens(cc, y2) ~ x, gaussian(), gaussian(), d4)
d5 <- data.frame(x = x, y = rbinom(n, 12, 0.5), tr = 12, cc = 0)
d5$y[1:3] <- 12; d5$cc[1:3] <- 1; d5$y2 <- d5$y
edge("binomial right-censored at its trials", y | trials(tr) + cens(cc) ~ x, binomial(), binomial(), d5)
