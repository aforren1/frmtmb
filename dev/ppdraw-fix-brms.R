# Fix round: brms 2.23.0 posterior_predict() against frmtmb predict() on
# the cases the side law changes: zib zi ~ x (3269 data), hurdle_poisson
# hu ~ x with no zeros, and se(s, sigma = TRUE) with sigma at 0.
# Usage: Rscript dev/ppdraw-fix-brms.R <lib>
args <- commandArgs(TRUE)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(brms); library(frmtmb)})
cat("brms", as.character(packageVersion("brms")), " frmtmb",
    find.package("frmtmb"), "\n")
summ <- function(label, p) {
  w <- function(a, b) mean(apply(p, 2, stats::quantile, b) -
                             apply(p, 2, stats::quantile, a))
  cat(sprintf(paste("%-26s mean %.4f  P(y = 0) %.4f  90%% width %.4f",
                    " 98%% width %.4f\n"), label, mean(p), mean(p == 0),
              w(0.05, 0.95), w(0.01, 0.99)))
}
frm_pred <- function(fit, nd = NULL) {
  set.seed(1)
  withCallingHandlers(
    predict(fit, newdata = nd, ndraws = 4000, summary = FALSE),
    message = function(m) invokeRestart("muffleMessage"),
    warning = function(w) invokeRestart("muffleWarning"))
}
set.seed(7)
dz <- data.frame(x = rnorm(200), N = sample(5:15, 200, TRUE))
dz$y <- rbinom(200, dz$N, plogis(0.3 + 0.5 * dz$x))
cat("\n== A1 zib, zi ~ x\n")
b <- brm(brms::bf(y | trials(N) ~ x, zi ~ x), family =
           brms::zero_inflated_binomial(), data = dz, chains = 4,
         iter = 2000, seed = 11, refresh = 0)
set.seed(1)
summ("brms posterior_predict", posterior_predict(b))
cat(sprintf("brms posterior_epred mean %.4f\n", mean(posterior_epred(b))))
f <- suppressWarnings(frm(bf(y | trials(N) ~ x, zi ~ x),
                          family = zero_inflated_binomial(), data = dz))
summ("frmtmb predict", frm_pred(f))
cat(sprintf("frmtmb fitted mean %.4f\n", mean(fitted(f)[, "Estimate"])))

cat("\n== A4 hurdle_poisson, hu ~ x, no zeros\n")
set.seed(12)
dh <- data.frame(x = rnorm(300))
dh$y <- rpois(300, exp(2.5 + 0.2 * dh$x))
b <- brm(brms::bf(y ~ x, hu ~ x), family = brms::hurdle_poisson(),
         data = dh, chains = 4, iter = 2000, seed = 11, refresh = 0)
set.seed(1)
summ("brms posterior_predict", posterior_predict(b))
cat(sprintf("brms posterior_epred mean %.4f\n", mean(posterior_epred(b))))
f <- suppressWarnings(frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(),
                          data = dh))
summ("frmtmb predict", frm_pred(f))
cat(sprintf("frmtmb fitted mean %.4f\n", mean(fitted(f)[, "Estimate"])))

cat("\n== se(s, sigma = TRUE), sigma at 0; rows s = 0.05 and 5\n")
set.seed(7)
dg <- data.frame(x = rnorm(200), s = runif(200, 0.2, 1))
dg$y <- 1 + dg$x + rnorm(200, 0, dg$s)
nd <- data.frame(x = c(0, 0), s = c(0.05, 5))
b <- brm(brms::bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(),
         data = dg, chains = 4, iter = 2000, seed = 11, refresh = 0)
cat("brms sigma quantiles (0.5, 0.9, 0.99):",
    format(stats::quantile(as_draws_df(b)$sigma, c(0.5, 0.9, 0.99)),
           digits = 4), "\n")
set.seed(1)
pb <- posterior_predict(b, newdata = nd)
cat("brms posterior_predict column sd:", format(apply(pb, 2, sd), digits = 4),
    "\n")
f <- suppressWarnings(frm(bf(y | se(s, sigma = TRUE) ~ x),
                          family = gaussian(), data = dg))
pf <- frm_pred(f, nd)
cat("frmtmb predict column sd:", format(apply(pf, 2, sd), digits = 4), "\n")
