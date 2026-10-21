# Which brms parameter carries the 0.18 gradient at frmtmb's
# hurdle_gamma estimate, and is it Stan's AD or the function? Finite
# differences of rstan::log_prob against rstan::grad_log_prob.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib", "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = "C:/Users/adf44/source/r/frmtmb-wt-cdfs/dev/stan-cache")
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
for (h in list.files("tests/testthat", "^helper.*[.]R$", full.names = TRUE)) sys.source(h, envir = globalenv())
run <- function(fam_f, fam_b, d) {
  fit <- frm(y | cens(cc, y2) ~ x, family = fam_f, data = d)
  r <- brms_lp_check(brms::bf(y | cens(cc, y2) ~ x), fam_b, d, fit, tol_grad = Inf)
  up <- rstan::unconstrain_pars(r$sf, r$pars)
  g <- rstan::grad_log_prob(r$sf, up, adjust_transform = FALSE)
  fd <- vapply(seq_along(up), function(j) {
    h <- 1e-6; e <- replace(numeric(length(up)), j, h)
    (rstan::log_prob(r$sf, up + e, adjust_transform = FALSE) -
       rstan::log_prob(r$sf, up - e, adjust_transform = FALSE)) / (2 * h)
  }, 0)
  print(rbind(stan_grad = g, fd_of_stan_lp = fd))
}
set.seed(104); n <- 200; keep <- rbinom(n, 1, 0.75)
y1 <- rpois(n, 3) * keep; y2 <- rnbinom(n, mu = 3, size = 2) * keep
y3 <- rbinom(n, 12, 0.35) * keep; y4 <- brms::rbeta_binomial(n, 12, 0.35, 5) * keep
y5 <- pmax(rpois(n, 3), 1) * keep; y6 <- pmax(rnbinom(n, mu = 3, size = 2), 1) * keep
y <- rgamma(n, 3, 3 / 2) * keep
set.seed(12)
d <- data.frame(y = y, x = rnorm(n), cc = sample(c(-1, 0, 0, 0, 1, 2), n, TRUE))
d$y2 <- d$y + 0.1; d$cc[d$y == 0] <- 0
cat("== hurdle_gamma, censored\n"); run(hurdle_gamma(), brms::hurdle_gamma(), d)
# the same rows with only one censoring code at a time
for (k in c(-1, 1, 2)) {
  dk <- d; dk$cc[dk$cc != k] <- 0
  cat("== hurdle_gamma, code", k, "only\n"); run(hurdle_gamma(), brms::hurdle_gamma(), dk)
}
