# The post-fit outputs the hurdle refusals used to guard, on the lane:
# predict() (simulated proportions against fitted()), and frmtmb.sample's
# log_lik() and posterior_predict() at a draw against brms's density at
# the stored columns. Data seed 20261060, sampler seed 3.
# Output: dev/ordmix-log-hurdle-postfit.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
set.seed(20261060)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- rlogis(n) + 0.6 * d$x
d$y <- ifelse(runif(n) < 0.2, 0L,
              1L + (u > -1 + 0.1 * d$x) + (u > 0.3) + (u > 1.5 - 0.1 * d$x))
st <- set_prior("student_t(3, 0, 2.5)", class = "Intercept")
for (case in c("gr", "cs")) {
  f <- if (case == "gr") bf(y | thres(gr = g) ~ x, hu ~ z) else
    bf(y ~ cs(x), hu ~ z)
  fit <- frm(f, family = hurdle_cumulative("probit"), data = d, prior = st)
  P <- fitted(fit)[, "Estimate", ]
  set.seed(1)
  # the estimates alone, so the proportions scatter by the draws only
  pr <- predict(fit, ndraws = 4000, propagate_error = FALSE)
  # each cell's binomial standard error at 4000 draws; the largest of
  # 1500 standard normals is about 3.4
  zz <- (pr - P) / sqrt(pmax(P * (1 - P), 1e-12) / 4000)
  cat(sprintf(paste0("%s: predict() %d x %d, max |z| of the proportions ",
                     "against fitted() = %.2f over %d cells\n"),
              case, nrow(pr), ncol(pr), max(abs(zz)), length(zz)))
  ds <- suppressMessages(suppressWarnings(
    frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3)))
  ll <- log_lik(ds)
  e <- posterior_epred(ds)
  for (s in c(1L, nrow(ll))) {
    ref <- log(e[s, , ][cbind(seq_len(n), d$y + 1L)])
    cat(sprintf("%s draw %d: max |log_lik - log epred[y]| = %.3g\n", case, s,
                max(abs(ll[s, ] - ref), na.rm = TRUE)))
  }
  pp <- posterior_predict(ds)
  cat(sprintf("%s: posterior_predict %d x %d, codes %s\n", case, nrow(pp),
              ncol(pp), paste(sort(unique(as.vector(pp))), collapse = " ")))
}
