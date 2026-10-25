# Reviewer, fix round 1, F2: the sampler's marginal of a bounded seam
# parameter whose likelihood ignores it, under the default prior (want
# uniform(lo, hi)) and under normal(0.3, 0.2) on the natural scale (want
# that normal truncated to (lo, hi)). Seeds below.
# Rscript dev/review-sarfcor-fr1-sample.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb.spatial); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.sample"), "\n")
frmtmb_local_registries(environment())
lo <- -0.4; hi <- 0.9
frmtmb_register_autocor(
  "flatt", args = "M",
  build = function(spec, rows, n_all, lookup) list(n = length(rows)),
  whiten = function(y, mu, sigma, theta, aux) {
    list(u = (y - mu) / sigma + 0 * theta[1], logjac = -length(y) * log(sigma))
  },
  unwhiten = function(u, mu, sigma, theta, aux) mu + sigma * u,
  npar = 1L, sigma_varies = FALSE,
  natural = function(theta, aux) c(qq = lo + (hi - lo) * plogis(theta[1])),
  bounds = function(aux) cbind(lower = lo, upper = hi))
set.seed(12)
d <- data.frame(x = rnorm(40)); d$y <- 1 + d$x + rnorm(40)
fit <- frm(y ~ x + flatt(1), data = d)
print(frmtmb_autocor_pars(fit))
ks <- function(lab, x, cdf, ess) {
  # thin to roughly independent draws
  th <- x[seq(1, length(x), by = max(1, floor(length(x) / ess)))]
  k <- suppressWarnings(ks.test(th, cdf))
  cat(sprintf("%-26s draws %d thinned %d mean %.4f  KS D %.4f p %.3f\n", lab,
              length(x), length(th), mean(x), k$statistic, k$p.value))
}
msgs <- character(0)
s1 <- withCallingHandlers(
  frm_sample(fit, chains = 4, iter = 3000, warmup = 1000, seed = 21, refresh = 0),
  message = function(m) {msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")},
  warning = function(w) invokeRestart("muffleWarning"))
cat(grep("qq", unlist(strsplit(msgs, "\n")), value = TRUE), sep = "\n")
dr <- posterior::as_draws_df(as.matrix(s1))
e1 <- posterior::ess_bulk(dr$qq)
cat("ESS qq (default):", round(e1), "\n")
ks("default, want U(lo, hi)", dr$qq, function(q) punif(q, lo, hi), e1)
s2 <- suppressWarnings(suppressMessages(frm_sample(
  fit, prior = set_prior("normal(0.3, 0.2)", class = "qq"),
  chains = 4, iter = 3000, warmup = 1000, seed = 22, refresh = 0)))
dr2 <- posterior::as_draws_df(as.matrix(s2))
e2 <- posterior::ess_bulk(dr2$qq)
cat("ESS qq (normal):", round(e2), "\n")
Z <- pnorm(hi, 0.3, 0.2) - pnorm(lo, 0.3, 0.2)
ks("normal(0.3,0.2), want TN", dr2$qq,
   function(q) (pnorm(q, 0.3, 0.2) - pnorm(lo, 0.3, 0.2)) / Z, e2)
# wrong reference, to show the test has power: the untruncated normal
# without Jacobian on theta would look like this
ks("  (same draws vs U(lo, hi))", dr2$qq, function(q) punif(q, lo, hi), e2)
