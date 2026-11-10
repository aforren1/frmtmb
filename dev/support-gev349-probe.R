# Lane support, part A: brms's row log_lik:349 (tests.log_lik.R, "log_lik
# for count and survival models works correctly", seed 233). The block's
# prep is rebuilt with the same random draws (its lines 236 to 254; no
# draw is taken between 254 and the GEV rows 342 to 348), and the GEV
# row is read the way the port's adapter reads it, row_lpdf() on
# doubles with exact = TRUE where the build has it, beside brms's
# dgen_extreme_value().
#
#   Rscript dev/support-gev349-probe.R   (FRMTMB_PORT_LIB picks the build)
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "\n")
inv_logit <- function(x) 1 / (1 + exp(-x))
set.seed(233L)
ns <- 25
nobs <- 10
trials <- sample(10:30, nobs, replace = TRUE)
dp <- list(eta = matrix(rnorm(ns * nobs), ncol = nobs),
           shape = rgamma(ns, 4), xi = runif(ns, -1, 0.5),
           phi = rgamma(ns, 1))
dp$sigma <- 1 / dp$shape
Y <- rbinom(nobs, size = trials, prob = rbeta(nobs, 1, 1))
i <- sample(nobs, 1)
dp$mu <- inv_logit(dp$eta)
# log_lik:273 (before the count rows): mu is exp(eta) from there on
dp$mu <- exp(dp$eta)
# log_lik:342 and :343: the identity link, so mu is read as it stands,
# and y = 0
y <- 0
mu <- dp$mu[, i]
fam <- gen_extreme_value()
ex <- if ("exact" %in% names(formals(frmtmb::row_lpdf))) list(exact = TRUE)
ll_f <- as.numeric(do.call(frmtmb::row_lpdf, c(list(fam, rep(y, ns),
  rep(y, ns), list(mu = mu, sigma = dp$sigma, xi = dp$xi), list(), NULL),
  ex)))
ll_b <- suppressWarnings(brms::dgen_extreme_value(
  y, mu = mu, sigma = dp$sigma, xi = dp$xi, log = TRUE))
t <- 1 + dp$xi * (y - mu) / dp$sigma
out <- t <= 0
cat(sprintf("draws %d, outside the support (t <= 0) %d\n", ns, sum(out)))
cat(sprintf("outside: frmtmb -Inf %d, finite %d, NaN %d; brms NaN %d\n",
            sum(ll_f[out] == -Inf), sum(is.finite(ll_f[out])),
            sum(is.nan(ll_f[out])), sum(is.nan(ll_b[out]))))
cat("outside, frmtmb:", format(ll_f[out], digits = 4), "\n")
ins <- !out
cat(sprintf("inside: %d draws, max rel diff frmtmb vs brms %.3g\n",
            sum(ins), max(abs(ll_f[ins] - ll_b[ins]) / abs(ll_b[ins]))))
# the probe of 0.71.0 (dev/rel071-gev-probe.R), at the point it names
p1 <- do.call(frmtmb::row_lpdf, c(list(fam, c(0, 0), c(0, 0),
  list(mu = c(1, -1), sigma = c(0.1, 1), xi = c(0.3, 0.3)), list(), NULL),
  ex))
cat("y = 0 at (mu 1, sigma 0.1, xi 0.3) and (mu -1, sigma 1, xi 0.3):",
    format(p1, digits = 6), "\n")
