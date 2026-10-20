# Reviewer of lane sampleprior, fix round 1: an independent reference
# for the negbinomial case. brms's priors read off brms's own Stan code
# (not default_prior()), the posterior of u = log(shape) by NESTED
# trapezoid quadrature in coordinates standardized at
# each u's conditional mode, not the lane's fixed grid.
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-agree-common.R")
cs <- agree_cases()$negbinomial
d <- cs$data
sc <- strsplit(as.character(brms::stancode(y ~ x, data = d,
                                           family = brms::negbinomial())), "\n")[[1]]
cat("brms lprior:", trimws(grep("lprior \\+=", sc, value = TRUE)), sep = "\n  ")
lpl <- grep("student_t_lpdf\\(Intercept", sc, value = TRUE)
nums <- as.numeric(strsplit(sub(".*Intercept \\| (.*)\\);.*", "\\1", lpl), ",")[[1]])
cat("Intercept prior read from the Stan code: student_t(", nums, ")\n")
xc <- d$x - mean(d$x)
lt <- function(v, df, mu, s) stats::dt((v - mu) / s, df, log = TRUE) - log(s)
# joint log density in (Ic, b) at u = log(shape), Jacobian of u included
ljoint <- function(a, b, u) {
  sh <- exp(u)
  mu <- exp(a + b * xc)
  sum(stats::dnbinom(d$y, size = sh, mu = mu, log = TRUE)) +
    lt(a, nums[1], nums[2], nums[3]) +
    0.4 * log(0.3) - lgamma(0.4) - 1.4 * u - 0.3 * exp(-u) + u
}
logm <- function(u) {
  f <- function(p) -ljoint(p[1], p[2], u)
  o <- stats::optim(c(log(mean(d$y)), 0), f, method = "BFGS", hessian = TRUE)
  L <- t(chol(solve(o$hessian)))
  top <- -o$value
  # trapezoid over a 101 x 101 grid of standardized coordinates, +-10
  zg <- seq(-10, 10, by = 0.2)
  Z <- as.matrix(expand.grid(zg, zg))
  P <- sweep(Z %*% t(L), 2, o$par, "+")
  eta <- P[, 1] + outer(P[, 2], xc)
  ll <- rowSums(matrix(stats::dnbinom(rep(d$y, each = nrow(P)), size = exp(u),
                                      mu = exp(eta), log = TRUE), nrow(P)))
  lj <- ll + lt(P[, 1], nums[1], nums[2], nums[3]) +
    0.4 * log(0.3) - lgamma(0.4) - 1.4 * u - 0.3 * exp(-u) + u
  mx <- max(lj)
  mx + log(sum(exp(lj - mx)) * 0.2^2) + log(abs(det(L)))
}
ug <- seq(-8, 40, by = 0.05)
t0 <- proc.time()[["elapsed"]]
lp <- vapply(ug, logm, 1)
w <- exp(lp - max(lp))
w <- w / sum(w)
m <- sum(w * ug)
s <- sqrt(sum(w * (ug - m)^2))
cdf <- cumsum(w)
q <- vapply(c(.05, .5, .95), function(p) stats::approx(cdf, ug, p, ties = "ordered")$y, 1)
cat(sprintf("reviewer reference: mean %.4f sd %.4f q05/50/95 %s | edge weights %.2g %.2g | %.0f s\n",
            m, s, paste(sprintf("%.3f", q), collapse = " "), w[1], w[length(w)],
            proc.time()[["elapsed"]] - t0))
cat("lane reference (nbref-lane.txt): mean 0.5116 sd 1.4362 q05/50/95 -1.010 0.218 2.819\n")
