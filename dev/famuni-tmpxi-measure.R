# Does brms's tmp_xi / scale_xi() parameterization remove the
# divergences frm_sample() shows on a scalar-xi GEV (review fix 5)? A
# standalone tmbstan comparison on the review's data (seeds 1401, 1403,
# 1402; n = 40; y ~ 1), both arms with frmtmb's own GEV density and
# brms's priors on mu and sigma:
#   arm log1p: xi = expm1(eta), normal(0, 2.5) on xi with its Jacobian,
#              the frm_sample(prior = set_prior("normal(0, 2.5)",
#              class = "xi")) arm of the review;
#   arm tmp:   xi = scale_xi(tmp_xi, y, mu, sigma), normal(0, 2.5) on
#              tmp_xi, brms's parameterization and default.
# 4 chains x 2000 (1000 warmup), seed 1. Output: dev/famuni-out/tmpxi.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(pkgbuild::has_build_tools())
suppressMessages({
  library(frmtmb)
  library(tmbstan)
})
ns <- asNamespace("frmtmb")
fam <- gen_extreme_value()
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
divs <- function(sf) {
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  sum(vapply(sp, function(s) sum(s[, "divergent__"]), 0))
}
# brms's priors: Intercept student_t(3, median, 2.5); sigma
# student_t(3, 0, 2.5) on sigma > 0 (log scale sampled, Jacobian added)
lt <- function(x, m, s) RTMB::dt((x - m) / s, 3, log = TRUE) - log(s)
amax <- function(v) {
  m <- v[1]
  for (i in seq_along(v)[-1]) m <- ns$ad_max(m, v[i])
  m
}
amin <- function(v) {
  m <- v[1]
  for (i in seq_along(v)[-1]) m <- ns$ad_min(m, v[i])
  m
}
scale_xi <- function(tmp, y, mu, sigma) {
  x <- (y - mu) / sigma
  b1 <- -1 / amin(x)
  b2 <- -1 / amax(x)
  lb <- ns$ad_min(b1, b2)
  ub <- ns$ad_max(b1, b2)
  lb + (ub - lb) / (1 + exp(-tmp))
}
for (cfg in list(c(-0.45, 1401), c(0.2, 1403), c(-0.85, 1402))) {
  set.seed(cfg[2])
  y <- rg(40, 0, 1, cfg[1])
  loc <- round(median(y), 1)
  nll_log1p <- function(p) {
    mu <- p[1]; ls <- p[2]; e <- p[3]
    xi <- expm1(e)
    d <- list(mu = mu, sigma = exp(ls), xi = xi, .eta_sigma = ls)
    -(sum(fam$lpdf(y, d, list())) + lt(mu, loc, 2.5) +
        lt(exp(ls), 0, 2.5) + ls + RTMB::dnorm(xi, 0, 2.5, log = TRUE) + e)
  }
  nll_tmp <- function(p) {
    mu <- p[1]; ls <- p[2]; tmp <- p[3]
    xi <- scale_xi(tmp, y, mu, exp(ls))
    d <- list(mu = mu, sigma = exp(ls), xi = xi, .eta_sigma = ls)
    -(sum(fam$lpdf(y, d, list())) + lt(mu, loc, 2.5) +
        lt(exp(ls), 0, 2.5) + ls + RTMB::dnorm(tmp, 0, 2.5, log = TRUE))
  }
  cat(sprintf("\n=== xi true %.2f, seed %d\n", cfg[1], cfg[2]))
  for (arm in c("log1p", "tmp")) {
    f <- if (arm == "log1p") nll_log1p else nll_tmp
    obj <- RTMB::MakeADFun(function(pr) f(pr$p), list(p = c(mean(y), 0, 0)),
                           silent = TRUE)
    sf <- suppressWarnings(tmbstan(obj, chains = 4, iter = 2000, seed = 1,
                                   refresh = 0, init = "par"))
    sm <- rstan::summary(sf)$summary
    dr <- as.matrix(sf)[, 1:3]
    xi <- if (arm == "log1p") expm1(dr[, 3]) else {
      apply(dr, 1, function(r) scale_xi(r[3], y, r[1], exp(r[2])))
    }
    cat(sprintf(paste0("  %-5s divergences %4d of 4000 | max Rhat %.3f, ",
                       "min bulk n_eff %.0f | xi median %.3f (2.5%% %.3f, ",
                       "97.5%% %.3f)\n"),
                arm, divs(sf), max(sm[1:3, "Rhat"]), min(sm[1:3, "n_eff"]),
                median(xi), quantile(xi, 0.025), quantile(xi, 0.975)))
  }
}
