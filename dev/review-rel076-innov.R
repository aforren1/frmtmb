# Review, cross-lane: under a first-order cov = TRUE ar(), ma() or
# arma() (acprior's innovation-sd reading of sigma) is frmtmb.sample's
# log_lik() (loglik's conditional) on the same sigma the sampler uses?
# For each structure, with the default priors and with a user prior on
# class sigma:
#   1. log_lik() against the dense leave-one-out conditional of each
#      group, with the draws' sigma as the MARGINAL sd and brms's
#      correlation of the structure;
#   2. the sampler's own log density (lp__, tmbstan has no bounds, so it
#      is the objective) minus the dense JOINT log-likelihood minus
#      brms's prior on the innovation sd sigma * c, minus the log
#      Jacobians of the chart (u -> natural) and of sigma_innov in
#      log sigma_marg: constant across draws if, and only if, the
#      sampler's density is the dense likelihood plus brms's prior at
#      the innovation sd;
#   3. loo() runs.
#   Rscript dev/review-rel076-innov.R <struct> > dev/review-rel076-out/innov-<struct>.txt
st <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), "\n")
set.seed(20261028)
ng <- 10; nt <- 6
d <- expand.grid(t = 1:nt, g = 1:ng)
d$x <- stats::rnorm(nrow(d))
mod <- switch(st, ar = list(ar = 0.7), ma = list(ma = 0.5),
              arma = list(ar = 0.6, ma = 0.3))
e <- unlist(lapply(1:ng, function(g) {
  as.numeric(stats::arima.sim(mod, nt, sd = 0.5))
}))
d$y <- 1 + 0.5 * d$x + e
fo <- switch(st, ar = y ~ x + ar(t, g, cov = TRUE),
             ma = y ~ x + ma(t, g, cov = TRUE),
             arma = y ~ x + arma(t, g, cov = TRUE))
corf <- function(L, ar, ma) {
  switch(st,
    ar = ar^L,
    ma = ifelse(L == 0, 1, ifelse(L == 1, ma / (1 + ma^2), 0)),
    arma = {
      r1 <- (1 + ar * ma) * (ar + ma) / (1 + 2 * ar * ma + ma^2)
      ifelse(L == 0, 1, r1 * ar^(L - 1))
    })
}
cfac <- function(ar, ma) {
  switch(st, ar = sqrt(1 - ar^2), ma = 1 / sqrt(1 + ma^2),
         arma = sqrt((1 - ar^2) / (1 + ma^2 + 2 * ar * ma)))
}
cond <- function(y, m, S) {
  vapply(seq_along(y), function(i) {
    o <- seq_along(y)[-i]
    A <- solve(S[o, o, drop = FALSE])
    cm <- m[i] + as.numeric(S[i, o, drop = FALSE] %*% A %*% (y[o] - m[o]))
    cv <- S[i, i] - as.numeric(S[i, o, drop = FALSE] %*% A %*% S[o, i])
    stats::dnorm(y[i], cm, sqrt(cv), log = TRUE)
  }, 0)
}
joint <- function(y, m, S) {
  L <- chol(S)
  z <- backsolve(L, y - m, transpose = TRUE)
  -0.5 * sum(z^2) - sum(log(diag(L))) - 0.5 * length(y) * log(2 * pi)
}
# brms's density for sigma: a half student_t(3, 0, sc) or a half
# normal(0, s), with brms's normalization (the half is a constant)
lst <- function(x, nu, mu, s) {
  lgamma((nu + 1) / 2) - lgamma(nu / 2) - 0.5 * log(nu * pi) - log(s) -
    (nu + 1) / 2 * log1p(((x - mu) / s)^2 / nu)
}
dp <- default_prior(fo, data = d, route = "sample")
print(dp[, c("prior", "class", "coef", "lb", "ub")])
ipr <- as.numeric(strsplit(gsub("student_t[(]|[)]", "",
                                dp$prior[dp$class == "Intercept"]), ",")[[1]])
run <- function(lab, prior, dens) {
  cat("\n==", st, lab, "\n")
  args <- list(fo, data = d, chains = 2, iter = 800, seed = 5, refresh = 0)
  if (!is.null(prior)) args$prior <- prior
  ds <- suppressWarnings(suppressMessages(do.call(frm_sample, args)))
  m <- as.matrix(ds$draws)
  cat("variables:", paste(colnames(m), collapse = " "), "\n")
  raw <- rstan::extract(ds$stanfit, permuted = FALSE)
  raw <- do.call(rbind, lapply(seq_len(dim(raw)[2]), function(ch) raw[, ch, ]))
  cat("raw columns:", paste(colnames(raw), collapse = " "), "\n")
  ar <- if ("ar[1]" %in% colnames(m)) m[, "ar[1]"] else rep(0, nrow(m))
  ma <- if ("ma[1]" %in% colnames(m)) m[, "ma[1]"] else rep(0, nrow(m))
  cat("ar in (-1, 1):", all(abs(ar) < 1), " ma in (-1, 1):", all(abs(ma) < 1),
      "\n")
  ll <- log_lik(ds)
  wl <- 0; res <- numeric(nrow(m))
  # the chart's raw coordinates: the last columns before lp__ are the
  # thetaac block, in the order ar then ma
  nac <- (st != "ma") + (st != "ar")
  pc <- ncol(raw) - 1L
  uac <- raw[, (pc - nac + 1):pc, drop = FALSE]
  nat <- cbind(if (st != "ma") ar, if (st != "ar") ma)
  chart <- isTRUE(all.equal(unname(nat), unname(-1 + 2 * stats::plogis(uac)),
                            tolerance = 1e-10))
  pacf <- isTRUE(all.equal(unname(nat), unname(uac / sqrt(1 + uac^2)),
                           tolerance = 1e-10))
  cat("raw coordinate of the coefficients: chart", chart, " internal",
      pacf, "\n")
  for (s in seq_len(nrow(m))) {
    mu <- m[s, "b_Intercept"] + m[s, "b_x"] * d$x
    sg <- m[s, "sigma"]
    ref <- numeric(nrow(d)); jl <- 0
    for (g in unique(d$g)) {
      r <- which(d$g == g)
      S <- sg^2 * corf(abs(outer(d$t[r], d$t[r], "-")), ar[s], ma[s])
      ref[r] <- cond(d$y[r], mu[r], S)
      jl <- jl + joint(d$y[r], mu[r], S)
    }
    wl <- max(wl, max(abs(ll[s, ] - ref)) / max(abs(ref)))
    si <- sg * cfac(ar[s], ma[s])
    # the raw coordinate is a chart v (natural = -1 + 2 plogis(v)) or,
    # without one, the internal th (natural = th / sqrt(1 + th^2)); which
    # one is read off the draws above
    p <- stats::plogis(uac[s, ])
    jac <- if (chart) sum(log(2) + log(p) + log1p(-p)) else
      sum(-1.5 * log1p(uac[s, ]^2))
    # brms's default on the Intercept, at the centered intercept
    ic <- m[s, "b_Intercept"] + m[s, "b_x"] * mean(d$x)
    res[s] <- raw[s, "lp__"] - jl - dens(si) - log(si) - jac -
      lst(ic, ipr[1], ipr[2], ipr[3])
  }
  cat("log_lik against the dense conditional (marginal sd): max relative",
      format(wl), "\n")
  cat("lp__ - (dense joint + prior at innovation sd + Jacobians): sd",
      format(stats::sd(res)), " range", format(diff(range(res))), "\n")
  # the same with the prior read at the MARGINAL sd, as 0.75.0 read it:
  # must NOT be constant
  res2 <- res + vapply(seq_len(nrow(m)), function(s) {
    si <- m[s, "sigma"] * cfac(ar[s], ma[s])
    dens(si) + log(si) - dens(m[s, "sigma"]) - log(m[s, "sigma"])
  }, 0)
  cat("control, prior read at the marginal sd instead: sd",
      format(stats::sd(res2)), "\n")
  lo <- suppressWarnings(loo(ds))
  cat("loo elpd", format(lo$estimates["elpd_loo", "Estimate"]),
      " pareto k > 0.7:", sum(lo$diagnostics$pareto_k > 0.7), "\n")
  cat("posterior means: sigma", format(mean(m[, "sigma"])), " innovation",
      format(mean(m[, "sigma"] * cfac(ar, ma))), "\n")
}
run("default priors", NULL, function(x) lst(x, 3, 0, 2.5))
run("user prior normal(0, 0.5) on class sigma",
    set_prior("normal(0, 0.5)", class = "sigma"),
    function(x) stats::dnorm(x, 0, 0.5, log = TRUE))
