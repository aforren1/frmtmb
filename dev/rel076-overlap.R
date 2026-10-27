# The semantic overlaps of the round's three lanes, which no lane could
# test with the others' code: run on the merged build (rellib-r14).
#   1. cosy() + mo(): frm_sample() (acprior's cosy chart and default,
#      mocode's coding), then log_lik() (loglik's conditional) and loo().
#      The reference is the dense conditional of each group at each
#      draw, from the draws' own columns, with the mo() mean written by
#      hand from brms's definition (D * sum(simo[1:(k - min)])).
#   2. ar(cov = TRUE): frm_sample() with the default priors (acprior's
#      sigma read at the innovation sd), then log_lik() against the
#      dense AR(1) conditional with the marginal sd, and loo().
#   3. sar(): the prior tables (acprior's one slot and the lb/ub at the
#      interval's ends) passed back as prior = to frm() and to
#      frm_sample(), then log_lik() of those draws (loglik's Jacobian
#      slot) against the dense conditional of the lag SAR.
#   Rscript dev/rel076-overlap.R <lib> [arms]   (arms: 1 2 3, default all)
a <- commandArgs(TRUE)
lib <- a[1]
arms <- if (length(a) > 1) a[-1] else c("1", "2", "3")
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    as.character(packageVersion("frmtmb")), "sample",
    as.character(packageVersion("frmtmb.sample")), "\n")

cond <- function(y, m, S) {
  n <- length(y)
  vapply(seq_len(n), function(i) {
    o <- seq_len(n)[-i]
    A <- solve(S[o, o, drop = FALSE])
    cm <- m[i] + as.numeric(S[i, o, drop = FALSE] %*% A %*% (y[o] - m[o]))
    cv <- S[i, i] - as.numeric(S[i, o, drop = FALSE] %*% A %*% S[o, i])
    stats::dnorm(y[i], cm, sqrt(cv), log = TRUE)
  }, 0)
}
grouped <- function(d, m, s, corf) {
  out <- numeric(nrow(d))
  for (g in unique(d$g)) {
    r <- which(d$g == g)
    S <- s^2 * corf(abs(outer(d$t[r], d$t[r], "-")))
    out[r] <- cond(d$y[r], m[r], S)
  }
  out
}
rel <- function(a, b) max(abs(a - b)) / max(abs(b))

if ("1" %in% arms) {
  cat("\n== 1. cosy() + mo(): frm_sample(), log_lik(), loo()\n")
  set.seed(20261027)
  d <- expand.grid(t = 1:4, g = 1:30)
  d$x <- stats::rnorm(nrow(d))
  d$k <- sample(1:4, nrow(d), replace = TRUE)
  eff <- c(0, 0.6, 0.9, 1.5)[d$k]
  d$y <- 1 + 0.5 * d$x + eff + stats::rnorm(30, 0, 0.4)[d$g] +
    stats::rnorm(nrow(d), 0, 0.9)
  t0 <- proc.time()[["elapsed"]]
  ds <- suppressMessages(frm_sample(y ~ x + mo(k) + cosy(t, g), data = d,
                                    chains = 2, iter = 600, seed = 3,
                                    refresh = 0))
  cat("frm_sample seconds:", round(proc.time()[["elapsed"]] - t0, 1), "\n")
  cat("variables:", paste(variables(ds), collapse = " "), "\n")
  m <- as.matrix(ds$draws)
  sx <- grep("^simo_mok1\\[", colnames(m))
  cat("simo weights:", length(sx), "(values", min(d$k), "to", max(d$k),
      ")\n")
  cat("cosy draws: min", format(min(m[, "cosy"])), "max",
      format(max(m[, "cosy"])), "(brms's [0, 1])\n")
  ll <- log_lik(ds)
  cat("log_lik dim:", dim(ll), " finite:", all(is.finite(ll)), "\n")
  D <- length(sx)
  worst <- 0
  for (s in seq_len(nrow(m))) {
    cs <- c(0, cumsum(m[s, sx]))
    mu <- m[s, "b_Intercept"] + m[s, "b_x"] * d$x +
      m[s, "bsp_mok"] * D * cs[d$k - min(d$k) + 1]
    rho <- m[s, "cosy"]
    ref <- grouped(d, mu, m[s, "sigma"], function(L) ifelse(L == 0, 1, rho))
    worst <- max(worst, rel(ll[s, ], ref))
  }
  cat("log_lik against the dense conditional, every draw: max relative",
      format(worst), "\n")
  lo <- suppressWarnings(loo::loo(ll))
  cat("loo elpd", format(lo$estimates["elpd_loo", "Estimate"]),
      "pareto k > 0.7:", sum(lo$diagnostics$pareto_k > 0.7), "\n")
  lo2 <- suppressWarnings(loo(ds))
  cat("loo(ds) elpd", format(lo2$estimates["elpd_loo", "Estimate"]),
      " identical to loo::loo(log_lik(ds)):",
      identical(lo2$estimates, lo$estimates), "\n")
}

if ("2" %in% arms) {
  cat("\n== 2. ar(cov = TRUE): default priors, log_lik(), loo()\n")
  set.seed(20261027)
  d <- expand.grid(t = 1:6, g = 1:8)
  d$x <- stats::rnorm(nrow(d))
  e <- unlist(lapply(1:8, function(g) {
    as.numeric(stats::arima.sim(list(ar = 0.9), 6, sd = 0.45))
  }))
  d$y <- 1 + 0.5 * d$x + e
  dp <- default_prior(y ~ x + ar(t, g, cov = TRUE), data = d,
                      route = "sample")
  print(dp[dp$class %in% c("ar", "sigma"), c("prior", "class", "lb", "ub")])
  ds <- suppressMessages(frm_sample(y ~ x + ar(t, g, cov = TRUE), data = d,
                                    chains = 2, iter = 600, seed = 3,
                                    refresh = 0))
  m <- as.matrix(ds$draws)
  cat("posterior means: sigma (marginal sd)", format(mean(m[, "sigma"])),
      " ar[1]", format(mean(m[, "ar[1]"])), " innovation sd",
      format(mean(m[, "sigma"] * sqrt(1 - m[, "ar[1]"]^2))), "\n")
  ll <- log_lik(ds)
  worst <- 0
  for (s in seq_len(nrow(m))) {
    mu <- m[s, "b_Intercept"] + m[s, "b_x"] * d$x
    a1 <- m[s, "ar[1]"]
    ref <- grouped(d, mu, m[s, "sigma"], function(L) a1^L)
    worst <- max(worst, rel(ll[s, ], ref))
  }
  cat("log_lik dim:", dim(ll), " against the dense AR(1) conditional",
      "with the marginal sd: max relative", format(worst), "\n")
  lo <- suppressWarnings(loo(ds))
  cat("loo elpd", format(lo$estimates["elpd_loo", "Estimate"]), "\n")
  # the same prior typed by the user: the table passed back
  ds2 <- suppressMessages(frm_sample(y ~ x + ar(t, g, cov = TRUE),
                                     data = d, prior = dp, chains = 2,
                                     iter = 600, seed = 3, refresh = 0))
  cat("table passed back as prior =: draws identical to the default's:",
      identical(unname(as.matrix(ds2$draws)), unname(m)), "\n")
}

if ("3" %in% arms) {
  cat("\n== 3. sar(): the prior tables, frm(), frm_sample(), log_lik()\n")
  suppressMessages(library(frmtmb.spatial))
  cat("frmtmb.spatial", as.character(packageVersion("frmtmb.spatial")),
      "\n")
  k <- 6
  A <- matrix(0, k^2, k^2)
  for (i in 1:k) for (j in 1:k) {
    v <- (i - 1) * k + j
    if (i < k) A[v, v + k] <- A[v + k, v] <- 1
    if (j < k) A[v, v + 1] <- A[v + 1, v] <- 1
  }
  W <- A / rowSums(A)
  n <- nrow(W)
  set.seed(2)
  d <- data.frame(x = stats::rnorm(n))
  IA <- diag(n) - 0.5 * W
  d$y <- as.vector(solve(IA, 1 + 0.5 * d$x + stats::rnorm(n, 0, 0.7)))
  d2 <- list(W = W)
  for (route in c("fit", "sample")) {
    dp <- default_prior(y ~ x + sar(W), data = d, data2 = d2,
                        route = route)
    cat("route", route, ": classes", paste(dp$class, collapse = " "),
        "| theta rows:", sum(dp$class == "theta"), "| lagsar lb, ub:",
        dp$lb[dp$class == "lagsar"], dp$ub[dp$class == "lagsar"], "\n")
  }
  dpf <- default_prior(y ~ x + sar(W), data = d, data2 = d2)
  f0 <- frm(y ~ x + sar(W), data = d, data2 = d2)
  f1 <- frm(y ~ x + sar(W), data = d, data2 = d2, prior = dpf)
  cat("frm(prior = table) logLik equal:",
      isTRUE(all.equal(as.numeric(logLik(f0)), as.numeric(logLik(f1)))),
      "\n")
  dps <- default_prior(y ~ x + sar(W), data = d, data2 = d2,
                       route = "sample")
  ds0 <- suppressMessages(frm_sample(y ~ x + sar(W), data = d, data2 = d2,
                                     chains = 2, iter = 500, seed = 3,
                                     refresh = 0))
  ds <- suppressMessages(frm_sample(y ~ x + sar(W), data = d, data2 = d2,
                                    prior = dps, chains = 2, iter = 500,
                                    seed = 3, refresh = 0))
  cat("frm_sample(prior = sample table): draws identical to the",
      "default's:", identical(unname(as.matrix(ds$draws)),
                              unname(as.matrix(ds0$draws))), "\n")
  m <- as.matrix(ds$draws)
  cat("variables:", paste(variables(ds), collapse = " "), "\n")
  ll <- log_lik(ds)
  worst <- 0
  for (s in seq_len(nrow(m))) {
    IAs <- diag(n) - m[s, "lagsar"] * W
    mu <- solve(IAs, m[s, "b_Intercept"] + m[s, "b_x"] * d$x)
    S <- m[s, "sigma"]^2 * solve(crossprod(IAs))
    ref <- cond(d$y, as.vector(mu), S)
    worst <- max(worst, rel(ll[s, ], ref))
  }
  cat("log_lik dim:", dim(ll), " against the dense conditional of the",
      "lag SAR: max relative", format(worst), "\n")
  lo <- suppressWarnings(loo(ds))
  cat("loo elpd", format(lo$estimates["elpd_loo", "Estimate"]), "\n")
}
cat("OVERLAP DONE\n")
