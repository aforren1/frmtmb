# Release review of 0.81.0: does frm_sample() of an ordinary GEV fit with
# a predicted xi (xi ~ x, xi ~ 1) give draws that put a fitted row outside
# its support, so that log_lik() of the draws has -Inf cells and loo()
# stops? The record says no such draw was constructed.
#
#   Rscript dev/review-rel081-gevxi.R <lib> <case> <seed>
#   case: xix (xi ~ x), xi1 (xi ~ 1), sc (scalar xi, the tmp_xi chart)
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1]
case <- args[2]
seed <- as.integer(args[3])
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", dirname(find.package("frmtmb")), "frmtmb",
    format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), "case", case, "seed", seed, "\n")
# GEV draws with an upper end (xi < 0): the largest rows sit near it
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
set.seed(seed)
n <- if (length(args) >= 5) as.integer(args[5]) else 30
xi0 <- if (length(args) >= 4) as.numeric(args[4]) else -0.6
d <- data.frame(x = stats::rnorm(n))
xi_true <- if (case == "xix") xi0 + 0.2 * d$x else rep(xi0, n)
d$y <- rg(n, 1, 1, xi_true)
# cens: the lowest fifth left-censored at their 20% quantile, so a scalar
# xi is blocked from the tmp_xi chart and sampled on the fitted objective
if (case == "cens") {
  qc <- stats::quantile(d$y, 0.2)
  d$cc <- ifelse(d$y < qc, "left", "none")
  d$y[d$y < qc] <- qc
}
f <- switch(case,
  xix = bf(y ~ 1, xi ~ x),
  xi1 = bf(y ~ 1, xi ~ 1),
  sc = bf(y ~ 1),
  cens = bf(y | cens(cc) ~ 1))
fit <- withCallingHandlers(frm(f, family = gen_extreme_value(), data = d),
  warning = function(w) {
    cat("fit warning:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
cat("fit logLik:", format(as.numeric(logLik(fit)), digits = 10), "\n")
t0 <- proc.time()[["elapsed"]]
ds <- withCallingHandlers(
  frm_sample(fit, chains = 4, iter = 1000, seed = seed, refresh = 0,
             cores = 1),
  warning = function(w) {
    cat("sample warning:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
cat("sample seconds:", round(proc.time()[["elapsed"]] - t0, 1), "\n")
ll <- log_lik(ds)
cat("log_lik cells:", length(ll), " -Inf:", sum(ll == -Inf),
    " +Inf:", sum(ll == Inf), " NaN:", sum(is.nan(ll)),
    " draws with a -Inf cell:", sum(apply(ll == -Inf, 1, any)), "of",
    nrow(ll), "\n")
ll1 <- log_lik(ds, newdata = ds$data)
cat("identical(log_lik(ds), log_lik(ds, newdata = ds$data)):",
    identical(ll, ll1), "\n")
if (any(ll == -Inf)) {
  ij <- which(ll == -Inf, arr.ind = TRUE)
  cat("rows with a -Inf cell:", paste(sort(unique(ij[, 2])), collapse = " "),
      "\n")
}
for (fn in c("loo", "psis", "waic")) {
  r <- tryCatch(suppressWarnings(switch(fn,
         loo = loo(ds), psis = psis(ds), waic = waic(ds))),
       error = function(e) e)
  if (inherits(r, "error")) {
    cat(fn, "(draws): ERROR:", conditionMessage(r), "\n")
  } else if (fn == "psis") {
    cat(fn, "(draws): ok\n")
  } else {
    est <- r$estimates
    cat(fn, "(draws): ok; elpd", format(est[1, "Estimate"], digits = 8), "\n")
  }
}
cat("GEVXI DONE\n")
