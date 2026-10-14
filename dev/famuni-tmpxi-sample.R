# frm_sample() on a scalar-xi GEV after the tmp_xi chart: the review's
# data (seeds 1401, 1403, 1402; n = 40; y ~ 1), 4 chains x 2000, seed 1,
# with brms's defaults, divergences counted. Also a check that every
# draw's xi keeps every row inside the support at that draw's mu and
# sigma, and the default prior table.
# Output: dev/famuni-out/tmpxi-sample.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
stopifnot(pkgbuild::has_build_tools())
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("frmtmb from", find.package("frmtmb"), "; sample from",
    find.package("frmtmb.sample"), "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
divs <- function(sf) {
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  sum(vapply(sp, function(s) sum(s[, "divergent__"]), 0))
}
for (cfg in list(c(-0.45, 1401), c(0.2, 1403), c(-0.85, 1402))) {
  set.seed(cfg[2])
  d <- data.frame(y = rg(40, 0, 1, cfg[1]))
  fit <- suppressWarnings(frm(y ~ 1, family = gen_extreme_value(), data = d))
  cat(sprintf("\n=== xi true %.2f, seed %d\n", cfg[1], cfg[2]))
  gp <- get_prior(y ~ 1, data = d, family = gen_extreme_value(),
                  route = "sample")
  print(as.data.frame(gp)[gp$prior != "(flat)", c("prior", "class", "dpar")],
        row.names = FALSE)
  msgs <- character(0)
  t0 <- proc.time()[["elapsed"]]
  dr <- withCallingHandlers(
    frm_sample(fit, chains = 4, iter = 2000, seed = 1, refresh = 0),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      cat("  WARN:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  cat("  frm_sample messages:", paste(msgs, collapse = " | "), "\n")
  cat(sprintf("  divergences %d of 4000 (%.0f s)\n", divs(dr$stanfit),
              proc.time()[["elapsed"]] - t0))
  print(summary(dr))
  dm <- as.data.frame(dr)
  mu <- dm[["b_Intercept"]]
  sg <- dm[["sigma"]]
  xi <- dm[["xi"]]
  tmin <- vapply(seq_along(xi), function(i) {
    min(1 + xi[i] * (d$y - mu[i]) / sg[i])
  }, 0)
  cat(sprintf("  draws with a row outside the support: %d of %d; min t %.3g\n",
              sum(tmin <= 0), length(tmin), min(tmin)))
  cat(sprintf("  xi median %.3f (2.5%% %.3f, 97.5%% %.3f)\n", median(xi),
              quantile(xi, 0.025), quantile(xi, 0.975)))
}
