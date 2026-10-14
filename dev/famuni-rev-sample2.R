# Reviewer probe: GEV divergences under frm_sample() against brms on the
# same data (brms samples a scalar xi through scale_xi()), with a
# positive-xi control. Seeds 1401 (xi -0.45), 1403 (xi 0.2), n = 40.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
divs <- function(sf) {
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  sum(vapply(sp, function(s) sum(s[, "divergent__"]), 0))
}
for (cfg in list(c(-0.45, 1401), c(0.2, 1403), c(-0.85, 1402))) {
  set.seed(cfg[2])
  n <- 40
  d <- data.frame(y = rg(n, 0, 1, cfg[1]))
  ws <- character(0)
  fit <- withCallingHandlers(frm(y ~ 1, family = gen_extreme_value(), data = d),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleMessage") })
  sp <- summary(fit)$spec_pars
  cat(sprintf("\n=== xi true %.2f: ML conv %d logLik %.4f xi %.7f (se %.4g)\n",
              cfg[1], fit$opt$convergence, as.numeric(logLik(fit)),
              sp["xi", 1], sp["xi", 2]))
  if (length(ws)) cat("   ML says:", substr(ws, 1, 250), sep = "\n   ")
  ev <- suppressWarnings(evd::fgev(d$y, std.err = FALSE))
  cat(sprintf("   evd::fgev: ll %.4f xi %.5f\n", -ev$deviance / 2,
              ev$estimate[["shape"]]))
  dr <- suppressMessages(frm_sample(fit, chains = 4, iter = 2000, seed = 1,
          prior = set_prior("normal(0, 2.5)", class = "xi"), refresh = 0))
  cat("   frm_sample: divergences", divs(dr$stanfit), "of 4000\n")
  print(summary(dr))
  bf <- suppressWarnings(suppressMessages(brms::brm(y ~ 1, data = d,
          family = brms::gen_extreme_value(), chains = 4, iter = 2000,
          seed = 1, backend = "rstan", refresh = 0, silent = 2)))
  cat("   brms: divergences", divs(bf$fit), "of 4000\n")
  print(summary(bf)$spec_pars)
}
