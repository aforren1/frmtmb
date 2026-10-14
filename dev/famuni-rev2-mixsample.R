# Review round 2: frm_sample() on a two-GEV mixture (scalar xi1, xi2,
# no chart and no warning), and on `xi ~ 1` (log1p by design, as brms).
# Seed 35 construction of dev/famuni-mixprobe.R at n = 200; seed 1402.
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
run <- function(lab, fit) {
  ws <- character(0)
  dr <- withCallingHandlers(suppressMessages(frm_sample(fit, chains = 4,
          iter = 2000, seed = 1, refresh = 0,
          prior = c(set_prior("normal(0, 2.5)", class = "xi1"),
                    set_prior("normal(0, 2.5)", class = "xi2")))),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") },
    error = function(e) stop(e))
  cat(lab, ": divergences", divs(dr$stanfit), "of 4000; warnings:",
      if (length(ws)) unique(substr(ws, 1, 80)) else "none", "\n")
}
set.seed(35)
d <- data.frame(y = c(rg(100, 0, 1, -0.3), rg(100, 6, 1, -0.3)))
fm <- suppressWarnings(frm(y ~ 1, family = mixture(gen_extreme_value(),
                                                   gen_extreme_value()), data = d))
print(summary(fm)$spec_pars)
r <- tryCatch(run("mixture(gev, gev), xi -0.3", fm), error = function(e)
  cat("mixture run ERROR:", conditionMessage(e), "\n"))
