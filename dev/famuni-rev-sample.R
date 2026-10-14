# Reviewer probe: does frm_sample() of a GEV with xi near -0.5 put draws
# where an observation lies outside the support (the soft penalty is a
# finite density there)? Seed 1401, n = 40.
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb from", find.package("frmtmb"), "; sample from",
    find.package("frmtmb.sample"), "\n")
rg <- function(n, m, s, k) m + s * (stats::rexp(n)^(-k) - 1) / k
for (cfg in list(c(-0.45, 1401), c(-0.85, 1402))) {
  set.seed(cfg[2])
  n <- 40
  d <- data.frame(y = rg(n, 0, 1, cfg[1]))
  fit <- frm(y ~ 1, family = gen_extreme_value(), data = d)
  cat("\n=== xi true", cfg[1], ": ML xi", summary(fit)$spec_pars["xi", 1],
      " max y", max(d$y), "\n")
  dr <- suppressMessages(frm_sample(fit, chains = 4, iter = 2000, seed = 1,
                   prior = set_prior("normal(0, 2.5)", class = "xi")))
  ll <- log_lik(dr)
  rowmin <- apply(ll, 1, min)
  cat("draws:", nrow(ll), "; draws whose smallest row log density is below -50:",
      sum(rowmin < -50), "; non-finite:", sum(!is.finite(rowmin)), "\n")
  cat("smallest row log density over all draws:", min(rowmin), "\n")
  sp <- tryCatch(rstan::get_sampler_params(dr$stanfit, inc_warmup = FALSE),
                 error = function(e) NULL)
  if (!is.null(sp)) {
    cat("divergences:", sum(vapply(sp, function(s) sum(s[, "divergent__"]), 0)),
        "\n")
  } else cat("sampler params not found; names(dr):", names(dr), "\n")
}
