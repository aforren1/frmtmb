# Release reviewer, 0.73.0: the loo_*() generics now owned by core and
# re-exported by frmtmb.sample, under each attach order of brms, loo and
# rstantools, with frmtmb.sample's weights refusal (core's
# check_weights_unbounded()) in the same session.
#   Rscript dev/rel073-rev-owners.R <lib> <order>
# order: none | brms_first | sample_first | loo_first | rstantools_first
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ord <- a[2]
q <- function(e) suppressWarnings(suppressMessages(e))
att <- function(p) q(library(p, character.only = TRUE))
switch(ord,
       none = att("frmtmb.sample"),
       brms_first = { att("brms"); att("frmtmb.sample") },
       sample_first = { att("frmtmb.sample"); att("brms") },
       loo_first = { att("loo"); att("frmtmb.sample") },
       rstantools_first = { att("rstantools"); att("frmtmb.sample") })
cat("order", ord, "| search:", paste(grep("^package:(brms|loo|rstantools|frmtmb)",
                                          search(), value = TRUE),
                                     collapse = " "), "\n")
set.seed(1)
n <- 60
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)), w = 1.5)
d$y <- 1 + 0.5 * d$x + rnorm(n)
fit <- q(frmtmb::frm(y ~ x, data = d))
mk <- function(fit, nd = 300L, sd = 0.03, seed = 5L) {
  est <- unlist(fit$estimates[names(fit$frame$par_template)])
  set.seed(seed)
  M <- matrix(rep(est, each = nd) + rnorm(nd * length(est), 0, sd), nd,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}
dr <- mk(fit)
for (g in c("loo_epred", "loo_linpred", "loo_predict",
            "loo_predictive_interval", "loo_R2", "loo_model_weights",
            "kfold", "reloo", "loo_moment_match", "loo_subsample")) {
  f <- get(g)  # what the user reaches from the search path
  owner <- environmentName(environment(f))
  on_draws <- tryCatch({
    r <- q(if (g == "loo_model_weights") f(dr, dr) else f(dr))
    paste("ok", paste(dim(as.array(r)), collapse = "x"))
  }, error = function(e) paste("ERROR", substr(conditionMessage(e), 1, 90)))
  on_fit <- tryCatch({q(f(fit)); "ran (should refuse)"},
                     error = function(e) substr(conditionMessage(e), 1, 70))
  cat(sprintf("%-24s reached %-12s | draws: %-40s | fit: %s\n", g, owner,
              on_draws, gsub("\\s+", " ", on_fit)))
}
r <- tryCatch(q(frm_sample(frmtmb::bf(y | weights(w) ~ x + (1 | obs)),
                           data = d, chains = 1, iter = 20, refresh = 0)),
              error = function(e) conditionMessage(e))
cat("frm_sample weights spike:", if (is.character(r)) substr(r, 1, 80) else
  "SAMPLED (should refuse)", "\n")
cat("END\n")
