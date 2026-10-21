# brms 2.23.0's loo_epred(ndraws =): which draws weight which. Traces
# brms:::validate_draw_ids() to record the subsample each half draws.
#   Rscript dev/loodraws-brms-ndraws.R  -> dev/loodraws-log/brms-ndraws.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(brms))
set.seed(1)
n <- 60
dd <- data.frame(x = rnorm(n))
dd$y <- rnorm(n, 1 + 0.5 * dd$x)
fit <- suppressWarnings(suppressMessages(brm(y ~ x, data = dd, chains = 4,
  iter = 2000, seed = 1, refresh = 0)))
ids <- list()
trace(brms:::validate_draw_ids, print = FALSE,
      exit = quote(ids[[length(ids) + 1L]] <<- returnValue()))
set.seed(42)
got <- suppressWarnings(suppressMessages(loo_epred(fit, ndraws = 100)))
untrace(brms:::validate_draw_ids)
cat("validate_draw_ids calls:", length(ids), "\n")
cat("subsamples identical:", identical(ids[[1]], ids[[2]]), "\n")
cat("draws shared by the two subsamples:", length(intersect(ids[[1]], ids[[2]])),
    "of", length(ids[[1]]), "\n")
ll <- log_lik(fit)[ids[[1]], ]
ps <- suppressWarnings(loo::psis(-ll, r_eff = loo::relative_eff(
  exp(ll), chain_id = rep(1L, 100))))
ep <- posterior_epred(fit)
other <- loo::E_loo(ep[ids[[2]], ], ps)$value
same <- loo::E_loo(ep[ids[[1]], ], ps)$value
cat("brms result == weights of subsample 1 on prediction of subsample 2:",
    identical(unname(got[, 1]), other), "\n")
cat("max |brms - same-draw answer| / max |same-draw answer|:",
    signif(max(abs(got[, 1] - same)) / max(abs(same)), 3), "\n")
# which pairing reproduces brms's number exactly
for (a in 1:2) for (b in 1:2) for (re in c("chain1", "none")) {
  lla <- log_lik(fit, draw_ids = ids[[a]])
  r_eff <- if (re == "chain1") loo::relative_eff(exp(lla),
                                                 chain_id = rep(1L, 100)) else 1
  psa <- suppressWarnings(loo::psis(-lla, r_eff = r_eff))
  v <- loo::E_loo(posterior_epred(fit, draw_ids = ids[[b]]), psa)$value
  cat(sprintf("weights from subsample %d, prediction %d, r_eff %s: max abs diff %s\n",
              a, b, re, signif(max(abs(got[, 1] - v)), 3)))
}
