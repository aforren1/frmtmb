# Review of lane gpgap, item 4: brms's absolute nugget and y's units. Does
# brms itself fit gp(x) when y is in larger units? The non-centered option
# "at brms's 1e-12" inherits whatever this shows. Seed 7701.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-brmsunits.R <log> <A>
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
A <- as.numeric(commandArgs(trailingOnly = TRUE)[1])
set.seed(7701)
n <- 60
x <- round(runif(n, 0, 6), 3)
d <- data.frame(x = x, y = A * (sin(2 * x) + 0.1 * rnorm(n)))
fr <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
bk <- Filter(function(k) identical(k$covstruct, "gp"), fr$frame$re_blocks)[[1]]
th <- fr$estimates$theta[bk$theta_idx]
cat(sprintf("frm: sigma %.4g sdgp %.4g rho %.4g\n", as.numeric(sigma(fr)),
            exp(th[[1]]), exp(th[[2]])))
u <- sort(unique(x))
K <- exp(2 * th[[1]]) * exp(-outer(u, u, "-")^2 / (2 * exp(2 * th[[2]])))
cat("kappa of the kernel without nugget:", signif(kappa(K, exact = TRUE), 3), "\n")
for (nug in c(1e-12, 1e-6 * exp(2 * th[[1]]))) {
  KK <- K; diag(KK) <- diag(KK) + nug
  cat(sprintf("chol at frm's estimates, nugget %.3g: %s\n", nug,
              if (inherits(try(chol(KK), silent = TRUE), "try-error")) "FAILS" else "ok"))
}
t0 <- proc.time()[["elapsed"]]
msgs <- character()
bm <- withCallingHandlers(
  tryCatch(brms::brm(brms::bf(y ~ gp(x)), data = d, family = gaussian(),
                     chains = 2, iter = 1000, warmup = 500, seed = 77, refresh = 0,
                     cores = 2, silent = 2),
           error = function(e) {
             cat("brm ERROR:", conditionMessage(e), "\n")
             NULL
           }),
  warning = function(w) {
    msgs <<- c(msgs, substr(conditionMessage(w), 1, 120))
    invokeRestart("muffleWarning")
  })
cat(sprintf("brm wall %.0f s\n", proc.time()[["elapsed"]] - t0))
print(unique(msgs))
if (!is.null(bm)) {
  print(summary(bm)$gp)
  print(summary(bm)$spec_pars)
  ds <- posterior::as_draws_df(bm)
  cat("draws:", nrow(ds), "sdgp median", signif(median(ds$sdgp_gpx), 4), "\n")
}
