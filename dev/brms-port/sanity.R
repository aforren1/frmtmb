# Loose plausibility check: point estimates from the CLEAN fits next to
# the posterior means the vignettes print. Not an agreement suite - the
# question is only whether signs and magnitudes line up.
HERE <- "C:/Users/adf44/source/r/frmtmb-wt-audit/dev/brms-port"
VIGS <- c("brms_overview", "brms_multilevel", "brms_distreg", "brms_nonlinear",
          "brms_phylogenetics", "brms_monotonic", "brms_multivariate",
          "brms_missings", "brms_customfamilies")
res <- list()
for (v in VIGS) {
  f <- file.path(HERE, "results-spell", paste0(v, ".rds"))
  if (file.exists(f)) res <- c(res, readRDS(f))
}
for (r in res) {
  if (!identical(r$kind, "model") || !identical(r$status, "OK")) next
  if (is.null(r$fit) || is.null(r$fit$fixef)) next
  cat("\n==", r$id, "\n   ", gsub("\n *", " ", substr(r$src_run %||% r$src, 1, 110)), "\n")
  fx <- r$fit$fixef; if (is.numeric(fx)) print(round(fx, 3)) else print(fx)
  cat("    logLik ", round(r$fit$loglik, 2),
      "  sigma ", round(r$fit$sigma, 3), "\n")
}
`%||%` <- function(a, b) if (is.null(a)) b else a
