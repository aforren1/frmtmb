# Reviewer of lane ordmix, re-check of B2: the margin of the
# degenerate-component warning (latent distance > 50) on SOUND fits:
# larger n, strong predictors (slopes scaled by `c`), and the cloglog and
# cauchit links with latent noise drawn from the link's own distribution.
# Two classes, class 1 slope 1.5 c with location +1, class 2 slope
# -0.8 c with location -1, as dev/ordmix-ident.R with c = 1. Seeds 1..R.
# Per fit: the largest reach over components, whether the degenerate
# warning fired, whether all standard errors are finite, and the
# relative error of the two slopes against the truth (best labeling).
# Usage: Rscript dev/ordmix-rev2-margin.R <link> <c> <n> [R]
args <- commandArgs(TRUE)
link <- args[1]
cc <- as.numeric(args[2])
n <- as.integer(args[3])
R <- if (length(args) >= 4) as.integer(args[4]) else 20L
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
noise <- switch(link,
  logit = function(n) rlogis(n),
  probit = function(n) rnorm(n),
  cloglog = function(n) log(-log(1 - runif(n))),
  cauchit = function(n) rcauchy(n))
truth <- c(1.5, -0.8) * cc
cut_at <- c(-1.5, 0, 1.5) * max(1, cc / 2)
for (s in seq_len(R)) {
  set.seed(s)
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, truth[1] * x + 1, truth[2] * x - 1) + noise(n)
  d <- data.frame(x = x, y = 1L + (lat > cut_at[1]) + (lat > cut_at[2]) +
                    (lat > cut_at[3]))
  w <- character(0)
  f <- tryCatch(withCallingHandlers(
    frm(bf(y ~ x), family = mixture(cumulative(link), cumulative(link)),
        data = d),
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
  if (is.character(f)) {
    cat(sprintf("REP link=%s c=%g n=%d seed=%d ERROR %s\n", link, cc, n, s,
                substr(f, 1, 150)))
    next
  }
  dg <- frmtmb:::mixture_ord_degeneracy(f, names(f$spec$responses)[1])
  fx <- suppressWarnings(fixef(f))
  b <- fx[c("mu1_x", "mu2_x"), "Estimate"]
  rel <- min(max(abs(b - truth) / abs(truth)),
             max(abs(rev(b) - truth) / abs(truth)))
  cat(sprintf("REP link=%s c=%g n=%d seed=%d reach=%.1f degwarn=%s se_finite=%s slope_relerr=%.3f maxabs_est=%.1f\n",
              link, cc, n, s, max(dg$reach),
              any(grepl("degenerate boundary", w)),
              all(is.finite(fx[, "Est.Error"])), rel,
              max(abs(fx[, "Estimate"]))))
}
