# Reviewer of lane ordmix: how often a two-component ordinal mixture's
# default fit is degenerate (a component running off to a step function:
# any |estimate| > 30 or a non-finite standard error), whether the fit
# warns at fit time, and whether the degenerate point has the higher
# log likelihood than the best of 6 other starts (seed 500 + s).
# The lane's own data-generating process (dev/ordmix-ident.R, "none").
# Usage: Rscript dev/ordmix-rev-degen.R <cum_cum|cum_sr> <n> [R]
args <- commandArgs(TRUE)
model <- args[1]
n <- as.integer(args[2])
R <- if (length(args) >= 3) as.integer(args[3]) else 40L
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
fam <- switch(model, cum_cum = mixture(cumulative(), cumulative()),
              cum_sr = mixture(cumulative(), sratio()))
cat("model", model, "n", n, "R", R, find.package("frmtmb"), "\n")
degen_of <- function(f) {
  fx <- suppressWarnings(fixef(f))
  any(abs(fx[, "Estimate"]) > 30) || any(!is.finite(fx[, "Est.Error"]))
}
for (s in seq_len(R)) {
  set.seed(s)
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  d <- data.frame(x = x, y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5))
  w <- character(0)
  f <- withCallingHandlers(frm(bf(y ~ x), family = fam, data = d),
                           warning = function(e) {
                             w <<- c(w, conditionMessage(e))
                             invokeRestart("muffleWarning")
                           })
  dg <- degen_of(f)
  ll <- as.numeric(logLik(f))
  set.seed(500 + s)
  best_nd <- -Inf
  best_any <- ll
  for (j in 1:6) {
    st <- lapply(f$frame$par_template, function(v) v + rnorm(length(v), 0, 1))
    g <- tryCatch(suppressWarnings(frm(bf(y ~ x), family = fam, data = d,
                                       start = st)),
                  error = function(e) NULL)
    if (is.null(g)) next
    lg <- as.numeric(logLik(g))
    best_any <- max(best_any, lg)
    if (!degen_of(g)) best_nd <- max(best_nd, lg)
  }
  cat(sprintf("REP model=%s n=%d seed=%d degenerate=%s ll=%.4f best_nondegen=%.4f best_any=%.4f fit_warnings=%d\n",
              model, n, s, dg, ll, best_nd, best_any, length(w)))
}
