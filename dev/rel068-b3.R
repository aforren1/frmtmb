# Copied unchanged from dev/ordmix-rev3-b3.R at the 0.68.0 release, except that it loads
# the release library rellib-r6 (dev/round-20261005.md).
# Reviewer of lane ordmix, final check of B3: is the step-function test
# "doubling the latent distances costs < 0.1" scale-free? Designs that
# move the log-likelihood's scale or its sensitivity: weights that sum
# to 1 or to n/10, small n, many categories, a binary predictor with a
# strong effect under the logit and the cauchit. The lane's two-class
# process (slopes 1.5 c and -0.8 c, locations +1 and -1), seeds 1..R.
# "Sound": every standard error finite and both slopes within 50% of
# the truth (best labeling).
# Usage: Rscript dev/ordmix-rev3-b3.R <design> [R]
args <- commandArgs(TRUE)
design <- args[1]
R <- if (length(args) >= 2) as.integer(args[2]) else 20L
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cfg <- switch(design,
  w_sum1 = list(n = 500, c = 1, link = "logit", w = "sum1"),
  w_tenth = list(n = 500, c = 1, link = "logit", w = "tenth"),
  n60 = list(n = 60, c = 1, link = "logit"),
  n100 = list(n = 100, c = 1, link = "logit"),
  n100_c5 = list(n = 100, c = 5, link = "logit"),
  n60_c5 = list(n = 60, c = 5, link = "logit"),
  cat8 = list(n = 200, c = 1, link = "logit", ncut = 7),
  bin_logit = list(n = 500, c = 10, link = "logit", binx = TRUE),
  bin_cauchit = list(n = 500, c = 10, link = "cauchit", binx = TRUE),
  cauchit_c10_n2000 = list(n = 2000, c = 10, link = "cauchit"))
noise <- switch(cfg$link, logit = rlogis, cauchit = rcauchy)
truth <- c(1.5, -0.8) * cfg$c
cuts <- if (is.null(cfg$ncut)) c(-1.5, 0, 1.5) * max(1, cfg$c / 2) else
  seq(-3, 3, length.out = cfg$ncut)
for (s in seq_len(R)) {
  set.seed(s)
  n <- cfg$n
  x <- if (isTRUE(cfg$binx)) rbinom(n, 1, 0.5) else rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, truth[1] * x + 1, truth[2] * x - 1) + noise(n)
  d <- data.frame(x = x, y = 1L + rowSums(outer(lat, cuts, ">")))
  d$w <- switch(cfg$w %||% "none", sum1 = rep(1 / n, n),
                tenth = rep(0.1, n), none = rep(1, n))
  f0 <- if (is.null(cfg$w)) bf(y ~ x) else bf(y | weights(w) ~ x)
  ww <- character(0)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(
    frm(f0, family = mixture(cumulative(cfg$link), cumulative(cfg$link)),
        data = d),
    warning = function(e) {
      ww <<- c(ww, conditionMessage(e))
      invokeRestart("muffleWarning")
    }), error = function(e) conditionMessage(e))
  if (is.character(f)) {
    cat(sprintf("REP design=%s seed=%d ERROR %s\n", design, s,
                substr(f, 1, 120)))
    next
  }
  dg <- frmtmb:::mixture_ord_degeneracy(f, names(f$spec$responses)[1])
  fx <- suppressWarnings(fixef(f))
  b <- fx[c("mu1_x", "mu2_x"), "Estimate"]
  rel <- min(max(abs(b - truth) / abs(truth)),
             max(abs(rev(b) - truth) / abs(truth)))
  sound <- all(is.finite(fx[, "Est.Error"])) && rel < 0.5
  cat(sprintf("REP design=%s seed=%d sound=%s degwarn=%s sharpen=%s near=%.1f maxse=%.3g relerr=%.2f\n",
              design, s, sound, any(grepl("degenerate boundary", ww)),
              paste(signif(dg$sharpen, 3), collapse = ","), max(dg$near),
              max(fx[, "Est.Error"]), rel))
}
