# Reviewer of lane ordmix, final check: the review's degenerate fits
# (dev/ordmix-rev-degen.R, cum + cum, n = 300, seeds 19, 33, 35, 38, 40)
# refitted with every weight 1, 10 and 1/n. The estimates are the same
# for constant weights; does the degenerate warning's verdict change?
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
n <- 300
for (s in c(19, 33, 35, 38, 40)) {
  set.seed(s)
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  d <- data.frame(x = x, y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5))
  out <- character(0)
  for (wv in c(1, 10, 1 / n)) {
    d$w <- wv
    w <- character(0)
    f <- withCallingHandlers(
      frm(bf(y | weights(w) ~ x), family = mixture(cumulative(),
                                                   cumulative()), data = d),
      warning = function(e) {
        w <<- c(w, conditionMessage(e))
        invokeRestart("muffleWarning")
      })
    dg <- frmtmb:::mixture_ord_degeneracy(f, "y")
    out <- c(out, sprintf("w=%.4g warn=%s sharpen=%s ll=%.4f", wv,
                          any(grepl("degenerate boundary", w)),
                          paste(signif(dg$sharpen, 3), collapse = ","),
                          as.numeric(logLik(f)) / wv))
  }
  cat("seed", s, ":", paste(out, collapse = " | "), "\n")
}
