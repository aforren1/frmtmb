# Reviewer short check, item 2: the ranking reads the lost list on H /
# mean weight, the fit-end standard-error check on H. On the sound fit of
# test-ordinal-mixture.R (cum_cum seed 3, n 500) with weights 1, 1/500,
# 1/5000: the warnings at the fit and at summary(), and the lost flags
# both ways.
#   Rscript ordstarts-short-p11.R > ordstarts-short-p11.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
wdata <- function(seed, n, w) {
  set.seed(seed)
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  data.frame(x = x, y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5), w = w)
}
grabw <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
for (seed in c(3, 5, 8)) {
  for (wv in c(1, 1 / 500, 1 / 5000)) {
    for (st in c(0, 10)) {
      r <- grabw(frm(bf(y | weights(w) ~ x),
                     family = mixture(cumulative(), cumulative()),
                     data = wdata(seed, 500, wv),
                     control = frmtmb_control(starts = st)))
      s <- grabw(summary(r$v))
      f <- r$v
      f$cache <- new.env(parent = emptyenv())
      raw <- ns$mixture_ord_thres_se(f, "y", force = TRUE, lost_inf = TRUE)
      f$cache <- new.env(parent = emptyenv())
      unit <- ns$mixture_ord_thres_se(f, "y", force = TRUE, lost_inf = TRUE,
                                      unit_lost = TRUE)
      cl <- ns$mixture_ord_start_class(r$v, "y")
      ms <- r$v$opt$mixture_starts
      cat(sprintf(paste0("seed %d w %-7g starts %2d | ll/w %.5f | ran %s |",
                         " lost raw %s unit %s | class %s | fit warnings: %s |",
                         " summary warnings: %s\n"),
                  seed, wv, st, as.numeric(logLik(f)) / wv,
                  if (is.null(ms)) "-" else ms[["starts"]],
                  paste(attr(raw, "lost"), collapse = ","),
                  paste(attr(unit, "lost"), collapse = ","),
                  if (length(cl$step)) "step" else if
                  (length(c(cl$boundary, cl$unplaced))) "flat" else "proper",
                  paste(substr(r$w, 1, 90), collapse = " || "),
                  paste(substr(s$w, 1, 90), collapse = " || ")))
    }
  }
}
cat("DONE p11\n")
