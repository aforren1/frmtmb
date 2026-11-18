# rel084-kt copy of dev/nlband-common.R: arm after = rellib-r23 then user lib (lane lib and rellib-r22 removed); FRMTMB_STAN_CACHE = the rel084 tree's dev/stan-cache. No other change.
# Common setup of lane nlband's scripts (round of 2026-11-16).
# NLBAND_ARM = before reads rellib-r22 (frmtmb 0.83.0 at c6383527); after
# (the default) puts the lane's library first.
arm <- Sys.getenv("NLBAND_ARM", "after")
libs <- c(if (arm == "after") "C:/Users/adf44/source/r/rellib-r23" else "C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-rel084/dev/stan-cache")
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "arm:", arm, "\n")
# run expr; keep its value, and the classes and texts of its conditions
tt <- function(expr) {
  w <- character(0); wc <- character(0); m <- character(0)
  mc <- character(0)
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); wc <<- c(wc, class(c)[1])
      invokeRestart("muffleWarning")
    },
    message = function(c) {
      m <<- c(m, conditionMessage(c)); mc <<- c(mc, class(c)[1])
      invokeRestart("muffleMessage")
    }),
    error = function(e) structure(list(msg = conditionMessage(e)),
                                  class = "tterr"))
  if (is.null(v)) v <- list()
  attr(v, "w") <- w; attr(v, "wc") <- wc
  attr(v, "m") <- m; attr(v, "mc") <- mc
  v
}
iserr <- function(x) inherits(x, "tterr")
llv <- function(r) {
  if (iserr(r) || is.null(r$opt)) return(NA_real_)
  v <- suppressWarnings(as.numeric(logLik(r)))
  if (is.finite(v)) v else NA_real_
}
wsum <- function(r) {
  wc <- attr(r, "wc")
  if (length(wc)) paste(unique(wc), collapse = "+") else "-"
}
msum <- function(r) {
  mc <- attr(r, "mc")
  if (length(mc)) paste(unique(mc), collapse = "+") else "-"
}
# the review's band data (dev/reviews/nlstop2-review-pole.R section 8):
# y = 1 / (x - 2) + N(0, 0.2), x ~ U(0, 4) without |x - 2| < hole
band_data <- function(seed, hole = 0.3, n = 400, pole = 2, sd = 0.2) {
  set.seed(seed)
  x <- runif(n, 0, 4)
  x <- x[abs(x - pole) > hole]
  d <- data.frame(x = x)
  d$y <- 1 / (d$x - pole) + rnorm(nrow(d), 0, sd)
  d
}
band_form <- bf(y ~ a / (x + b), a + b ~ 1, nl = TRUE)
