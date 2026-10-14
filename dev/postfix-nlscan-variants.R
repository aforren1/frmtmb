# Lane postfix: variants of the nonlinear-intercept scan, swapped into
# the installed lane build, on the models of dev/postfix-nlstart.R.
#   Rscript dev/postfix-nlscan-variants.R <lib> <variant>
# variant "common": every intercept at one common power-of-ten value
# (the first version of the rule); "coord10": the common value, then two
# coordinate passes over each intercept on the powers of ten with zero;
# "built": the installed rule (common value, then two coordinate passes
# on the 1-2-5 grid with zero).
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
variant <- a[2]
ns <- asNamespace("frmtmb")
make_scan <- function(coord_grid) {
  function(obj, frame, start, prior_entries, placed) {
    if (!"beta" %in% names(obj$par)) return(NULL)
    idx <- integer(0)
    for (lp in frame[["linpreds"]]) {
      if (!identical(lp[["par"]], "beta")) next
      if (!is.null(lp[["nl_body"]]) || !is.null(lp[["constant"]])) next
      nl <- frame[["spec"]]$responses[[lp[["resp"]]]]$nlpars
      if (!lp[["dpar"]] %in% nl) next
      icpt <- match("(Intercept)", colnames(lp[["X"]]))
      if (!is.na(icpt)) idx <- c(idx, as.integer(lp[["idx"]][icpt]))
    }
    if (!length(idx)) return(NULL)
    f <- function(p) tryCatch(as.numeric(obj$fn(p)), error = function(e) NA)
    tpl <- frmtmb:::make_start(frame, start, prior_entries)
    best <- NULL
    bv <- Inf
    for (cv in as.vector(t(c(1, -1) %o% 10^(-2:3)))) {
      tpl[["beta"]][idx] <- cv
      p <- frmtmb:::outer_from_template(tpl, obj, frame)
      v <- f(p)
      if (is.finite(v) && v < bv) {
        best <- p
        bv <- v
      }
    }
    if (is.null(best) || is.null(coord_grid)) return(best)
    pos <- which(names(obj$par) == "beta")[idx]
    for (pass in 1:2) {
      for (k in pos) {
        for (cv in coord_grid) {
          p <- best
          p[k] <- cv
          v <- f(p)
          if (is.finite(v) && v < bv) {
            best <- p
            bv <- v
          }
        }
      }
    }
    best
  }
}
scan <- switch(variant,
  common = make_scan(NULL),
  coord10 = make_scan(c(0, as.vector(t(c(1, -1) %o% 10^(-2:3))))),
  built = NULL)
if (!is.null(scan)) {
  # the closure keeps its grid; the internals it calls are named
  unlockBinding("nl_intercept_scan", ns)
  assign("nl_intercept_scan", scan, envir = ns)
}
cat("variant", variant, "\n")
sys.source("dev/postfix-nlstart.R", envir = new.env())
