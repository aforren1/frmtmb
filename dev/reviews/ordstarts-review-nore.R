# Review of lane ordstarts: does the MEDIUM (the start ranking does not
# read the threshold standard-error test) reach the shipped default,
# the 10 starts of an ordinal mixture WITHOUT random effects? The lane
# says "code reading, not measured". Data mk(seed, n, mix = TRUE) of
# dev/optima-csmix.R, y ~ x + z, as dev/ordstart-starts.R. Per seed:
# starts = 0 (A0) and the default (D, 10 starts on a cold fit), with
# their warnings; the class and the forced threshold standard error.
#   Rscript dev/reviews/ordstarts-review-nore.R <family> <n> <seeds> <out>
args <- commandArgs(trailingOnly = TRUE)
famname <- args[1]
nobs <- as.integer(args[2])
seeds <- eval(parse(text = args[3]))
out <- args[4]
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
fam <- switch(famname,
  cum_cum = function() mixture(cumulative(), cumulative()),
  cum_sratio = function() mixture(cumulative(), sratio()))
q <- function(e) suppressWarnings(suppressMessages(e))
grab <- function(expr) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    },
    message = function(x) invokeRestart("muffleMessage")),
    error = function(e) e)
  list(v = v, w = w, el = proc.time()[["elapsed"]] - t0)
}
info <- function(r, tag) {
  f <- r$v
  if (inherits(f, "error")) {
    return(stats::setNames(data.frame(NA, NA, NA, NA, NA, NA, r$el),
                           paste0(tag, c("ll", "code", "deg", "nwarn",
                                         "cls", "tse", "el"))))
  }
  resp <- names(f$spec$responses)[1]
  cl <- ns$mixture_ord_start_class(f, resp)
  tse <- suppressWarnings(max(q(ns$mixture_ord_thres_se(f, resp)),
                              na.rm = TRUE))
  stats::setNames(data.frame(as.numeric(logLik(f)), f$opt$convergence,
    any(grepl("ended at a degenerate boundary", r$w)), length(r$w),
    if (length(cl$step)) "step" else if (length(cl$boundary)) "flat" else
      "proper", tse, r$el),
    paste0(tag, c("ll", "code", "deg", "nwarn", "cls", "tse", "el")))
}
first <- TRUE
for (s in seeds) {
  d <- mk(s, n = nobs, mix = TRUE)
  a <- grab(frm(bf(y ~ x + z), family = fam(), data = d,
                control = frmtmb_control(starts = 0)))
  b <- grab(frm(bf(y ~ x + z), family = fam(), data = d))
  row <- cbind(data.frame(family = famname, n = nobs, seed = s),
               info(a, "A0_"), info(b, "D_"))
  utils::write.table(row, out, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = first, append = !first)
  first <- FALSE
  cat(famname, nobs, s, "A0", row$A0_ll, row$A0_deg, "D", row$D_ll,
      row$D_deg, "\n")
}
cat("DONE\n")
