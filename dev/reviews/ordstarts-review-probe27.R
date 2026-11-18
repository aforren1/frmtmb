# Review of lane ordstarts: what is the point 1.52 above the default fit
# that 3 of 8 jittered refits (sd 0.5, set.seed(9100 + k)) of cum_cum
# int n = 300 seed 27 reach from A, above B (starts = 10 and 30)?
# Same jitter as dev/reviews/ordstarts-review-rerun.R.
#   Rscript dev/reviews/ordstarts-review-probe27.R <family> <n> <G> <seed>
args <- commandArgs(trailingOnly = TRUE)
famname <- args[1]
nobs <- as.integer(args[2])
G <- as.integer(args[3])
seed <- as.integer(args[4])
re <- "int"
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
src <- readLines("dev/optima-csmix.R")
eval(parse(text = src[grep("^mk <- ", src):(grep("^one <- ", src) - 1L)]))
lane <- readLines("dev/ordstarts-run.R")
eval(parse(text = lane[grep("^mkg <- ", lane):(grep("^form <- ", lane) -
                                                  1L)]))
fam <- switch(famname,
  cum_cum = function() mixture(cumulative(), cumulative()),
  cum_sratio = function() mixture(cumulative(), sratio()))
q <- function(e) suppressWarnings(suppressMessages(e))
d <- mkg(seed)
A <- q(frm(bf(y ~ x + z + (1 | g)), family = fam(), data = d))
A0 <- A
A0$control$starts <- 0L
llA <- as.numeric(logLik(A))
cat(famname, nobs, seed, "A", llA, "\n")
resp <- names(A$spec$responses)[1]
for (k in 1:8) {
  set.seed(9100 + k)
  st <- A$estimates
  for (nm in setdiff(names(st), "b")) {
    st[[nm]] <- st[[nm]] + stats::rnorm(length(st[[nm]]), 0, 0.5)
  }
  r <- tryCatch(q(refit(A0, d$y, start = st)), error = function(e) e)
  if (inherits(r, "error")) {
    cat("k", k, "error\n")
    next
  }
  ll <- as.numeric(logLik(r))
  cl <- ns$mixture_ord_start_class(r, resp)
  q(ns$se_fit_analysis(r, force = TRUE))
  tse <- q(ns$mixture_ord_thres_se(r, resp))
  dg <- q(ns$mixture_ord_degeneracy(r, resp))
  w <- character()
  withCallingHandlers(ns$mixture_ord_degenerate_check(r, resp),
                      warning = function(x) {
                        w <<- c(w, conditionMessage(x))
                        invokeRestart("muffleWarning")
                      })
  cat(sprintf("k %d ll - A %+.4f code %d cls %s tse %s sharpen %s near %s",
              k, ll - llA, r$opt$convergence,
              if (length(cl$step)) "step" else if (length(cl$boundary))
                "flat" else "proper",
              paste(signif(tse, 3), collapse = ","),
              paste(signif(dg$sharpen, 3), collapse = ","),
              paste(signif(dg$near, 3), collapse = ",")),
      "| max|par|", signif(max(abs(r$opt$par)), 3), "\n")
  if (length(w)) cat("   degenerate check:", substr(w, 1, 200), "\n")
}
