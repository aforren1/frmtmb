# Lane wcheck, item 2972: the fuzz specs with weights(w) and ar1() that
# reached objectives near -1e15 on rellib-r8. Since 0.73.0 the fuzz
# design halves w under ar1() (sum(w) < n always), so the tier no longer
# draws the runaway; this restores the pre-0.73.0 draw (w * 2) and fits
# every spec of the 300-row plan with re = ar1 and aterm = weights, then
# the same with the halved weights.
#   Rscript dev/wcheck-2972-fuzz.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages({library(frmtmb); library(testthat)})
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
source("tests/testthat/helper-fuzz.R")
plan <- fuzz_plan(seed = 20260901L, size = 300L)
sel <- plan[plan$re == "ar1" & plan$aterm == "weights", , drop = FALSE]
cat("plan rows:", nrow(plan), "| ar1 + weights:", nrow(sel), "\n")
one <- function(sp, d, tag) {
  r <- fuzz_fit_one(sp, d)
  ok <- is.null(r$error) && !is.null(r$value)
  if (!ok) {
    m <- gsub("[[:space:]]+", " ", conditionMessage(r$error %||% simpleError("?")))
    t <- if (grepl("has no maximum", m)) "REFUSED" else "ERROR"
    cat(sprintf("  %-7s %s | %s\n", tag, t, substr(m, 1, 110)))
  } else {
    f <- r$value
    cat(sprintf("  %-7s ok logLik %.6g | betad %s | conv %s | warnings %d\n",
                tag, as.numeric(logLik(f)),
                paste(sprintf("%.4g", f$estimates$betad), collapse = " "),
                f$opt$convergence, length(r$warnings)))
  }
}
`%||%` <- function(a, b) if (is.null(a)) b else a
for (i in seq_len(nrow(sel))) {
  sp <- as.list(sel[i, ])
  cat(sprintf("SPEC %s seed=%s\n", fuzz_spec_key(sp), sp$seed))
  d <- fuzz_data(sp)
  if (is.null(d)) { cat("  no data\n"); next }
  cat(sprintf("  sum(w) - n: halved %.3f, pre-0.73.0 %.3f\n",
              sum(d$w) - nrow(d), sum(2 * d$w) - nrow(d)))
  one(sp, d, "halved")
  d2 <- d; d2$w <- 2 * d$w
  one(sp, d2, "w * 2")
}
cat("DONE\n")
