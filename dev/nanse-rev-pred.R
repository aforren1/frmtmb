# Reviewer: does a prediction along a lost direction get an honest
# (infinite / NaN) standard error, or a finite one from the
# pseudo-inverse? Spread ridge y ~ a + b, a ~ 0 + f, b ~ 1.
#   Rscript dev/nanse-rev-pred.R base|lane|merge
args <- commandArgs(trailingOnly = TRUE)
arm <- if (length(args)) args[1] else "lane"
libs <- switch(arm,
  base = "C:/Users/adf44/source/r/rellib-r5",
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
spread <- function(k, m = 5, seed = 1) {
  set.seed(seed)
  f <- factor(rep(seq_len(k), each = m))
  data.frame(f = f, y = rnorm(k)[f] + rnorm(k * m, 0, 0.5))
}
cap <- function(expr) {
  w <- character()
  v <- withCallingHandlers(tryCatch(expr, error = function(e) {
    paste("ERROR:", conditionMessage(e))
  }), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  list(v = v, w = w)
}
for (k in c(10L, 60L)) {
  dd <- spread(k)
  r <- cap(frm(bf(y ~ a + b, a ~ 0 + f, b ~ 1, nl = TRUE), data = dd,
               family = gaussian()))
  fit <- r$v
  cat("\n== k =", k, "\n  frm() warnings:", length(r$w), "\n")
  for (x in r$w) cat("   ", substr(x, 1, 160), "\n")
  nd <- dd[c(1, 6), , drop = FALSE]
  # mu = a_k + b is identified; a alone and b alone are not
  for (dp in list(NULL, "a", "b")) {
    p <- cap(fitted(fit, newdata = nd, dpar = dp, summary = TRUE))
    lab <- if (is.null(dp)) "mu" else dp
    if (is.character(p$v)) {
      cat("  fitted(dpar =", lab, "):", p$v, "\n")
    } else {
      cat(sprintf("  fitted(dpar = %s): est %s  Est.Error %s  warnings %d\n",
                  lab, paste(signif(p$v[, 1], 4), collapse = " "),
                  paste(signif(p$v[, 2], 4), collapse = " "), length(p$w)))
    }
    for (x in p$w) cat("     ", substr(x, 1, 160), "\n")
  }
  s <- cap(summary(fit))
  cat("  summary warnings:", length(s$w), "\n")
  out <- capture.output(print(s$v))
  cat(paste0("   | ", head(out, 14)), sep = "\n")
  cat(paste0("   | ", tail(out, 8)), sep = "\n")
  ci <- cap(confint(fit))
  cat("  confint rows finite:", sum(is.finite(ci$v[, 1])), "of",
      nrow(ci$v), " warnings", length(ci$w), "\n")
}
