# Lane wt-arcovsample: RENDER the Rd topics this lane changed and grep
# the rendered text, because `%` starts a comment in Rd even inside
# verbatim macros and reading the source does not show that
# (lane-rules.md).
#   Rscript dev/arcovsample-rd.R
.libPaths(c("C:/Users/adf44/source/r/wt-arcovsample-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
W <- "C:/Users/adf44/source/r/frmtmb-wt-arcovsample"
files <- c(file.path(W, "man/frmtmb-sampling-api.Rd"),
           file.path(W, "man/frmtmb-autocor.Rd"),
           file.path(W, "extensions/frmtmb.sample/man/sample-log_lik.Rd"),
           file.path(W, "extensions/frmtmb.sample/man/sample-loo.Rd"))
want <- c("arma_cond_resp", "arma_cond_dpars", "one-step",
          "CONDITIONAL density", "cov = FALSE", "cov = TRUE",
          "max(p, q)")
for (f in files) {
  out <- tempfile(fileext = ".txt")
  tools::Rd2txt(f, out = out)
  txt <- paste(readLines(out, warn = FALSE), collapse = "\n")
  cat("\n==== ", basename(f), " (", nchar(txt), " chars rendered)\n",
      sep = "")
  for (w in want) {
    n <- length(gregexpr(w, txt, fixed = TRUE)[[1L]])
    if (regmatches(txt, regexpr(w, txt, fixed = TRUE))[1L] %in% w ||
          grepl(w, txt, fixed = TRUE)) {
      cat("  present x", n, ": ", w, "\n", sep = "")
    }
  }
  # any bare % would have eaten the rest of its line
  if (grepl("%", txt, fixed = TRUE)) cat("  NOTE: a % survived\n")
}
