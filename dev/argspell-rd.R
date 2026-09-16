## Render every Rd file this lane touched and grep the OUTPUT, not the
## source: `%` starts a comment in Rd even inside \preformatted{}, and a
## backtick span beginning with `r ` is executed as inline R by roxygen
## markdown. Both fail silently in the source.
##
## Run: Rscript dev/argspell-rd.R
.libPaths(c("C:/Users/adf44/source/r/argspell-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
files <- c(
  "man/fitted.frmtmb_fit.Rd", "man/predict.frmtmb_fit.Rd",
  "man/residuals.frmtmb_fit.Rd", "man/simulate.frmtmb_fit.Rd",
  "man/conditional_effects.Rd", "man/pp_check.Rd",
  "man/frm_bootstrap.Rd", "man/dharma_residuals.Rd",
  "man/confint.frmtmb_fit.Rd", "man/hypothesis.Rd",
  "man/frmtmb-sampling-api.Rd",
  "extensions/frmtmb.sample/man/posterior_epred.Rd",
  "extensions/frmtmb.sample/man/sample-posterior_summary.Rd")
bad <- 0L
for (f in files) {
  if (!file.exists(f)) { cat("MISSING", f, "\n"); bad <- bad + 1L; next }
  out <- tempfile(fileext = ".txt")
  ok <- tryCatch({ tools::Rd2txt(f, out = out); TRUE },
                 error = function(e) { cat("RENDER FAILED", f, ":",
                                           conditionMessage(e), "\n"); FALSE })
  if (!ok) { bad <- bad + 1L; next }
  txt <- readLines(out, warn = FALSE)
  cat(sprintf("%-52s %4d lines\n", f, length(txt)))
  # an unescaped % eats the rest of its line, so a line that ends where
  # a % was is the symptom; look for the literal instead
  pct <- grep("%", txt, value = TRUE, fixed = TRUE)
  if (length(pct)) cat("   literal %:", length(pct), "line(s)\n")
  rr <- grep("`r ", txt, value = TRUE, fixed = TRUE)
  if (length(rr)) {
    cat("   INLINE-R SPAN LEFT IN OUTPUT:\n")
    cat(paste0("     ", rr), sep = "\n")
    bad <- bad + 1L
  }
}
cat(sprintf("\nRd files with a problem: %d of %d\n", bad, length(files)))

cat("\n---- fitted() usage, as rendered ----\n")
out <- tempfile(fileext = ".txt")
tools::Rd2txt("man/fitted.frmtmb_fit.Rd", out = out)
txt <- readLines(out, warn = FALSE)
i <- grep("^Usage:", txt)
if (length(i)) writeLines(txt[i[1L]:min(i[1L] + 10L, length(txt))])
cat("\n---- predict() usage, as rendered ----\n")
tools::Rd2txt("man/predict.frmtmb_fit.Rd", out = out)
txt <- readLines(out, warn = FALSE)
i <- grep("^Usage:", txt)
if (length(i)) writeLines(txt[i[1L]:min(i[1L] + 10L, length(txt))])
