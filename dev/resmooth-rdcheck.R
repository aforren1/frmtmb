# Lane wt-resmooth, nits round. Render the Rd files this round touched
# and grep the RENDERED text, because `%` starts a comment in Rd even
# inside verbatim macros and reading the source does not show that.
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (f in c("man/frm_linpred.Rd", "man/simulate.frmtmb_fit.Rd",
            "man/pp_check.Rd", "man/conditional_effects.Rd",
            "man/frm_bootstrap.Rd", "man/frm_lp_basis.Rd")) {
  tt <- tempfile()
  tools::Rd2txt(f, out = tt)
  x <- readLines(tt, warn = FALSE)
  cat("==", f, "|", length(x), "lines |",
      sum(nchar(x) > 0), "non-empty\n")
  hit <- grep("by = f|invert|percent|8 percent", x, value = TRUE)
  if (length(hit)) cat(paste0("   ", hit), sep = "\n")
  # a stray unescaped percent sign eats the rest of its Rd line, so the
  # rendered text losing one is the symptom to look for
  src <- readLines(f, warn = FALSE)
  bad <- grep("(^|[^\\\\])%", src, value = TRUE)
  if (length(bad)) {
    cat("   UNESCAPED % in the source:\n")
    cat(paste0("     ", substr(bad, 1, 70)), sep = "\n")
  }
}
cat("DONE\n")
