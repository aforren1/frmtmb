# Verify the two edited manual pages by RENDERING them, not by reading
# the source: `%` starts a comment in Rd even inside verbatim macros, and
# a dropped character is invisible in the .Rd file.

.libPaths(c("C:/Users/adf44/source/r/wt-gradcheck-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (f in c("man/frmtmb_control.Rd", "man/diagnose.Rd")) {
  tf <- tempfile()
  tools::Rd2txt(f, out = tf)
  x <- readLines(tf, warn = FALSE)
  cat("=== ", f, " rendered lines:", length(x), "\n")
  i <- grep(paste(c("DOES NOT COVER", "TWO SCALES",
                    "IMPORTANCE-CORRECTED", "still warned", "297",
                    "four ways", "grad_proj_par", "0.45 s",
                    "Newton step", "bound holds", "held in place"),
                  collapse = "|"), x)
  for (k in sort(unique(c(i, i + 1L)))) {
    if (k <= length(x)) cat("  |", x[[k]], "\n")
  }
  cat("  stray percent signs:", sum(grepl("%", x)), "\n\n")
}
