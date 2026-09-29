# Render the two Rd pages this lane touched and print the new sections.
# `%` starts a comment in Rd even inside verbatim macros, so an Rd is
# verified by RENDERING it, not by reading the source (lane-rules.md).
setwd("C:/Users/adf44/source/r/frmtmb-wt-csfactor")
render <- function(rd) {
  o <- tempfile()
  tools::Rd2txt(rd, out = o)
  readLines(o, warn = FALSE)
}
x <- render("man/frm.Rd")
cat("frm.Rd rendered lines:", length(x), "\n")
for (key in c("Monotonic effects", "Ordinal thresholds",
              "Category-specific effects", "treatment-contrast",
              "category boundary", "2.7e5")) {
  cat(sprintf("  %-28s at %s\n", key,
              paste(grep(key, x, fixed = TRUE), collapse = ", ")))
}
i <- grep("Category-specific effects", x, fixed = TRUE)
if (length(i)) {
  cat("\n", paste(x[i[1L]:min(length(x), i[1L] + 50L)], collapse = "\n"),
      "\n")
}
v <- render("man/variables.Rd")
cat("\nvariables.Rd bcs_ at",
    paste(grep("bcs_", v, fixed = TRUE), collapse = ", "), "\n")
