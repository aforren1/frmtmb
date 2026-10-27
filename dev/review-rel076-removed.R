# Review: names in frmtmb 0.75.0's namespace (rellib-r13) that 0.76.0
# (rellib-r14) lacks, and whether any extension's R/ or tests still
# names one; also frmtmb.* namespaces' frmtmb:: or frmtmb::: references
# that 0.76.0 does not define.
#   Rscript dev/review-rel076-removed.R
rel <- "C:/Users/adf44/source/r/frmtmb-wt-rel076"
nm <- function(lib) {
  callr_out <- system2("C:/Program Files/R/R-4.6.1/bin/Rscript.exe",
    c("-e", shQuote(sprintf(
      ".libPaths(c('%s','C:/Users/adf44/AppData/Local/R/win-library/4.6'));cat(ls(asNamespace('frmtmb'),all.names=TRUE),sep='\\n')",
      lib))), stdout = TRUE)
  callr_out
}
old <- nm("C:/Users/adf44/source/r/rellib-r13")
new <- nm("C:/Users/adf44/source/r/rellib-r14")
gone <- setdiff(old, new)
cat("0.75.0 names:", length(old), " 0.76.0 names:", length(new),
    " removed:", length(gone), "\n")
cat("  ", gone, "\n")
cat("added:", setdiff(new, old), "\n")
files <- c(list.files(file.path(rel, "extensions"), "[.][Rr]$",
                      recursive = TRUE, full.names = TRUE),
           list.files(file.path(rel, "tests"), "[.][Rr]$", recursive = TRUE,
                      full.names = TRUE),
           list.files(file.path(rel, "R"), "[.][Rr]$", full.names = TRUE),
           list.files(file.path(rel, "vignettes"), "[.][Rr]md$",
                      full.names = TRUE),
           list.files(file.path(rel, "extensions"), "[.][Rr]md$",
                      recursive = TRUE, full.names = TRUE))
for (g in gone) {
  hit <- character()
  for (f in files) {
    l <- readLines(f, warn = FALSE)
    w <- grep(paste0("\\b", gsub(".", "[.]", g, fixed = TRUE), "\\b"), l)
    if (length(w)) hit <- c(hit, paste0(sub(rel, "", f), ":", w))
  }
  cat("removed", g, "named at:", if (length(hit)) hit else "none", "\n")
}
# every frmtmb::name / frmtmb:::name in the extensions and the tests
refs <- character()
for (f in files) {
  l <- readLines(f, warn = FALSE)
  m <- regmatches(l, gregexpr("frmtmb:::?[A-Za-z._][A-Za-z0-9._]*", l))
  refs <- c(refs, unlist(m))
}
refs <- unique(sub("^frmtmb:::?", "", refs))
cat("frmtmb:: names referenced:", length(refs), " undefined in 0.76.0:",
    setdiff(refs, new), "\n")
