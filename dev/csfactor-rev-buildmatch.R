# Does C:/Users/adf44/source/r/wt-csfactor-lib match the worktree source?
# Compare the deparsed body of every function the diff touched against a
# fresh parse of the R/ files, so a stale install is caught by content
# rather than by a timestamp.
.libPaths(c("C:/Users/adf44/source/r/wt-csfactor-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-csfactor"
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
cat("sample ", as.character(packageVersion("frmtmb.sample")), "from",
    dirname(system.file(package = "frmtmb.sample")), "\n")

src <- new.env()
files <- c("frame.R", "predict.R", "compat.R", "fit.R", "confint.R",
           "influence.R")
for (f in files) sys.source(file.path(WT, "R", f), envir = src,
                            keep.source = FALSE)
ns <- asNamespace("frmtmb")
nm <- ls(src, all.names = TRUE)
bad <- character()
miss <- character()
for (n in nm) {
  if (!is.function(get(n, src))) next
  if (!exists(n, envir = ns, inherits = FALSE)) { miss <- c(miss, n); next }
  a <- paste(deparse(get(n, src)), collapse = "\n")
  b <- paste(deparse(get(n, envir = ns)), collapse = "\n")
  if (!identical(a, b)) bad <- c(bad, n)
}
cat("functions compared:", length(nm), "\n")
cat("MISSING from installed namespace:", length(miss), "\n")
if (length(miss)) cat("  ", paste(miss, collapse = " "), "\n")
cat("DIFFERENT body:", length(bad), "\n")
if (length(bad)) cat("  ", paste(bad, collapse = " "), "\n")

# the three new functions must exist
for (n in c("cs_term_design", "cs_term_newdata", "check_cs_identified",
            "cs_newdata_columns")) {
  cat(n, ":", exists(n, envir = ns, inherits = FALSE), "\n")
}
# installed test file vs worktree test file
for (tf in c("test-cs-factor.R", "test-v17.R", "test-review-v29.R",
             "test-brms-likelihood.R")) {
  p1 <- file.path(WT, "tests", "testthat", tf)
  p2 <- system.file("", package = "frmtmb")
  cat("testfile", tf, "exists in worktree:", file.exists(p1), "\n")
}
cat("done\n")
