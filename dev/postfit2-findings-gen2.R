# Fills the GENERATED-CELEVEL-* blocks of dev/postfit2-findings.md
# (section 7b) from the logs, and the gated-tier block of section 7.
#   Rscript dev/postfit2-findings-gen2.R
wt <- "C:/Users/adf44/source/r/frmtmb-wt-postfit2"
lg <- file.path(wt, "dev/postfit2-log")
rd <- function(f) readLines(file.path(lg, f), warn = FALSE)
res <- function(files) {
  unlist(lapply(files, function(f) {
    p <- file.path(lg, f)
    if (!file.exists(p)) return(paste("MISSING", f))
    r <- grep("^RESULT", readLines(p, warn = FALSE), value = TRUE)
    if (!length(r)) paste("NO RESULT", f) else paste0(r, "   [", f, "]")
  }))
}
f <- file.path(wt, "dev/postfit2-findings.md")
x <- readLines(f, warn = FALSE)
put <- function(x, tag, val) {
  i <- which(x == tag)
  if (length(i) != 1L) stop("marker ", tag, " found ", length(i), " times")
  c(x[seq_len(i - 1L)], val, x[-seq_len(i)])
}
x <- put(x, "GENERATED-CELEVEL-BASE", rd("celevel-base.txt"))
x <- put(x, "GENERATED-CELEVEL-LANE", rd("celevel-lane.txt"))
x <- put(x, "GENERATED-CELEVEL-BRMS", rd("celevel-brms.txt"))
x <- put(x, "GENERATED-CELEVEL-TESTS", c(
  "base build:",
  res(c("before-core-test-ce-options.R.txt",
        "before-sample-test-postfit-draws.R.txt",
        "before-core-test-brms-methods.R.txt")),
  "this lane, after the fix, one process per file:",
  res(c("after3-core-test-ce-options.R.txt",
        "after3-core-test-ce-bands.R.txt", "after3-core-test-nlf.R.txt",
        "after3-sample-test-postfit-draws.R.txt",
        "after3-sample-test-conditional-effects-draws.R.txt",
        "after3-sample-test-draws-methods.R.txt",
        "after3-sample-test-re-formula-draws.R.txt",
        "after3-sample-test-predfix-new-levels.R.txt")),
  "the whole ungated suites again, after the fix (dev/postfit2-suite.sh):",
  rd("suite-summary.txt"),
  "the rest of the gated tier (dev/postfit2-gated.sh):",
  rd("gated-tier.txt"),
  "R CMD check --as-cran, final (dev/postfit2-check.ps1):",
  vapply(c("frmtmb", "frmtmb.sample"), function(p) {
    f <- file.path(wt, "dev/postfit2-check", p, paste0(p, ".Rcheck"),
                   "00check.log")
    if (!file.exists(f)) return(paste(p, "NO CHECK LOG"))
    s <- grep("^Status:", readLines(f, warn = FALSE), value = TRUE)
    paste(p, if (length(s)) s else "NO STATUS LINE")
  }, "")))
writeLines(x, f)
cat("filled\n")
