# genrev round 3: run frmtmb.sample's test-draws-methods.R against a
# library stack, counting failures AND errors, and print the message of
# every failing expectation so the REASON for a failure is on record.
a <- commandArgs(trailingOnly = TRUE)
LIBS <- strsplit(a[1], ";", fixed = TRUE)[[1]]
FILE <- a[2]
.libPaths(c(LIBS, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
library(testthat)
suppressMessages(library(frmtmb.sample))
td <- dirname(FILE)
for (h in list.files(td, pattern = "^(helper|setup).*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = globalenv())
}
rep <- ListReporter$new()
test_file(FILE, reporter = rep, package = "frmtmb.sample")
res <- rep$get_results()
np <- nf <- ne <- ns <- 0L
fails <- character()
for (t in res) {
  for (r in t$results) {
    cl <- class(r)[1]
    if (cl == "expectation_success") np <- np + 1L
    else if (cl == "expectation_skip") ns <- ns + 1L
    else {
      if (cl == "expectation_error") ne <- ne + 1L else nf <- nf + 1L
      fails <- c(fails, sprintf("  [%s] %s :: %s", cl, t$test,
        gsub("\\s+", " ", substr(conditionMessage(r), 1, 160))))
    }
  }
}
cat(sprintf("LIBS %s\n", paste(basename(LIBS), collapse = " > ")))
cat(sprintf("frmtmb.sample source: %s\n",
            dirname(system.file(package = "frmtmb.sample"))))
cat(sprintf("ASSERT %d PASS %d FAIL %d ERROR %d SKIP %d\n",
            np + nf + ne, np, nf, ne, ns))
if (length(fails)) cat(fails, sep = "\n")
cat("GENREVDONE\n")
