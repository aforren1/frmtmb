# Run ONE test file in ONE process against a named library, counting
# failures AND errors, and naming every failing block. A runner that
# sums `failed` and not `error` prints a clean line for a file that
# aborted (dev/lane-rules.md).
#   Rscript dev/brmsnames-runtest.R <LIB> <file> [<LIB2> ...]
av <- commandArgs(trailingOnly = TRUE)
libs <- av[-2]
.libPaths(c(libs, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
file <- av[2]
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = normalizePath("dev/brmsnames-log/emptycache", mustWork = FALSE))
library(testthat)
pkg <- regmatches(file, regexpr("frmtmb[.][a-z]+", file))
if (!length(pkg)) pkg <- "frmtmb"
library(pkg, character.only = TRUE)
cat("libs: ", paste(libs, collapse = " | "), "\n")
cat(pkg, " from: ", dirname(system.file(package = pkg)), "\n", sep = "")
cat("frmtmb from: ", dirname(system.file(package = "frmtmb")), "\n")
res <- test_file(file, reporter = "silent", package = pkg)
df <- as.data.frame(res)
np <- sum(df$passed); nf <- sum(df$failed); ne <- sum(df$error)
ns <- sum(df$skipped)
cat("FILE   ", file, "\n")
cat("BLOCKS ", nrow(df), "\n")
cat("PASS   ", np, "\n")
cat("FAIL   ", nf, "\n")
cat("ERROR  ", ne, "\n")
cat("SKIP   ", ns, "\n")
for (i in seq_len(nrow(df))) {
  if (df$failed[i] > 0 || df$error[i]) {
    cat("  failing: ", df$test[i], "\n", sep = "")
    for (r in res[[i]]$results) {
      if (inherits(r, "expectation_failure") ||
          inherits(r, "expectation_error"))
        cat("      ", substr(gsub("[\r\n]+", " ", conditionMessage(r)),
                             1, 160), "\n")
    }
  }
}
cat("DONE\n")
