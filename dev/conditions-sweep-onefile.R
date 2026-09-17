# Lane wt-conditions: run ONE test file in ONE process, counting
# failures AND errors, and record the class of every condition the
# file's own expect_error() calls caught.
#   Rscript dev/conditions-sweep-onefile.R lane|base <file> <outdir>
# Writes <outdir>/<file with / as _>.log and .tsv.
av <- commandArgs(trailingOnly = TRUE)
arm <- av[1L]; file <- av[2L]; outdir <- av[3L]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                             mustWork = FALSE))
suppressMessages(library(testthat))
pkg <- regmatches(file, regexpr("frmtmb[.][a-z]+", file))
if (!length(pkg)) pkg <- "frmtmb"
suppressMessages(library(pkg, character.only = TRUE))
cat("frmtmb from:", dirname(system.file(package = "frmtmb")), "\n")
cat("StanHeaders:", format(tryCatch(packageVersion("StanHeaders"),
                                    error = function(e) "absent")), "\n")

# expect_error() returns the condition it caught; a trace on both the
# namespace and the attached binding sees every call a test makes
# (dev/conditions-traceprobe.R established that one binding is not
# enough)
sweep_env <- new.env()
sweep_env$rows <- list()
sweep_record <- function(v) {
  if (inherits(v, "condition")) {
    sweep_env$rows[[length(sweep_env$rows) + 1L]] <- data.frame(
      file = file, class = paste(class(v), collapse = "/"),
      message = substr(gsub("[\r\n\t]+", " ", conditionMessage(v)), 1, 300),
      stringsAsFactors = FALSE)
  }
}
for (w in list(asNamespace("testthat"),
               as.environment("package:testthat"))) {
  suppressMessages(trace("expect_error", where = w, print = FALSE,
    exit = bquote(.(sweep_record)(returnValue(default = NULL)))))
}

res <- test_file(file, reporter = "silent", package = pkg)
df <- as.data.frame(res)
line <- sprintf("BLOCKS %d PASS %d FAIL %d ERROR %d SKIP %d  %s",
                nrow(df), sum(df$passed), sum(df$failed), sum(df$error),
                sum(df$skipped), file)
cat(line, "\n")
for (i in seq_len(nrow(df))) {
  if (df$failed[i] > 0 || df$error[i]) {
    cat("  failing: ", df$test[i], "\n", sep = "")
    for (r in res[[i]]$results) {
      if (inherits(r, c("expectation_failure", "expectation_error")))
        cat("      ", substr(gsub("[\r\n]+", " ", conditionMessage(r)),
                             1, 600), "\n")
    }
  }
}
rows <- do.call(rbind, sweep_env$rows)
if (is.null(rows)) {
  rows <- data.frame(file = character(0), class = character(0),
                     message = character(0))
}
cat("CAUGHT", nrow(rows), "\n")
utils::write.table(rows, file.path(outdir, paste0(gsub("/", "_", file),
                                                  ".tsv")),
                   sep = "\t", row.names = FALSE, quote = TRUE)
cat("DONE\n")
