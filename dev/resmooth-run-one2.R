# Lane wt-resmooth. Same as dev/resmooth-run-one.R, against the SECOND
# private library. It exists because the first library was being read by
# a long gated run when the last two fixes landed, and installing into a
# library another R process is reading is what destroyed this machine's
# shared library twice (dev/lane-rules.md).
#   Rscript dev/resmooth-run-one2.R <file>
args <- commandArgs(trailingOnly = TRUE)
f <- args[[1L]]
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib2",
            "C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages(library(frmtmb))
stopifnot(identical(normalizePath(find.package("frmtmb")),
                    normalizePath("C:/Users/adf44/source/r/wt-resmooth-lib2/frmtmb")))
pkg <- if (grepl("extensions/", f, fixed = TRUE)) {
  sub(".*extensions/([^/]+)/.*", "\\1", f)
} else "frmtmb"
if (pkg != "frmtmb") suppressMessages(library(pkg, character.only = TRUE))
res <- testthat::test_file(f, package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = "silent")
df <- as.data.frame(res)
tot <- function(nm) sum(df[[nm]])
cat(sprintf("RESULT %s lib=lib2 pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), tot("passed"), tot("failed"), tot("error"),
            tot("skipped"), tot("warning")))
bad <- df[df$failed > 0 | df$error > 0, , drop = FALSE]
if (nrow(bad)) {
  for (i in seq_len(nrow(bad))) {
    cat("  BAD:", bad$test[[i]], "\n")
    for (r in res[[which(vapply(res, function(x) identical(x$test,
                                                           bad$test[[i]]),
                                NA))[1L]]]$results) {
      if (inherits(r, c("expectation_failure", "expectation_error"))) {
        cat("    ", gsub("\n", " | ", substr(conditionMessage(r), 1, 400)),
            "\n")
      }
    }
  }
}
