# Lane remeasure, bin 2P G3: the case where the adapter's weighting is
# ABSENT. A scratch copy of dev/brmsport-bin2-adapters.R with the
# `ll * w[i]` line removed is staged beside core's test helpers exactly
# as dev/brmsport-bin2-run.R stages the real one, and the generated
# log_lik file is recorded against it. brms's log_lik:32 asserts the
# weighted value, so it must stop holding while log_lik:19 (the same
# density unweighted) keeps holding.
#
#   Rscript dev/remeasure-bin2-absent.R \
#     > dev/remeasure-out/probe/bin2-absent.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = "frmtmb")
stage <- "dev/remeasure-out/absent-stage"
unlink(stage, recursive = TRUE)
dir.create(stage, recursive = TRUE)
file.copy(list.files("tests/testthat", "^helper.*[.]R$", full.names = TRUE),
          stage)
ad <- readLines("dev/brmsport-bin2-adapters.R")
hit <- grep("ll <- ll * w[i]", ad, fixed = TRUE)
stopifnot(length(hit) == 1L)
cat("removed adapter line", hit, ":", trimws(ad[hit]), "\n")
writeLines(ad[-hit], file.path(stage, "helper-brms-bin2-adapters.R"))
f <- "dev/brmsport-bin2/frmtmb/test-brms-suite-bin2-log_lik.R"
file.copy(f, stage)
rec <- normalizePath(file.path(stage, "rec.tsv"), winslash = "/",
                     mustWork = FALSE)
Sys.setenv(FRMTMB_BRMSPORT_RECORD = rec)
suppressMessages({library(testthat); library(frmtmb)})
r <- as.data.frame(test_file(file.path(stage, basename(f)),
                             package = "frmtmb",
                             env = testthat::test_env("frmtmb"),
                             reporter = "silent", stop_on_failure = FALSE))
cat("RESULT pass=", sum(r$passed), " fail=", sum(r$failed), " err=",
    sum(r$error), "\n", sep = "")
x <- utils::read.delim(rec, header = FALSE, quote = "",
                       colClasses = "character")
x <- x[x[[1]] == "assert" & x[[3]] %in% c("log_lik:19", "log_lik:32"), ]
for (k in seq_len(nrow(x))) {
  cat(x[k, 3], " held=", x[k, 5], " msg=", substr(x[k, 7], 1, 300), "\n",
      sep = "")
}
