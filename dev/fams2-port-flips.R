# Every ported brms assertion that now HOLDS on the lane build where the
# committed verdict says otherwise, for one package, uncapped.
#   Rscript dev/fams2-port-flips.R <pkg> <file>
args <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
root <- "C:/Users/adf44/source/r/frmtmb-wt-fams2"
tdir <- if (args[1] == "frmtmb") file.path(root, "tests/testthat") else
  file.path(root, "extensions", args[1], "tests/testthat")
suppressPackageStartupMessages(library(testthat))
rep <- ListReporter$new()
testthat::test_file(file.path(tdir, args[2]), package = args[1],
                    env = testthat::test_env(args[1]), reporter = rep,
                    load_package = "installed")
msgs <- unlist(lapply(rep$get_results(), function(t) {
  vapply(t$results, function(r) {
    if (inherits(r, "expectation_failure") || inherits(r, "expectation_error"))
      conditionMessage(r) else NA_character_
  }, "")
}))
msgs <- msgs[!is.na(msgs)]
ids <- regmatches(msgs, regexpr("[a-z_-]+:[0-9]+(?= now HOLDS)", msgs, perl = TRUE))
cat(args[1], args[2], ": failures", length(msgs), "; now HOLDS", length(ids), "\n")
cat(sort(ids), sep = "\n")
other <- msgs[!grepl("now HOLDS", msgs)]
if (length(other)) { cat("OTHER FAILURES:\n"); cat(other, sep = "\n") }
