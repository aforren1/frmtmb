# Lane prilabel: one package's whole suite in ONE R session, as covr and
# R CMD check run it (lane-env.md, "Test isolation"). Arg 1: the package
# (frmtmb or frmtmb.sample).
local({
  p <- commandArgs(TRUE)[1]
  .libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
              "C:/Users/adf44/source/r/rellib-r9",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
  wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
  d <- if (p == "frmtmb") file.path(wt, "tests/testthat") else
    file.path(wt, "extensions", p, "tests/testthat")
  suppressMessages(library(testthat))
  suppressMessages(library(p, character.only = TRUE))
  cat("lib:", p, dirname(find.package(p)), "; frmtmb",
      dirname(find.package("frmtmb")), "\n")
  res <- as.data.frame(testthat::test_dir(d, package = p,
                                          load_package = "none",
                                          stop_on_failure = FALSE,
                                          reporter = "silent"))
  cat(sprintf("ONESESSION %s files=%d pass=%d fail=%d err=%d skip=%d\n", p,
              length(unique(res$file)), sum(res$passed), sum(res$failed),
              sum(res$error), sum(res$skipped)))
  bad <- res[res$failed > 0 | res$error, c("file", "test")]
  if (nrow(bad)) print(bad, row.names = FALSE)
})
