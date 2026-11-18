# Line coverage of core R/ by core's own test suite, in ONE R process:
# covr::environment_coverage() traces the namespace of a keep.source
# build of the worktree (installed in the lane library with
# R CMD INSTALL --with-keep.source) and runs a driver that calls
# testthat::test_file() on each test file in turn. One process, not
# package_coverage()'s parent plus child, to stay within the lane's
# two-process limit while the static scripts run.
#   Rscript dev/codesize-covr-env.R [pattern]
# Writes dev/codesize-cov/cov.rds (covr object) and
# dev/codesize-cov/results.tsv (pass/fail/err/skip/seconds per file).
# Settings as .github/workflows/test-coverage.yaml: NOT_CRAN=true, gated
# brms fit tests off. Unlike CI, the extensions are on the library path
# (rellib-r22), so tests that use an installed extension run.
a <- commandArgs(trailingOnly = TRUE)
pat <- if (length(a)) a[1] else "^test-.*[.]R$"
wt <- "C:/Users/adf44/source/r/frmtmb-wt-codesize"
.libPaths(c("C:/Users/adf44/source/r/wt-codesize-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           FRMTMB_STAN_CACHE = file.path(wt, "dev/stan-cache"),
           R_PARALLEL_PORT = "11517")
out <- file.path(wt, "dev/codesize-cov")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
library(testthat)
library(frmtmb)
cat("frmtmb from", find.package("frmtmb"), "\n")
fs <- list.files(file.path(wt, "tests/testthat"), pat)
res <- file.path(out, if (length(a)) "results-probe.tsv" else "results.tsv")
drv <- file.path(out, "driver.R")
writeLines(c(
  sprintf("fs <- %s", paste(deparse(fs), collapse = "")),
  sprintf("res <- %s", deparse(res)),
  sprintf("td <- %s", deparse(file.path(wt, "tests/testthat"))),
  sprintf("snap <- %s", deparse(sub("[.]tsv$", "-snap.rds", res))),
  "cat('file\\tpass\\tfail\\terr\\tskip\\tsec\\tmb\\n', file = res)",
  "i <- 0",
  "for (f in fs) {",
  "  t0 <- proc.time()[[3]]",
  "  r <- tryCatch(as.data.frame(testthat::test_file(file.path(td, f),",
  "    package = 'frmtmb', reporter = 'silent',",
  "    env = testthat::test_env('frmtmb'))), error = function(e) NULL)",
  "  s <- if (is.null(r)) c(NA, NA, 1, NA) else",
  "    c(sum(r$passed), sum(r$failed), sum(r$error), sum(r$skipped))",
  "  g <- gc(); i <- i + 1",
  "  cat(f, s, round(proc.time()[[3]] - t0, 1), round(sum(g[, 2])),",
  "      sep = '\\t', file = res, append = TRUE)",
  "  cat('\\n', file = res, append = TRUE)",
  "  # a snapshot of the counters, so a run stopped early still reports",
  "  if (i %% 15 == 0) saveRDS(structure(as.list(covr:::.counters),",
  "    class = 'coverage'), snap)",
  "}"), drv)
t0 <- proc.time()[[3]]
cov <- covr::environment_coverage(asNamespace("frmtmb"), test_files = drv)
saveRDS(cov, file.path(out, if (length(a)) "cov-probe.rds" else "cov.rds"))
cat("done in", round(proc.time()[[3]] - t0), "s\n")
df <- as.data.frame(cov)
cat("records", nrow(df), "files", length(unique(df$filename)), "\n")
print(head(sort(table(df$filename), decreasing = TRUE)))
print(covr::percent_coverage(cov))
