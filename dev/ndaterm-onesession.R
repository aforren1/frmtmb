# Lane ndaterm: one package's whole suite in ONE session, in file order,
# as covr and R CMD check run it. Usage:
#   Rscript dev/ndaterm-onesession.R <pkg>
# Log: the caller's redirection, dev/ndaterm-out/onesession-<pkg>.txt
pkg <- commandArgs(TRUE)[1]
base <- "C:/Users/adf44/source/r"
.libPaths(c(file.path(base, "wt-ndaterm-lib"), file.path(base, "rellib-r9"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ndaterm/dev/stan-cache")
suppressMessages(library(testthat))
options(testthat.progress.max_fails = Inf)
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package("frmtmb"), " ", pkg, "from", find.package(pkg),
    "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-ndaterm"
dir <- if (pkg == "frmtmb") file.path(wt, "tests/testthat") else
  file.path(wt, "extensions", pkg, "tests/testthat")
res <- testthat::test_dir(dir, package = pkg, load_package = "none",
                          stop_on_failure = FALSE,
                          reporter = testthat::SummaryReporter$new(
                            max_reports = Inf, show_praise = FALSE))
df <- as.data.frame(res)
cat(sprintf("ONESESSION %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            pkg, sum(df$passed), sum(df$failed), sum(df$error),
            sum(df$skipped), sum(df$warning)))
bad <- df[df$failed > 0 | df$error, c("file", "test")]
if (nrow(bad)) print(bad, row.names = FALSE)
