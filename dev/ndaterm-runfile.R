# Lane ndaterm: run ONE test file in this process and print one RESULT
# line. Usage: Rscript dev/ndaterm-runfile.R <pkg> <file> <arm>
# arm "mine" puts the lane library first; "r9" is the base build alone;
# "coreonly" is the lane core (a scratch library) over rellib-r9.
# pkg is frmtmb or an extension name; file is a path to a test file.
a <- commandArgs(TRUE)
pkg <- a[1]; f <- a[2]; arm <- if (length(a) > 2) a[3] else "mine"
base <- "C:/Users/adf44/source/r"
libs <- c(if (arm == "mine") file.path(base, "wt-ndaterm-lib"),
          if (arm == "coreonly") "C:/Users/adf44/AppData/Local/Temp/3/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/ndaterm-corelib",
          file.path(base, "rellib-r9"),
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ndaterm/dev/stan-cache")
suppressMessages(library(testthat))
options(testthat.progress.max_fails = Inf)
suppressMessages(library(pkg, character.only = TRUE))
cat("lib:", find.package("frmtmb"), " ", pkg, "from", find.package(pkg),
    "\n")
res <- testthat::test_file(f, package = pkg,
                           env = testthat::test_env(pkg),
                           reporter = testthat::ProgressReporter$new(
                             show_praise = FALSE),
                           stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf("RESULT %s pass=%d fail=%d error=%d skip=%d warn=%d\n",
            basename(f), sum(df$passed), sum(df$failed),
            sum(df$error), sum(df$skipped), sum(df$warning)))
