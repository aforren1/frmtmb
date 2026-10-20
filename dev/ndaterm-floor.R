# Lane ndaterm: frmtmb.sample's fix on the 0.71.0 core. The lane's
# frmtmb.sample alone in a scratch library, then rellib-r9.
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/3/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/ndaterm-libsample", "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-ndaterm/dev/stan-cache")
suppressMessages({library(testthat); library(frmtmb.sample)})
cat("frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\nfrmtmb.sample", find.package("frmtmb.sample"), "\n")
f <- "extensions/frmtmb.sample/tests/testthat/test-newdata-aterms-draws.R"
df <- as.data.frame(test_file(f, package = "frmtmb.sample",
                              env = test_env("frmtmb.sample"),
                              reporter = "silent", stop_on_failure = FALSE))
cat(sprintf("RESULT pass=%d fail=%d error=%d skip=%d\n", sum(df$passed),
            sum(df$failed), sum(df$error), sum(df$skipped)))
