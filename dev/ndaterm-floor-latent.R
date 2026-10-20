# Lane ndaterm, fix round 1: frmtmb.latent's hmm fix on the 0.71.0
# core. The lane's frmtmb.latent alone in a scratch library, then
# rellib-r9. Log dev/ndaterm-out/floor-latent.txt
.libPaths(c("C:/Users/adf44/AppData/Local/Temp/3/claude/c--Users-adf44-source-r-frmtmb/7e21965d-23fa-4f24-9801-8cad9bc17105/scratchpad/ndaterm-latlib", "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages({library(testthat); library(frmtmb.latent)})
cat("frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "\nfrmtmb.latent", find.package("frmtmb.latent"), "\n")
f <- "extensions/frmtmb.latent/tests/testthat/test-hmm-aterms.R"
df <- as.data.frame(test_file(f, package = "frmtmb.latent",
                              env = test_env("frmtmb.latent"),
                              reporter = "silent", stop_on_failure = FALSE))
cat(sprintf("RESULT pass=%d fail=%d error=%d skip=%d\n", sum(df$passed),
            sum(df$failed), sum(df$error), sum(df$skipped)))
