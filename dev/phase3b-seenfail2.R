# Punch round 2: run the two changed test files against a given build,
# so the new blocks are seen failing on the released build (rellib-r3)
# and on the round-1 build (phase3b-lib, before its reinstall).
# Usage: Rscript dev/phase3b-seenfail2.R <lib> <label>
# Output: dev/phase3b-log/seenfail-punch2-<label>.txt
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(testthat))
out_file <- sprintf("dev/phase3b-log/seenfail-punch2-%s.txt", a[2])
td <- "extensions/frmtmb.eam/tests/testthat"
out <- c(sprintf("build: %s %s", find.package("frmtmb.eam"),
                 packageVersion("frmtmb.eam")),
         sprintf("run at %s", format(Sys.time())))
for (f in c("test-wiener-cdf.R", "test-contaminant.R")) {
  r <- as.data.frame(test_file(file.path(td, f), reporter = "silent",
                               package = "frmtmb.eam",
                               load_package = "installed"))
  bad <- r[r$failed > 0 | r$error, c("test", "nb", "failed", "error")]
  out <- c(out, "", sprintf("%s: %d expectations, %d failed, %d blocks errored",
                            f, sum(r$nb), sum(r$failed), sum(r$error)),
           capture.output(print(bad, right = FALSE)))
}
unlink(file.path(td, "_problems"), recursive = TRUE)
writeLines(out, out_file)
cat(out, sep = "\n")
