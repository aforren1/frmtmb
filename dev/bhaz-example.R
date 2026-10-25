# Run the help example of frmtmb_register_aterm() on the lane build and
# check that it leaves the addition-term registry as it found it.
#   Rscript dev/bhaz-example.R > dev/bhaz-log/example.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
before <- ns$frmtmb_aterm_registry$reg
feat0 <- frm_compat_features()$name
example("frmtmb_register_aterm", package = "frmtmb", run.dontrun = FALSE)
cat("registry identical after the example:",
    identical(ns$frmtmb_aterm_registry$reg, before), "\n")
cat("knownsd registered after the example:",
    "knownsd" %in% names(ns$frmtmb_aterm_registry$reg), "\n")
cat("compat features identical:",
    identical(frm_compat_features()$name, feat0), "\n")
