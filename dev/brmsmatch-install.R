## Roxygenise and install frmtmb.sample into this lane's private
## library. Core is taken unchanged from the round's shared reference
## build, which is read-only.
.libPaths(c("C:/Users/adf44/source/r/brmsmatch-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkg <- "C:/Users/adf44/source/r/frmtmb-wt-brmsmatch/extensions/frmtmb.sample"
roxygen2::roxygenise(pkg, load_code = roxygen2::load_pkgload)
cat("ROXYGENISED\n")
