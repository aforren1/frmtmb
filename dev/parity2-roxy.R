# Lane parity2: roxygenise core and print every message, to read the
# "Caused by error" line the install script's filter shows alone.
.libPaths(c("C:/Users/adf44/source/r/wt-parity2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
arg <- commandArgs(TRUE)[1]
if (is.na(arg)) arg <- "C:/Users/adf44/source/r/frmtmb-wt-parity2"
roxygen2::roxygenise(arg)
