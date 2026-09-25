# Library paths for lane wt-predfix. PREDFIX_ARM=base reads the released
# 0.62.0 build alone; anything else puts the lane's private library first.
PREDFIX_ARM <- Sys.getenv("PREDFIX_ARM", "lane")
.libPaths(c(if (PREDFIX_ARM != "base") "C:/Users/adf44/source/r/predfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm", PREDFIX_ARM, "frmtmb", format(packageVersion("frmtmb")),
    "from", dirname(find.package("frmtmb")), "\n")
