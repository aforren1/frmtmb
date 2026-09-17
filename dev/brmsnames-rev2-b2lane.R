## Reviewer recheck: B2 without gr(by =), which frmtmb does not support.
## Lane names against brms's B2 fit (dev/brmsnames-rev2-log/brms-d.txt).
##   Rscript dev/brmsnames-rev2-b2lane.R
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb)})
source("dev/brmsnames-rev2-data.R")
d <- rev2_data()
fit <- suppressWarnings(frm(bf(y ~ ns(x, df = 3) + bs(z, df = 4) + mo(xo) +
                                 (1 | g)), family = gaussian(), data = d))
print(variables(fit))
print(frmtmb:::brms_par_labels(fit))
