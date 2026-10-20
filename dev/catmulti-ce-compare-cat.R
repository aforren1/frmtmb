# conditional_effects(categorical = TRUE) on the categorical(refcat = "b") fit: frmtmb
# (lane) against brms 2.23.0's on the same data, grid matched on x.
# Output: dev/catmulti-out/ce-compare-cat.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib", "C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("dev/catmulti-data.R")
dm <- catmulti_cat_data(11)
fm <- frm(bf(f ~ x), family = categorical(refcat = "b"), data = dm)
cf <- suppressMessages(conditional_effects(fm, categorical = TRUE))[[1]]
bm <- readRDS("dev/catmulti-out/brms-cat-refb.rds")
cb <- suppressMessages(brms::conditional_effects(bm, categorical = TRUE))[[1]]
cat("rows frmtmb", nrow(cf), "brms", nrow(cb), "\n")
cat("x grid identical:", isTRUE(all.equal(sort(unique(cf$x)), sort(unique(cb$x)))), "\n")
k <- match(paste(cb$x, cb$cats__), paste(cf$x, cf$cats__))
cat("matched rows:", sum(!is.na(k)), "\n")
d <- abs(cf$estimate__[k] - cb$estimate__) / cb$se__
cat(sprintf("|frmtmb - brms| / brms se: median %.3f max %.3f\n", median(d), max(d)))
r <- cf$se__[k] / cb$se__
cat(sprintf("se ratio frmtmb / brms: median %.3f range %.3f to %.3f\n", median(r), min(r), max(r)))
