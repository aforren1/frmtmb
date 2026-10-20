# brms-B correction: stancode() of a categorical category with its own
# formula, for the first and the second predicted category.
# Output: dev/catmulti-out/brms-b.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
source("dev/catmulti-data.R")
dd <- catmulti_cat_data(42)
tr <- function(nm, e) cat(nm, ":", tryCatch({e; "OK"}, error = function(x) conditionMessage(x)), "\n")
tr("refcat = c, mua ~ z (first predicted)",
   stancode(bf(f ~ x, mua ~ z), data = dd, family = categorical(refcat = "c")))
tr("refcat = c, mub ~ z (second predicted)",
   stancode(bf(f ~ x, mub ~ z), data = dd, family = categorical(refcat = "c")))
tr("default, muc ~ z (second predicted)",
   stancode(bf(f ~ x, muc ~ z), data = dd, family = categorical()))
