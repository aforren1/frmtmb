# Reviewer: a brms multinomial(refcat =) object as family =.
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib", "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
suppressPackageStartupMessages(library(frmtmb))
source("dev/catmulti-rev-data.R")
dm <- rev_mult_data(31)
f <- frm(bf(Yc | trials(N) ~ x), family = brms::multinomial(refcat = "y3"), data = dm)
cat("brms::multinomial(refcat = 'y3') dpars:", names(fixef_by_dpar(f)), "refcat", family(f)$refcat, "\n")
f <- frm(bf(Yc | trials(N) ~ x), family = brms::multinomial(refcat = NA), data = dm)
cat("brms::multinomial(refcat = NA) dpars:", names(fixef_by_dpar(f)), "\n")
f <- frm(bf(Yc | trials(N) ~ x), family = brmsfamily("multinomial", refcat = "y2"), data = dm)
cat("brmsfamily('multinomial', refcat = 'y2') dpars:", names(fixef_by_dpar(f)), "\n")
