# Lane prilabel: the full variables() order of the brms fit of
# dev/prilabel-brms-fit.R, without the r_ and z_ entries.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
fit <- readRDS(file.path(wt, "dev/prilabel-log/brms-fit.rds"))
v <- variables(fit)
print(v[!grepl("^r_|^z_", v)])
