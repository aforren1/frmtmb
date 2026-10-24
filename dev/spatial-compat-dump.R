# Dump the resolved compatibility table, frm_compat() over every pair,
# and the vocabulary, with the packages in SPATIAL_PKGS loaded.
# Usage: Rscript dev/spatial-compat-dump.R <out.rds>
libs <- strsplit(Sys.getenv("SPATIAL_LIBS"), ";", fixed = TRUE)[[1]]
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
pk <- strsplit(Sys.getenv("SPATIAL_PKGS"), ",", fixed = TRUE)[[1]]
for (p in pk) suppressMessages(loadNamespace(p))
tab <- frm_compat()
ft <- frm_compat_features()
key <- function(a, b) ifelse(a < b, paste(a, b, sep = " x "),
                             paste(b, a, sep = " x "))
tab$pair <- key(tab$feature_a, tab$feature_b)
out <- commandArgs(trailingOnly = TRUE)[1]
saveRDS(list(tab = as.data.frame(tab), features = ft,
             lib = dirname(find.package("frmtmb")),
             version = as.character(packageVersion("frmtmb")),
             pkgs = pk), out)
cat("lib", dirname(find.package("frmtmb")), "pairs", nrow(tab),
    "features", nrow(ft), "\n")
