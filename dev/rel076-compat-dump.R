# Dump frm_compat() over every pair and frm_compat_features(), with the
# packages named loaded, from the library stack named. Adapted from
# dev/spatial-compat-dump.R.
#   Rscript dev/rel076-compat-dump.R <out.rds> <lib;lib> <pkg,pkg,...>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(a[2], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
pk <- strsplit(a[3], ",", fixed = TRUE)[[1]]
for (p in pk) suppressMessages(loadNamespace(p))
tab <- as.data.frame(frm_compat())
ft <- frm_compat_features()
libs <- vapply(c("frmtmb", pk), function(p) {
  paste(p, format(packageVersion(p)), dirname(find.package(p)))
}, "")
saveRDS(list(tab = tab, features = ft, libs = libs), a[1])
cat(paste0("lib: ", libs, "\n"), sep = "")
cat("pairs", nrow(tab), "features", nrow(ft), "\n")
