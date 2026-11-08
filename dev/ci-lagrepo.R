# A file:// repository whose index is the installed library itself, so
# that dev/release/cran-lag.R can be seen to exit 0 when nothing lags
# (the real user library has never been fully current).
#
#   Rscript dev/ci-lagrepo.R <out dir> <lib1;lib2;...>
args <- commandArgs(trailingOnly = TRUE)
out <- args[[1L]]
libs <- strsplit(args[[2L]], ";", fixed = TRUE)[[1L]]
ip <- installed.packages(lib.loc = c(libs, .Library),
                         fields = c("Depends", "Imports", "LinkingTo",
                                    "Suggests"))
ip <- ip[!duplicated(ip[, "Package"]), , drop = FALSE]
cols <- c("Package", "Version", "Depends", "Imports", "LinkingTo",
          "Suggests")
dir.create(file.path(out, "src", "contrib"), recursive = TRUE,
           showWarnings = FALSE)
db <- ip[, cols, drop = FALSE]
db[is.na(db)] <- ""
write.dcf(db, file.path(out, "src", "contrib", "PACKAGES"))
cat("wrote", nrow(db), "entries\n")
