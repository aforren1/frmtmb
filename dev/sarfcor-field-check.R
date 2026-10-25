# Which structures frmtmb_covstruct_info() calls fields, against the
# list frmtmb.sample's draws_cov_fast() kept through 0.74.0; and the
# full_rank property of every structure.
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb.spatial))
nm <- frmtmb:::cs_names()
info <- lapply(nm, frmtmb_covstruct_info)
field <- nm[vapply(info, `[[`, NA, "field")]
cat("structures:", length(nm), "\n")
cat("fields by the registry:", sort(field), "\n")
cat("old list:", sort(c("smooth", "gp", "hsgp", "car", "spde")), "\n")
cat("identical:", identical(sort(field), sort(c("smooth", "gp", "hsgp", "car", "spde"))), "\n")
fr <- vapply(info, function(i) as.character(i$full_rank), "")
print(table(fr))
cat("full_rank TRUE:", sort(nm[fr == "TRUE"]), "\n")
