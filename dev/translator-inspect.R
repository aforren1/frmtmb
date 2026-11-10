# Print the L1 fields of measured pairs' records. Diagnostic only.
#
#   Rscript dev/translator-inspect.R <outdir> <id> [<id> ...]
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-translator-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
for (id in args[-1]) {
  r <- readRDS(file.path(args[1], paste0(id, ".rds")))
  cat("==", id, "\n")
  str(r[grep("^(l1|l4|frm_warn|verdict|reason|level|decided)", names(r))])
}
