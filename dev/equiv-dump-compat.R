# Dump frm_compat() and its features on the build under test, so the
# sampling frame is a file and not a recollection. Run once per build:
#   Rscript dev/equiv-dump-compat.R [rellib] [outdir]
# Loading every extension matters: the pair table grows with each one.
args <- commandArgs(trailingOnly = TRUE)
rellib <- if (length(args) >= 1) args[1] else
  "C:/Users/adf44/source/r/rellib-r15"
# the 3 MB pair table is machine output (gitignored *-out/); the
# tracked record of the frame is sample.tsv and strata.tsv
out <- if (length(args) >= 2) args[2] else "dev/equiv-frame-out"
.libPaths(c(rellib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressPackageStartupMessages(library(frmtmb))
t0 <- frm_compat()
ext <- c("frmtmb.survival", "frmtmb.spatial", "frmtmb.latent",
         "frmtmb.sample", "frmtmb.learn", "frmtmb.eam",
         "frmtmb.coupling", "frmtmb.ode")
for (p in ext) {
  if (requireNamespace(p, quietly = TRUE)) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}
cp <- frm_compat()
ft <- frm_compat_features()
cat("core-only pairs:", nrow(t0), "\n")
cat("all-loaded pairs:", nrow(cp), "\n")
print(table(cp$status))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
cp$in_core_only <- paste(cp$feature_a, cp$feature_b) %in%
  c(paste(t0$feature_a, t0$feature_b), paste(t0$feature_b, t0$feature_a))
utils::write.table(as.data.frame(cp), file.path(out, "compat.tsv"),
                   sep = "\t", row.names = FALSE, quote = TRUE,
                   qmethod = "double")
utils::write.table(as.data.frame(ft), file.path(out, "features.tsv"),
                   sep = "\t", row.names = FALSE, quote = TRUE,
                   qmethod = "double")
print(table(ft$kind))
