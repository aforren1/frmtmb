# Dump frmtmb's compatibility table and vocabulary on one build, with a
# chosen set of extensions loaded, for the before/after comparison of
# rule 4 of the move (the table reports the same rows with the
# extension loaded).
#
#   Rscript dev/mvnlatent-compat-dump.R <lib or "base"> <out.rds> [pkg ...]
a <- commandArgs(trailingOnly = TRUE)
base <- c("C:/Users/adf44/source/r/rellib-r11",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(a[1], "base")) base else c(a[1], base))
suppressMessages(library(frmtmb))
pk <- a[-(1:2)]
for (p in pk) suppressMessages(loadNamespace(p))
tb <- frm_compat()
ft <- frm_compat_features()
cat("frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
for (p in pk) {
  cat(p, as.character(packageVersion(p)), dirname(find.package(p)), "\n")
}
cat("rows", nrow(tb), "features", nrow(ft), "\n")
saveRDS(list(tb = as.data.frame(tb), ft = as.data.frame(ft)), a[2])
