# Reviewer, re-check: frmtmb.eam's test-family.R fails on the lane build.
# Which compat cells of the eam families are not "refused" for an
# addition term the family does not accept, and does subset() actually
# run with an eam family? Per arm:
#   Rscript dev/aterms2-rev-26-eamcompat.R <base|lane>
# Log: dev/aterms2-rev-log-26-eamcompat.txt
arm <- commandArgs(TRUE)[[1]]
libs <- c("C:/Users/adf44/source/r/wt-aterms2-lib",
          "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
ft <- frm_compat_features()
ats <- ft$name[ft$kind == "aterm"]
cat(arm, "aterm features:", ats, "\n")
for (fm in c("wiener", "gddm", "lba", "rdm", "wiener_gng")) {
  tb <- frm_compat(fm)
  other <- ifelse(tb$feature_a == fm, tb$feature_b, tb$feature_a)
  keep <- other %in% ats
  x <- tb[keep & tb$status != "refused", ]
  oth <- ifelse(x$feature_a == fm, x$feature_b, x$feature_a)
  cat(arm, fm, ": not refused:", paste0(oth, "=", x$status, collapse = " "),
      "\n")
}
