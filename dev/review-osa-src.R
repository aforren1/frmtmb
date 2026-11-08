# Review lane osa: dump RTMB/TMB oneStepPredict code paths the lane relies on.
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
library(frmtmb)
cat("lib:", "frmtmb", as.character(packageVersion("frmtmb")),
    dirname(dirname(system.file(package = "frmtmb"))), "\n")
for (p in c("RTMB", "TMB", "RTMBdist")) {
  cat(p, as.character(packageVersion(p)), system.file(package = p), "\n")
}
cat("==== RTMB::oneStepPredict\n")
print(RTMB::oneStepPredict)
cat("==== TMB::oneStepPredict\n")
print(TMB::oneStepPredict)
ns <- asNamespace("RTMB")
for (f in grep("OSA|osa", ls(ns, all.names = TRUE), value = TRUE)) {
  cat("==== RTMB:::", f, "\n")
  print(get(f, ns))
}
