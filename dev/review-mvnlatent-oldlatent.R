# Does the released frmtmb.latent (rellib-r11) load over the lane's
# core? And sample 0.21.0's mixture_mvn path over the lane's core.
lane <- "C:/Users/adf44/source/r/wt-mvnlatent-lib"
rel <- "C:/Users/adf44/source/r/rellib-r11"
.libPaths(c(rel, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
loadNamespace("frmtmb", lib.loc = lane)
cat("frmtmb", format(packageVersion("frmtmb", lib.loc = lane)),
    dirname(getNamespaceInfo("frmtmb", "path")), "\n")
cat("latent candidate", format(packageVersion("frmtmb.latent", lib.loc = rel)),
    "\n")
r <- tryCatch({
  loadNamespace("frmtmb.latent", lib.loc = rel)
  "loaded"
}, error = function(e) conditionMessage(e))
cat("load:", r, "\n")
cat("hmm usable:", isNamespaceLoaded("frmtmb.latent"), "\n")
