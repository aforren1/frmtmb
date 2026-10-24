# Does frmtmb.latent 0.6.2 (rellib-r11) load over the lane's core, and
# what does it then offer for mixture_mvn()?
.libPaths(c("C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
loadNamespace("frmtmb", lib.loc = "C:/Users/adf44/source/r/wt-mvnlatent-lib")
cat("frmtmb", format(packageVersion("frmtmb", lib.loc = "C:/Users/adf44/source/r/wt-mvnlatent-lib")),
    "us_sigma exported:", "us_sigma" %in% getNamespaceExports("frmtmb"), "\n")
r <- tryCatch({
  loadNamespace("frmtmb.latent", lib.loc = "C:/Users/adf44/source/r/rellib-r11")
  paste("loaded frmtmb.latent", format(packageVersion("frmtmb.latent",
        lib.loc = "C:/Users/adf44/source/r/rellib-r11")),
        "; mixture_mvn exported:",
        "mixture_mvn" %in% getNamespaceExports("frmtmb.latent"))
}, error = function(e) paste("LOAD ERROR:", conditionMessage(e)))
cat(r, "\n")
