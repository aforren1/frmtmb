# Reviewer, priority 3: which methods did insight, marginaleffects and
# emmeans have to write FOR brmsfit, and does frmtmb_fit have the same
# one? A method those packages wrote for brmsfit is usually there
# BECAUSE brms's shapes are unusual: insight::get_residuals.brmsfit
# exists only to pull the Estimate column out of brms's matrix. Item
# 2.6f gave frmtmb those shapes, so every such method is a place the
# item can have left a hole.
#
#   Rscript dev/shapes-rev-methodgap.R

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))

own <- unique(c(grep("[.]frmtmb_fit$", ls(asNamespace("frmtmb")),
                     value = TRUE),
                as.character(utils::methods(class = "frmtmb_fit"))))
own <- sub("[.]frmtmb_fit$", "", own)

for (pk in c("insight", "marginaleffects", "emmeans")) {
  suppressPackageStartupMessages(
    library(pk, character.only = TRUE))
  nm <- ls(asNamespace(pk), all.names = TRUE)
  brmsm <- sub("[.]brmsfit$", "",
               grep("[.]brmsfit$", nm, value = TRUE))
  miss <- setdiff(brmsm, own)
  cat("==== ", pk, ": ", length(brmsm), " methods for brmsfit; ",
      length(brmsm) - length(miss), " of those generics also have a ",
      "frmtmb_fit method\n", sep = "")
  cat("  no frmtmb_fit method for: ", paste(sort(miss), collapse = ", "),
      "\n\n")
}
