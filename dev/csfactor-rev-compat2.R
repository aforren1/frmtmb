# Does the edited compat registry note reach frm_compat()'s table?
#   Rscript dev/csfactor-rev-compat2.R
.libPaths(c("C:/Users/adf44/source/r/wt-csfactor-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
tb <- as.data.frame(frm_compat(feature_a = "cs_pred()"))
cat("rows for cs_pred():", nrow(tb), "\n")
print(table(tb$status))
w <- grepl("model.matrix", tb$note, fixed = TRUE)
cat("rows carrying the new note:", sum(w), "\n")
if (any(w)) {
  print(tb[w, c("feature_a", "kind_a", "feature_b", "kind_b", "status")])
  cat("\nnote text:\n")
  cat(strwrap(unique(tb$note[w]), 74), sep = "\n")
  cat("\n")
}
cat("\nordinal_cs group members among feature_b:\n")
print(sort(unique(tb$feature_b[tb$status == "works"]))[1:min(40L,
      length(unique(tb$feature_b[tb$status == "works"])))])
cat("\nrows whose note mentions 'both sides':",
    sum(grepl("both sides", tb$note, fixed = TRUE)), "\n")
cat("done\n")
