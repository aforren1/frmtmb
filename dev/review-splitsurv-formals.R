arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-splitsurv-lib","C:/Users/adf44/source/r/rellib-r11","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "lane") { ns1 <- asNamespace("frmtmb"); ns2 <- asNamespace("frmtmb.survival") } else { ns1 <- asNamespace("frmtmb.spline"); ns2 <- asNamespace("frmtmb") }
out <- list()
for (f in c("frm_curve","frm_curve_deriv","frm_curve_feature")) out[[f]] <- formals(get(f, ns1))
for (f in c("royston_parmar","rp_floored")) out[[f]] <- formals(get(f, if (arm=="lane") ns2 else ns1))
for (f in c("cox","cox_baseline")) out[[f]] <- formals(get(f, ns2))
saveRDS(out, sprintf("dev/review-splitsurv-out/formals-%s.rds", arm))
cat(arm, "exports:", if (arm=="lane") paste(sort(getNamespaceExports("frmtmb.survival")), collapse=" ") else paste(sort(getNamespaceExports("frmtmb.spline")), collapse=" "), "\n")
