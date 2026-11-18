# Lane squant, fix round 2: quantile_elf's compat status for every
# covariance structure and smooth/gp special, core only.
# Run: Rscript dev/squant-compat4.R > dev/squant-out/compat4-<arm>.txt
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- frm_compat("quantile_elf")
o <- ifelse(q$feature_a == "quantile_elf", q$feature_b, q$feature_a)
k <- ifelse(q$feature_a == "quantile_elf", q$kind_b, q$kind_a)
sel <- k == "covstruct" | o %in% c("s()", "t2()", "gp_pred()", "mi()")
for (i in which(sel)) cat(sprintf("%-10s %-11s %s\n", o[i], q$status[i],
                                  substr(q$note[i], 1, 90)))
print(table(q$status))
