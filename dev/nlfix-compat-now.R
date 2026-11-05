# Rows of the nl structure in frm_compat(), base build and lane build.
lib <- Sys.getenv("NLFIX_LIB",
                  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
cp <- frm_compat("nl")
cp <- cp[cp$feature_a == "nl" | cp$feature_b == "nl", ]
other <- ifelse(cp$feature_a == "nl", cp$feature_b, cp$feature_a)
kind <- ifelse(cp$feature_a == "nl", cp$kind_b, cp$kind_a)
tab <- data.frame(other, kind, status = cp$status,
                  note = substr(cp$note, 1, 70))
tab <- tab[order(tab$kind, tab$other), ]
print(tab[tab$kind != "family" & tab$kind != "method", ], row.names = FALSE)
print(table(tab$kind, tab$status))
print(tab[tab$kind == "family", c("other", "status")], row.names = FALSE)
