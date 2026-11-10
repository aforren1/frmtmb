# Reviewer scratch: does the new test's core_aterms survive extensions
# loaded earlier in the session?
.libPaths(c("C:/Users/adf44/source/r/wt-compat-lib",
            "C:/Users/adf44/source/r/rellib-r19",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (p in c("frmtmb.eam", "frmtmb.survival", "frmtmb.learn")) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
ns <- asNamespace("frmtmb")
str(utils::head(ns$frmtmb_compat_contrib$features, 3))
ft <- frm_compat_features()
contrib <- names(unlist(ns$frmtmb_compat_contrib$features))
print(contrib[grepl("dec|bhaz", contrib)])
core_aterms <- setdiff(ft$name[ft$kind == "aterm"], contrib)
print(core_aterms)
