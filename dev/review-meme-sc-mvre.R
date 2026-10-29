# Review of lane meme, short check: the mvre cor group row on base.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r15", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-meme"
source(file.path(wt, "dev/review-meme-shapes.R")); source(file.path(wt, "dev/review-meme-sc-shapes.R"))
d <- sc_data(); sh <- sc_shapes(bf, mvbf, set_rescor, cumulative, gp)
t <- as.data.frame(default_prior(sh$mvre$f, data = d, route = "sample"))
cat(arm, "\n"); print(t[t$class %in% c("cor", "sd"), c("prior", "class", "coef", "group", "resp")])
