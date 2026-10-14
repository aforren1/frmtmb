# Lane dblbar: frmtmb.sample's names on a (1 + f || g) fit. Data seed
# 20261013; sampler seed 1. Usage: Rscript dev/dblbar-smoke2.R base|lane
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
dd <- dblbar_data()
sink(file.path(wt, sprintf("dev/dblbar-log/smoke2-%s.txt", arm)), split = TRUE)
fit <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
dr <- frm_sample(fit, chains = 2, iter = 600, seed = 1, refresh = 0)
print(grep("^r_g.g01|^theta|^sd_|^cor_", variables(dr), value = TRUE))
s <- summary(dr)
print(class(s)); print(names(s))
print(head(rownames(s), 8))
ad <- as_draws_df(dr)
print(grep("^sd_|^cor_", names(ad), value = TRUE))
ps <- posterior_summary(dr)
print(grep("^sd_|^cor_", rownames(ps), value = TRUE))
print(VarCorr(dr)$g$sd)
print(head(names(suppressWarnings(posterior_samples(dr))), 20))
sink()
