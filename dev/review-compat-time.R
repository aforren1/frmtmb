# Reviewer scratch (round 2): frm_compat() cost, all extensions loaded.
#   Rscript dev/review-compat-time.R <lib;lib>
args <- commandArgs(trailingOnly = TRUE)
.libPaths(c(strsplit(args[1], ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
for (p in c("frmtmb.survival", "frmtmb.spatial", "frmtmb.latent",
            "frmtmb.sample", "frmtmb.learn", "frmtmb.eam",
            "frmtmb.coupling", "frmtmb.ode")) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}
cat("frmtmb", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
invisible(frm_compat("dec()", "thres()"))
tm <- function(expr, k) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(k)) eval(expr)
  (proc.time()[["elapsed"]] - t0) / k
}
cat("frm_compat(a, b), mean of 20:",
    round(tm(quote(frm_compat("dec()", "thres()")), 20), 3), "s\n")
cat("frmtmb_compat_rules_tbl(), mean of 20:",
    round(tm(quote(ns$frmtmb_compat_rules_tbl()), 20), 3), "s\n")
cat("frm_compat(), mean of 3:", round(tm(quote(frm_compat()), 3), 3), "s\n")
