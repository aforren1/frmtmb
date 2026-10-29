# Reviewer: does the exact gp() nugget explain the L1 gap of p094
# (lognormal x gp, -94.4 nats)? frmtmb's nugget is set to brms's size
# in this process only (assignInNamespace), and L1 is read again.
# Usage: Rscript gp-nugget.R <nugget> [key]
args <- commandArgs(trailingOnly = TRUE)
nug <- as.numeric(args[1])
key <- if (length(args) >= 2) args[2] else "lognormal x gp"
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/review-equiv-cache",
                                             winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
utils::assignInNamespace("gp_nugget", nug, "frmtmb")
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
smp <- read.delim("dev/equiv-data/sample.tsv")
row <- smp[paste(smp$feature_a, smp$feature_b, sep = " x ") == key, ]
rec <- E$eq_measure(key, row)
cat(sprintf("NUGGET %g %s: l1_const=%s extra=%s grad=%s l1_err=%s frm_err=%s\n",
            nug, key, signif(rec$l1_const, 4),
            paste(signif(rec$l1_extra_const, 4), collapse = ","),
            signif(rec$l1_grad, 3), substr(rec$l1_error %||% "", 1, 100),
            substr(rec$frm_error %||% "", 1, 100)))
