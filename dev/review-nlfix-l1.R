# One spec of dev/review-nlfix-specs.R through the equiv instrument (unchanged).
args <- commandArgs(trailingOnly = TRUE)
nm <- args[1]
outdir <- "dev/review-nlfix-l1-out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-nlfix-lib2", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win", NOT_CRAN = "true")
suppressPackageStartupMessages({ library(frmtmb); library(testthat); library(frmtmb.survival); library(frmtmb.eam) })
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) sys.source(f, envir = E)
specs <- new.env(parent = E)
sys.source("dev/review-nlfix-specs.R", envir = specs)
sp <- specs$l1_specs[[nm]]
if (is.null(sp)) stop("no spec ", nm)
E$eq_specs[["variant"]] <- sp
rec <- E$eq_measure("variant", list(status = "untested"))
g <- function(x) if (is.null(rec[[x]])) NA else paste(format(rec[[x]], digits = 4), collapse = " ")
cat(sprintf("RESULT %s verdict=%s l1_ok=%s const=%s grad=%s extra=%s | frm_err=%s | brms_err=%s | l1_err=%s | reason=%s\n", nm,
            g("verdict"), g("l1_ok"), g("l1_const"), g("l1_grad"), g("l1_extra_const"),
            substr(g("frm_error"), 1, 200), substr(g("brms_error"), 1, 200), substr(g("l1_error"), 1, 200), substr(g("reason"), 1, 200)),
    file = file.path(outdir, paste0(nm, ".txt")))
str(rec[["l1"]], max.level = 2, nchar.max = 300)
cat("L4 record:\n"); str(rec[c("l4", "l4_ok", "l4_z", "l4_route", "l2", "frm_logLik")], max.level = 2, nchar.max = 200)
