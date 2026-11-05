# L1 runs of the lane's nl models through the equiv instrument
# (dev/equiv-lib.R eq_measure), one spec per R process:
#   Rscript dev/nlfix-l1.R <spec-name> [outdir]
# Writes <outdir>/<spec-name>.txt with the record's verdict and level
# numbers. The instrument is reused unchanged; the specs are below.
args <- commandArgs(trailingOnly = TRUE)
nm <- args[1]
outdir <- if (length(args) >= 2) args[2] else "dev/nlfix-l1-out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(strsplit(Sys.getenv("FRMTMB_EQUIV_LIB",
  "C:/Users/adf44/source/r/wt-nlfix-lib2;C:/Users/adf44/source/r/rellib-r16"),
  ";", fixed = TRUE)[[1]], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
specs <- new.env(parent = E)
sys.source("dev/nlfix-l1-specs.R", envir = specs)
sp <- specs$l1_specs[[nm]]
if (is.null(sp)) stop("no spec ", nm)
E$eq_specs[["variant"]] <- sp
rec <- E$eq_measure("variant", list(status = "untested"))
keep <- c("verdict", "reason", "method_reached", "frm_error", "brms_error", "l1_const",
          "l1_grad", "l1_points", "l1_ok", "l1_extra_const", "l1_error", "l4", "l2",
          "l4_error")
sink(file.path(outdir, paste0(nm, ".txt")))
cat("spec:", nm, "\nlib:", find.package("frmtmb"), "\n")
cat("bf:", sp()$bf_frm %||% sp()$bf, "\nfamily:", sp()$family %||% "gaussian()",
    "\n")
str(rec[intersect(names(rec), keep)], max.level = 3, nchar.max = 400)
cat("l1:\n")
str(rec[["l1"]], max.level = 2, nchar.max = 300)
sink()
