# Measure one sampled pair as dev/equiv-measure.R does, but score it
# against the LIVE frm_compat() status of the build under test instead of
# the status frozen into dev/equiv-data/sample.tsv when the sample was
# drawn. The flow (dev/equiv-flow.R) reads `row$status` to tell "compat
# row wrong" from "refusal correct" and "equivalent", so with the frozen
# status a pair whose row was corrected after the draw keeps "compat row
# wrong" on every build. Nothing else differs from equiv-measure.R.
#
#   Rscript dev/compat-equiv-live.R <id> <outdir>
#
# The library is read from FRMTMB_EQUIV_LIB, a ";"-separated list.
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
outdir <- args[2]
lib <- Sys.getenv("FRMTMB_EQUIV_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!nzchar(Sys.getenv("FRMTMB_STAN_CACHE"))) {
  Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache",
                                               winslash = "/"))
}
if (!nzchar(Sys.getenv("R_MAKEVARS_USER"))) {
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
}
t0 <- proc.time()[["elapsed"]]
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
  requireNamespace("brms")
  requireNamespace("rstan")
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
eq_env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = eq_env)
}
smp <- read.delim("dev/equiv-data/sample.tsv")
row <- smp[smp$id == id, ]
if (nrow(row) != 1L) stop("no sampled pair ", id)
live <- frm_compat(row$feature_a, row$feature_b)
live <- live$status[(live$feature_a == row$feature_a &
                       live$feature_b == row$feature_b) |
                      (live$feature_a == row$feature_b &
                         live$feature_b == row$feature_a)]
if (length(live) != 1L) stop("no live compat row for ", id)
cat("pair", id, ":", row$feature_a, "x", row$feature_b,
    "[ frozen", row$status, "| live", live, "]\n")
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
row$status_frozen <- row$status
row$status <- live
key <- paste(row$feature_a, row$feature_b, sep = " x ")
rec <- eq_env$eq_measure(key, row)
rec$seconds <- round(proc.time()[["elapsed"]] - t0, 1)
rec$frmtmb_version <- as.character(packageVersion("frmtmb"))
saveRDS(rec, file.path(outdir, paste0(id, ".rds")))
flat <- eq_env$eq_flatten(rec, row)
flat$status_frozen <- row$status_frozen
utils::write.table(flat, file.path(outdir, paste0(id, ".tsv")), sep = "\t",
                   row.names = FALSE, quote = TRUE, qmethod = "double")
cat("verdict", id, ":", flat$verdict, "level", flat$level, "-",
    flat$reason, "\n")
cat("frm_error:", if (is.null(rec$frm_error)) "<none>" else rec$frm_error,
    "\n")
cat("brms_error:", if (is.null(rec$brms_error)) "<none>" else
  rec$brms_error, "\n")
