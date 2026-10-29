# Measure one sampled pair against brms 2.23.0.
#
#   Rscript dev/equiv-measure.R <id> [outdir] [--force]
#
# Writes <outdir>/<id>.rds (the whole record) and <outdir>/<id>.tsv (one
# row). A pair whose .tsv exists is skipped, which is what makes a run
# resumable; --force remeasures it. The library is read from
# FRMTMB_EQUIV_LIB (default rellib-r15), so a later build is measured by
# pointing that variable at its library.
args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
outdir <- if (length(args) >= 2 && !startsWith(args[2], "--")) args[2] else
  "dev/equiv-out"
force <- "--force" %in% args
lib <- Sys.getenv("FRMTMB_EQUIV_LIB", "C:/Users/adf44/source/r/rellib-r15")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
tsv <- file.path(outdir, paste0(id, ".tsv"))
if (file.exists(tsv) && !force) {
  cat("skip", id, "(measured)\n")
  quit(save = "no")
}
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
# brms reads me() from the search path (helper-brms.R with_brms_me());
# frmtmb exports no me(), so this masks nothing of frmtmb's
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
# testthat runs the helpers in an environment whose parent is the
# package namespace, and the translator calls package internals
# (cs_entry()); the instrument is built the same way
eq_env <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = eq_env)
}

smp <- read.delim("dev/equiv-data/sample.tsv")
row <- smp[smp$id == id, ]
if (nrow(row) != 1L) stop("no sampled pair ", id)
key <- paste(row$feature_a, row$feature_b, sep = " x ")
cat("pair", id, ":", key, "[", row$status, "]\n")
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")

rec <- eq_env$eq_measure(key, row)
eq_flatten <- eq_env$eq_flatten
rec$seconds <- round(proc.time()[["elapsed"]] - t0, 1)
rec$frmtmb_version <- as.character(packageVersion("frmtmb"))
saveRDS(rec, file.path(outdir, paste0(id, ".rds")))
flat <- eq_flatten(rec, row)
utils::write.table(flat, tsv, sep = "\t", row.names = FALSE,
                   quote = TRUE, qmethod = "double")
cat("verdict", id, ":", flat$verdict, "level", flat$level, "-",
    flat$reason, "\n")
