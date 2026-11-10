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
# sample.tsv's status is the pair's STRATUM, frozen at the draw; it sets
# the weights and the per-stratum rates and stays in the `status` column.
# The verdict is scored against the build's own row: scored against the
# frozen status, a pair whose row was corrected after the draw read
# "compat row wrong" on every later build (dev/compat-findings.md). The
# row is read with the extensions the frozen statuses were read with,
# in a process of its own, so that none of them loads into this one
# (dev/equiv-live-row.R).
Sys.setenv(FRMTMB_EQUIV_LIB = lib)
# A child that died with no word (2 of 42 in a parallel batch,
# dev/compat-findings.md section 11) was once reported as a missing row.
# Its exit status and output are kept apart from the row's absence, and
# a failed child is retried once: the row is deterministic, the process
# is not. The script's full path puts the worktree on the child's command
# line, so a process listing shows whose child it is.
live_row <- function(id) {
  tag <- paste0("live-", id, "-", Sys.getpid())
  why <- character()
  for (attempt in 1:2) {
    out <- file.path(tempdir(), paste0(tag, "-", attempt, ".tsv"))
    log <- file.path(tempdir(), paste0(tag, "-", attempt, ".log"))
    rc <- system2(file.path(R.home("bin"), "Rscript"),
                  c(shQuote(normalizePath("dev/equiv-live-row.R",
                                          winslash = "/")),
                    id, shQuote(out)), stdout = log, stderr = log)
    if (rc == 0L && file.exists(out)) {
      return(read.delim(out, stringsAsFactors = FALSE)$live_status)
    }
    tail <- if (file.exists(log)) {
      utils::tail(readLines(log, warn = FALSE), 8)
    }
    why <- sprintf("exit status %s%s; output: %s", rc,
                   if (rc == 0L) ", but no row file" else "",
                   if (length(tail)) paste(tail, collapse = " | ") else
                     "(none)")
    cat("live-row child, attempt", attempt, "failed:", why, "\n")
  }
  stop("dev/equiv-live-row.R failed twice for ", id, ": ", why)
}
live <- live_row(id)
if (length(live) != 1L || is.na(live)) {
  stop("no frm_compat() row for ", key)
}
cat("pair", id, ":", key, "[ stratum", row$status, "| live", live, "]\n")
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")

row_live <- row
row_live$status <- live
rec <- eq_env$eq_measure(key, row_live)
eq_flatten <- eq_env$eq_flatten
rec$seconds <- round(proc.time()[["elapsed"]] - t0, 1)
rec$frmtmb_version <- as.character(packageVersion("frmtmb"))
rec$stratum_status <- row$status
rec$live_status <- live
saveRDS(rec, file.path(outdir, paste0(id, ".rds")))
flat <- eq_flatten(rec, row)
flat$stratum_status <- row$status
flat$live_status <- live
utils::write.table(flat, tsv, sep = "\t", row.names = FALSE,
                   quote = TRUE, qmethod = "double")
cat("verdict", id, ":", flat$verdict, "level", flat$level, "-",
    flat$reason, "\n")
