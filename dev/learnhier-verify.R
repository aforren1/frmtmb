# Lane `learnhier`: is every result file a COMPLETE record?
#
#   Rscript dev/learnhier-verify.R <dir> [rec|imp] ["YYYY-MM-DD HH:MM"]
#
# QUARANTINE BY TIME, NOT BY COUNT, AND DO IT FIRST. This machine's R
# library has been destroyed five times (`dev/machine-library.md`), and
# the fifth took the round's shared reference build with it at 17:24 on
# 2026-09-10. A fit whose namespace vanishes mid-run can still reach
# `saveRDS()` and write a record that passes every structural check
# below, so completeness is necessary and NOT sufficient: a run that
# overlaps a library loss has to be discarded on its timestamp. The
# third argument is that cutoff; every file modified after it is listed
# and the script exits non-zero.
#
# The timestamp half needs no R and should be run WITHOUT it while a
# library is being restored, because loading from a half-restored
# library is the thing that must not happen:
#
#   find <dir> -name '*.rds' -newermt "2026-09-10 17:24:00"
#
# A process killed part way through `saveRDS()` can leave a file that
# still opens. What it cannot leave is the full field set, because the
# record is written in ONE call at the end of the run and every field is
# put there by the line before it. So the test is structural rather than
# a checksum, and it names the missing field rather than saying "bad".
#
# It also reads each file twice and compares, which catches a truncated
# stream that decompresses to a short but parseable object.
#
# AND it checks that the field set is CONSTANT WITHIN A DESIGN, which is
# the check that earned its keep: two optional fields were added to the
# record while the bandit arm was already running, so 90 of the 120
# bandit files predate them. That is a schema drift and not a truncation
# (`oracle_beta` and `oracle_betad` are the oracle's INPUTS, which only
# the rlddm summary reads, and `oracle` itself is in every file), but a
# verifier that could not tell the two apart would be no use.
source("dev/learnhier-env.R")
a <- commandArgs(trailingOnly = TRUE)
DIR <- a[[1L]]
KIND <- if (length(a) >= 2L) a[[2L]] else "rec"
CUTOFF <- if (length(a) >= 3L) as.POSIXct(a[[3L]]) else NULL

# Written by the line before the save, so absence means the save did not
# complete.
need_rec <- c("design", "seed", "ns", "nt", "rows", "sim_s",
              "truth_drawn", "truth_fitted", "fit_s", "ok", "total_s",
              "peak_r_mb")
need_rec_ok <- c("conv", "max_grad", "pdHess", "n_bad_se", "logLik",
                 "n_par", "fixed", "vc", "oracle")
need_imp <- c("design", "seed", "draws", "ns", "nt", "lap_s", "lap_ok",
              "imp_s", "imp_ok", "outcome", "warns")
need_imp_ok <- c("moves", "rounds", "capped", "spread", "stalled",
                 "ess_min", "ess_median", "mcse", "imp", "imp_ll")
# Added to the record after part of the bandit arm had run. Reported,
# never required.
optional <- c("oracle_beta", "oracle_betad")

fs <- list.files(DIR, pattern = "[.]rds$", full.names = TRUE)
cat("files ", length(fs), " in ", DIR, "\n", sep = "")

# The timestamp gate, before anything is read.
mt <- file.info(fs)$mtime
cat("modified from ", format(min(mt)), " to ", format(max(mt)), "\n",
    sep = "")
suspect <- character(0)
if (!is.null(CUTOFF)) {
  suspect <- basename(fs[mt > CUTOFF])
  cat("cutoff ", format(CUTOFF), ": ", length(suspect),
      " file(s) written after it\n", sep = "")
  if (length(suspect)) {
    cat(paste0("  SUSPECT ", suspect, collapse = "\n"), "\n")
  }
}

bad <- character(0)
rows <- list()
for (f in fs) {
  r <- try(readRDS(f), silent = TRUE)
  if (inherits(r, "try-error")) {
    bad <- c(bad, paste0(basename(f), ": unreadable"))
    next
  }
  r2 <- try(readRDS(f), silent = TRUE)
  if (inherits(r2, "try-error") || !identical(names(r), names(r2))) {
    bad <- c(bad, paste0(basename(f), ": two reads disagree"))
    next
  }
  if (KIND == "imp") {
    miss <- setdiff(need_imp, names(r))
    if (isTRUE(r$imp_ok)) miss <- c(miss, setdiff(need_imp_ok, names(r)))
    rows[[length(rows) + 1L]] <- data.frame(
      file = basename(f), design = paste0("imp-", r$draws),
      seed = r$seed, n_a = length(r$moves), n_b = length(r$imp),
      shape = "", opt = 0L,
      finite = all(is.finite(r$imp)), stringsAsFactors = FALSE)
  } else {
    miss <- setdiff(need_rec, names(r))
    if (isTRUE(r$ok)) miss <- c(miss, setdiff(need_rec_ok, names(r)))
    if (grepl("^rlddm", basename(f)) && isTRUE(r$ok) &&
        !("ndt" %in% names(r))) {
      miss <- c(miss, "ndt")
    }
    rows[[length(rows) + 1L]] <- data.frame(
      file = basename(f), design = r$design, seed = r$seed,
      n_a = length(r$fixed), n_b = length(r$vc),
      shape = paste0(r$ns, "x", r$nt),
      opt = sum(optional %in% names(r)),
      finite = if (isTRUE(r$ok)) {
        all(is.finite(r$fixed[grepl("[.]est$", names(r$fixed))]))
      } else NA, stringsAsFactors = FALSE)
  }
  if (length(miss)) {
    bad <- c(bad, paste0(basename(f), ": missing ",
                         paste(miss, collapse = ",")))
  }
}
tab <- do.call(rbind, rows)
cat("\nper design. A field count that VARIES inside a design is the",
    " shape a partial write would take.\n", sep = "")
drift <- 0L
for (dz in sort(unique(tab$design))) {
  s <- tab[tab$design == dz, , drop = FALSE]
  va <- sort(unique(s$n_a))
  vb <- sort(unique(s$n_b))
  if (length(va) > 1L || length(vb) > 1L) drift <- drift + 1L
  cat("  ", dz, ": ", nrow(s), " files, ",
      length(unique(s$seed)), " distinct seeds, fields ",
      paste(va, collapse = "/"), " and ", paste(vb, collapse = "/"),
      if (nzchar(s$shape[[1L]])) {
        paste0(", design ", paste(unique(s$shape), collapse = "/"))
      } else "",
      ", optional fields present on ", sum(s$opt > 0L), " of ", nrow(s),
      ", all estimates finite: ", all(s$finite), "\n", sep = "")
}
if (length(bad)) {
  cat("\nINCOMPLETE FILES ", length(bad), ":\n", sep = "")
  cat(paste0("  ", utils::head(bad, 20L), collapse = "\n"), "\n")
  quit(save = "no", status = 1L)
}
cat("\nall ", nrow(tab), " files are complete records; ", drift,
    " designs show a varying field count\n", sep = "")
if (length(suspect)) {
  cat("but ", length(suspect),
      " of them were written after the cutoff and are NOT usable\n",
      sep = "")
  quit(save = "no", status = 1L)
}
