# What the survey's rates would read with this lane's remeasured pairs,
# WITHOUT touching the record. Builds a scratch root that holds a copy of
# dev/equiv-data/ (with the adjudication rows this lane proposes applied
# to the copy) and an outdir of per-pair rows: each pair's row from the
# record (dev/equiv-data/results.tsv, measured verdict restored), replaced
# by this lane's measurement where it has one. Then runs
# dev/equiv-summary.R there without --record.
#
#   Rscript dev/translator-project.R <scratch root> <adjudication edits>
#                                    <outdir> [<outdir> ...]
#
# <adjudication edits> is a TSV with the adjudication columns; a row with
# verdict "DROP" removes that id's row from the copy, any other row
# replaces or adds it.
args <- commandArgs(trailingOnly = TRUE)
root <- args[1]
edits_p <- args[2]
outs <- args[-(1:2)]
rd <- function(p) read.delim(p, stringsAsFactors = FALSE,
                             na.strings = c("NA"), colClasses = "character")
unlink(root, recursive = TRUE)
dir.create(file.path(root, "dev", "equiv-data"), recursive = TRUE)
dir.create(file.path(root, "dev", "proj-out"))
file.copy(list.files("dev/equiv-data", full.names = TRUE),
          file.path(root, "dev", "equiv-data"))
file.copy("dev/equiv-summary.R", file.path(root, "dev"))

adj <- rd("dev/equiv-data/adjudication.tsv")
ed <- rd(edits_p)
adj <- adj[!adj$id %in% ed$id, , drop = FALSE]
adj <- rbind(adj, ed[ed$verdict != "DROP", names(adj), drop = FALSE])
utils::write.table(adj, file.path(root, "dev/equiv-data/adjudication.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)

res <- rd("dev/equiv-data/results.tsv")
added <- c("verdict_measured", "adjudicated", "nodefect", "measured",
           "top_level", "weak_evidence")
mine <- character()
for (i in seq_len(nrow(res))) {
  id <- res$id[i]
  src <- NULL
  for (o in outs) {
    p <- file.path(o, paste0(id, ".tsv"))
    if (file.exists(p)) src <- p
  }
  dest <- file.path(root, "dev", "proj-out", paste0(id, ".tsv"))
  if (!is.null(src)) {
    file.copy(src, dest)
    mine <- c(mine, id)
  } else {
    r <- res[i, , drop = FALSE]
    r$verdict <- r$verdict_measured
    # the level the record shows for an adjudicated row is the
    # adjudication's; the measured one is not kept, and the copied
    # adjudication puts it back
    r <- r[, setdiff(names(r), added), drop = FALSE]
    utils::write.table(r, dest, sep = "\t", row.names = FALSE, quote = TRUE,
                       qmethod = "double", na = "NA")
  }
}
cat("pairs from this lane's outdirs:", length(mine), "\n")
old <- setwd(root)
on.exit(setwd(old))
rc <- system2(file.path(R.home("bin"), "Rscript.exe"),
              c("dev/equiv-summary.R", "dev/proj-out"))
if (rc != 0) stop("equiv-summary.R failed")
