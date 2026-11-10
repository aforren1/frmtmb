# The live frm_compat() status of sampled pairs, read from the same table
# the frozen statuses came from: frmtmb with every extension that
# dev/equiv-dump-compat.R loads (its `ext` vector is read from that file,
# so the two lists cannot drift). Run in a process of its own, so that
# the extensions never load into the process that measures a pair.
#
#   Rscript dev/equiv-live-row.R <id,id,...|all> <out.tsv>
#
# The library is read from FRMTMB_EQUIV_LIB, a ";"-separated list.
args <- commandArgs(trailingOnly = TRUE)
lib <- Sys.getenv("FRMTMB_EQUIV_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ext <- NULL
for (e in as.list(parse("dev/equiv-dump-compat.R"))) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
        identical(e[[2L]], as.name("ext"))) {
    ext <- eval(e[[3L]], baseenv())
  }
}
if (!length(ext)) stop("no ext vector in dev/equiv-dump-compat.R")
loaded <- character()
for (p in ext) {
  # An absent extension is skipped, as the dump skips it; one that is
  # installed but fails to load stops, because a quiet skip would read a
  # different table from the one the frozen statuses came from
  if (!nzchar(system.file(package = p))) next
  suppressPackageStartupMessages(library(p, character.only = TRUE))
  loaded <- c(loaded, p)
}
smp <- utils::read.delim("dev/equiv-data/sample.tsv",
                         stringsAsFactors = FALSE)
if (!identical(args[1], "all")) {
  smp <- smp[smp$id %in% strsplit(args[1], ",", fixed = TRUE)[[1]], ]
}
cp <- as.data.frame(frm_compat())
smp$live_status <- vapply(seq_len(nrow(smp)), function(i) {
  a <- smp$feature_a[i]
  b <- smp$feature_b[i]
  s <- cp$status[(cp$feature_a == a & cp$feature_b == b) |
                   (cp$feature_a == b & cp$feature_b == a)]
  if (length(s) == 1L) s else NA_character_
}, "")
smp$extensions <- paste(loaded, collapse = ",")
smp$frmtmb_lib <- find.package("frmtmb")
# Written beside the target and renamed, so a child stopped mid-write
# leaves no file that the caller reads as complete
part <- paste0(args[2], ".part")
utils::write.table(smp, part, sep = "\t", row.names = FALSE,
                   quote = TRUE, qmethod = "double")
if (!file.rename(part, args[2])) stop("cannot rename ", part)
