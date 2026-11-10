# Tabulate dev/compat-m-race.sh output: per child, its exit status, the
# extensions it loaded, and how many pairs have no live row.
#   Rscript dev/compat-m-race-read.R <outdir>
out <- commandArgs(trailingOnly = TRUE)[1]
errs <- sort(list.files(out, "[.]err$", full.names = TRUE))
for (e in errs) {
  tag <- sub("[.]err$", "", basename(e))
  el <- readLines(e, warn = FALSE)
  rc <- sub("^rc=", "", grep("^rc=", el, value = TRUE))
  tsv <- file.path(out, paste0(tag, ".tsv"))
  if (file.exists(tsv)) {
    d <- utils::read.delim(tsv, stringsAsFactors = FALSE)
    cat(sprintf("%-7s rc=%s pairs=%d na=%d ext=%d %s\n", tag, rc, nrow(d),
                sum(is.na(d$live_status)),
                length(strsplit(d$extensions[1], ",")[[1]]),
                paste(d$id[is.na(d$live_status)], collapse = " ")))
  } else {
    cat(sprintf("%-7s rc=%s NO TSV\n", tag, rc))
  }
  other <- setdiff(el, grep("^rc=", el, value = TRUE))
  if (length(other)) cat("   stderr:", utils::tail(other, 5), sep = "\n   ")
}
