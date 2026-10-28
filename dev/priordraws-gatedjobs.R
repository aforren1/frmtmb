# The gated job list: every file of the ungated run that skipped
# something, to be run again with the gates open.
L <- readLines("dev/release/priordraws-ungated.log")
pkg <- NA_character_
keep <- character(0)
for (l in L) {
  if (grepl("^== .* files ==$", l)) pkg <- sub("^== ([^ ]+) .*", "\\1", l)
  if (grepl("^RESULT ", l) && grepl("skip=[1-9]", l)) {
    keep <- c(keep, paste(pkg, sub("^RESULT ([^ ]+) .*", "\\1", l)))
  }
}
j <- readLines("dev/priordraws-all.jobs")
key <- vapply(strsplit(j, " "), function(p) paste(p[1], basename(p[2])), "")
writeLines(j[key %in% keep], "dev/priordraws-gated.jobs")
cat(sum(key %in% keep), "of", length(keep), "\n")
