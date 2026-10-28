# The rest of the gated tier: files of the ungated run that ran nothing
# at all (a skip at file level reports no skip count).
L <- readLines("dev/release/priordraws-ungated.log")
pkg <- NA_character_
keep <- character(0)
for (l in L) {
  if (grepl("^== .* files ==$", l)) pkg <- sub("^== ([^ ]+) .*", "\\1", l)
  if (grepl("^RESULT .* pass=0 fail=0 err=0 skip=0", l)) {
    keep <- c(keep, paste(pkg, sub("^RESULT ([^ ]+) .*", "\\1", l)))
  }
}
j <- readLines("dev/priordraws-all.jobs")
key <- vapply(strsplit(j, " "), function(p) paste(p[1], basename(p[2])), "")
con <- file("dev/priordraws-gated2.jobs", "wb")
writeLines(j[key %in% keep], con)
close(con)
cat(sum(key %in% keep), "of", length(keep), "\n")
