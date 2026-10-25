# The hunks of `git diff <file>` that belong to one item, for a file two
# items change: S1's hunks are the ones that add keep_missing() or
# require_fit_packages() (fix round 1), and the rest are the other
# item's. Writes the file's header and the chosen hunks to stdout.
#   Rscript dev/smallfix-split-hunks.R <file> <S1|other>
a <- commandArgs(trailingOnly = TRUE)
d <- system2("git", c("diff", "-U3", "--", a[1]), stdout = TRUE,
             stderr = FALSE)
at <- grep("^@@", d)
if (!length(at)) quit(save = "no")
head <- d[seq_len(at[1L] - 1L)]
ends <- c(at[-1L] - 1L, length(d))
keep <- character(0)
for (k in seq_along(at)) {
  h <- d[at[k]:ends[k]]
  s1 <- any(grepl("^[+].*(keep_missing|require_fit_packages)", h))
  if ((a[2] == "S1") == s1) keep <- c(keep, h)
}
if (length(keep)) writeLines(c(head, keep))
