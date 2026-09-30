# Lane ceplot punch 2: the multivariate "old_levels" item leaves the
# divergence list of dev/ceplot-findings.md; it was a defect, now fixed.
f <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-findings.md"
s <- paste(readLines(f), collapse = "\n")
old <- "- With `sample_new_levels = \"old_levels\"` on a multivariate `fitted()`,
  each response chooses its seen levels on its own; brms chooses once
  per call, so a grouping factor two responses share reads one group in
  brms and possibly two here. (Punch round 1 corrected this item: it
  said \"differs only for a block shared by two responses\", and the
  same defect within ONE response, per block rather than per factor,
  was the review's B1, now fixed.)
"
stopifnot(lengths(regmatches(s, gregexpr(old, s, fixed = TRUE))) == 1L)
s <- sub(old, "", s, fixed = TRUE)
writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], f)
cat("edited\n")
