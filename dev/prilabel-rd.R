# Lane prilabel: render the Rd pages this lane changed and show the new
# passages, so a stray % or markup is seen in the output, not the source.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
pat <- paste("Group-level standard deviations|One u_j is also shared",
             "named as brms names it|row per coefficient|take over",
             "parameters it reaches", sep = "|")
for (f in c("man/set_prior.Rd", "man/frmtmb-student-re.Rd",
            "man/frmtmb-multimembership.Rd", "man/default_prior.Rd",
            "extensions/frmtmb.sample/man/frm_sample.Rd")) {
  o <- capture.output(tools::Rd2txt(file.path(wt, f)))
  cat("==", f, "\n")
  for (k in grep(pat, o)) cat(o[k:min(length(o), k + 8L)], sep = "\n")
}
