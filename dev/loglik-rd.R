# Lane loglik: render the help sections this lane wrote, to check the
# math and the markup by reading the output (lane-rules: verify an Rd by
# rendering it). Output: dev/loglik-rd.txt
files <- c("man/frmtmb-autocor.Rd", "extensions/frmtmb.spatial/man/sar.Rd",
           "man/frmtmb_register_autocor.Rd",
           "extensions/frmtmb.sample/man/sample-log_lik.Rd",
           "man/frmtmb-sampling-api.Rd")
pat <- "pointwise log-density|The Jacobian|Residual covariance terms|autocor_row_loglik"
for (f in files) {
  o <- capture.output(tools::Rd2txt(f, options = list(underline_titles = FALSE)))
  i <- grep(pat, o)
  cat("=====", f, "\n")
  keep <- sort(unique(unlist(lapply(i, function(k) k:min(length(o), k + 12)))))
  cat(o[keep], sep = "\n")
}
