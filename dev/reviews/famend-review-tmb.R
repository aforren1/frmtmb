for (fn in list(TMB::MakeADFun, RTMB::MakeADFun)) {
  f <- capture.output(print(fn))
  i <- grep("last.par.best|value.best", f)
  cat("matches", length(i), "\n")
  if (length(i)) print(f[sort(unique(unlist(lapply(i, function(k) max(1, k - 5):(k + 1)))))])
}
