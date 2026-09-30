f <- "C:/Users/adf44/source/r/frmtmb-wt-sampfix/dev/sampfix-findings.md"
s <- paste(gsub("\r", "", readLines(f)), collapse = "\n")
old <- "  without it (`dev/sampfix-log/nits-laplace-noisna.txt`). The residual
  risk stays: a probe whose own code errors at both fills passes the
  call on, and only the per-draw check stands behind it, which cannot
  see a result that is non-finite at every draw."
new <- "  without it (`dev/sampfix-log/nits-laplace-noisna.txt`). The probe
  itself was the hole: an error at both fills let the call through. It
  now re-raises that error, which the call's first draw would raise
  anyway, so a probe that fails on its own account cannot wave a call
  through again."
stopifnot(grepl(old, s, fixed = TRUE))
writeLines(sub(old, new, s, fixed = TRUE), f)
cat("ok\n")
