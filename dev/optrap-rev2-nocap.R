# Reviewer of lane optrap, re-check: predictor_resid_sd() without its cap
# of 5e6 design cells (the least squares runs on every design), sourced
# after library(frmtmb); this process only.
local({
  ns <- asNamespace("frmtmb")
  f <- get("predictor_resid_sd", ns)
  b <- deparse(body(f))
  i <- grep("<= 5e+06", b, fixed = TRUE)
  stopifnot(length(i) == 1L)
  b[i] <- sub("<= 5e+06", "<= Inf", b[i], fixed = TRUE)
  body(f) <- parse(text = b)[[1]]
  utils::assignInNamespace("predictor_resid_sd", f, "frmtmb")
  cat("predictor_resid_sd patched: no cap\n")
})
