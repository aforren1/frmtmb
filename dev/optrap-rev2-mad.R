# Reviewer of lane optrap, re-check: a candidate change to
# predictor_resid_sd(): the spread of the least-squares residuals by
# their MAD (stats::mad(), consistent for the sd under normality)
# instead of their sd, so that outliers in heavy-tailed data do not set
# the unit. Sourced after library(frmtmb); this process only.
local({
  ns <- asNamespace("frmtmb")
  f <- get("predictor_resid_sd", ns)
  b <- deparse(body(f))
  i <- grep("s <- stats::sd(z)", b, fixed = TRUE)
  stopifnot(length(i) == 1L)
  b[i] <- sub("stats::sd(z)", "stats::mad(z)", b[i], fixed = TRUE)
  body(f) <- parse(text = b)[[1]]
  utils::assignInNamespace("predictor_resid_sd", f, "frmtmb")
  cat("predictor_resid_sd patched: mad of the residuals\n")
})
