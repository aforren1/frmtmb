p <- "extensions/frmtmb.sample/NEWS.md"
x <- paste(readLines(p), collapse = "\n")
old <- "  `tmp_xi`. On its `log1p` scale the support's edge was a wall in the
  sampled space. On n = 40 data sets at `xi` -0.45, 0.20 and -0.85
  (4 chains x 2000), 957, 227 and 3388 of 4000 transitions diverged that
  way, with Rhat 1.97 on the last; through `tmp_xi`, none did, and the
  posterior of `xi` is brms's to 0.01 at the first two."
new <- "  `tmp_xi`. On its `log1p` scale the support's edge was a wall in the
  sampled space. On the review's n = 40 data sets at `xi` -0.45, 0.20
  and -0.85 (4 chains x 2000), 1187, 263 and 3856 of 4000 transitions
  diverged, with Rhat 2.30 on the last, against 0, 1 and 0 in brms.
  Through `tmp_xi` none diverge, and the posterior median of `xi` is
  within 0.005 of brms's at the first two."
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], p)
