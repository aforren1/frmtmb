setwd("C:/Users/adf44/source/r/frmtmb-wt-famuni")
ed <- function(p, old, new) {
  x <- paste(readLines(p), collapse = "\n")
  if (!grepl(old, x, fixed = TRUE)) stop("not found in ", p)
  writeLines(strsplit(sub(old, new, x, fixed = TRUE), "\n")[[1]], p)
}
ed("extensions/frmtmb.sample/NEWS.md",
"  within 0.005 of brms's at the first two. One difference
  from brms remains: frmtmb holds `xi` above -1 on every route, where
  brms's range goes on (`dev/famuni-out/tmpxi.txt`, `tmpxi-sample.txt`).
  A prior on class `xi` is then a prior on `tmp_xi`, as in brms. A
  `cens()` or `trunc()` response, or a second GEV response, keeps the
  `log1p` scale, and `frm_sample()` warns that it does.",
"  within 0.005 of brms's at the first two. One difference
  from brms remains: frmtmb holds `xi` above -1 on every route, where
  brms's range goes on (`dev/famuni-out/tmpxi.txt`, `tmpxi-sample.txt`).
  The hold is material where the data put mass there: on the xi -0.85
  data brms puts 29 percent of its posterior below -1, and
  `frm_sample()`'s posterior is brms's conditioned on `xi > -1` (median
  of `xi` -0.876 against -0.877; the review's round 2). A prior on
  class `xi` is then a prior on `tmp_xi`, as in brms. A `cens()` or
  `trunc()` response, a second GEV response, or GEV components of a
  `mixture()` keep the `log1p` scale, and `frm_sample()` warns that they
  do.")
ed("extensions/frmtmb.eam/NEWS.md",
"  Under an identity link the density, and the fitted objective, is `NaN`
  for a `bs` at or below zero or a `bias` outside (0, 1), where the
  series alone returned a finite number for a bias of 1.2.",
"  Under an identity link the density, and the fitted objective, is `NaN`
  for a `bs` at or below zero, a `bias` outside (0, 1), or an `ndt` at
  or below zero (`link_ndt = \"identity\"`), where the series alone
  returned a finite number for a bias of 1.2.")
cat("ok\n")
