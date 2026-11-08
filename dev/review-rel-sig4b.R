# Short check of the 0.80.0 release fixes: bernoulli sig4 55, default and
# with options(frmtmb.opt_basis = FALSE), every condition, the SEs, and
# whether the returned tape carries the body's basis.
#   Rscript dev/review-rel-sig4b.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
Sys.setenv(REV_ARM = strsplit(lib, ";")[[1]][1])
src <- readLines("dev/review-nlstart-fresh2.R")
eval(parse(text = src[1:(grep("^cat\\(\"ROWHEAD", src) - 1L)]))
d <- gen("bernoulli", "sig4", 55 * 1000 + match("sig4", names(B)))
fm <- form("bernoulli", "sig4")
for (o in list(list(), list(frmtmb.opt_basis = FALSE))) {
  r <- withr::with_options(o, tt(frm(fm, family = bernoulli(), data = d)))
  ob <- r$obj$env$.frmtmb_opt_basis
  se <- tryCatch(suppressWarnings(sqrt(diag(vcov(r)))), error = function(e) NA)
  cat(if (length(o)) "opt_basis off" else "default", sprintf("ll %.6f", llv(r)),
      "basis", if (is.null(ob)) "none" else paste("nl", isTRUE(ob$nl)),
      "se", format(se, digits = 3), "\n")
  for (i in seq_along(attr(r, "w"))) {
    cat("  W[", attr(r, "wc")[i], "]", substr(attr(r, "w")[i], 1, 160), "\n")
  }
  for (i in seq_along(attr(r, "m"))) {
    cat("  M[", attr(r, "mc")[i], "]", substr(attr(r, "m")[i], 1, 120), "\n")
  }
}
cat("SIG4B DONE\n")
