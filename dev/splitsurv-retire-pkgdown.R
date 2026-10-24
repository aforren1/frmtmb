# The retirement half of core's _pkgdown.yml: frmtmb.spline leaves the
# menu, and its subsite's pages redirect to the pages that now carry the
# same functions. Separate from dev/splitsurv-edit-pkgdown.R so that the
# retirement is its own hunk. Run once.
f <- "_pkgdown.yml"
x <- readLines(f)
one <- function(s) {
  i <- which(x == s)
  if (length(i) != 1L) stop("not once: ", s)
  i
}
i <- one("      - text: frmtmb.spline (curve inference and spline survival)")
stopifnot(grepl("frmtmb.spline/$", x[i + 1L]))
x <- x[-c(i, i + 1L)]
i <- one("redirects:")
red <- c(
  "frmtmb.spline/index.html" = "frmtmb.survival/index.html",
  "frmtmb.spline/reference/index.html" = "reference/index.html",
  "frmtmb.spline/reference/frm_curve.html" = "reference/frm_curve.html",
  "frmtmb.spline/reference/frm_curve_deriv.html" =
    "reference/frm_curve_deriv.html",
  "frmtmb.spline/reference/frm_curve_feature.html" =
    "reference/frm_curve_feature.html",
  "frmtmb.spline/reference/royston_parmar.html" =
    "frmtmb.survival/reference/royston_parmar.html",
  "frmtmb.spline/reference/rp_floored.html" =
    "frmtmb.survival/reference/rp_floored.html",
  "frmtmb.spline/reference/frmtmb.spline-package.html" =
    "frmtmb.survival/reference/frmtmb.survival-package.html",
  "frmtmb.spline/articles/curve-inference.html" =
    "articles/curve-inference.html",
  "frmtmb.spline/articles/royston-parmar.html" =
    "frmtmb.survival/articles/royston-parmar.html",
  "frmtmb.spline/news/index.html" = "frmtmb.survival/news/index.html")
lines <- sprintf("  - [\"%s\", \"%s\"]", names(red), unname(red))
x <- append(x, lines, after = i)
writeLines(x, f)
cat(length(lines), "redirects\n")
