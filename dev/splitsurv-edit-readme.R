f <- "README.md"
x <- readLines(f)
r1 <- function(old, new) {
  i <- which(x == old)
  if (length(i) != 1L) stop("not once: ", old)
  x[i] <<- new
}
r1("| [frmtmb.spline](https://aforren1.github.io/frmtmb/frmtmb.spline/) | curve inference on any fitted smooth, and the flexible parametric survival family `royston_parmar()` |",
   "| [frmtmb.survival](https://aforren1.github.io/frmtmb/frmtmb.survival/) | survival families: brms's `cox()` with an M-spline baseline, and the flexible parametric `royston_parmar()` |")
r1("remotes::install_github(\"aforren1/frmtmb\", subdir = \"extensions/frmtmb.spline\")",
   "remotes::install_github(\"aforren1/frmtmb\", subdir = \"extensions/frmtmb.survival\")")
r1("  one), `von_mises()` circular responses, and `cox()` proportional",
   "  one), and `von_mises()` circular responses. brms's `cox()`")
r1("  hazards with a flexible baseline and Laplace frailties.",
   "  proportional hazards model is in frmtmb.survival.")
writeLines(x, f)
