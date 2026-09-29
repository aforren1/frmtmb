# Lane splinecurve: render each changed Rd to text so it can be read.
#   Rscript dev/splinecurve-rd.R > dev/splinecurve-rd.log
man <- "C:/Users/adf44/source/r/frmtmb-wt-release/extensions/frmtmb.spline/man"
for (f in c("frm_curve.Rd", "frm_curve_deriv.Rd", "frm_curve_feature.Rd")) {
  cat("\n=================== ", f, "\n", sep = "")
  tools::Rd2txt(file.path(man, f), options = list(underline_titles = FALSE))
}
