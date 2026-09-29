# Reviewer, lane splinecurve: render the three changed Rd files.
#   Rscript splinecurve-rev-06-rd.R > splinecurve-rev-log/06-rd.txt
man <- "C:/Users/adf44/source/r/frmtmb-wt-release/extensions/frmtmb.spline/man"
for (f in c("frm_curve.Rd", "frm_curve_deriv.Rd", "frm_curve_feature.Rd")) {
  cat("########", f, "\n")
  tools::Rd2txt(file.path(man, f), out = stdout(),
                options = list(underline_titles = FALSE, width = 80))
}
