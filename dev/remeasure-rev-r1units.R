# Reviewer probe for R1: par_units of fit_imp2 under autoscale = FALSE,
# and 0.69.1's par_units, to show why the switch does not restore SEs.
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
data("nhanes", package = "mice")
b <- bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) + set_rescor(FALSE)
for (a in list(NULL, FALSE)) {
  f <- suppressWarnings(frm(b, data = nhanes,
                            control = frmtmb_control(autoscale = a)))
  cat(format(packageVersion("frmtmb")), "autoscale", format(a),
      "par_units:", format(f$par_units %||% "NULL"), "\n")
}
