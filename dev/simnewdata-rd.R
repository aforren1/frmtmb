# Render the two Rd files this lane changed, as a reader sees them.
#   Rscript dev/simnewdata-rd.R > dev/simnewdata-log/rd.txt
for (f in c("man/simulate.frmtmb_fit.Rd", "man/pp_check.Rd")) {
  cat("=====", f, "\n")
  tools::Rd2txt(f, options = list(underline_titles = FALSE))
}
