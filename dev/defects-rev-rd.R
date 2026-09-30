# Reviewer of lane defects: render the changed Rd files and check them.
for (f in c("man/print.frmtmb_family.Rd", "man/residuals.frmtmb_fit.Rd",
            "man/ranef.Rd", "man/mixture.Rd", "man/fitted.frmtmb_fit.Rd",
            "man/predict.frmtmb_fit.Rd",
            "extensions/frmtmb.sample/man/sample-as_draws.Rd")) {
  out <- capture.output(tools::Rd2txt(f, options = list(underline_titles = FALSE)))
  msgs <- tools::checkRd(f)
  cat("==", f, "lines", length(out), "checkRd", length(msgs), "\n")
  if (length(msgs)) print(msgs)
}
out <- capture.output(tools::Rd2txt("man/print.frmtmb_family.Rd"))
cat(out, sep = "\n")
