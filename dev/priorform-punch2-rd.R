# Render frm_sample's Rd and confirm the replace rule reads as written.
rd <- "extensions/frmtmb.sample/man/frm_sample.Rd"
print(tools::checkRd(rd))
out <- capture.output(tools::Rd2txt(rd))
i <- grep("REPLACES", out)
writeLines(out[(i - 3):(i + 6)])
