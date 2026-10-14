# Render the two Rd pages the consolidation edited (com_poisson) and
# print the paragraphs that name it, to check the rendered text.
for (f in c("man/brmsfamily.Rd", "man/frmtmb-families.Rd")) {
  x <- capture.output(tools::Rd2txt(f, options = list(underline_titles = FALSE)))
  i <- grep("com_poisson|Conway", x)
  cat("==", f, "\n")
  for (j in i) cat(x[j:min(length(x), j + 9)], sep = "\n")
}
