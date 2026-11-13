# Lists the DHARMa help pages and vignette lines that mention censoring.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
cat("DHARMa", format(packageVersion("DHARMa")), "\n")
h <- tools::Rd_db("DHARMa")
for (n in names(h)) {
  t <- capture.output(tools::Rd2txt(h[[n]]))
  hit <- grep("censor|truncat|integerResponse|randomi", t, ignore.case = TRUE)
  if (length(hit)) {
    cat("==", n, "\n")
    for (i in hit) cat(t[max(1, i - 2):min(length(t), i + 3)], sep = "\n")
  }
}
