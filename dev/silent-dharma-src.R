# Prints DHARMa's getQuantile() source and its documentation details.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
print(DHARMa::getQuantile)
h <- tools::Rd_db("DHARMa")
t <- capture.output(tools::Rd2txt(h[["getQuantile.Rd"]]))
cat(t, sep = "\n")
