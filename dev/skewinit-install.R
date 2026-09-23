LIB <- "C:/Users/adf44/source/r/skewinit-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
pkg <- "C:/Users/adf44/source/r/frmtmb-wt-skewinit"
for (f in list.files(file.path(pkg, "R"), "[.]R$", full.names = TRUE)) {
  invisible(parse(f))
}
cat("parse ok\n")
roxygen2::roxygenise(pkg, roclets = c("namespace", "rd"))
cat("roxygen ok\n")
