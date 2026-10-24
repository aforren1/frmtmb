# What the internal us_sigma() of the base build does with a theta of
# the wrong length (it recycles with a warning). Log:
# dev/mvnlatent-out/us-sigma-base-length.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
f <- get("us_sigma", asNamespace("frmtmb"))
cat("frmtmb", as.character(packageVersion("frmtmb")), "\n")
for (k in list(list(numeric(5), 3L), list(numeric(4), 2L),
               list(numeric(7), 3L))) {
  r <- tryCatch(f(k[[1]], k[[2]]), error = function(e) conditionMessage(e),
                warning = function(w) paste("WARNING:", conditionMessage(w)))
  cat("length", length(k[[1]]), "d", k[[2]], "->\n")
  print(r)
}
