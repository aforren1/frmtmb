# Reviewer of lane defects: print the deparsed text of one function, to
# find a mutation site.
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(trailingOnly = TRUE)
x <- deparse(get(a[1], envir = asNamespace("frmtmb")), width.cutoff = 500L)
cat(grep(a[2], x, value = TRUE, fixed = TRUE), sep = "\n")
