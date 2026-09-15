# Lane eamhier, punch round 1: does the recovery table RENDER?
#
# An unescaped `%` starts a comment in Rd even inside \preformatted{},
# so the coverage column lost its percent sign and anything after it.
# Reading the .Rd source does not show that; rendering does, which is
# why this renders.
#
# Run: Rscript --vanilla dev/eamhier-scripts/eamhier-rdcheck.R
rd <- "extensions/frmtmb.eam/man/wiener.Rd"
out <- tempfile(fileext = ".txt")
tools::Rd2txt(rd, out = out)
x <- readLines(out, warn = FALSE)
i <- grep("Wald coverage|mu intercept|mu condition effect|^  bs  ",
          x)
i <- c(i, grep("sd(mu | s) ", x, fixed = TRUE),
       grep("sd(log bs | s) ", x, fixed = TRUE))
cat("== the rendered recovery table ==\n")
writeLines(x[sort(unique(i))])
cat("\n== every percent sign that survived rendering ==\n")
writeLines(grep("%", x, value = TRUE))
