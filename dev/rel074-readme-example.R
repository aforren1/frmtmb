# The README's "Example" block, as the release fix writes it (the release
# review, finding 2), run on rellib-r12. The block is read from README.md
# itself, so this runs what a reader copies.
#   Rscript dev/rel074-readme-example.R > dev/rel074-log/readme-example.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-release")
x <- readLines("README.md")
a <- which(x == "## Example")
s <- which(x == "```r")
s <- s[s > a][1]
e <- which(x == "```")
e <- e[e > s][1]
code <- x[(s + 1L):(e - 1L)]
cat("the block, lines", s + 1L, "to", e - 1L, "of README.md:\n")
cat(code, sep = "\n")
cat("---- run\n")
set.seed(20261023)
for (ex in parse(text = code)) {
  cat("> ", deparse1(ex), "\n", sep = "")
  r <- withVisible(eval(ex, globalenv()))
  if (r$visible) print(r$value)
}
cat("README EXAMPLE DONE\n")
