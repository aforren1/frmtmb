# Reviewer, fix round 1: brms 2.23.0's Stan code for cens() + trunc()
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
d <- data.frame(y = c(2, 3, 5, 9, 4), x = rnorm(5), cc = c(-1, 0, 1, 2, 0), y2 = c(2, 3, 5, 9, 6))
for (fam in list(poisson(), gaussian())) {
  sc <- stancode(y | cens(cc, y2) + trunc(lb = 2, ub = 9) ~ x, data = d, family = fam)
  ln <- strsplit(sc, "\n")[[1]]
  i <- grep("cens[n] == 0", ln, fixed = TRUE)[1]
  cat("----", fam$family, "\n"); cat(ln[(i - 2):(i + 22)], sep = "\n")
}
