.libPaths(c("C:/Users/adf44/source/r/rellib-r16","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
g <- factor(rep(seq_len(15L), each = 10)); z <- rnorm(150); u <- rnorm(15, 0, 0.7)[g]
for (s in c(1e5, 1e6)) {
  d <- data.frame(y = (0.3 + 0.5 * z + u + 0.5 * rt(150, 5)) * s, z, g)
  m <- tryCatch({suppressWarnings(frm(bf(y ~ z + (1 | g), nu = 4), family = student(), data = d, importance = 300L)); "fit"}, error = function(e) substr(conditionMessage(e), 1, 140))
  cat("rellib-r16 student y*", s, ":", m, "\n")
}
