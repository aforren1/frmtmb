# Lane prilabel, fix round 1: brms 2.23.0's id order of group-level
# terms (terms_re() sorts a predictor's terms by group), on formulas
# whose written order is not sorted. No compile.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(41)
n <- 240
dd <- data.frame(x = rnorm(n), gs = factor(rep(1:8, length.out = n)),
                 g = factor(rep(1:6, length.out = n)),
                 h = factor(rep(c("p", "q"), each = 3, length.out = n)),
                 k = factor(rep(1:4, length.out = n)))
dd$y <- rnorm(n)
for (f in list(bf(y ~ x + (1 | gs) + (1 | g:h)),
               bf(y ~ x + (1 | k) + (1 + x | g), sigma ~ (1 | k) + (1 | h)),
               bf(y ~ x + (1 | gs) + (0 + x | g) + (1 | g)))) {
  re <- as.data.frame(brms:::frame_re(brms:::brmsterms(f), data = dd))
  cat("\n==", deparse1(f$formula), "| sigma:",
      deparse1(f$pforms$sigma), "\n")
  print(re[, c("id", "group", "coef", "dpar")], row.names = FALSE)
}
