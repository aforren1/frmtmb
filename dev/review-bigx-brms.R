.libPaths(c("C:/Users/adf44/source/r/rellib-r18", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), find.package("brms"), "\n")
set.seed(501); n <- 200
x1 <- rnorm(n); d <- data.frame(x1 = x1, x2a = x1 + 1e-9 * rnorm(n))
d$y <- 1 + 0.5 * x1 + rnorm(n)
sd1 <- standata(bf(y ~ x1 + x2a), data = d)
cat("x2 = x1 + 1e-9 noise: brms K =", sd1$K, " columns:", colnames(sd1$X), "\n")
xx <- runif(n, 1000, 1010); d$xx <- xx; d$yy <- sin((xx - 1000) / 2) + rnorm(n, 0, 0.3)
for (deg in c(5, 8)) {
  sd2 <- standata(bf(as.formula(sprintf("yy ~ poly(xx, %d, raw = TRUE)", deg))), data = d)
  cat("raw poly", deg, "at [1000,1010]: brms K =", sd2$K, "\n")
}
sc <- stancode(bf(y ~ x1 + x2a), data = d)
cat(grep("Xc|means_X", strsplit(sc, "\n")[[1]], value = TRUE)[1:4], sep = "\n")
