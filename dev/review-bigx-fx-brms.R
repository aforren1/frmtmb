.libPaths(c("C:/Users/adf44/source/r/rellib-r18", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(11); d <- data.frame(x = rnorm(50) + 10); d$y <- 1 + 0.4 * d$x + rnorm(50)
sc <- stancode(bf(y ~ x), data = d, prior = set_prior("normal(0, 100)", class = "Intercept", lb = 3))
l <- strsplit(sc, "\n")[[1]]
cat(grep("Intercept", l, value = TRUE), sep = "\n")
