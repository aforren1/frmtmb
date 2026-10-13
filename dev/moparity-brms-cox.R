# brms 2.23.0's cox(): the baseline simplex's declaration, prior and
# default_prior row. No compile.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(2); n <- 300; x <- rnorm(n)
tt <- rexp(n, exp(-0.5 + 0.7 * x)); ct <- rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
fo <- bf(time | cens(cens) ~ x)
sc <- stancode(fo, data = dd, family = cox())
cat(grep("bhaz", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
dp <- default_prior(fo, data = dd, family = cox())
print(dp[, c("prior", "class", "coef", "source")])
sd <- standata(fo, data = dd, family = cox())
cat("Kbhaz", sd$Kbhaz, "con_sbhaz", sd$con_sbhaz, "\n")
