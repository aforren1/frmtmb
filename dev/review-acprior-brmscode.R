# Reviewer: brms 2.23.0's full Stan code for cov = TRUE ar(), to see
# whether sigma is the innovation or the marginal sd.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
d <- data.frame(y = rnorm(40), x = rnorm(40), t = rep(1:10, 4),
                g = rep(1:4, each = 10))
sc <- stancode(y ~ x + ar(time = t, gr = g, p = 1, cov = TRUE), data = d)
cat(sc)
cat("\n\n===== ma cov\n")
sc2 <- stancode(y ~ x + ma(time = t, gr = g, cov = TRUE), data = d)
ln <- strsplit(sc2, "\n")[[1]]
i <- grep("cholesky_cor_ma1", ln)[1]
cat(ln[i:(i + 30)], sep = "\n")
