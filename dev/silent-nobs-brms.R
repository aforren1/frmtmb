# Item 1562: brms 2.23.0's nobs() on the same mvbf with subset(). An
# empty brmsfit holds the data, which is all nobs.brmsfit() reads.
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(5)
d6 <- data.frame(x = rnorm(80))
d6$y1 <- d6$x + rnorm(80)
d6$y2 <- d6$x + rnorm(80)
d6$s1 <- seq_len(80) <= 72
d6$s2 <- seq_len(80) <= 72
d6$s3 <- seq_len(80) <= 60
f <- brm(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s2) ~ x) +
           set_rescor(FALSE), data = d6, empty = TRUE)
cat("A nobs()", nobs(f), " nobs(resp = 'y1')", nobs(f, resp = "y1"), "\n")
g <- brm(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s3) ~ x) +
           set_rescor(FALSE), data = d6, empty = TRUE)
cat("B nobs()", nobs(g), " nobs(resp = 'y2')", nobs(g, resp = "y2"), "\n")
print(body(brms:::nobs.brmsfit))
