# What brms 2.23.0 treats as one gp(gr = TRUE) position: rows that
# differ in the last bit of a coordinate.
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
print(body(brms:::match_rows))
set.seed(1)
x <- c(1/3, 1/3 * (1 + 2^-52), seq(0.5, 3, length.out = 8))
d <- data.frame(x = x, y = rnorm(10))
sd <- brms::standata(bf(y ~ gp(x)), data = d)
cat("brms NBgp_1 / Xgp rows:", nrow(sd$Xgp_1), " Jgp_1:", sd$Jgp_1, "\n")
cat("identical(x[1], x[2]):", identical(x[1], x[2]), "\n")
