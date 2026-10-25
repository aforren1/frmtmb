# brms 2.23.0, interval-censored cox(): the Stan code scores the
# interval [Y, rcens] with the cumulative baseline hazard at Y for both
# ends, so the interval's probability is exactly 0. No compile:
# stancode() and standata().
#   Rscript dev/bhaz-brms-icens.R > dev/bhaz-log/brms-icens.txt 2>&1
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
d <- data.frame(t = rexp(30) + 0.1, c2 = rep(c(0, 2), 15), x = rnorm(30))
d$u <- d$t * 1.5
fo <- bf(t | cens(c2, u) + bhaz(df = 6) ~ x)
sc <- stancode(fo, data = d, family = cox())
x <- strsplit(as.character(sc), "\n")[[1]]
cat("--- Stan lines for the censoring and the baseline ---\n")
cat(grep("Jicens|rcens|cbhaz\\[|Zcbhaz", x, value = TRUE), sep = "\n")
sd <- standata(fo, data = d, family = cox())
cat("--- standata ---\n")
cat("names with cbhaz:", grep("cbhaz", names(sd), value = TRUE), "\n")
cat("rows of Zcbhaz", nrow(sd$Zcbhaz), "= N", sd$N,
    "; no basis at rcens is passed\n")
# the interval term as the generated code computes it, at equal weights
# and mu = 1: log_diff_exp(lcdf(rcens | cbhaz[n]), lcdf(Y | cbhaz[n]))
s <- rep(1 / sd$Kbhaz, sd$Kbhaz)
cb <- as.vector(sd$Zcbhaz %*% s)
i <- which(d$c2 == 2)[1L]
lcdf <- function(cbhaz) log1p(-exp(-cbhaz))
cat(sprintf("row %d: Y %.4f rcens %.4f; lcdf at Y %.6f, at rcens as coded %.6f",
            i, d$t[i], d$u[i], lcdf(cb[i]), lcdf(cb[i])),
    "\n interval log probability as coded:",
    log(exp(lcdf(cb[i])) - exp(lcdf(cb[i]))), "\n")
