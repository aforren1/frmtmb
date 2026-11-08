# Review lane osa, short check: seed = NULL with row blocks. The uniforms
# frmtmb draws again must be the ones TMB drew block by block from the
# session stream, and the stream must end where n draws leave it.
#   Rscript dev/review-osa-rng2.R
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-rng2\n")
set.seed(901)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(rep(1:30, each = 10)))
d$y <- rpois(n, exp(0.8 + 0.5 * d$x + rnorm(30, 0, 0.5)[d$g]))
f0 <- frm(y ~ x, family = poisson(), data = d)
f1 <- frm(y ~ x + (1 | g), family = poisson(), data = d)
b <- fixef(f0)[, "Estimate"]
mu <- exp(b[[1]] + b[[2]] * d$x)
set.seed(5)
r1 <- residuals(f0, type = "osa", seed = NULL)[, "Estimate"]
s_after <- .Random.seed
set.seed(5)
U <- runif(n)
s_ref <- .Random.seed
cat(sprintf(paste("seed = NULL, %d rows in blocks: max|r - qnorm(F - U p)|",
                  "with U = set.seed(5); runif(n): %.3g; stream after equals",
                  "n draws: %s\n"), n,
            max(abs(r1 - qnorm(ppois(d$y, mu) - U * dpois(d$y, mu)))),
            identical(s_after, s_ref)))
set.seed(5)
r1b <- residuals(f0, type = "osa", seed = NULL)[, "Estimate"]
cat("seed = NULL twice from the same state agree:", identical(r1, r1b), "\n")
set.seed(42)
s0 <- .Random.seed
ra <- residuals(f1, type = "osa")[, "Estimate"]
cat("(1 | g), default seed: session stream unchanged:", identical(s0, .Random.seed), "\n")
rb <- residuals(f1, type = "osa")[, "Estimate"]
cat("(1 | g), two calls identical:", identical(ra, rb), "\n")
# a row's position in its step, read with the documented uniforms
op <- options(frmtmb.osa_block = Inf)
rc <- residuals(f1, type = "osa")[, "Estimate"]
options(op)
cat(sprintf("(1 | g): blocked vs one call max|diff| %.3g\n", max(abs(ra - rc))))
cat("END review-osa-rng2\n")
