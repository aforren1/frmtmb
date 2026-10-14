# Reviewer: does frmtmb.sample on rellib-r8 name a correlated block's
# draws theta_<k> too (a categorical (1 | ID | g) fit)?
# Output: dev/famsimplex-rev-sample-base.log
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), "\n")
cat(La_library(), "\n")
set.seed(3)
n <- 300
d <- data.frame(x = rnorm(n), g = factor(rep(1:15, length.out = n)))
u <- matrix(rnorm(30, 0, 0.5), 15)
E <- cbind(0, 0.3 + 0.5 * d$x + u[d$g, 1], -0.2 + u[d$g, 2])
P <- exp(E) / rowSums(exp(E))
d$f <- factor(c("a", "b", "c")[apply(P, 1, function(p) sample(3, 1, prob = p))])
fit <- frm(bf(f ~ x + (1 | ID | g)), family = categorical(), data = d)
ds <- suppressWarnings(suppressMessages(frm_sample(fit, chains = 1,
                                                    iter = 200, refresh = 0,
                                                    seed = 1)))
v <- variables(ds)
cat("non-r_ variables:", v[!grepl("^r_", v)], "\n")
