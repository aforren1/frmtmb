# Lane sampfix, script 09: what the one non-finite cell of hypothesis()
# and the two of loo_compare(ds, ds) in dev/sampfix-log/03-lane.txt are.
.libPaths(c("C:/Users/adf44/source/r/wt-sampfix-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb)); suppressMessages(library(frmtmb.sample))
set.seed(1212L)
dd <- data.frame(g = factor(rep(1:6, each = 5L)), t = rep(1:5, 6L))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(6L, 0, 0.5)[as.integer(dd$g)] +
  rnorm(nrow(dd), 0, 0.7)
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
ds <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 2, iter = 300, refresh = 0, seed = 3)))
h <- hypothesis(ds, "x > 0")$hypothesis
print(h)
cat("draws of b_x above 0: ", sum(ds$draws[, "b_x"] > 0), " of ",
    nrow(ds$draws), "\n", sep = "")
lc <- suppressWarnings(loo_compare(ds, ds))
print(unclass(lc))
