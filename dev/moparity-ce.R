# conditional_effects() and mo_search on mo() crossed with a factor and
# on a group-level mo() term.
# Usage: Rscript dev/moparity-ce.R
.libPaths(c("C:/Users/adf44/source/r/wt-moparity-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1)
n <- 300
dd <- data.frame(x = sample(0:3, n, TRUE),
                 f = factor(sample(c("a", "b", "c"), n, TRUE)),
                 g = factor(sample(paste0("g", 1:20), n, TRUE)))
bg <- rnorm(20, 0, 0.5)
dd$y <- 1 + (1 + (dd$f == "b") - (dd$f == "c") + bg[dd$g]) *
  c(0, 1, 1.5, 1.7)[dd$x + 1] + rnorm(n, 0, 0.5)
fit <- frm(bf(y ~ mo(x) * f), family = gaussian(), data = dd)
ce <- conditional_effects(fit, effects = "x:f")
d <- ce[[1]]
cat("ce x:f rows", nrow(d), "| cols", paste(names(d)[1:4], collapse = " "), "\n")
# the curve at each level of f against fitted() at the same rows
nd <- data.frame(x = d$x, f = d$f)
cat("max |ce - fitted(newdata)|:",
    max(abs(d$estimate__ - fitted(fit, newdata = nd)[, 1])), "\n")
fit2 <- suppressMessages(frm(bf(y ~ mo(x) + (mo(x) | g)),
                             family = gaussian(), data = dd))
ce2 <- conditional_effects(fit2, effects = "x")
cat("ce x on (mo(x) | g): rows", nrow(ce2[[1]]), "estimate",
    paste(format(ce2[[1]]$estimate__, digits = 4), collapse = " "), "\n")
ce3 <- conditional_effects(fit2, effects = "x", re_formula = NULL)
cat("ce x, re_formula = NULL: rows", nrow(ce3[[1]]), "\n")
print(summary(fit2))
cat("mo_search on fit (runs, gain):", fit$opt$mo_search, "\n")
