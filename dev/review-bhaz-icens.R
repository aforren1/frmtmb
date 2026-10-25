# Reviewer: brms 2.23.0 Stan code for an interval-censored cox() model.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- data.frame(t = rexp(30) + 0.1, c2 = rep(c(0, 2), 15), x = rnorm(30))
d$u <- d$t * 1.5
r <- tryCatch(stancode(bf(t | cens(c2, u) + bhaz(df = 6) ~ x), data = d,
                       family = cox()), error = function(e) conditionMessage(e))
x <- strsplit(as.character(r), "\n")[[1]]
cat(grep("cox_l|rcens|Zcbhaz|cbhaz", x, value = TRUE), sep = "\n")
