# Reviewer, punch round 1: brms on `y ~ . - w` refit on data without w.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(6)
d <- data.frame(y = rnorm(30), x1 = rnorm(30), w = rnorm(30))
fit <- brm(y ~ . - w, data = d, empty = TRUE)
cat("brms stored:", deparse1(fit$formula$formula), "\n")
r <- tryCatch(colnames(standata(fit$formula, data = d[, c("y", "x1")])$X),
              error = function(e) paste("ERR:", conditionMessage(e)))
cat("standata(stored formula, data without w):", r, "\n")
cat("stats:", deparse1(formula(terms(y ~ . - w, data = d))), "\n")
