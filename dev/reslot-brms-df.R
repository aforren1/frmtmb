# Lane reslot, fix round 1: does brms 2.23.0 take a plain data frame as
# its prior? validate_prior() only.
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
set.seed(1)
dd <- data.frame(y = rnorm(40), x = rnorm(40), g = factor(rep(1:8, 5)))
df1 <- data.frame(prior = "normal(0, 1)", class = "b", coef = "",
                  group = "", resp = "", dpar = "", nlpar = "", lb = "",
                  ub = "", source = "user")
r <- tryCatch(validate_prior(df1, y ~ x + (1 | g), data = dd),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)
df2 <- data.frame(prior = "normal(0, 1)", class = "b")
r <- tryCatch(validate_prior(df2, y ~ x + (1 | g), data = dd),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)
cat("as.brmsprior(df2):\n")
print(tryCatch(as.brmsprior(df2), error = function(e) conditionMessage(e)))
