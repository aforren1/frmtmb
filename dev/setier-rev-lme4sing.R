# Reviewer of lane setier: what lme4 does with a singular fit (level,
# wording, control), read from the installed lme4.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
library(lme4)
cat("lme4", format(packageVersion("lme4")), "\n")
print(formals(lmerControl)$check.conv.singular)
print(formals(glmerControl)$check.conv.singular)
f <- deparse(lme4:::checkConv)
i <- grep("singular|isSingular", f)
print(f[sort(unique(c(i - 1, i, i + 1)))])
print(lme4:::.makeCC)
set.seed(30)
d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
d$y <- 1 + 0.5 * d$x + rnorm(100)
cnd <- list()
m <- withCallingHandlers(lmer(y ~ x + (1 | g), data = d, REML = FALSE),
  condition = function(x) cnd[[length(cnd) + 1L]] <<- x)
for (x in cnd) cat(class(x)[1:2], ":", conditionMessage(x), "\n")
