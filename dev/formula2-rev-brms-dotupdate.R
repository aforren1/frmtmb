# Reviewer: what brms stores as the formula of a `y ~ .` fit, and what
# update(newdata = ) then refits (brm(empty = TRUE) compiles nothing).
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(6)
d <- data.frame(y = rnorm(30), x1 = rnorm(30), x2 = rnorm(30))
fit <- brm(y ~ ., data = d, empty = TRUE)
cat("stored:", deparse1(fit$formula$formula), "\n")
d2 <- d; d2$extra <- rnorm(30)
up <- tryCatch(update(fit, newdata = d2, testmode = TRUE),
               error = function(e) conditionMessage(e))
if (is.character(up)) cat("update error:", up, "\n") else
  cat("after update(newdata = d2):", deparse1(up$formula$formula), "\n")
