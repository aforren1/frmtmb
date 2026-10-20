# Does the 0.71.0 advice, hypothesis(method = "profile"), work on the
# names the refusal gives it? On rellib-r9. Output: dev/catmulti-out/profile-advice.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(6)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:20, 10)))
b <- matrix(rnorm(60, 0, 0.6), 20, 3)
d$y <- 1 + b[d$g, 1] + (0.5 + b[d$g, 2]) * d$x + (0.2 + b[d$g, 3]) * d$z + rnorm(n)
fit <- suppressWarnings(frm(bf(y ~ x + z + (1 + x + z | g)), data = d))
r <- tryCatch(hypothesis(fit, "cor_g__Intercept__x = 0", method = "profile"),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(r)) r else r$hypothesis)
r <- tryCatch(hypothesis(fit, "cor_g__Intercept__x = 0", class = NULL, method = "profile"),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(r)) r else r$hypothesis)
r <- tryCatch(hypothesis(fit, "sigma = 1", class = NULL, method = "profile"),
              error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(r)) r else r$hypothesis)
