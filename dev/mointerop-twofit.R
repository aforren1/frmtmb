.libPaths(c("C:/Users/adf44/source/r/wt-mointerop-lib", "C:/Users/adf44/source/r/rellib-r14", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(5); n <- 200
d <- data.frame(x = sample(1:2, n, TRUE), z = rnorm(n))
d$y <- 1 + 0.8 * d$x + 0.5 * d$z + rnorm(n); d$x01 <- d$x - 1L
f1 <- frm(y ~ mo(x) + z, data = d); f2 <- frm(y ~ x01 + z, data = d); m <- lm(y ~ x01 + z, data = d)
print(rbind(mo = fixef(f1)[, 1], lin = fixef(f2)[c(1, 3, 2), 1], lm = coef(m)[c(1, 3, 2)]), digits = 12)
print(f1$opt$message); print(f1$opt$iterations); print(max(abs(f1$opt$gradient %||% NA)))
print(c(logLik(f1), logLik(f2), logLik(m)), digits = 15)
str(f1$opt[c("polish", "mo_search", "convergence", "evaluations")])
f3 <- frm(y ~ mo(x) + z, data = d, verbose = TRUE)
d$x3 <- sample(1:3, n, TRUE); d$y3 <- 1 + c(0, .5, 1.2)[d$x3] + 0.5 * d$z + rnorm(n)
f4 <- frm(y3 ~ mo(x3) + z, data = d); m4 <- lm(y3 ~ factor(x3) + z, data = d)
print(rbind(mo = fitted(f4, newdata = data.frame(x3 = 1:3, z = 0))[, 1], lm = predict(m4, newdata = data.frame(x3 = 1:3, z = 0))), digits = 12)
