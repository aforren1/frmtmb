# Re-check item 4: do influence() refits of a PLAIN predictor re-read an
# environment constant changed after the fit, as trials(k) does? Seed 93.
.libPaths(c("C:/Users/adf44/source/r/wt-formrobust-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(93)
n <- 60
d <- data.frame(x = rnorm(n), g = factor(rep(1:10, each = 6)))
d$y <- rnorm(n, 0.5 * d$x)
ck <- 0.5
f1 <- frm(bf(y ~ x + offset(x * ck) + (1 | g)), data = d)
a <- influence(f1, groups = "g")$fixed
ck <- 2
b <- influence(f1, groups = "g")$fixed
cat("offset(x * ck): influence refits unchanged after ck change:",
    identical(a, b), "\n")
cat("ratio of slope refits after/before:", range(b[, 2] / a[, 2]), "\n")
nd <- d[1:3, ]
ck <- 0.5; p1 <- fitted(f1, newdata = nd)[, 1]
ck <- 2; p2 <- fitted(f1, newdata = nd)[, 1]
cat("fitted(newdata) moves with ck:", !identical(p1, p2), "\n")
cat("-- I(x * ck) on the base build and the lane: influence works? --\n")
ck <- 0.5
f2 <- frm(bf(y ~ I(x * ck) + (1 | g)), data = d)
r <- tryCatch({influence(f2, groups = "g"); "ok"}, error = function(e) conditionMessage(e))
cat("influence(I(x * ck)):", substr(r, 1, 150), "\n")
