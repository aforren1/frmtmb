# Item 1: fitted(mv) on the lane against brms 2.23.0's printed values
# for the same three rows (dev/predfix-log/brms.txt).
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
set.seed(20260921)
n <- 150
d <- data.frame(x = rnorm(n), g = factor(rep(1:15, each = 10)))
u <- rnorm(15, 0, 0.8)
e <- matrix(rnorm(2 * n), n, 2) %*% chol(matrix(c(1, 0.7, 0.7, 1), 2))
d$y1 <- 1 + 0.5 * d$x + e[, 1]
d$y2 <- -0.3 * d$x + e[, 2]
fr <- frm(bf(mvbind(y1, y2) ~ x) + gaussian() + set_rescor(TRUE), data = d)
print(fitted(fr, newdata = d[1:3, ]))
