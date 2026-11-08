# Cost of one frm_linpred() call on new rows, the unit of a node of the
# draw's integral. Seed 4101, the data of readers-item1.R.
source("dev/readers-setup.R")
set.seed(4101)
J <- 40; n <- 400
g <- factor(sample(seq_len(J), n, TRUE))
u <- rnorm(J, 0, 1)
x <- runif(n, -1, 1)
dd <- data.frame(x = x, g = g, zo = rnorm(n), sz = 0.4)
dd$y <- rpois(n, exp(0.5 + 0.3 * x + u[g]))
fit <- frm(bf(y ~ x + me(zo, sz) + (1 | g), hu ~ me(zo, sz)),
           family = hurdle_poisson(), data = dd)
for (nr in c(3, 100)) {
  nd <- data.frame(x = seq(-1, 1, length.out = nr), g = "1", zo = 0, sz = 0.4)
  t <- system.time(for (i in 1:200) frm_linpred(fit, newdata = nd,
                                                 type = "response"))
  cat("rows", nr, ": ms per call", 1000 * t[["elapsed"]] / 200, "\n")
  t2 <- system.time(for (i in 1:5) fitted(fit, newdata = nd))
  cat("rows", nr, ": ms per fitted()", 1000 * t2[["elapsed"]] / 5, "\n")
}
