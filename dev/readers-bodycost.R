# The body Jacobian's cost after the fix: one shift input per part. A copy of
# dev/review-readers-bodycost.R that resets the memory peak per size.
source("dev/readers-setup.R")
set.seed(8202)
n <- 600
d <- data.frame(l = rnorm(n), sx = 0.3)
d$xo <- d$l + rnorm(n, 0, 0.3)
d$y <- exp(0.2 + 0.8 * d$l) + rnorm(n, 0, 0.3)
fB <- frm(bf(y ~ exp(a), a ~ me(xo, sx), nl = TRUE), family = gaussian(),
          data = d)
for (m in c(1000, 4000, 20000)) {
  nd <- data.frame(xo = seq(-2, 2, length.out = m), sx = 0.3)
  nd0 <- transform(nd, sx = 1e-12)
  invisible(gc(reset = TRUE))
  t1 <- system.time(p1 <- frm_linpred(fB, newdata = nd, se.fit = TRUE))
  t0 <- system.time(p0 <- frm_linpred(fB, newdata = nd0, se.fit = TRUE))
  cat(sprintf("rows %5d: with draw %6.2f s, sx = 1e-12 %6.2f s, max memory %s MB\n",
              m, t1[["elapsed"]], t0[["elapsed"]],
              format(round(sum(gc()[, 6])))))
}
