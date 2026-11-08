# Review of lane readers: lp_basis_nl() tapes one input per row for each
# body part with a me() draw (nl_basis_extra()), and takes the full Jacobian,
# rows x (coefficients + rows). Cost of frm_linpred(se.fit = TRUE) on new
# rows of the body fit of dev/review-readers-body.R part B (data seed 8202),
# against the same rows with sx = 1e-12 (the same tape, a draw of no size).
source("dev/review-readers-setup.R")
set.seed(8202)
n <- 600
d <- data.frame(l = rnorm(n), sx = 0.3)
d$xo <- d$l + rnorm(n, 0, 0.3)
d$y <- exp(0.2 + 0.8 * d$l) + rnorm(n, 0, 0.3)
fB <- frm(bf(y ~ exp(a), a ~ me(xo, sx), nl = TRUE), family = gaussian(),
          data = d)
for (m in c(250, 500, 1000, 2000, 4000)) {
  nd <- data.frame(xo = seq(-2, 2, length.out = m), sx = 0.3)
  nd0 <- transform(nd, sx = 1e-12)
  gc()
  t1 <- system.time(p1 <- frm_linpred(fB, newdata = nd, se.fit = TRUE))
  t0 <- system.time(p0 <- frm_linpred(fB, newdata = nd0, se.fit = TRUE))
  cat(sprintf("rows %5d: with draw %6.2f s, sx = 1e-12 %6.2f s, max memory %s MB\n",
              m, t1[["elapsed"]], t0[["elapsed"]],
              format(round(sum(gc()[, 6])))))
}
