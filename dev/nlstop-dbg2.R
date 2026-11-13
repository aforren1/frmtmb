# Lane nlstop: debug of item 2584, the boundary probe with a Newton step
# in the other parameters, on the Hessian at the fit.
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
set.seed(8)
d <- suppressMessages(mgcv::gamSim(eg = 6, n = 200, verbose = FALSE))
d$z <- runif(nrow(d)); d$g <- factor(sample(c("a", "b", "c"), nrow(d), TRUE))
d$g8 <- factor(sample(letters[1:8], nrow(d), TRUE))
d$gr <- factor(rep(1:20, length.out = nrow(d)))
d$y <- d$y + 0.6 * rnorm(20)[d$gr]
f <- frm(bf(y ~ s(x1), sigma ~ s(x2)), data = d)
p <- f$opt$par
H <- ns$fit_outer_hessian(f)$H
obj <- f$obj
for (j in 5:6) for (down in c(3, 5, 10)) {
  q <- p; q[j] <- q[j] - down
  fq <- obj$fn(q)
  g <- drop(obj$gr(q))
  o <- setdiff(seq_along(p), j)
  st <- -solve(H[o, o], g[o])
  q1 <- q; q1[o] <- q1[o] + st
  f1 <- obj$fn(q1)
  g1 <- drop(obj$gr(q1))
  st2 <- -solve(H[o, o], g1[o])
  q2 <- q1; q2[o] <- q2[o] + st2
  f2 <- obj$fn(q2)
  cat(sprintf("j %d down %g: f(p) %.5f f(q) %.5f one step %.5f two %.5f\n",
              j, down, f$opt$objective, fq, f1, f2))
}
