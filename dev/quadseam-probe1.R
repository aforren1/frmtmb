# frozen GK tape against a freshly anchored one, at the boundary probe
# points of the item-3477 fit (log sd moved by -2, and the up grid)
.libPaths(c("C:/Users/adf44/source/r/rellib-r20","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(3)
n <- 120
d1 <- data.frame(g = factor(rep(1:20, each = 6)), x = rnorm(n))
d1$y <- rbinom(n, 1, plogis(0.2 + 0.5 * d1$x))
f <- frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1, quadrature = TRUE)
cat("env names:", head(ls(f$obj$env), 60), "\n")
gk <- list(b = structure(list(dim = 1, adaptive = FALSE, debug = FALSE,
           method = "marginal_gk"), class = "GK"))
nll <- frmtmb:::build_objective(f$frame)
fresh <- function(q) {
  tpl <- f$estimates
  tpl$beta[] <- q[1:2]; tpl$theta[] <- q[3]
  t0 <- proc.time()[[3]]
  o <- RTMB::MakeADFun(nll, tpl, random = "b", integrate = gk, silent = TRUE)
  v <- o$fn(o$par)
  attr(v, "secs") <- proc.time()[[3]] - t0
  v
}
p <- f$opt$par
lap <- frm(bf(y ~ x + (1 | g)), family = bernoulli(), data = d1)
for (th in c(p[3] - 2, p[3], log(c(0.01, 0.03, 0.1, 0.3, 1)), p[3] + 10)) {
  q <- p; q[3] <- th
  a <- f$obj$fn(q); b <- fresh(q)
  ql <- lap$opt$par; ql[3] <- th
  cat(sprintf("theta %8.3f frozen %.10f fresh %.10f d %.3e laplace %.10f secs %.2f\n",
              th, a, b, a - b, lap$obj$fn(ql), attr(b, "secs")))
}
