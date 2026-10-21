# hurdle_gamma: brms's gradient is 0.18 at frmtmb's censored estimate
# while the values agree. frmtmb's own gradient there, and against
# finite differences of its objective. Data as test-brms-lik-cdfs.R's.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib", "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(104); n <- 200; keep <- rbinom(n, 1, 0.75)
y1 <- rpois(n, 3) * keep; y2 <- rnbinom(n, mu = 3, size = 2) * keep
y3 <- rbinom(n, 12, 0.35) * keep; y4 <- brms::rbeta_binomial(n, 12, 0.35, 5) * keep
y5 <- pmax(rpois(n, 3), 1) * keep; y6 <- pmax(rnbinom(n, mu = 3, size = 2), 1) * keep
y <- rgamma(n, 3, 3 / 2) * keep
set.seed(12)
d <- data.frame(y = y, x = rnorm(n), cc = sample(c(-1, 0, 0, 0, 1, 2), n, TRUE))
d$y2 <- d$y + 0.1
d$cc[d$y == 0] <- 0
fit <- frm(y | cens(cc, y2) ~ x, family = hurdle_gamma(), data = d)
p <- fit$opt$par
g <- fit$obj$gr(p)
fd <- numDeriv::grad(fit$obj$fn, p)
print(rbind(par = p, ad = as.vector(g), fd = fd))
cat("conv", fit$opt$convergence, fit$opt$message, "\n")
o2 <- nlminb(p, fit$obj$fn, fit$obj$gr, control = list(rel.tol = 1e-14, iter.max = 500))
cat("polished: fn", fit$obj$fn(p), "->", o2$objective, " par change", max(abs(o2$par - p)), "\n")
