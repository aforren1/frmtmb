.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
gen <- function(seed, n = 500) {
  set.seed(seed)
  x <- rnorm(n); cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 2 * x, -1.5 * x) + rlogis(n)
  data.frame(x = x, cls = cls, y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5))
}
fam <- mixture(cumulative(), cumulative(), order = "mu")
d <- gen(11)
fit <- frm(bf(y ~ x), family = fam, data = d)
obj <- fit$obj
o <- nlminb(fit$opt$par, obj$fn, obj$gr, control = list(rel.tol = 1e-15, x.tol = 0, iter.max = 5000, eval.max = 5000))
cat("continued nlminb:", -o$objective, o$message, max(abs(obj$gr(o$par))), "\n"); print(o$par)
o2 <- optim(fit$opt$par, obj$fn, obj$gr, method = "BFGS", control = list(maxit = 5000, reltol = 1e-15))
cat("BFGS:", -o2$value, o2$convergence, max(abs(obj$gr(o2$par))), "\n"); print(o2$par)
H <- optimHess(o2$par, obj$fn, obj$gr); print(eigen(H)$values)
# a profile along theta: maximize over the rest at fixed theta
for (th in c(-1.5, -1, -0.5, -0.24, 0, 0.5)) {
  f <- function(q) obj$fn(append(q, th, after = 2))
  g <- function(q) obj$gr(append(q, th, after = 2))[-3]
  oo <- nlminb(fit$opt$par[-3], f, g, control = list(rel.tol = 1e-14))
  cat("theta1 logodds", th, "profile ll", -oo$objective, "b", round(oo$par[1:2], 3), "\n")
}
