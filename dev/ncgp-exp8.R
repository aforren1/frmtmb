# y = exp(0.8 x) + 0.3 e, seed 5701: why does frm() stop 0.68 nats short? Reads
# the tape's value and gradient against the closed form at the closed-form
# maximum and along a line through it, then what a restart does.
#   Rscript dev/ncgp-exp8.R
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(5701)
n <- 100
x <- sort(runif(n, 0, 10))
e0 <- rnorm(n)
y <- exp(0.8 * x) + 0.3 * e0
d <- data.frame(x = x, y = y)
w <- NULL
fit <- withCallingHandlers(frm(bf(y ~ gp(x)), data = d, family = gaussian(), control = frmtmb_control(check_se = "ignore")),
  warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
cat("conv", fit$opt$convergence, fit$opt$message, "ll", -fit$opt$objective, "\n")
print(w)
obj <- fit$obj
print(names(obj$par))
print(obj$par)
print(fit$opt$par)
print(max(abs(obj$gr(fit$opt$par))))
print(obj$gr(fit$opt$par))
D2 <- outer(x, x, "-")^2
marg <- function(p) {
  s2 <- exp(2 * p[2]); K <- s2 * exp(-D2 / (2 * exp(2 * p[3])))
  diag(K) <- diag(K) * (1 + 1e-12)
  diag(K) <- diag(K) + exp(2 * p[4])
  R <- chol(K)
  r <- backsolve(R, y - p[1], transpose = TRUE)
  0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
}
# p in frm order: beta, betad(log sigma), theta1, theta2
order_p <- function(p) c(p[["beta"]], p[names(p) == "theta"][1], p[names(p) == "theta"][2], p[["betad"]])
pf <- fit$opt$par
cat("fn at opt", obj$fn(pf), "marg", marg(order_p(pf)), "\n")
# gradient of the marginal by central differences at frm's point
ng <- function(p, h = 1e-4) vapply(seq_along(p), function(i) {
  e <- replace(numeric(length(p)), i, h)
  (marg(p + e) - marg(p - e)) / (2 * h)
}, 0)
cat("num grad marg (beta, th1, th2, logsigma):", ng(order_p(pf)), "\n")
gf <- obj$gr(pf)
cat("tape grad same order:", gf[names(pf) == "beta"], gf[names(pf) == "theta"], gf[names(pf) == "betad"], "\n")
# noise of the tape's value along lscale
cf <- nlminb(order_p(pf), marg, control = list(iter.max = 1000, eval.max = 2000, rel.tol = 1e-12))
cat("cf from frm's point: ll", -cf$objective, "par", cf$par, "\n")
pc <- pf
pc[names(pf) == "beta"] <- cf$par[1]
pc[names(pf) == "theta"] <- cf$par[2:3]
pc[names(pf) == "betad"] <- cf$par[4]
cat("tape at cf max", obj$fn(pc), "gr", obj$gr(pc), "\n")
for (eps in c(1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 1e-6)) {
  p2 <- pc; p2[names(pf) == "theta"][2] <- p2[names(pf) == "theta"][2] + eps
  m2 <- order_p(p2)
  cat(sprintf("th2 + %.0e: tape %.8f marg %.8f diff %.2e\n", eps, obj$fn(p2), marg(m2), obj$fn(p2) - marg(m2)))
}
st <- nlminb(pc, obj$fn, obj$gr, control = list(iter.max = 300, eval.max = 400))
cat("nlminb from cf max: conv", st$convergence, st$message, "ll", -st$objective, "\n")
st2 <- nlminb(pf, obj$fn, obj$gr, control = list(iter.max = 300, eval.max = 400))
cat("nlminb from frm's stop: conv", st2$convergence, st2$message, "ll", -st2$objective, "\n")
b <- fit$obj$env$last.par.best
Hd <- fit$obj$env$spHess(b, random = TRUE)
ev <- eigen(as.matrix(Hd), only.values = TRUE)$values
cat("inner H: max eig", max(ev), "min", min(ev), "\n")
