# 2939: (1 + f || gr(g, cov = A)) against its numeric spelling
source("dev/refits-lib.R")
set.seed(11)
ng <- 12
A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
d8 <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                            levels = paste0("g", 1:ng)),
                 f = factor(rep(c("a", "b", "c"), 40)))
La <- t(chol(A))
ra <- 0.8 * La %*% rnorm(ng)
rb <- 0.6 * La %*% rnorm(ng)
rc <- 1.2 * La %*% rnorm(ng)
gi <- as.integer(d8$g)
d8$y <- 1 + ifelse(d8$f == "a", ra[gi], ifelse(d8$f == "b", rb[gi],
                                                 rc[gi])) + rnorm(120, 0, .5)
d8$fb <- as.numeric(d8$f == "b")
d8$fc <- as.numeric(d8$f == "c")
h1 <- frm(y ~ 1 + (1 + f || gr(g, cov = A)), data = d8, data2 = list(A = A))
h2 <- frm(y ~ 1 + (1 + fb + fc || gr(g, cov = A)), data = d8,
          data2 = list(A = A))
print(VarCorr(h1)$g$sd)
print(VarCorr(h2)$g$sd)
print(c(logLik(h1), logLik(h2)), digits = 15)
cat("h1 at h2's estimates:\n")
print(h1$obj$fn(h2$opt$par), digits = 15)
print(h2$obj$fn(h1$opt$par), digits = 15)
print(rbind(h1 = h1$opt$par, h2 = h2$opt$par))
print(h1$opt$convergence); print(h2$opt$convergence)
for (b in h1$frame$re_blocks) print(c(b$covstruct, b$label, b$cnms))
for (b in h2$frame$re_blocks) print(c(b$covstruct, b$label, b$cnms))
# the exact marginal likelihood at each fit's estimate
Z <- function(v) { M <- matrix(0, 120, ng); M[cbind(1:120, gi)] <- v; M }
Zs <- list(Z(rep(1, 120)), Z(d8$fb), Z(d8$fc))
ZAZ <- lapply(Zs, function(M) M %*% A %*% t(M))
nll <- function(p) {
  V <- diag(exp(2 * p[5]), 120)
  for (k in 1:3) V <- V + exp(2 * p[k + 1]) * ZAZ[[k]]
  -mvtnorm::dmvnorm(d8$y, rep(p[1], 120), V, log = TRUE)
}
o <- optim(c(1, 0, 0, 0, 0), nll, method = "BFGS",
           control = list(reltol = 1e-14, maxit = 2000))
o <- optim(o$par, nll, method = "Nelder-Mead",
           control = list(reltol = 1e-15, maxit = 5000))
cat("exact ML:", format(-o$value, digits = 15), " sds:",
    format(exp(o$par[2:5]), digits = 6), "\n")
