# 2939: a factor in || over gr(cov = A) or gr(prec = Q), against two
# references: the same columns written as numeric variables (lme4's split,
# which frmtmb fits already), and the exact gaussian marginal likelihood.
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
d8$fa <- as.numeric(d8$f == "a")
d8$fb <- as.numeric(d8$f == "b")
d8$fc <- as.numeric(d8$f == "c")
f1 <- go("factor ||", frm(y ~ 1 + (0 + f || gr(g, cov = A)), data = d8,
                          data2 = list(A = A)))
f2 <- go("numeric ||", frm(y ~ 1 + (0 + fa + fb + fc || gr(g, cov = A)),
                           data = d8, data2 = list(A = A)))
f3 <- go("cor = FALSE", frm(y ~ 1 + (0 + f | gr(g, cov = A, cor = FALSE)),
                            data = d8, data2 = list(A = A)))
if (inherits(f1, "frmtmb_fit")) {
  cat("\nlogLik factor, numeric, cor = FALSE:\n")
  print(c(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
          as.numeric(logLik(f3))), digits = 15)
  cat("VarCorr factor:\n"); print(VarCorr(f1))
  cat("VarCorr numeric:\n"); print(VarCorr(f2))
  # the exact marginal likelihood of y ~ N(b0, sum_k s_k^2 Z_k A Z_k' +
  # sigma^2 I), maximized directly
  Z <- lapply(c("a", "b", "c"), function(l) {
    M <- matrix(0, 120, ng); M[cbind(1:120, gi)] <- d8$f == l; M
  })
  ZAZ <- lapply(Z, function(M) M %*% A %*% t(M))
  nll <- function(p) {
    V <- diag(exp(2 * p[5]), 120)
    for (k in 1:3) V <- V + exp(2 * p[k + 1]) * ZAZ[[k]]
    -mvtnorm::dmvnorm(d8$y, rep(p[1], 120), V, log = TRUE)
  }
  o <- optim(c(1, 0, 0, 0, 0), nll, method = "BFGS",
             control = list(reltol = 1e-14, maxit = 2000))
  o <- optim(o$par, nll, method = "Nelder-Mead",
             control = list(reltol = 1e-15, maxit = 5000))
  cat("exact ML logLik:", format(-o$value, digits = 15), "\n")
  cat("exact sds a b c, sigma:", format(exp(o$par[2:5]), digits = 8), "\n")
  cat("frmtmb - exact:", format(as.numeric(logLik(f1)) + o$value,
                                digits = 4), "\n")
  nd <- data.frame(g = factor(c("g1", "g5"), levels = levels(d8$g)),
                   f = factor(c("b", "c"), levels = c("a", "b", "c")),
                   fa = c(0, 0), fb = c(1, 0), fc = c(0, 1))
  go("predict newdata factor vs numeric", {
    p1 <- fitted(f1, newdata = nd)[, 1]
    p2 <- fitted(f2, newdata = nd)[, 1]
    max(abs(p1 - p2))
  })
  go("ranef dimnames", dimnames(ranef(f1)$g))
  go("variables sd", grep("^sd_", variables(f1), value = TRUE))
  go("summary", { print(summary(f1)); "ok" })
  go("simulate", dim(as.matrix(simulate(f1, nsim = 2, seed = 1))))
}
## gr(prec = Q)
Q <- solve(A)
g1 <- go("factor || prec", frm(y ~ 1 + (0 + f || gr(g, prec = Q)),
                               data = d8, data2 = list(Q = Q)))
g2 <- go("numeric || prec", frm(y ~ 1 + (0 + fa + fb + fc ||
                                           gr(g, prec = Q)),
                                data = d8, data2 = list(Q = Q)))
if (inherits(g1, "frmtmb_fit")) {
  print(c(as.numeric(logLik(g1)), as.numeric(logLik(g2)),
          as.numeric(logLik(f1))), digits = 15)
}
## an intercept and a factor: (1 + f || gr(g, cov = A))
h1 <- go("1 + f ||", frm(y ~ 1 + (1 + f || gr(g, cov = A)), data = d8,
                         data2 = list(A = A)))
h2 <- go("1 + fb + fc numeric ||", frm(y ~ 1 + (1 + fb + fc ||
                                                  gr(g, cov = A)),
                                       data = d8, data2 = list(A = A)))
if (inherits(h1, "frmtmb_fit")) {
  print(c(as.numeric(logLik(h1)), as.numeric(logLik(h2))), digits = 15)
}
cat("\nDONE\n")
