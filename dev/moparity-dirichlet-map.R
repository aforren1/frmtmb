# A MAP fit under a simo prior against a hand-written log posterior:
# the objective at the fit's optimum, and the optimum of the hand
# function found by optim() from another start.
# Usage: Rscript dev/moparity-dirichlet-map.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args) && nzchar(args[[1]])) args[[1]] else
  "C:/Users/adf44/source/r/wt-moparity-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(3)
n <- 120
dm <- data.frame(inc = sample(0:3, n, TRUE), z = rnorm(n))
dm$y <- 1 + 0.6 * c(0, 1, 1.4, 1.5)[dm$inc + 1] + 0.3 * dm$z + rnorm(n)
a <- c(2, 1.5, 3)
pr <- set_prior("dirichlet(c(2, 1.5, 3))", class = "simo", coef = "moinc1")
fit <- frm(bf(y ~ mo(inc) + z), family = gaussian(), data = dm, prior = pr)
est <- fit$estimates
cat("estimates:", format(unlist(est), digits = 8), "\n")
# the hand log posterior, in the fit's parameters
hand <- function(beta, logsig, zeta) {
  w <- frmtmb:::mo_simplex(zeta)
  mo <- 3 * c(0, cumsum(w))[dm$inc + 1]
  mu <- beta[1] + beta[2] * dm$z + beta[3] * mo
  ll <- sum(dnorm(dm$y, mu, exp(logsig), log = TRUE))
  # stereographic sphere point, as mo_simplex() builds it
  D <- 3
  P <- -1 / sqrt(D)
  r2 <- sum(zeta^2)
  u <- (2 * as.vector(frmtmb:::mo_sphere_basis(D) %*% zeta) +
          (r2 - 1) * P) / (r2 + 1)
  lp <- lgamma(sum(a)) - sum(lgamma(a)) + sum((a - 1) * log(w)) +
    sum(log(abs(u))) + (D - 1) * log(4) - (D - 1) * log(1 + r2) -
    D * log(2)
  ll + lp
}
pv <- function(p) hand(p[1:3], p[4], p[5:6])
p_fit <- c(est$beta, est$betad, est$zeta1)
h_at <- pv(p_fit)
obj_at <- -fit$obj$fn(fit$opt$par)
cat(sprintf("MAP objective -fn %.12f | hand log posterior %.12f | gap %.3g\n",
            obj_at, h_at, obj_at - h_at))
cat("max |gradient| at the fit:", max(abs(fit$obj$gr(fit$opt$par))), "\n")
set.seed(9)
best <- -Inf
best_p <- NULL
for (k in 1:10) {
  st <- p_fit + rnorm(6, 0, 0.3)
  o1 <- optim(st, function(p) -pv(p), control = list(maxit = 20000,
                                                     reltol = 1e-14))
  o2 <- optim(o1$par, function(p) -pv(p), method = "BFGS",
              control = list(reltol = 1e-15, maxit = 5000))
  if (-o2$value > best) {
    best <- -o2$value
    best_p <- o2$par
  }
}
cat(sprintf("optim of the hand function, best of 10 starts: value %.12f (fit %.12f, fit - best %.3g), max |param - fit| %.3g, weights fit %s | optim %s\n",
            best, h_at, h_at - best, max(abs(best_p - p_fit)),
            paste(format(frmtmb:::mo_simplex(p_fit[5:6]), digits = 6),
                  collapse = " "),
            paste(format(frmtmb:::mo_simplex(best_p[5:6]), digits = 6),
                  collapse = " ")))
# the ML fit is not this one: the prior moves the weights
ml <- frm(bf(y ~ mo(inc) + z), family = gaussian(), data = dm)
cat("ML weights:", format(frmtmb:::mo_simplex(ml$estimates$zeta1),
                           digits = 6), "\n")
# a concentration below 1/2 has no MAP
r <- tryCatch(frm(bf(y ~ mo(inc) + z), family = gaussian(), data = dm,
                  prior = set_prior("dirichlet(0.4)", class = "simo",
                                    coef = "moinc1")),
              error = function(e) conditionMessage(e))
cat("alpha 0.4:", substr(r, 1, 90), "\n")
