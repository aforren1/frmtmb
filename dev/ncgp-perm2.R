# How reliably does the factor's order move eta's standard error? dev/ncgp-perm.R read
# one design per size. Here five seeds at n 40, 100, 200 and two noises, new positions,
# sorted against farthest-first, one process per arm.
#   Rscript dev/ncgp-perm2.R <sorted|perm>
arm <- commandArgs(TRUE)[1]
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "arm", arm, "\n")
if (arm == "perm") source("dev/ncgp-perm-patch.R")
ex <- parse("dev/ncgp-truth.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]])[1] %in% c("marginal", "closed_max", "REL")) eval(e)
}
for (n in c(40L, 100L, 200L)) for (seed in 1:5) for (s in c(1e-3, 3e-4)) {
  set.seed(seed)
  x <- round(stats::runif(n, 0, 6), 3)
  y <- sin(2 * x) + s * stats::rnorm(n)
  d <- data.frame(x = x, y = y)
  fit <- suppressWarnings(frm(bf(y ~ gp(x)), data = d, family = gaussian()))
  u <- sort(unique(x))
  J <- match(x, u)
  D2 <- outer(u, u, "-")^2
  th <- fit$estimates$theta
  pf <- c(unname(fit$estimates$beta[1]), th[1], th[2], unname(fit$estimates$betad[1]))
  S <- solve(numDeriv::hessian(function(p) marginal(p, y, D2, J), pf))
  nd <- data.frame(x = seq(min(x) + 0.13, max(x) - 0.17, length.out = 25L))
  cf_eta <- function(p) {
    K <- exp(2 * p[2]) * exp(-D2 / (2 * exp(2 * p[3])))
    ks <- exp(2 * p[2]) * exp(-outer(nd$x, u, "-")^2 / (2 * exp(2 * p[3])))
    V <- K[J, J, drop = FALSE]
    diag(V) <- diag(V) * (1 + 1e-12) + exp(2 * p[4])
    R <- chol(V)
    W <- backsolve(R, t(ks[, J, drop = FALSE]), transpose = TRUE)
    list(m = p[1] + as.vector(t(W) %*% backsolve(R, y - p[1], transpose = TRUE)),
         v = exp(2 * p[2]) - colSums(W^2))
  }
  ce <- cf_eta(pf)
  G <- numDeriv::jacobian(function(p) cf_eta(p)$m, pf)
  ref <- sqrt(pmax(ce$v, 0) + rowSums((G %*% S) * G))
  fv <- suppressWarnings(stats::fitted(fit, newdata = nd))
  r <- fv[, 2L] / ref
  cat(sprintf("NCGPQ %-6s n %3d seed %d s %.0e conv %d se/cf min %.4f median %.4f max %.4f | est err/cf se max %.2e\n",
              arm, n, seed, s, fit$opt$convergence, min(r), stats::median(r), max(r),
              max(abs(fv[, 1L] - ce$m) / ref)))
}
