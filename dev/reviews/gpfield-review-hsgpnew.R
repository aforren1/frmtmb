# Review of lane gpfield, 1919: an approximate field at new positions
# inside and outside its basis's boundary. fitted() at new x against the
# basis written here from brms's formula (brms:::eigen_fun_laplacian:
# phi_j(x) = L^-1/2 sin(sqrt(lambda_j) (x + L)), x centered on brms's
# scaled training positions), and predict()'s spread there.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
set.seed(9200)
n <- 30
d <- data.frame(x = sort(stats::runif(n, 0, 10)))
d$y <- 2 * sin(d$x) + stats::rnorm(n, 0, 0.2)
q <- gf_fit(bf(y ~ gp(x, k = 15)), data = d, family = gaussian())
fit <- q$fit
cat("warn:", unique(sub(":.*", "", q$warn)), "\n")
bk <- Filter(function(b) b$covstruct == "hsgp", fit$frame$re_blocks)[[1]]
dmax <- bk$gp_dmax
ux <- unique(d$x)
ctr <- mean(ux / dmax)
L <- pi * bk$aux_omega[1, 1]^-1 / 2
lam <- bk$aux_omega[, 1]
phi <- function(x) {
  xs <- x / dmax - ctr
  outer(xs, lam, function(a, w) sin(w * (a + L)) / sqrt(L))
}
lp <- fit$frame$linpreds[[asNamespace("frmtmb")$gp_block_lp(bk, fit$frame)]]; Ztr <- as.matrix(asNamespace("frmtmb")$lp_Z(fit, lp))[, bk$c_idx]
cat("basis at training rows, max |mine - frmtmb|:",
    max(abs(phi(d$x) - Ztr)), " L", L, "boundary in x units: [",
    (ctr - L) * dmax, ",", (ctr + L) * dmax, "]\n")
xs <- c(5, 10.5, 14, 16.5, 20, 30)
b <- fit$estimates$b[bk$b_idx]
mine <- fit$estimates$beta[1] + as.vector(phi(xs) %*% b)
ft <- fitted(fit, newdata = data.frame(x = xs))[, 1]
set.seed(1)
p <- predict(fit, newdata = data.frame(x = xs), summary = FALSE, ndraws = 2000)
for (j in seq_along(xs)) {
  cat(sprintf("HSGPNEW x %5.1f | fitted %.5f mine %.5f | predict mean %.3f sd %.3f\n",
              xs[j], ft[j], mine[j], mean(p[, j]), sd(p[, j])))
}
