# The new numeric assertions of test-refit-plumbing.R over 20 ulp-scale
# perturbations of the response, y * (1 + 2e-15 * k), k = 1..20.
source("dev/refits-lib.R")
## 2939: factor || gr(cov) against the numeric spelling and the exact
## marginal density
set.seed(11)
ng <- 12
A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
d <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                           levels = paste0("g", 1:ng)),
                f = factor(rep(c("a", "b", "c"), 40)))
L <- t(chol(A))
r <- cbind(0.8 * L %*% rnorm(ng), 0.6 * L %*% rnorm(ng),
           1.2 * L %*% rnorm(ng))
gi <- as.integer(d$g)
y0 <- 1 + r[cbind(gi, as.integer(d$f))] + rnorm(120, 0, 0.5)
d$fa <- as.numeric(d$f == "a"); d$fb <- as.numeric(d$f == "b")
d$fc <- as.numeric(d$f == "c")
out <- NULL
for (k in 1:20) {
  d$y <- y0 * (1 + 2e-15 * k)
  f1 <- frm(y ~ 1 + (0 + f || gr(g, cov = A)), data = d, data2 = list(A = A))
  f2 <- frm(y ~ 1 + (0 + fa + fb + fc || gr(g, cov = A)), data = d,
            data2 = list(A = A))
  h1 <- frm(y ~ 1 + (1 + f || gr(g, cov = A)), data = d, data2 = list(A = A))
  h2 <- frm(y ~ 1 + (1 + fb + fc || gr(g, cov = A)), data = d,
            data2 = list(A = A))
  g1 <- frm(y ~ 1 + (0 + f || gr(g, prec = Q)), data = d,
            data2 = list(Q = solve(A)))
  vc <- VarCorr(f1); sds <- vc$g$sd[, "Estimate"]
  V <- diag(vc$residual__$sd[, "Estimate"]^2, 120)
  for (j in 1:3) {
    Z <- matrix(0, 120, ng); Z[cbind(1:120, gi)] <- as.integer(d$f) == j
    V <- V + sds[j]^2 * Z %*% A %*% t(Z)
  }
  R <- chol(V)
  e <- backsolve(R, d$y - fixef(f1)[1, 1], transpose = TRUE)
  ex <- -sum(log(diag(R))) - 0.5 * sum(e^2) - 60 * log(2 * pi)
  rel <- function(a, b) abs(as.numeric(a) - as.numeric(b)) / abs(as.numeric(b))
  out <- rbind(out, c(k = k, f12 = rel(logLik(f1), logLik(f2)),
                      exact = rel(logLik(f1), ex),
                      h12 = rel(logLik(h1), logLik(h2)),
                      prec = rel(logLik(g1), logLik(f1))))
}
print(signif(out, 3))
cat("max rel: f12", max(out[, "f12"]), " exact", max(out[, "exact"]),
    " h12", max(out[, "h12"]), " prec", max(out[, "prec"]), "\n")
cat("test tolerances: 1e-10, 1e-9, 1e-10, 1e-8\n")

## 3434: influence() deletion against the written-out basis
set.seed(8)
dl <- data.frame(x = runif(60, 0, 3), z = rnorm(60),
                 blk = factor(rep(1:5, 12)))
yl <- sin(dl$x) + 0.3 * dl$z + rnorm(60, 0, 0.3)
rat <- numeric(20)
for (k in 1:20) {
  dl$y <- yl * (1 + 2e-15 * k)
  fl <- frm(y ~ z + s(x, k = 5), data = dl)
  inf <- influence(fl, groups = "blk", data = dl)
  si <- fl$frame$linpreds[[1L]]$smooths[[1L]]
  keep <- dl[dl$blk != "1", ]
  M <- frmtmb:::smooth_basis_at(si, keep)
  nr <- sum(si$nr)
  for (j in seq_len(nr)) keep[[paste0("r", j)]] <- M[, j]
  keep$fx <- M[, nr + 1L]
  keep$id <- factor(1)
  fo <- as.formula(paste0("y ~ z + fx + homdiag(0 + ",
                          paste0("r", seq_len(nr), collapse = " + "),
                          " | id)"))
  ref <- fixef(frm(fo, data = keep))[, 1]
  got <- inf$fixed["1", c("z", "s(x).fx1")]
  moved <- abs(got - inf$fixed_full[c("z", "s(x).fx1")])
  rat[k] <- max(abs(got - ref[c("z", "fx")]) / moved)
}
print(signif(rat, 3))
cat("influence ratio max", max(rat), "(test bound 1e-3)\n")
