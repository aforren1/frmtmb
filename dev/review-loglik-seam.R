# Reviewer of lane loglik: wrong jacobian slots the frame-build check
# may pass, and the cost of the dense fallback by n. Seed 81.
# Output: dev/review-loglik-log/seam.txt
LIB <- Sys.getenv("RV_LIB", "C:/Users/adf44/source/r/wt-loglik-lib")
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")

ll_dense_cond <- function(y, m, S) {
  n <- length(y)
  vapply(seq_len(n), function(i) {
    o <- seq_len(n)[-i]
    A <- solve(S[o, o, drop = FALSE])
    r <- y[o] - m[o]
    cm <- m[i] + as.numeric(S[i, o, drop = FALSE] %*% A %*% r)
    cv <- S[i, i] - as.numeric(S[i, o, drop = FALSE] %*% A %*% S[o, i])
    stats::dnorm(y[i], cm, sqrt(cv), log = TRUE)
  }, 0)
}
rel <- function(a, b) max(abs(a - b)) / max(abs(b))
report <- function(nm, expr) {
  r <- tryCatch(expr, error = function(e) paste("REFUSED:",
                                                conditionMessage(e)))
  cat(sprintf("%-46s %s\n", nm, if (is.character(r)) substr(r, 1, 160)
              else sprintf("ACCEPTED, rel err vs dense = %.3e", r)))
}

set.seed(81)
n <- 40
d <- data.frame(x = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(n)
V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)

fixb <- function(spec, rows, n_all, lookup) {
  V <- lookup(spec$V, "x()")[rows, rows, drop = FALSE]
  L <- t(chol(V))
  list(L = L, Linv = forwardsolve(L, diag(nrow(L))),
       logdet = sum(log(diag(L))))
}
# sigma_varies = TRUE, as fcor(): sigma is one value per row
fixw <- function(y, mu, sigma, theta, aux) {
  list(u = as.vector(aux$Linv %*% ((y - mu) / sigma)),
       logjac = -aux$logdet - sum(log(rep_len(sigma, length(y)))))
}
fixu <- function(u, mu, sigma, theta, aux) mu + sigma * (aux$L %*% u)

# (a) forgets to divide by sigma: right at the check's sigma = 1
frmtmb_register_autocor("rvnosig", args = "V", build = fixb, whiten = fixw,
                        unwhiten = fixu, sigma_varies = TRUE,
                        jacobian = function(sigma, theta, aux) aux$Linv)
# (b) the help page's example `Linv / sigma` with sigma_varies = TRUE:
# divides ROWS by sigma, where K = Linv D^-1 divides COLUMNS
frmtmb_register_autocor("rvrows", args = "V", build = fixb, whiten = fixw,
                        unwhiten = fixu, sigma_varies = TRUE,
                        jacobian = function(sigma, theta, aux) {
                          aux$Linv / sigma
                        })
# (c) a theta-dependent map, rho = tanh(theta) with start 0, whose
# jacobian ignores rho: right at the start only
P <- matrix(0, n, n)
P[cbind(1:(n - 1), 2:n)] <- 1
P <- P + t(P)
Wr <- P / rowSums(P)
sb <- function(spec, rows, n_all, lookup) {
  W <- lookup(spec$W, "x()")[rows, rows, drop = FALSE]
  list(W = W)
}
sw <- function(y, mu, sigma, theta, aux) {
  rho <- tanh(theta[1])
  A <- diag(length(y)) - rho * aux$W
  list(u = as.vector(A %*% (y - mu)) / sigma,
       logjac = as.numeric(determinant(A)$modulus) - length(y) * log(sigma))
}
su <- function(u, mu, sigma, theta, aux) {
  rho <- tanh(theta[1])
  A <- diag(nrow(aux$W)) - rho * aux$W
  mu + sigma * solve(A, u)
}
reg_s <- function(name, jac) {
  frmtmb_register_autocor(name, args = "W", build = sb, whiten = sw,
                          unwhiten = su, npar = 1L,
                          start = function(aux) 0, sigma_varies = FALSE,
                          natural = function(theta, aux) tanh(theta),
                          jacobian = jac)
}
reg_s("rvnorho", function(sigma, theta, aux) diag(nrow(aux$W)) / sigma)
# (d) the transpose: K' gives a different Q = K K' on an asymmetric W
reg_s("rvtrans", function(sigma, theta, aux) {
  t(diag(nrow(aux$W)) - tanh(theta[1]) * aux$W) / sigma
})
# (e) the right one, the guard's absent case
reg_s("rvright", function(sigma, theta, aux) {
  (diag(nrow(aux$W)) - tanh(theta[1]) * aux$W) / sigma
})

# reference for the fixed-matrix terms: D V D with sigma ~ x
for (nm in c("rvnosig", "rvrows")) {
  report(paste(nm, "sigma ~ x"), {
    fit <- frm(bf(stats::as.formula(paste0("y ~ x + ", nm, "(V)")),
                  sigma ~ x), data = d, data2 = list(V = V),
               dry_run = "objective")
    class(fit) <- setdiff(class(fit), "frmtmb_unfitted")
    fit$estimates$beta[] <- c(1, 0.5)
    fit$estimates$betad[] <- c(log(0.7), 0.6)
    dp <- eval_dpars(fit)$y
    s <- as.numeric(dp$sigma)
    ref <- ll_dense_cond(d$y, as.numeric(dp$mu), outer(s, s) * V)
    rel(autocor_row_loglik(fit, eval_dpars(fit), "y"), ref)
  })
  report(paste(nm, "constant sigma 1.7"), {
    fit <- frm(stats::as.formula(paste0("y ~ x + ", nm, "(V)")), data = d,
               data2 = list(V = V), dry_run = "objective")
    class(fit) <- setdiff(class(fit), "frmtmb_unfitted")
    fit$estimates$beta[] <- c(1, 0.5)
    fit$estimates$betad[] <- log(1.7)
    ref <- ll_dense_cond(d$y, 1 + 0.5 * d$x, 1.7^2 * V)
    rel(autocor_row_loglik(fit, eval_dpars(fit), "y"), ref)
  })
}
for (nm in c("rvnorho", "rvtrans", "rvright")) {
  report(paste(nm, "rho = tanh(0.8)"), {
    fit <- frm(stats::as.formula(paste0("y ~ x + ", nm, "(W)")), data = d,
               data2 = list(W = Wr), dry_run = "objective")
    class(fit) <- setdiff(class(fit), "frmtmb_unfitted")
    fit$estimates$beta[] <- c(1, 0.5)
    fit$estimates$betad[] <- log(0.9)
    fit$estimates$thetaac[] <- 0.8
    A <- diag(n) - tanh(0.8) * Wr
    S <- 0.9^2 * solve(crossprod(A))
    ref <- ll_dense_cond(d$y, 1 + 0.5 * d$x, S)
    rel(autocor_row_loglik(fit, eval_dpars(fit), "y"), ref)
  })
}

# the dense fallback by n: time and peak allocation of one call
cat("\nfallback cost (one autocor_row_loglik() call, no jacobian slot)\n")
frmtmb_register_autocor("rvfb", args = "W", build = sb, whiten = sw,
                        unwhiten = su, npar = 1L, start = function(aux) 0,
                        natural = function(theta, aux) tanh(theta),
                        sigma_varies = FALSE)
for (nn in c(250, 500, 1000, 2000)) {
  set.seed(nn)
  dd <- data.frame(x = rnorm(nn))
  dd$y <- 1 + rnorm(nn)
  Pn <- matrix(0, nn, nn)
  Pn[cbind(1:(nn - 1), 2:nn)] <- 1
  Pn <- Pn + t(Pn)
  Wn <- Pn / rowSums(Pn)
  t_build <- system.time(
    fit <- frm(y ~ x + rvfb(W), data = dd, data2 = list(W = Wn),
               dry_run = "objective"))[["elapsed"]]
  class(fit) <- setdiff(class(fit), "frmtmb_unfitted")
  fit$estimates$thetaac[] <- 0.5
  dp <- eval_dpars(fit)
  gc(reset = TRUE)
  m0 <- sum(gc()[, 2])
  tt <- system.time(v <- autocor_row_loglik(fit, dp, "y"))[["elapsed"]]
  m1 <- sum(gc()[, 6])
  cat(sprintf("n = %5d  build %.2f s  call %.3f s  peak Mb %.0f  (8 n^2 = %.0f Mb)  finite %s\n",
              nn, t_build, tt, m1 - m0, 8 * nn^2 / 2^20, all(is.finite(v))))
}
cat("DONE\n")
