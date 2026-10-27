# Reviewer of lane loglik, fix round 1: slips of a jacobian slot that the
# lane's argument may not cover, how often the check runs, and its cost
# on a 1600-row sar(). Seed 101. Output: dev/review-loglik-log/fr1-seam2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
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
unfit <- function(fit) {
  class(fit) <- setdiff(class(fit), "frmtmb_unfitted")
  fit
}
report <- function(nm, expect, expr) {
  r <- tryCatch(expr, error = function(e) paste("REFUSED:",
                                                conditionMessage(e)))
  got <- if (is.character(r)) "REFUSED" else "ACCEPTED"
  cat(sprintf("%-44s expected %-8s got %-8s %s\n", nm, expect, got,
              if (is.character(r)) substr(r, 10, 90) else
                sprintf("rel err vs dense = %.3e", r)))
}

set.seed(101)
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
fixw <- function(y, mu, sigma, theta, aux) {
  list(u = as.vector(aux$Linv %*% ((y - mu) / sigma)),
       logjac = -aux$logdet - sum(log(rep_len(sigma, length(y)))))
}
fixu <- function(u, mu, sigma, theta, aux) mu + sigma * (aux$L %*% u)
right_fix <- function(sigma, theta, aux) {
  aux$Linv * rep(1 / rep_len(sigma, nrow(aux$Linv)), each = nrow(aux$Linv))
}
reg_fix <- function(name, jac) {
  frmtmb_register_autocor(name, args = "V", build = fixb, whiten = fixw,
                          unwhiten = fixu, sigma_varies = TRUE,
                          jacobian = jac)
}
# S1: right up to a per-row sign: K'K, so every log_lik() value, is
# unchanged
sgn <- rep(c(1, -1), length.out = n)
reg_fix("rvsign", function(sigma, theta, aux) sgn * right_fix(sigma, theta, aux))
# S2: the sigma of the first call (the start's, sigma = 1) kept for later
reg_fix("rvcache", local({
  s0 <- NULL
  function(sigma, theta, aux) {
    if (is.null(s0)) s0 <<- sigma
    right_fix(s0, theta, aux)
  }
}))
# S2b: the first row's sigma for every row
reg_fix("rvfirst", function(sigma, theta, aux) right_fix(sigma[1], theta, aux))
fix_case <- function(nm, expect) {
  report(paste(nm, "sigma ~ x"), expect, {
    fit <- unfit(frm(bf(stats::as.formula(paste0("y ~ x + ", nm, "(V)")),
                        sigma ~ x), data = d, data2 = list(V = V),
                     dry_run = "objective"))
    fit$estimates$beta[] <- c(1, 0.5)
    fit$estimates$betad[] <- c(log(0.7), 0.6)
    dp <- eval_dpars(fit)$y
    s <- as.numeric(dp$sigma)
    ref <- ll_dense_cond(d$y, as.numeric(dp$mu), outer(s, s) * V)
    rel(autocor_row_loglik(fit, eval_dpars(fit), "y"), ref)
  })
}
fix_case("rvsign", "either")
fix_case("rvcache", "REFUSED")
fix_case("rvfirst", "REFUSED")

# error-SAR-like maps over W (A (y - mu)), mean eta, rho = lo + (hi - lo) plogis(theta)
P <- matrix(0, n, n)
P[cbind(1:(n - 1), 2:n)] <- 1
P <- P + t(P)
Wr <- P / rowSums(P)                   # asymmetric, real eigenvalues
dg <- 1 / sqrt(rowSums(P))
Ws <- P * outer(dg, dg)                # symmetric
sb <- function(spec, rows, n_all, lookup) {
  W <- lookup(spec$W, "x()")[rows, rows, drop = FALSE]
  ev <- Re(eigen(W, only.values = TRUE)$values)
  list(W = W, ev = ev, lo = 1 / min(ev), hi = 1 / max(ev))
}
rho_of <- function(theta, aux) aux$lo + (aux$hi - aux$lo) / (1 + exp(-theta[1]))
sw <- function(y, mu, sigma, theta, aux) {
  A <- diag(length(y)) - rho_of(theta, aux) * aux$W
  list(u = as.vector(A %*% (y - mu)) / sigma,
       logjac = sum(log1p(-rho_of(theta, aux) * aux$ev)) -
         length(y) * log(sigma))
}
su <- function(u, mu, sigma, theta, aux) {
  A <- diag(nrow(aux$W)) - rho_of(theta, aux) * aux$W
  mu + sigma * solve(A, u)
}
reg_s <- function(name, jac, st = function(aux) 0) {
  frmtmb_register_autocor(
    name, args = "W", build = sb, whiten = sw, unwhiten = su, npar = 1L,
    start = st, natural = function(theta, aux) c(rho = rho_of(theta, aux)),
    bounds = function(aux) cbind(aux$lo, aux$hi), sigma_varies = FALSE,
    jacobian = jac)
}
right_s <- function(sigma, theta, aux) {
  (diag(nrow(aux$W)) - rho_of(theta, aux) * aux$W) / sigma
}
# S3: right for a symmetric W only (symmetrizes W)
reg_s("rvsym", function(sigma, theta, aux) {
  (diag(nrow(aux$W)) - rho_of(theta, aux) * (aux$W + t(aux$W)) / 2) / sigma
})
# S4: the transpose, which on a symmetric W is the same matrix
reg_s("rvtsym", function(sigma, theta, aux) t(right_s(sigma, theta, aux)))
reg_s("rvok", right_s)
# S5: rho frozen at the start's value, with the start at logit 0, 8, 14
for (s0 in c(0, 8, 14)) {
  reg_s(paste0("rvfrozen", s0), local({
    t0 <- s0
    function(sigma, theta, aux) right_s(sigma, t0, aux)
  }), st = local({
    t0 <- s0
    function(aux) t0
  }))
}
sar_case <- function(nm, W, expect, rho = 0.3) {
  report(paste(nm, if (isSymmetric(W)) "symmetric W" else "asymmetric W"),
         expect, {
    fit <- unfit(frm(stats::as.formula(paste0("y ~ x + ", nm, "(W)")),
                     data = d, data2 = list(W = W), dry_run = "objective"))
    aux <- fit$frame$autocor$y$aux
    fit$estimates$beta[] <- c(1, 0.5)
    fit$estimates$betad[] <- log(0.9)
    fit$estimates$thetaac[] <- qlogis((rho - aux$lo) / (aux$hi - aux$lo))
    A <- diag(n) - rho * W
    eta <- 1 + 0.5 * d$x
    ref <- ll_dense_cond(d$y, eta,
                         0.81 * solve(crossprod(A)))
    rel(autocor_row_loglik(fit, eval_dpars(fit), "y"), ref)
  })
}
sar_case("rvsym", Wr, "REFUSED")
sar_case("rvsym", Ws, "ACCEPTED")
sar_case("rvtsym", Ws, "ACCEPTED")
sar_case("rvtsym", Wr, "REFUSED")
sar_case("rvok", Wr, "ACCEPTED")
for (s0 in c(0, 8, 14)) sar_case(paste0("rvfrozen", s0), Wr, "REFUSED")

# S6: a one-sided parameter, rho = theta on (-1, 0], start 0, a map that
# stops for theta > 0, and a slot that ignores rho
rho1 <- function(theta) {
  if (!inherits(theta, "advector") && theta[1] > 0) stop("rho must be <= 0")
  theta[1]
}
frmtmb_register_autocor(
  "rvoneside", args = "W", build = sb,
  whiten = function(y, mu, sigma, theta, aux) {
    A <- diag(length(y)) - rho1(theta) * aux$W
    list(u = as.vector(A %*% (y - mu)) / sigma,
         logjac = sum(log1p(-rho1(theta) * aux$ev)) - length(y) * log(sigma))
  },
  unwhiten = function(u, mu, sigma, theta, aux) {
    A <- diag(nrow(aux$W)) - rho1(theta) * aux$W
    mu + sigma * solve(A, u)
  },
  npar = 1L, start = function(aux) 0,
  natural = function(theta, aux) c(rho = theta[1]), sigma_varies = FALSE,
  jacobian = function(sigma, theta, aux) diag(nrow(aux$W)) / sigma)
report("rvoneside (ignores rho) asymmetric W", "REFUSED", {
  fit <- unfit(frm(y ~ x + rvoneside(W), data = d, data2 = list(W = Wr),
                   dry_run = "objective"))
  fit$estimates$beta[] <- c(1, 0.5)
  fit$estimates$betad[] <- log(0.9)
  fit$estimates$thetaac[] <- -0.5
  A <- diag(n) + 0.5 * Wr
  ref <- ll_dense_cond(d$y, 1 + 0.5 * d$x,
                       0.81 * solve(crossprod(A)))
  rel(autocor_row_loglik(fit, eval_dpars(fit), "y"), ref)
})

# how often the slot is called: at frm() and per draw of log_lik()
calls <- 0L
reg_s("rvcount", function(sigma, theta, aux) {
  calls <<- calls + 1L
  right_s(sigma, theta, aux)
})
fit <- frm(y ~ x + rvcount(W), data = d, data2 = list(W = Wr),
           dry_run = "objective")
cat("slot calls during frm():", calls, "\n")
calls <- 0L
lab <- c(brms_par_labels(fit), "lp__")
K <- 5L
m <- matrix(0, K, length(lab), dimnames = list(NULL, lab))
m[, "b_Intercept"] <- 1
m[, "sigma"] <- 1
m[, "thetaac_1"] <- seq(-0.5, 0.5, length.out = K)
ds <- structure(list(stanfit = NULL, draws = m, fit = fit),
                class = "frmtmb_draws")
invisible(log_lik(ds))
cat("slot calls during log_lik() of", K, "draws:", calls, "\n")

# the check's cost on a 1600-row sar(), and frm()'s build time
k <- 40
nn <- as.integer(k * k)
id <- matrix(seq_len(nn), k, k)
Pn <- matrix(0, nn, nn)
Pn[cbind(as.vector(id[-k, ]), as.vector(id[-1, ]))] <- 1
Pn[cbind(as.vector(id[, -k]), as.vector(id[, -1]))] <- 1
Pn <- Pn + t(Pn)
W16 <- Pn / rowSums(Pn)
d16 <- data.frame(x = rnorm(nn))
d16$y <- 1 + rnorm(nn)
tb <- system.time(f16 <- frm(y ~ x + sar(W16), data = d16,
                             data2 = list(W16 = W16),
                             dry_run = "objective"))[["elapsed"]]
ac <- f16$frame$autocor$y
e <- frmtmb:::ac_block_entry(ac)
st <- e$start(ac$aux)
bad <- function(what) stop(what)
tc <- min(vapply(1:5, function(i) system.time(
  frmtmb:::autocor_seam_check_jacobian(e, ac, ac$aux, st, nn, bad))[["elapsed"]], 0))
tm <- min(vapply(1:5, function(i) system.time(
  frmtmb:::autocor_seam_check_map(e, ac, ac$aux, st, nn))[["elapsed"]], 0))
cat(sprintf("sar n = %d: frm(dry_run) %.2f s; check_jacobian %.3f s; whole check_map %.3f s\n",
            nn, tb, tc, tm))
cat("DONE\n")
