# Item 7: a reference where no posterior read exists from either start. Field integrated
# out by importance sampling from a t proposal at the exact Laplace mode, in 256-bit linear
# algebra on the kernel, n small. Gives (1) the exact marginal likelihood against frm()'s
# Laplace objective, (2) a grid-polished maximum, (3) a pseudo-marginal posterior at brms's
# default priors against frm()'s MAP. Seed 7100 + design index. Args: case
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
case <- args[1]
nm <- 8L
set.seed(7100 + match(case, c("pois16", "pois12", "binom7", "binom5", "pois4", "bern")))
x <- sort(runif(nm, 0, 6))
fam <- switch(case, pois16 = "pois", pois12 = "pois", pois4 = "pois", binom7 = "binom", binom5 = "binom", bern = "binom")
lm0 <- switch(case, pois16 = 16, pois12 = 12, pois4 = 4, binom7 = 0, binom5 = 0, bern = 0)
tr <- switch(case, pois16 = NA, pois12 = NA, pois4 = NA, binom7 = 1e7, binom5 = 1e5, bern = 20)
eta_true <- lm0 + sin(x)
y <- if (fam == "pois") rpois(nm, exp(eta_true)) else rbinom(nm, tr, plogis(eta_true))
if (fam == "pois") { d <- data.frame(y = y, x = x); fml <- y ~ gp(x, scale = FALSE); fm <- poisson() } else { d <- data.frame(y = y, x = x, tr = tr); fml <- y | trials(tr) ~ gp(x, scale = FALSE); fm <- binomial() }
r <- rev_fit(fml, d, family = fm)
g <- rev_gp_read(r$fit)
cat(sprintf("RNG %s y %s\n", case, paste(y, collapse = " ")))
cat(sprintf("RNG frm ll %.6f mu %.5f sd %.5f rho %.5f conv %d w %d\n", g$ll, g$beta[1], g$sd, g$rho, g$conv, length(r$warn)))
loglik <- function(eta) if (fam == "pois") sum(dpois(y, exp(eta), log = TRUE)) else sum(dbinom(y, tr, plogis(eta), log = TRUE))
w_of <- function(eta) if (fam == "pois") exp(eta) else tr * plogis(eta) * (1 - plogis(eta))
g_of <- function(eta) if (fam == "pois") y - exp(eta) else y - tr * plogis(eta)
D2 <- outer(x, x, "-")^2
# Exact reference at (mu, sd, rho): eigen of K at 256 bits is overkill at n 8 with nugget 1e-12; a
# Cholesky in double holds to 1e-4 of the pivot's scale, and everything below is in eta's coordinates.
is_ll <- function(mu, sd, rho, zs, tdf = 5) {
  K <- sd^2 * exp(-D2 / (2 * rho^2)); diag(K) <- diag(K) * (1 + 1e-12)
  Ki <- chol2inv(chol(K))
  eta <- rep(mu, nm) + 0
  # Newton from the observed values
  eta <- if (fam == "pois") log(pmax(y, 0.5)) else qlogis(pmin(pmax(y / tr, 1e-8), 1 - 1e-8))
  for (it in 1:60) {
    W <- w_of(eta)
    H <- Ki; diag(H) <- diag(H) + W
    gr <- g_of(eta) - as.vector(Ki %*% (eta - mu))
    st <- solve(H, gr)
    eta <- eta + st
    if (max(abs(st)) < 1e-12) break
  }
  W <- w_of(eta); H <- Ki; diag(H) <- diag(H) + W
  Rh <- chol(H)
  lpK <- function(e) -0.5 * sum((e - mu) * (Ki %*% (e - mu))) - sum(log(diag(chol(K)))) - 0.5 * nm * log(2 * pi)
  lap <- loglik(eta) + lpK(eta) + 0.5 * nm * log(2 * pi) - sum(log(diag(Rh)))
  # t proposal around the mode, shape H^-1
  draws <- backsolve(Rh, t(zs))  # nm x S
  S <- ncol(draws)
  chi <- attr(zs, "chi")
  scale <- sqrt(tdf / chi)
  e <- eta + sweep(draws, 2, scale, "*")
  lt <- lgamma((tdf + nm) / 2) - lgamma(tdf / 2) - nm / 2 * log(tdf * pi) + sum(log(diag(Rh))) -
    (tdf + nm) / 2 * log1p(rowSums(zs^2) / chi)
  lw <- vapply(seq_len(S), function(i) loglik(e[, i]) + lpK(e[, i]), 0) - lt
  a <- max(lw)
  ll <- a + log(mean(exp(lw - a)))
  ess <- sum(exp(lw - a))^2 / sum(exp(lw - a)^2)
  list(ll = ll, lap = lap, ess = ess, sdw = sd(exp(lw - a)) / mean(exp(lw - a)))
}
mkz <- function(S, tdf = 5) { z <- matrix(rnorm(S * nm), S, nm); attr(z, "chi") <- rchisq(S, tdf); z }
set.seed(1)
zs <- mkz(4000)
a <- is_ll(g$beta[1], g$sd, g$rho, zs)
cat(sprintf("RNG at frm's point: IS exact ll %.6f Laplace by hand %.6f frm %.6f | IS-frm %.2e hand-frm %.2e ESS %.0f of 4000\n", a$ll, a$lap, g$ll, a$ll - g$ll, a$lap - g$ll, a$ess))
# (2) the maximum of the exact marginal with common random numbers
fo <- function(p) { v <- tryCatch(is_ll(p[1], exp(p[2]), exp(p[3]), zs)$ll, error = function(e) NA); if (is.finite(v)) -v else 1e10 }
best <- NULL
for (st in list(c(g$beta[1], log(g$sd), log(g$rho)), c(mean(log(pmax(y, 1))), log(1), log(1)), c(g$beta[1], log(g$sd) + 1, log(g$rho) + 0.7), c(g$beta[1], log(g$sd) - 1, log(g$rho) - 0.7))) {
  o <- optim(st, fo, method = "Nelder-Mead", control = list(maxit = 250, reltol = 1e-8))
  if (is.null(best) || o$value < best$value) best <- o
}
cat(sprintf("RNG IS max ll %.5f mu %.4f sd %.4f rho %.4f | frm's point on the same draws %.5f | short by %.3e\n", -best$value, best$par[1], exp(best$par[2]), exp(best$par[3]), a$ll, a$ll + best$value))
# (3) pseudo-marginal Metropolis at frm_sample's default priors, read from default_prior at MAP
pri <- function(p) dnorm(p[1], 0, 10, log = TRUE) + (dt(exp(p[2]) / 2.5, 3, log = TRUE) + p[2]) + (-1.5 * log(1) + invgamma_lp(exp(p[3])) + p[3])
invgamma_lp <- function(z) { a0 <- 1.5; b0 <- 1; a0 * log(b0) - lgamma(a0) - (a0 + 1) * log(z) - b0 / z }
map <- rev_fit(fml, d, family = fm, prior = c(set_prior("normal(0, 10)", class = "Intercept"), set_prior("student_t(3, 0, 2.5)", class = "sdgp"), set_prior("inv_gamma(1.5, 1)", class = "lscale")))
if (!inherits(map$fit, "error")) {
  gm <- rev_gp_read(map$fit)
  cat(sprintf("RNG frm MAP mu %.5f sd %.5f rho %.5f conv %d\n", gm$beta[1], gm$sd, gm$rho, gm$conv))
}
S2 <- 600
cur <- c(g$beta[1], log(g$sd), log(g$rho))
Hs <- tryCatch(solve(optimHess(cur, fo)), error = function(e) diag(0.05, 3))
if (any(eigen(Hs)$values <= 0)) Hs <- diag(0.05, 3)
Lp <- t(chol(Hs * 2.38^2 / 3))
curv <- -fo(cur) + pri(cur)
set.seed(2)
chain <- matrix(NA, 4000, 3)
acc <- 0
for (i in 1:4000) {
  zz <- mkz(S2)
  if (i %% 25 == 1) curv <- { v <- tryCatch(is_ll(cur[1], exp(cur[2]), exp(cur[3]), zz)$ll, error = function(e) -Inf); v + pri(cur) }
  pr <- cur + as.vector(Lp %*% rnorm(3))
  v <- tryCatch(is_ll(pr[1], exp(pr[2]), exp(pr[3]), zz)$ll, error = function(e) -Inf)
  lp <- v + pri(pr)
  if (is.finite(lp) && log(runif(1)) < lp - curv) { cur <- pr; curv <- lp; acc <- acc + 1 }
  chain[i, ] <- cur
}
ch <- chain[-(1:1000), ]
q <- apply(ch, 2, quantile, c(.05, .5, .95))
cat(sprintf("RNG pseudo-marginal accept %.2f mu %.4f %.4f %.4f | sd %.4f %.4f %.4f | rho %.4f %.4f %.4f | post sd of (mu, log sd, log rho) %.4f %.4f %.4f\n", acc / 4000, q[1, 1], q[2, 1], q[3, 1],
            exp(q[1, 2]), exp(q[2, 2]), exp(q[3, 2]), exp(q[1, 3]), exp(q[2, 3]), exp(q[3, 3]), sd(ch[, 1]), sd(ch[, 2]), sd(ch[, 3])))
if (exists("gm")) cat(sprintf("RNG |z| MAP in posterior sd: mu %.2f log sd %.2f log rho %.2f\n", abs(gm$beta[1] - mean(ch[, 1])) / sd(ch[, 1]), abs(log(gm$sd) - mean(ch[, 2])) / sd(ch[, 2]), abs(log(gm$rho) - mean(ch[, 3])) / sd(ch[, 3])))
cat("DONE\n")
