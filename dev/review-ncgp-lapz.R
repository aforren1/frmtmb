# Item 7: frm()'s Laplace objective against a hand Laplace in z, and exact/Laplace by importance
# sampling with the Gaussian as a control, n 8, the same designs as review-ncgp-ngref.R.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
for (case in c("pois16", "pois12", "binom7", "binom5", "pois4", "bern")) {
  nm <- 8L
  set.seed(7100 + match(case, c("pois16", "pois12", "binom7", "binom5", "pois4", "bern")))
  x <- sort(runif(nm, 0, 6))
  fam <- switch(case, pois16 = "pois", pois12 = "pois", pois4 = "pois", "binom")
  lm0 <- switch(case, pois16 = 16, pois12 = 12, pois4 = 4, 0)
  tr <- switch(case, binom7 = 1e7, binom5 = 1e5, bern = 20, NA)
  y <- if (fam == "pois") rpois(nm, exp(lm0 + sin(x))) else rbinom(nm, tr, plogis(lm0 + sin(x)))
  if (fam == "pois") { d <- data.frame(y = y, x = x); fml <- y ~ gp(x, scale = FALSE); fm <- poisson() } else { d <- data.frame(y = y, x = x, tr = tr); fml <- y | trials(tr) ~ gp(x, scale = FALSE); fm <- binomial() }
  r <- rev_fit(fml, d, family = fm); f <- r$fit; g <- rev_gp_read(f)
  K <- exp(-outer(x, x, "-")^2 / (2 * g$rho^2)); diag(K) <- diag(K) * (1 + 1e-12)
  bk <- f$frame$re_blocks[[1]]; p <- bk$gp_order
  L <- matrix(0, nm, nm); L[p, ] <- t(chol(K[p, p]))
  J <- g$sd * L
  ll <- function(eta) if (fam == "pois") sum(dpois(y, exp(eta), log = TRUE)) else sum(dbinom(y, tr, plogis(eta), log = TRUE))
  w_of <- function(eta) if (fam == "pois") exp(eta) else tr * plogis(eta) * (1 - plogis(eta))
  g_of <- function(eta) if (fam == "pois") y - exp(eta) else y - tr * plogis(eta)
  z <- f$estimates$b
  for (it in 1:100) {
    eta <- g$beta[1] + as.vector(J %*% z)
    H <- crossprod(J * sqrt(w_of(eta))); diag(H) <- diag(H) + 1
    gr <- as.vector(crossprod(J, g_of(eta))) - z
    st <- solve(H, gr); z <- z + st
    if (max(abs(st)) < 1e-14) break
  }
  eta <- g$beta[1] + as.vector(J %*% z)
  H <- crossprod(J * sqrt(w_of(eta))); diag(H) <- diag(H) + 1
  Rz <- chol(H)
  lapz <- ll(eta) - 0.5 * sum(z^2) - sum(log(diag(Rz)))
  # exact over Laplace: importance sampling in z, a t proposal, Gaussian control on the same draws
  set.seed(11)
  S <- 20000L; tdf <- 7
  zs <- matrix(rnorm(S * nm), S, nm); chi <- rchisq(S, tdf); sc <- sqrt(tdf / chi)
  dr <- t(backsolve(Rz, t(zs))) * sc
  lt <- lgamma((tdf + nm) / 2) - lgamma(tdf / 2) - nm / 2 * log(tdf * pi) + sum(log(diag(Rz))) - (tdf + nm) / 2 * log1p(rowSums(zs^2) / chi)
  q <- rowSums(zs^2) * sc^2
  lg <- -0.5 * q + sum(log(diag(Rz))) - nm / 2 * log(2 * pi)
  lg_val <- lapz + sum(log(diag(Rz))) - 0.5 * q  # the Laplace Gaussian, integrating to exp(lapz)
  lw <- vapply(seq_len(S), function(i) { zi <- z + dr[i, ]; ll(g$beta[1] + as.vector(J %*% zi)) - 0.5 * sum(zi^2) - nm / 2 * log(2 * pi) }, 0) - lt
  lw0 <- lg_val - nm / 2 * log(2 * pi) - lt
  a <- max(lw0)
  ratio <- mean(exp(lw - a)) / mean(exp(lw0 - a))
  se <- sd(exp(lw - a) - ratio * exp(lw0 - a)) / sqrt(S) / mean(exp(lw0 - a))
  cat(sprintf("RLZ %-7s frm ll %.8f hand Laplace in z %.8f diff %.2e | exact - Laplace (IS, control) %.5f +- %.5f | conv %d w %d sd %.3f rho %.3f mu %.3f\n", case, g$ll, lapz, g$ll - lapz, log(ratio), se / ratio, g$conv, length(r$warn), g$sd, g$rho, g$beta[1]))
  flush.console()
}
cat("DONE\n")
