# Item 5: kriging mean and SE at low noise against a hand delta method on the closed
# form: var = cond var at theta_hat + J V_theta J', with V_theta from the marginal's
# numeric Hessian. The lane says its test is only a regression guard.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
arm <- commandArgs(TRUE)[1]
rev_libs(arm)
nug <- if (arm == "lane") 1e-12 else 1e-6
krig <- function(p, s, X0, ...) {
  # p: mu, log sigma, log sd, log rho. eta at new positions, y's marginal
  mu <- p[1]; sigma <- exp(p[2]); sd <- exp(p[3]); rho <- exp(p[4])
  C <- exp(-s$D2 / (2 * rho^2)); diag(C) <- diag(C) + s$nug
  r <- (sigma / sd)^2
  # kept at position level: n = m in these designs
  A <- C; diag(A) <- diag(A) + r
  R <- chol(A)
  Ks <- exp(-outer(X0, s$ux, function(a, b) (a - b)^2) / (2 * rho^2))
  al <- backsolve(R, backsolve(R, s$y - mu, transpose = TRUE))
  mean <- mu + as.vector(Ks %*% al)
  w <- backsolve(R, t(Ks), transpose = TRUE)
  # brms's and the fit's new-position variance has a nugget of its own? read without
  cv <- sd^2 * pmax(1 - colSums(w^2), 0)
  list(mean = mean, cv = cv)
}
for (cs in list(list(n = 40, sg = 0.3, om = 2), list(n = 40, sg = 1e-2, om = 2), list(n = 60, sg = 1e-3, om = 2), list(n = 60, sg = 1e-4, om = 1), list(n = 100, sg = 1e-3, om = 2), list(n = 30, sg = 3e-5, om = 1))) {
  set.seed(1500 + cs$n)
  x <- sort(runif(cs$n, 0, 8)); y <- sin(cs$om * x) + cs$sg * rnorm(cs$n)
  r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x))
  f <- r$fit; g <- rev_gp_read(f)
  s <- cf_setup(y, x, nug = nug); s$ux <- x
  P <- c(g$beta[[1]], log(g$sigma), log(g$sd), log(g$rho))
  nll <- function(p) -cf_ll_at(s, p[1], exp(p[2]), exp(p[3]), exp(p[4]))
  V <- solve(optimHess(P, nll, control = list(ndeps = rep(1e-3, 4))))
  x0 <- c(0.1, 1.234, 3.3, 4.01, 7.3, 8.4, 9)
  k0 <- krig(P, s, x0)
  J <- vapply(1:4, function(i) { h <- 1e-5; e <- replace(numeric(4), i, h); (krig(P + e, s, x0)$mean - krig(P - e, s, x0)$mean) / (2 * h) }, numeric(7))
  dv <- rowSums((J %*% V) * J)
  se0 <- sqrt(k0$cv + dv)
  pr <- fitted(f, newdata = data.frame(x = x0))
  # at the fitted positions: conditional var by the marginal
  cat(sprintf("RK %s n %d sg %g sd %.3f rho %.3f sigma %.3g conv %d w %d\n", arm, cs$n, cs$sg, g$sd, g$rho, g$sigma, g$conv, length(r$warn)))
  cat(sprintf("RK   mean  max |frm - cf| %.2e (cf cv sqrt %.2e .. %.2e)\n", max(abs(pr[, 1] - k0$mean)), min(sqrt(k0$cv)), max(sqrt(k0$cv))))
  cat(sprintf("RK   se ratio frm/cf at new x: %s\n", paste(sprintf("%.4f", pr[, 2] / se0), collapse = " ")))
  rownames(V) <- NULL
  zz <- sqrt(diag(as.matrix(vcov(f, full = TRUE))))
  flush.console()
}
cat("DONE\n")
