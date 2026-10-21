# Reviewer, fix round 1: the log-scale simplex_rdirichlet() and the
# .eta_shape removal, run on the reinstalled lane build.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
rd <- frmtmb:::simplex_rdirichlet
p <- exp(c(0, 0.3, -0.2, -2.5))
p <- p / sum(p)
K <- length(p)
z <- qnorm(1 - 1e-6 / 2)
n <- 1e5
cat("## moments against the analytic Dirichlet, n = 1e5, seed 7\n")
for (phi in c(1e-6, 1e-3, 0.01, 0.5, 20, 1e4, 1e8)) {
  set.seed(7)
  a <- matrix(rep(p * phi, each = n), n)
  x <- rd(a, letters[1:K])
  m_err <- (colMeans(x) - p) / (sqrt(p * (1 - p) / (phi + 1)) / sqrt(n))
  v <- p * (1 - p) / (phi + 1)
  sev <- apply(x, 2, function(c) sd((c - mean(c))^2)) / sqrt(n)
  v_err <- (apply(x, 2, var) - v) / sev
  cv <- -p[1] * p[2] / (phi + 1)
  pr <- (x[, 1] - mean(x[, 1])) * (x[, 2] - mean(x[, 2]))
  c_err <- (mean(pr) - cv) / (sd(pr) / sqrt(n))
  cat(sprintf(paste("phi %-6g NaN rows %d, rows sum to 1 (max |s-1| %.1e);",
                    "max |z| mean %.2f var %.2f cov %.2f (bound %.2f)\n"),
              phi, sum(!is.finite(rowSums(x))), max(abs(rowSums(x) - 1)),
              max(abs(m_err)), max(abs(v_err)), abs(c_err), z))
}
cat("\n## rows where every log-gamma is very negative: vertex shares\n")
for (phi in c(1e-6, 1e-9, 1e-12)) {
  set.seed(8)
  a <- matrix(rep(p * phi, each = n), n)
  x <- rd(a, letters[1:K])
  share <- tabulate(max.col(x, ties.method = "first"), K) / n
  vert <- mean(apply(x, 1, max) > 1 - 1e-12)
  zs <- (share - p) / sqrt(p * (1 - p) / n)
  cat(sprintf("phi %g: share at a vertex %.5f; vertex shares %s vs p %s; max |z| %.2f\n",
              phi, vert, paste(sprintf("%.4f", share), collapse = " "),
              paste(sprintf("%.4f", p), collapse = " "), max(abs(zs))))
}
cat("\n## extreme alphas\n")
probe("alpha 0 in one cell (phi * p underflow)", {
  a <- matrix(c(0, 1e-3, 1e-3, 1e-3), 1)
  rd(a[rep(1, 5), , drop = FALSE], letters[1:4])
})
probe("alpha 1e300", {
  a <- matrix(c(1e300, 1e300, 1e300, 1e300), 1)
  rd(a, letters[1:4])
})
cat("\n## predict() at phi 0.02 (data of dev/nb2dm-rev-smallphi.R)\n")
set.seed(41)
nn <- 200
d <- data.frame(x = rnorm(nn))
P <- exp(cbind(0, 0.3 + 0.4 * d$x, -0.2))
P <- P / rowSums(P)
d$N <- sample(5:30, nn, TRUE)
d$Y <- t(vapply(seq_len(nn), function(i) {
  lg <- log(stats::rgamma(3, P[i, ] * 0.02 + 1)) +
    log(stats::runif(3)) / (P[i, ] * 0.02)
  q <- exp(lg - max(lg))
  as.numeric(rmultinom(1, d$N[i], q / sum(q)))
}, numeric(3)))
colnames(d$Y) <- c("a", "b", "c")
fit <- frm(Y | trials(N) ~ x, family = dirichlet_multinomial(), data = d)
probe("phi", exp(fixef_by_dpar(fit)$phi))
probe("predict(ndraws = 2000) dims and NA", {
  set.seed(2)
  pr <- predict(fit, ndraws = 2000)
  list(dim = dim(pr), any_na = anyNA(pr),
       rowsum_ok = all(abs(rowSums(pr[, "Estimate", ]) - d$N) < 1e-8 * d$N))
})
probe("simulate(nsim = 200) rows sum to trials", {
  s <- simulate(fit, nsim = 200, seed = 3)
  all(vapply(s, function(m) all(rowSums(m) == d$N), NA))
})

cat("\n## negbinomial2 without .eta_shape: tape values at extreme sigma\n")
nb2 <- brmsfamily("negbinomial2")
nb <- negbinomial()
y <- c(0, 2, 9, 40)
tape2 <- function(ls) {
  p0 <- c(rep(log(3), 4), rep(ls, 4))
  tp <- RTMB::MakeTape(function(q) {
    frmtmb:::row_lpdf(nb2, y, y, list(mu = exp(q[1:4]), sigma = exp(q[5:8]),
                                      .eta_mu = q[1:4], .eta_sigma = q[5:8]),
                      list(), NULL)
  }, p0)
  list(v = as.numeric(tp(p0)), j = tp$jacobian(p0))
}
tape1 <- function(ls) {
  p0 <- c(rep(log(3), 4), rep(-ls, 4))
  tp <- RTMB::MakeTape(function(q) {
    frmtmb:::row_lpdf(nb, y, y, list(mu = exp(q[1:4]), shape = exp(q[5:8]),
                                     .eta_mu = q[1:4], .eta_shape = q[5:8]),
                      list(), NULL)
  }, p0)
  list(v = as.numeric(tp(p0)), j = tp$jacobian(p0))
}
for (ls in c(-700, -40, -20, -5, 0, 3, 20, 40, 700)) {
  a <- tape2(ls)
  b <- tape1(ls)
  ref <- if (abs(ls) < 30) dnbinom(y, size = exp(-ls), mu = 3, log = TRUE) else
    rep(NA, 4)
  cat(sprintf("log sigma %5g: nb2 %s | nb(shape=1/sigma) %s | dnbinom %s | max rel nb2/nb %.2e\n",
              ls, paste(format(a$v, digits = 10), collapse = " "),
              paste(format(b$v, digits = 10), collapse = " "),
              paste(format(ref, digits = 10), collapse = " "),
              max(abs(a$v / b$v - 1), na.rm = TRUE)))
}
probe("negbinomial2 non-count refusal", tryCatch(
  frm(I(y + 0.5) ~ 1, family = nb2,
      data = data.frame(y = rpois(30, 3))), error = conditionMessage))
probe("Poisson-limit fit, seed 4 of runoff", {
  set.seed(4)
  dd <- data.frame(x = rnorm(300))
  dd$y <- rpois(300, exp(1 + 0.3 * dd$x))
  m2 <- suppressWarnings(frm(y ~ x, family = nb2, data = dd))
  m1 <- suppressWarnings(frm(y ~ x, family = nb, data = dd))
  c(ll_nb2 = logLik(m2), ll_nb = logLik(m1), conv2 = m2$opt$convergence)
})
