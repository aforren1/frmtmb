# Lane adhess: step rules for the outer Hessian of a fit WITH random
# effects (central differences of the Laplace gradient). Each candidate
# is compared with a reference: Richardson extrapolation (numDeriv, six
# halvings) of the AD gradient in standard-error units around the
# estimate, at two base steps whose agreement is printed (refself).
#   ADHESS_ARM=lane|base Rscript dev/adhess-restep.R [model ...]
# Candidates:
#   abs   fd_hessian() as 0.71.0 builds it (absolute 1e-3)
#   des   1e-3 / max(1, rms of the design column) for design coefficients
#   curv  abs, then every column whose rise exceeds 1e-2 retaken at a
#         rise of 1e-4 (fd_hessian()'s retake, unconditional)
#   desR  des with Richardson (steps h and 2h)
source("C:/Users/adf44/source/r/frmtmb-wt-adhess/dev/adhess-common.R")
args <- commandArgs(TRUE)
fd_steps <- function(par, gr, h, richardson = FALSE) {
  n <- length(par)
  H <- matrix(0, n, n)
  col <- function(i, hi) {
    dp <- par; dp[i] <- par[i] + hi
    g1 <- as.numeric(gr(dp))
    dp[i] <- par[i] - hi
    g2 <- as.numeric(gr(dp))
    (g1 - g2) / (2 * hi)
  }
  for (i in seq_len(n)) {
    c1 <- col(i, h[i])
    if (richardson) c1 <- (4 * c1 - col(i, 2 * h[i])) / 3
    H[, i] <- c1
  }
  0.5 * (H + t(H))
}
design_scale <- function(fit) {
  frame <- fit$frame
  pn <- names(fit$opt$par)
  s <- rep(1, length(pn))
  for (lp in frame[["linpreds"]] %||% list()) {
    if (!is.null(lp[["constant"]]) || !is.null(lp[["nl_body"]]) ||
          length(lp[["mo"]] %||% list()) ||
          !lp[["par"]] %in% c("beta", "betad")) next
    X <- lp[["X"]]
    if (!ns$design_ok(X)) next
    pos <- ns$outer_index_of(lp, pn, frame)
    if (length(pos) != ncol(X)) next
    rms <- sqrt(as.numeric(Matrix::colSums(X^2)) / nrow(X))
    s[pos] <- pmax(s[pos], rms)
  }
  s
}
curv_retake <- function(par, gr, H0) {
  d <- diag(H0)
  redo <- which(is.finite(d) & d > 0 & 0.5 * d * 1e-6 > 1e-2)
  h <- rep(1e-3, length(par))
  h[redo] <- sqrt(2 * 1e-4 / d[redo])
  list(H = fd_steps(par, gr, h), n = length(redo))
}
# des, then a non-design column whose rise exceeds 1e-2 retaken at a
# rise of 1e-4
desc_retake <- function(par, gr, h, s) {
  H <- fd_steps(par, gr, h)
  d <- diag(H)
  redo <- which(s == 1 & is.finite(d) & d > 0 & 0.5 * d * h^2 > 1e-2)
  if (!length(redo)) return(H)
  h[redo] <- sqrt(2 * 1e-4 / d[redo])
  fd_steps(par, gr, h)
}
ref_hessian <- function(obj, p, D, d0) {
  g0 <- function(dz) as.numeric(obj$gr(p + dz / D)) / D
  J <- numDeriv::jacobian(g0, rep(0, length(p)), method = "Richardson",
                          method.args = list(eps = d0, d = d0, zero.tol = 1,
                                             r = 6, v = 2))
  J <- (J + t(J)) / 2
  J * outer(D, D)
}
set.seed(1)
mk <- list()
mk$pois_xbig_re <- function() {
  n <- 300
  d <- data.frame(x = rnorm(n, 500, 100), g = factor(rep(1:20, length.out = n)))
  d$y <- rpois(n, exp(0.3 + 0.004 * (d$x - 500) + rnorm(20, 0, 0.3)[d$g]))
  list(f = bf(y ~ x + (1 | g)), fam = poisson(), d = d,
       lme4 = quote(lme4::glmer(y ~ x + (1 | g), data = d, family = poisson)))
}
mk$pois_re <- function() {
  n <- 300
  d <- data.frame(x = rnorm(n), g = factor(rep(1:20, length.out = n)))
  d$y <- rpois(n, exp(0.3 + 0.4 * d$x + rnorm(20, 0, 0.3)[d$g]))
  list(f = bf(y ~ x + (1 | g)), fam = poisson(), d = d,
       lme4 = quote(lme4::glmer(y ~ x + (1 | g), data = d, family = poisson)))
}
mk$bern_xbig_re <- function() {
  n <- 600
  d <- data.frame(x = rnorm(n, 0, 100), g = factor(rep(1:40, length.out = n)))
  d$y <- rbinom(n, 1, plogis(0.2 + 0.01 * d$x + rnorm(40, 0, 0.5)[d$g]))
  list(f = bf(y ~ x + (1 | g)), fam = bernoulli(), d = d,
       lme4 = quote(lme4::glmer(y ~ x + (1 | g), data = d, family = binomial)))
}
mk$negbin_re <- function() {
  n <- 400
  d <- data.frame(x = rnorm(n, 50, 10), g = factor(rep(1:25, length.out = n)))
  d$y <- rnbinom(n, mu = exp(1 + 0.03 * (d$x - 50) + rnorm(25, 0, 0.4)[d$g]),
                 size = 3)
  list(f = bf(y ~ x + (1 | g)), fam = negbinomial(), d = d)
}
mk$gamma_re <- function() {
  n <- 400
  d <- data.frame(x = rnorm(n, 200, 50), g = factor(rep(1:25, length.out = n)))
  mu <- exp(0.5 + 0.004 * (d$x - 200) + rnorm(25, 0, 0.3)[d$g])
  d$y <- rgamma(n, shape = 4, rate = 4 / mu)
  list(f = bf(y ~ x + (1 | g)), fam = Gamma(link = "log"), d = d)
}
mk$pois_slope_re <- function() {
  n <- 600
  g <- factor(rep(1:30, length.out = n))
  d <- data.frame(x = rnorm(n), g = g)
  u0 <- rnorm(30, 0, 0.4); u1 <- rnorm(30, 0, 0.2)
  d$y <- rpois(n, exp(0.5 + 0.3 * d$x + u0[g] + u1[g] * d$x))
  list(f = bf(y ~ x + (x | g)), fam = poisson(), d = d,
       lme4 = quote(lme4::glmer(y ~ x + (x | g), data = d, family = poisson)))
}
mk$pois_n1e4 <- function() {
  n <- 1e4
  d <- data.frame(x = rnorm(n), g = factor(rep(1:100, length.out = n)))
  d$y <- rpois(n, exp(0.3 + 0.4 * d$x + rnorm(100, 0, 0.3)[d$g]))
  list(f = bf(y ~ x + (1 | g)), fam = poisson(), d = d)
}
mk$pois_xbig_n1e4 <- function() {
  n <- 1e4
  d <- data.frame(x = rnorm(n, 500, 100), g = factor(rep(1:100, length.out = n)))
  d$y <- rpois(n, exp(0.3 + 0.004 * (d$x - 500) + rnorm(100, 0, 0.3)[d$g]))
  list(f = bf(y ~ x + (1 | g)), fam = poisson(), d = d)
}
mk$sleep_lmm <- function() {
  data("sleepstudy", package = "lme4", envir = environment())
  list(f = bf(Reaction ~ Days + (Days | Subject)), fam = gaussian(),
       d = sleepstudy)
}
mk$gauss_sigma_re <- function() {
  n <- 400
  d <- data.frame(x = rnorm(n), z = rnorm(n, 100, 20),
                  g = factor(rep(1:20, length.out = n)))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 0.5)[d$g] +
    rnorm(n, 0, exp(-0.5 + 0.01 * (d$z - 100)))
  list(f = bf(y ~ x + (1 | g), sigma ~ z), fam = gaussian(), d = d)
}
mk$nl_re <- function() {
  n <- 300
  d <- data.frame(x = 5 + runif(n), g = factor(rep(1:15, length.out = n)))
  d$y <- (2 + rnorm(15, 0, 0.2)[d$g]) * exp(0.8 * (d$x - 5)) +
    rnorm(n, 0, 0.3)
  list(f = bf(y ~ a * exp(b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
       fam = gaussian(), d = d,
       start = list(beta = c(2 * exp(-0.8 * 5), 0.8)))
}
mk$pois_xsmall_re <- function() {
  n <- 300
  d <- data.frame(x = rnorm(n, 0, 0.01), g = factor(rep(1:20, length.out = n)))
  d$y <- rpois(n, exp(0.3 + 40 * d$x + rnorm(20, 0, 0.3)[d$g]))
  list(f = bf(y ~ x + (1 | g)), fam = poisson(), d = d)
}
mk$bern_x10_re <- function() {
  n <- 1000
  d <- data.frame(x = rnorm(n, 0, 10), g = factor(rep(1:50, length.out = n)))
  d$y <- rbinom(n, 1, plogis(-0.3 + 0.08 * d$x + rnorm(50, 0, 0.6)[d$g]))
  list(f = bf(y ~ x + (1 | g)), fam = bernoulli(), d = d,
       lme4 = quote(lme4::glmer(y ~ x + (1 | g), data = d, family = binomial)))
}
todo <- if (length(args)) args else names(mk)
for (nm in todo) {
  set.seed(match(nm, names(mk)))
  cs <- mk[[nm]]()
  r <- cap(if (is.null(cs$start)) frm(cs$f, family = cs$fam, data = cs$d) else
             frm(cs$f, family = cs$fam, data = cs$d, start = cs$start))
  f <- r$v
  if (is_err(f)) { cat(sprintf("%-15s ERROR %s\n", nm, f)); next }
  obj <- f$obj
  p <- ns$sdr_outer_point(obj)
  saved <- ns$obj_state_save(obj)
  gr <- function(q) obj$gr(q)
  t0 <- proc.time()[[3]]
  H0 <- ns$fd_hessian(p, gr)$H
  tabs <- proc.time()[[3]] - t0
  D <- sqrt(abs(diag(H0)))
  R1 <- ref_hessian(obj, p, D, 1e-1)
  R2 <- ref_hessian(obj, p, D, 1e-2)
  ns$obj_state_restore(obj, saved)
  V1 <- tryCatch(solve(R1), error = function(e) R1 * NaN)
  V2 <- tryCatch(solve(R2), error = function(e) R2 * NaN)
  if (!all(is.finite(V2))) {
    cat(sprintf("%-15s reference not invertible, lost %s, warn %s\n",
                nm, lost_of(f), wtag(r$w)))
    next
  }
  s <- design_scale(f)
  hdes <- 1e-3 / pmax(1, s)
  cand <- list(abs = H0,
               des = fd_steps(p, gr, hdes),
               curv = curv_retake(p, gr, H0)$H,
               desR = fd_steps(p, gr, hdes, richardson = TRUE),
               des2 = fd_steps(p, gr, 1e-3 / s),
               desC = desc_retake(p, gr, hdes, s))
  ns$obj_state_restore(obj, saved)
  ev <- eigen(cov2cor(V2), only.values = TRUE)$values
  cat(sprintf("%-15s %s | np %d | n %d | code %d | kappa %.1e | refself %s | scale %s | warn %s | lost %s\n",
              nm, arm, length(p), nrow(cs$d), f$opt$convergence,
              max(ev) / min(ev),
              paste(sprintf("%.1e", cov_err(V1, V2)), collapse = "/"),
              paste(signif(s, 3), collapse = ","), wtag(r$w), lost_of(f)))
  for (cn in names(cand)) {
    V <- tryCatch(solve(cand[[cn]]), error = function(e) NULL)
    e <- if (is.null(V)) c(NaN, NaN) else cov_err(V, V2)
    cat(sprintf("    %-5s SE err %.2e | cor err %.2e\n", cn, e[1], e[2]))
  }
  Vf <- suppressWarnings(vcov(f, full = TRUE))
  e <- cov_err(Vf, V2)
  cat(sprintf("    vcov  SE err %.2e | cor err %.2e (the fit's own)\n",
              e[1], e[2]))
  if (!is.null(cs$lme4)) {
    d <- cs$d
    l <- suppressWarnings(eval(cs$lme4))
    sel <- sqrt(diag(as.matrix(vcov(l))))
    k <- seq_along(sel)
    cat(sprintf("    lme4 fixed SE / ref - 1: %s\n",
                paste(sprintf("%+.2e", sel / sqrt(diag(V2))[k] - 1),
                      collapse = " ")))
  }
}
