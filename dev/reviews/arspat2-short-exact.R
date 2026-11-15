# Short check, item 1: the margin objective against the observed-rows fit
# (fcor, M cut to the observed rows) or a dense reference (sar: the t
# margin of the observed rows, Sigma_oo = [(A'A)^-1]_oo, b by RTMB's
# Laplace), at three outer points; AD gradient against numDeriv off the
# optimum. Lane only.
source("dev/reviews/arspat2-short-common.R")
n <- 60
M <- 0.5^abs(outer(1:n, 1:n, "-"))
mk <- function(seed = 4, df = 4) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 10)))
  d$y <- 1 + 0.5 * d$x + rnorm(6)[d$h] +
    as.vector(t(chol(M)) %*% rnorm(n)) * sqrt(df / rchisq(1, df))
  d
}
grad_chk <- function(fit) {
  p <- fit$opt$par + 0.1 * seq_along(fit$opt$par) / length(fit$opt$par)
  g_ad <- as.vector(fit$obj$gr(p))
  g_fd <- numDeriv::grad(function(q) fit$obj$fn(q), p)
  max(abs(g_ad - g_fd)) / max(abs(g_fd))
}
nu_idx <- function(fit) grep("^nu", ns$outer_par_labels(fit))
fcor_case <- function(lab, f, d, miss, profile = FALSE, nus = c(3, 6, 21)) {
  d$y[miss] <- NA
  ok <- !is.na(d$y)
  ctl <- frmtmb_control(profile = profile)
  fo <- f
  fo[[3]] <- do.call(substitute, list(fo[[3]], list(M = quote(Mo))))
  a <- quiet(frm(update(f, y | mi() ~ .), data = d, data2 = list(M = M),
                 family = student(), control = ctl))
  b <- quiet(frm(fo, data = d[ok, ], data2 = list(Mo = M[ok, ok]),
                 family = student(), control = ctl))
  stopifnot(identical(names(a$opt$par), names(b$opt$par)))
  inu <- nu_idx(a)
  r <- vapply(0:2, function(k) {
    p <- pts(a$opt$par, k)
    p[inu] <- log(nus[k + 1L] - 1)
    reld(a$obj$fn(p), b$obj$fn(p))
  }, 0)
  cat(sprintf("%-26s m %2d unread %2d latent %d | rel %s | grad rel %.2e | conv %d\n",
              lab, sum(!ok), length(a$frame$mi_map$y$unread),
              length(a$estimates$miss), paste(sprintf("%.1e", r), collapse = " "),
              grad_chk(a), a$opt$convergence))
}
d <- mk()
fcor_case("fcor (1|h)", y ~ x + (1 | h) + fcor(M), d,
          c(2, 8, 15, 22, 29, 36, 43, 50, 57))
fcor_case("fcor (1|h) profile", y ~ x + (1 | h) + fcor(M), d,
          c(2, 8, 15, 22, 29, 36, 43, 50, 57), profile = TRUE)
fcor_case("fcor nu 1.02/1.1/1.5", y ~ x + (1 | h) + fcor(M), mk(df = 1.2),
          c(3, 9, 17, 25, 31, 44, 52), nus = c(1.02, 1.1, 1.5))
set.seed(40)
fcor_case("fcor m = n/2 (1|h)", y ~ x + (1 | h) + fcor(M), d,
          sort(sample(n, n / 2)))

# sar against the dense reference
sar_case <- function(lab, W, ty, seed, nmiss) {
  n <- nrow(W)
  ev <- Re(eigen(W, only.values = TRUE)$values)
  set.seed(seed)
  d <- data.frame(x = rnorm(n), h = factor(rep(1:7, length.out = n)))
  A0 <- diag(n) - 0.8 / max(ev) * W
  d$y <- as.vector(solve(A0, 1 + 0.5 * d$x + rnorm(7)[d$h] +
                           sqrt(4 / rchisq(1, 4)) * rnorm(n)))
  miss <- sort(sample(n, nmiss))
  ok <- !(seq_len(n) %in% miss)
  d$y[miss] <- NA
  fo <- stats::as.formula(paste0("y | mi() ~ x + (1 | h) + sar(W, type = '",
                                 ty, "')"))
  fit <- quiet(frm(fo, family = student(), data = d, data2 = list(W = W)))
  ac <- fit$frame$autocor[[1L]]
  lab_p <- ns$outer_par_labels(fit)
  nm <- names(fit$opt$par)
  X <- cbind(1, d$x)
  Z <- stats::model.matrix(~ 0 + h, d)
  yo <- d$y[ok]
  ko <- sum(ok)
  ref <- function(p) {
    rho <- unname(ns$autocor_natural(p[nm == "thetaac"], ac))
    A <- diag(n) - rho * W
    B <- if (ty == "lag") solve(A) else diag(n)
    So <- solve(crossprod(A))[ok, ok]
    Soi <- solve(So)
    ldS <- as.numeric(determinant(So)$modulus)
    Bo <- B[ok, , drop = FALSE]
    beta <- unname(p[nm == "beta"])
    bd <- unname(p[nm == "betad"])
    lsd <- unname(p[nm == "theta"])
    f <- function(q) {
      sig <- exp(bd[1])
      nu <- 1 + exp(bd[2])
      r <- yo - as.vector(Bo %*% (X %*% beta + Z %*% q$b))
      Q <- sum(r * as.vector(Soi %*% r)) / sig^2
      -(lgamma((nu + ko) / 2) - lgamma(nu / 2) - 0.5 * ko * log(nu * pi) -
          0.5 * (ldS + 2 * ko * log(sig)) - 0.5 * (nu + ko) * log1p(Q / nu) +
          sum(RTMB::dnorm(q$b, 0, exp(lsd), log = TRUE)))
    }
    obj <- RTMB::MakeADFun(f, list(b = numeric(ncol(Z))), random = "b",
                           silent = TRUE)
    obj$fn(numeric(0))
  }
  inu <- nu_idx(fit)
  r <- vapply(0:2, function(k) {
    p <- pts(fit$opt$par, k)
    p[inu] <- log(c(2, 5, 20)[k + 1L] - 1)
    reld(fit$obj$fn(p), ref(p))
  }, 0)
  cat(sprintf("%-26s m %2d unread %2d latent %d | rel %s | grad rel %.2e | conv %d | labels %s\n",
              lab, nmiss, length(fit$frame$mi_map$y$unread),
              length(fit$estimates$miss), paste(sprintf("%.1e", r), collapse = " "),
              grad_chk(fit), fit$opt$convergence, paste(lab_p, collapse = ",")))
}
Wb <- rook(7)
sar_case("sar lag, binary W 7x7", Wb, "lag", 21, 15)
Wr <- rook(8)
Wr <- Wr / rowSums(Wr)
sar_case("sar error, row-std W 8x8", Wr, "error", 22, 16)
