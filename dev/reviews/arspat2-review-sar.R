# Review of lane arspat2, attack 1: sar() lag with a W that is NOT row
# standardized (binary rook adjacency, 7 x 7), student(), mi(), (1 | h);
# a multivariate model where one response is autocorrelated; two
# autocorrelated responses; and the AD gradient against finite differences
# on the margin route. Reference for sar: the dense t margin of the
# observed rows with b integrated by RTMB's Laplace (as test-sar.R).
source("dev/reviews/arspat2-review-common.R")
k <- 7
W <- rook(k)
n <- nrow(W)
ev <- eigen(W, only.values = TRUE)$values
cat("W binary rook 7x7: eigen range", format(range(ev), digits = 5), "\n")
set.seed(21)
d <- data.frame(x = rnorm(n), h = factor(rep(1:7, length.out = n)))
rho0 <- 0.8 / max(ev)
A0 <- diag(n) - rho0 * W
tsc <- sqrt(4 / rchisq(1, 4))
d$y <- as.vector(solve(A0, 1 + 0.5 * d$x + rnorm(7)[d$h] + tsc * rnorm(n)))
miss <- sort(sample(n, 15))
ok <- !(seq_len(n) %in% miss)
X <- cbind(1, d$x)
Z <- stats::model.matrix(~ 0 + h, d)
ref <- function(p, rho, y, lag) {
  A <- diag(n) - rho * W
  B <- if (lag) solve(A) else diag(n)
  Bo <- B[ok, , drop = FALSE]
  So <- solve(crossprod(A))[ok, ok]
  Soi <- solve(So)
  ldS <- as.numeric(determinant(So)$modulus)
  yo <- y[ok]
  kk <- sum(ok)
  f <- function(q) {
    sig <- exp(q$lsig)
    nu <- 1 + exp(q$lnu)
    r <- yo - as.vector(Bo %*% (X %*% q$beta + Z %*% q$b))
    Q <- sum(r * as.vector(Soi %*% r)) / sig^2
    -(lgamma((nu + kk) / 2) - lgamma(nu / 2) - 0.5 * kk * log(nu * pi) -
        0.5 * (ldS + 2 * kk * log(sig)) - 0.5 * (nu + kk) * log1p(Q / nu) +
        sum(RTMB::dnorm(q$b, 0, exp(q$lsd), log = TRUE)))
  }
  q0 <- list(beta = p[1:2], lsig = p[[3]], lnu = p[[4]], lsd = p[[5]],
             b = numeric(ncol(Z)))
  obj <- RTMB::MakeADFun(f, q0, random = "b", silent = TRUE)
  obj$fn(p[1:5])
}
fits <- list()
for (ty in c("lag", "error")) {
  dm <- d
  dm$y[miss] <- NA
  fo <- stats::as.formula(paste0("y | mi() ~ x + (1 | h) + sar(W, type = '",
                                 ty, "')"))
  fit <- quiet(frm(fo, family = student(), data = dm, data2 = list(W = W)))
  fits[[ty]] <- fit
  ac <- fit$frame[["autocor"]][[1L]]
  iac <- which(names(fit$opt$par) == "thetaac")
  inu <- which(names(fit$opt$par) == "betad")[2L]
  cat(sprintf("sar %s: lo %.5f hi %.5f, labels %s, route margin %s\n", ty,
              ac$aux$lo, ac$aux$hi,
              paste(names(fit$opt$par), collapse = " "),
              ns$autocor_seam_t_missing(fit$frame)$y$margin))
  for (kk in 0:2) {
    p <- fit$opt$par + 0.05 * kk * seq_along(fit$opt$par) / length(fit$opt$par)
    p[inu] <- log(c(2, 5, 20)[kk + 1L])
    rho <- unname(ns$autocor_natural(p[iac], ac))
    a <- fit$obj$fn(p)
    b <- ref(unname(p[-iac]), rho, dm$y, ty == "lag")
    cat(sprintf("  point %d rho %.4f: fn %.8f ref %.8f rel %.2e\n", kk, rho,
                a, b, abs(a - b) / abs(b)))
  }
}

# AD gradient against Richardson finite differences, off the optimum
grad_chk <- function(lab, fit) {
  p <- fit$opt$par + 0.1 * seq_along(fit$opt$par) / length(fit$opt$par)
  g_ad <- as.vector(fit$obj$gr(p))
  g_fd <- numDeriv::grad(function(q) fit$obj$fn(q), p)
  cat(sprintf("%-28s max |gr AD - FD| %.2e, max |gr| %.2e, rel %.2e\n", lab,
              max(abs(g_ad - g_fd)), max(abs(g_fd)),
              max(abs(g_ad - g_fd)) / max(abs(g_fd))))
}
grad_chk("sar lag margin route", fits$lag)
grad_chk("sar error margin route", fits$error)
nf <- 60
M <- 0.5^abs(outer(1:nf, 1:nf, "-"))
set.seed(4)
df <- data.frame(x = rnorm(nf), h = factor(rep(1:6, each = 10)),
                 z = rnorm(nf))
df$y <- 1 + 0.5 * df$x + rnorm(6)[df$h] +
  as.vector(t(chol(M)) %*% rnorm(nf)) * sqrt(4 / rchisq(1, 4))
df$y2 <- 0.5 * df$x + as.vector(t(chol(M)) %*% rnorm(nf)) *
  sqrt(4 / rchisq(1, 4))
df$z <- 0.3 * df$x + rnorm(6)[df$h] + rnorm(nf)
ms <- c(2, 8, 15, 22, 29, 36, 43, 50)
dm <- df
dm$y[ms] <- NA
ff <- quiet(frm(bf(y | mi() ~ x + (1 | h) + fcor(M)), family = student(),
                data = dm, data2 = list(M = M)))
grad_chk("fcor margin route", ff)
fg <- quiet(frm(bf(y | mi() ~ x + fcor(M), sigma ~ 1 + (1 | h)),
                family = gaussian(), data = dm, data2 = list(M = M)))
grad_chk("fcor gaussian sigma RE", fg)

# multivariate, one response autocorrelated: separable, so logLik is the
# sum of the observed-rows fit of y and the fit of z
okr <- !is.na(dm$y)
mv <- quiet(frm(bf(y | mi() ~ x + (1 | h) + fcor(M), family = student()) +
                  bf(z ~ x + (1 | h), family = gaussian()) + set_rescor(FALSE),
                data = dm, data2 = list(M = M)))
fy <- quiet(frm(bf(y ~ x + (1 | h) + fcor(Mo)), family = student(),
                data = dm[okr, ], data2 = list(Mo = M[okr, okr])))
fz <- quiet(frm(bf(z ~ x + (1 | h)), family = gaussian(), data = dm))
tm <- ns$autocor_seam_t_missing(mv$frame)
cat(sprintf("mv one autocorrelated: route %s, logLik mv %.8f, sum %.8f, diff %+.2e, conv %s %s %s\n",
            if (length(tm)) tm$y$margin else "none",
            as.numeric(logLik(mv)), as.numeric(logLik(fy)) + as.numeric(logLik(fz)),
            as.numeric(logLik(mv)) - as.numeric(logLik(fy)) - as.numeric(logLik(fz)),
            mv$opt$convergence, fy$opt$convergence, fz$opt$convergence))
# two autocorrelated responses, both with missing values, RE on the first
dm2 <- dm
dm2$y2[c(5, 12, 19, 40, 47)] <- NA
ok2 <- !is.na(dm2$y2)
mv2 <- quiet(frm(bf(y | mi() ~ x + (1 | h) + fcor(M), family = student()) +
                   bf(y2 | mi() ~ x + fcor(M), family = student()) +
                   set_rescor(FALSE),
                 data = dm2, data2 = list(M = M)))
fy2 <- quiet(frm(bf(y2 ~ x + fcor(Mo)), family = student(),
                 data = dm2[ok2, ], data2 = list(Mo = M[ok2, ok2])))
tm <- ns$autocor_seam_t_missing(mv2$frame)
cat(sprintf("mv two autocorrelated: routes %s, logLik mv %.8f, sum %.8f, diff %+.2e\n",
            paste(vapply(tm, function(t) as.character(t$margin), ""), collapse = ","),
            as.numeric(logLik(mv2)), as.numeric(logLik(fy)) + as.numeric(logLik(fy2)),
            as.numeric(logLik(mv2)) - as.numeric(logLik(fy)) - as.numeric(logLik(fy2))))
# weights / se() beside fcor(): refused?
w <- tryCatch(quiet(frm(bf(y | mi() + weights(z^2) ~ x + (1 | h) + fcor(M)),
                        family = student(), data = dm, data2 = list(M = M))),
              error = function(e) conditionMessage(e))
cat("weights beside fcor + mi:", if (is.character(w)) substr(w, 1, 120) else "FITTED", "\n")
