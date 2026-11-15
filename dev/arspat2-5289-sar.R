# Backlog 5289 for sar() (frmtmb.spatial): student() with sar(), mi() and
# (1 | h). frmtmb's ML objective against an independent reference written
# here: the dense multivariate t margin of the observed rows (mean
# (B (X beta + Z b))_o, scale sigma^2 [(A'A)^-1]_oo, B = A^-1 for the lag
# type and I for the error type, A = I - rho W), with b integrated by
# RTMB's Laplace approximation, which is what the margin route should
# equal. fcor(M) is the control (B = I, scale sigma^2 M_oo); it checks the
# reference's parameter map. 3 points per case, nu near 3, 6 and 21.
# Base build when ARSPAT2_LIB="". Data seed 7, 7 x 7 rook lattice
# (n = 49, row-standardized), 12 missing values, 7 groups.
lib <- Sys.getenv("ARSPAT2_LIB", "C:/Users/adf44/source/r/wt-arspat2-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(dirname(system.file(package = "frmtmb"))), "\n")
cat("frmtmb.spatial", format(packageVersion("frmtmb.spatial")), "from",
    dirname(dirname(system.file(package = "frmtmb.spatial"))), "\n")
ns <- asNamespace("frmtmb")
grid <- expand.grid(r = 1:7, c = 1:7)
n <- nrow(grid)
Aj <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(grid$r[i] - grid$r[j]) + abs(grid$c[i] - grid$c[j]) == 1)
})
W <- Aj / rowSums(Aj)
M <- 0.5^abs(outer(1:n, 1:n, "-"))
set.seed(7)
d <- data.frame(x = rnorm(n), h = factor(rep(1:7, length.out = n)))
eta <- 1 + 0.5 * d$x + rnorm(7, 0, 1)[d$h]
e <- rnorm(n) * sqrt(4 / rchisq(1, 4))
d$ye <- eta + as.vector(solve(diag(n) - 0.4 * W, e))
d$yl <- as.vector(solve(diag(n) - 0.4 * W, eta + e))
d$yf <- eta + as.vector(t(chol(M)) %*% e)
miss <- sort(sample(n, 12))
d$ye[miss] <- NA
d$yl[miss] <- NA
d$yf[miss] <- NA
ok <- !is.na(d$ye)
X <- cbind(1, d$x)
Z <- stats::model.matrix(~ 0 + h, d)

# the reference: -log of the Laplace integral over b, at outer p
ref_ll <- function(p, type, rho, y) {
  if (type == "fcor") {
    B <- diag(n)
    S <- M
  } else {
    A <- diag(n) - rho * W
    S <- solve(crossprod(A))
    B <- if (type == "lag") solve(A) else diag(n)
  }
  Bo <- B[ok, , drop = FALSE]
  So <- S[ok, ok]
  Soi <- solve(So)
  ldS <- as.numeric(determinant(So)$modulus)
  yo <- y[ok]
  k <- sum(ok)
  f <- function(q) {
    beta <- q$beta
    sig <- exp(q$lsig)
    nu <- 1 + exp(q$lnu)
    r <- yo - as.vector(Bo %*% (X %*% beta + Z %*% q$b))
    Q <- sum(r * as.vector(Soi %*% r)) / sig^2
    ll <- lgamma((nu + k) / 2) - lgamma(nu / 2) - 0.5 * k * log(nu * pi) -
      0.5 * (ldS + 2 * k * log(sig)) - 0.5 * (nu + k) * log1p(Q / nu)
    -(ll + sum(RTMB::dnorm(q$b, 0, exp(q$lsd), log = TRUE)))
  }
  obj <- RTMB::MakeADFun(f, list(beta = p[1:2], lsig = p[3], lnu = p[4],
                                 lsd = p[5], b = numeric(ncol(Z))),
                         random = "b", silent = TRUE)
  -obj$fn(p[1:5])
}

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
cases <- list(
  list(type = "fcor", f = yf | mi() ~ x + (1 | h) + fcor(M), y = d$yf),
  list(type = "error", f = ye | mi() ~ x + (1 | h) + sar(W, type = "error"),
       y = d$ye),
  list(type = "lag", f = yl | mi() ~ x + (1 | h) + sar(W, type = "lag"),
       y = d$yl)
)
worst <- 0
for (cs in cases) {
  fit <- quiet(frm(cs$f, data = d, data2 = list(M = M, W = W),
                   family = student()))
  lab <- ns$outer_par_labels(fit)
  cat(sprintf("\n=== %s: labels %s\n", cs$type, paste(lab, collapse = " ")))
  inu <- grep("^nu", lab)
  iac <- which(names(fit$opt$par) == "thetaac")
  ac <- fit$frame$autocor[[1L]]
  for (k in 0:2) {
    p <- fit$opt$par + 0.05 * k * seq_along(fit$opt$par) /
      length(fit$opt$par)
    p[inu] <- c(log(2), log(5), log(20))[k + 1L]
    rho <- if (length(iac)) {
      ac$aux$lo + (ac$aux$hi - ac$aux$lo) / (1 + exp(-p[iac]))
    } else NA
    a <- -fit$obj$fn(p)
    b <- ref_ll(p[setdiff(seq_along(p), iac)], cs$type, rho, cs$y)
    worst <- max(worst, abs(a - b))
    cat(sprintf("  point %d (rho %.4f): frmtmb %.8f  reference %.8f  diff %+.3e\n",
                k, rho, a, b, a - b))
  }
  cat("  mi_laplace_exact", ns$mi_laplace_exact(fit$frame), "\n")
}
cat(sprintf("\nlargest |diff| over every point: %.3e\n", worst))
