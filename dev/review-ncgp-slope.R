# Item 3: the compressed Hessian with a correlated block ahead of the field, b first as
# shipped. Gradient plain, compressed, central difference; the fit on each tape; and
# the closed form as truth. Seeds 600 + s, the same designs as review-ncgp-compress.R.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
orig_scope <- ns$gp_tape_scope
plain_scope <- function(frame, expr) expr
setscope <- function(f) { unlockBinding("gp_tape_scope", ns); assign("gp_tape_scope", f, envir = ns); lockBinding("gp_tape_scope", ns) }
# Dense marginal of y ~ gp(x) + (z | g) with X intercept + optional covariate, parameters
# (log sigma, log sd, log rho, log sg1, log sg2, atanh r).
cf_dense <- function(d, frm_pars) {
  n <- nrow(d)
  D2 <- outer(d$x, d$x, "-")^2
  gi <- as.integer(d$g)
  same <- outer(gi, gi, "==") * 1
  f <- function(p) {
    sigma <- exp(p[1]); sd <- exp(p[2]); rho <- exp(p[3]); s1 <- exp(p[4]); s2 <- exp(p[5]); r <- tanh(p[6])
    K <- sd^2 * (exp(-D2 / (2 * rho^2)) + diag(1e-12, n))
    G <- same * (s1^2 + r * s1 * s2 * outer(d$z, rep(1, n)) + r * s1 * s2 * outer(rep(1, n), d$z) + s2^2 * outer(d$z, d$z))
    V <- K + G; diag(V) <- diag(V) + sigma^2
    R <- tryCatch(chol(V), error = function(e) NULL)
    if (is.null(R)) return(1e300)
    X <- matrix(1, n, 1)
    Wy <- backsolve(R, d$y, transpose = TRUE); WX <- backsolve(R, X, transpose = TRUE)
    b <- solve(crossprod(WX), crossprod(WX, Wy))
    e <- Wy - WX %*% b
    0.5 * sum(e^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
  }
  f
}
for (seed in 2:3) {
  set.seed(600 + seed)
  n <- 70
  d <- data.frame(x = runif(n, 0, 6), z = runif(n, 0, 6), g = factor(sample(letters[1:6], n, TRUE)),
                  f = factor(sample(c("u", "v"), n, TRUE)), xg = seq(0, 6, length.out = n),
                  xr = rep(sort(runif(14, 0, 6)), each = 5))
  sg <- c(0.3, 0.01, 1e-3, 0.1)[seed]
  d$y <- sin(2 * d$x) + 0.3 * d$z + rnorm(6)[as.integer(d$g)] * 0.5 + sin(d$xg) + sin(d$xr) + rnorm(n, sd = sg)
  for (fml in list(y ~ (z | g) + gp(x), y ~ gp(x) + (z | g))) {
    setscope(orig_scope)
    r1 <- rev_fit(fml, d)
    setscope(plain_scope)
    r0 <- rev_fit(fml, d)
    setscope(orig_scope)
    f1 <- r1$fit; f0 <- r0$fit
    p <- f1$opt$par
    g1 <- f1$obj$gr(p); g0 <- f0$obj$gr(p)
    fd <- vapply(seq_along(p), function(i) { h <- 1e-4; e <- replace(numeric(length(p)), i, h); (f1$obj$fn(p + e) - f1$obj$fn(p - e)) / (2 * h) }, 0)
    cat(sprintf("RSL seed %d sg %g %s b order %s\n", seed, sg, deparse(fml), paste(vapply(f1$frame$re_blocks, function(k) k$covstruct, ""), collapse = ",")))
    cat(sprintf("RSL   compressed ll %.5f conv %d w %d msg %s %s | plain ll %.5f conv %d w %d msg %s\n", logLik(f1), f1$opt$convergence, length(r1$warn), f1$opt$message, paste(round(f1$opt$par[names(p) == "theta"], 3), collapse = ","),
                logLik(f0), f0$opt$convergence, length(r0$warn), f0$opt$message))
    cat("RSL   names ", names(p), "\n")
    cat("RSL   gr comp ", sprintf("%.5f", g1), "\n")
    cat("RSL   gr plain", sprintf("%.5f", g0), "\n")
    cat("RSL   central ", sprintf("%.5f", fd), "\n")
    flush.console()
  }
  obj <- cf_dense(d)
  best <- NULL
  th <- f1$estimates$theta
  bks <- f1$frame$re_blocks
  gpb <- Filter(function(k) identical(k$covstruct, "gp"), bks)[[1]]
  usb <- Filter(function(k) !identical(k$covstruct, "gp"), bks)[[1]]
  print(th); print(usb$covstruct)
  starts <- list(c(log(sg), 0.5, 0, 0, -1, 0), c(log(exp(f1$estimates$betad[[1]])), th[gpb$theta_idx], th[usb$theta_idx][1:2], 0),
                 c(log(sg), 1, 0.5, -1, -2, 0.5), c(log(0.5), 0, -1, 0, -2, 0))
  for (st in starts) {
    o <- nlminb(st, obj, control = list(rel.tol = 1e-12, iter.max = 800, eval.max = 1200))
    if (is.null(best) || o$objective < best$objective) best <- o
    cat(sprintf("RSL   cf from start: ll %.5f\n", -o$objective))
  }
  cat(sprintf("RSL cf max ll %.5f sigma %.4g sd %.4g rho %.4g sg1 %.4g sg2 %.4g r %.3f | shipped gap %.3e\n", -best$objective, exp(best$par[1]), exp(best$par[2]),
              exp(best$par[3]), exp(best$par[4]), exp(best$par[5]), tanh(best$par[6]), as.numeric(logLik(f1)) + best$objective))
}
cat("DONE\n")
