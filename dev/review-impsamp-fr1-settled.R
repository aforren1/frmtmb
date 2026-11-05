# Fix round 1, item 4: is "Use importance = 128 or more" true? bern_slope_sep, data
# seed 3 (the lane's natural firing that names a count). Fit at 8 to 4096 draws over
# draw seeds; read each against my exact marginal at its own estimates; and the
# Laplace fit, which the message tells the user to check.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
set.seed(3)
g <- factor(rep(1:20, each = 4L))
x <- rnorm(80)
y <- as.integer(x + rnorm(20, 0, 2)[g] + rnorm(80, 0, .1) > 0)
d <- data.frame(y, x, g)
f <- bf(y ~ x + (x | g))
idx <- split(seq_along(y), g)
g1 <- gh_rule(41L)
gr <- as.matrix(expand.grid(1:41, 1:41))
xn <- matrix(g1$x[gr], nrow(gr), 2L)
lw <- rowSums(matrix(log(g1$w)[gr], nrow(gr), 2L)) + rowSums(xn^2)
exact <- function(p) {
  S <- as.matrix(us_sigma(p[3:5], 2L))
  ch <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(ch)) return(NA_real_)
  Pi <- chol2inv(ch)
  ld <- 2 * sum(log(diag(ch)))
  lpr <- function(U) -0.5 * rowSums((U %*% Pi) * U) - log(2 * pi) - 0.5 * ld
  tot <- 0
  for (ii in idx) {
    ll <- function(U) {
      eta <- outer(p[1] + p[2] * x[ii], rep(1, nrow(U))) + cbind(1, x[ii]) %*% t(U)
      colSums(matrix(dbinom(rep(y[ii], nrow(U)), 1, plogis(as.vector(eta)), log = TRUE), length(ii), nrow(U)))
    }
    neg <- function(u) {
      U <- matrix(u, 1L, 2L)
      -(ll(U) + lpr(U))
    }
    o <- optim(c(0, 0), neg, method = "BFGS", control = list(reltol = 1e-10, maxit = 300))
    h <- optimHess(o$par, neg)
    cc <- tryCatch(chol(h), error = function(e) NULL)
    Lc <- if (is.null(cc)) t(ch) else backsolve(cc, diag(2))
    U <- matrix(o$par, nrow(xn), 2L, byrow = TRUE) + sqrt(2) * xn %*% t(Lc)
    a <- ll(U) + lpr(U) + lw
    tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(2) + sum(log(abs(diag(Lc)))))
  }
  tot
}
lap <- suppressWarnings(frm(f, family = bernoulli(), data = d))
v <- frmtmb:::varcorr_matrices(lap)[[1L]]
cat(sprintf("Laplace: conv %d ll %.3f sd %.3f %.3f cor %.3f exact at its estimates %.3f\n", lap$opt$convergence,
            as.numeric(logLik(lap)), sqrt(v[1, 1]), sqrt(v[2, 2]), v[1, 2] / sqrt(v[1, 1] * v[2, 2]),
            -exact(unname(lap$opt$par))))
for (N in c(8L, 128L, 256L, 512L, 1024L, 2048L, 4096L)) for (s in c(1L, 2L, 3L)) {
  w <- character(0)
  fi <- withCallingHandlers(
    tryCatch(frm(f, family = bernoulli(), data = d, importance = N, control = frmtmb_control(importance_seed = s)),
             error = function(e) e),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 60))
      invokeRestart("muffleWarning")
    })
  if (inherits(fi, "error")) {
    m <- conditionMessage(fi)
    cat(sprintf("N %d seed %d STOP: %s\n", N, s, sub(".*of its draws, which is [^.]*[.0-9]* of [0-9]+. ", "", m)))
    next
  }
  vv <- frmtmb:::varcorr_matrices(fi)[[1L]]
  ex <- -exact(unname(fi$opt$par))
  cat(sprintf("N %d seed %d ll %.3f exact-at-fit %.3f err %+.3f (%.1f MCSE) mcse %.3f ess %.4f rounds %d capped %s sd %.2f %.2f cor %.3f warns %d %s\n",
              N, s, as.numeric(logLik(fi)), ex, as.numeric(logLik(fi)) - ex,
              (as.numeric(logLik(fi)) - ex) / fi$importance$mcse, fi$importance$mcse, fi$importance$ess_min,
              fi$importance$rounds, fi$importance$capped, sqrt(vv[1, 1]), sqrt(vv[2, 2]),
              vv[1, 2] / sqrt(vv[1, 1] * vv[2, 2]), length(w), paste(w, collapse = ";")))
}
