# Fix round 2, item 3: are the capped fits' answers usable? bern_slope_sep at data
# seeds 3, 5, 7, 9, 2: 64 and 512 draws, draw seeds 1 and 2, against my exact marginal
# at each fit's estimates, the Laplace fit, and how far the estimates spread.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
g1 <- gh_rule(41L)
gr <- as.matrix(expand.grid(1:41, 1:41))
xn <- matrix(g1$x[gr], nrow(gr), 2L)
lw <- rowSums(matrix(log(g1$w)[gr], nrow(gr), 2L)) + rowSums(xn^2)
for (ds in c(3L, 5L, 7L, 9L, 2L)) {
  set.seed(ds)
  g <- factor(rep(1:20, each = 4L))
  x <- rnorm(80)
  y <- as.integer(x + rnorm(20, 0, 2)[g] + rnorm(80, 0, .1) > 0)
  d <- data.frame(y, x, g)
  f <- bf(y ~ x + (x | g))
  idx <- split(seq_along(y), g)
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
      cc <- tryCatch(chol(optimHess(o$par, neg)), error = function(e) NULL)
      Lc <- if (is.null(cc)) t(ch) else backsolve(cc, diag(2))
      U <- matrix(o$par, nrow(xn), 2L, byrow = TRUE) + sqrt(2) * xn %*% t(Lc)
      a <- ll(U) + lpr(U) + lw
      tot <- tot - (max(a) + log(sum(exp(a - max(a)))) + log(2) + sum(log(abs(diag(Lc)))))
    }
    tot
  }
  lap <- suppressWarnings(frm(f, family = bernoulli(), data = d))
  v <- frmtmb:::varcorr_matrices(lap)[[1L]]
  cat(sprintf("data seed %d Laplace: conv %d sd %.2f %.2f cor %.3f slope %.3f, exact ll at it %.3f\n", ds, lap$opt$convergence,
              sqrt(v[1, 1]), sqrt(v[2, 2]), v[1, 2] / sqrt(v[1, 1] * v[2, 2]), lap$opt$par[[2]], -exact(unname(lap$opt$par))))
  for (N in c(64L, 512L)) for (s in 1:2) {
    w <- character(0)
    fi <- withCallingHandlers(
      tryCatch(frm(f, family = bernoulli(), data = d, importance = N, control = frmtmb_control(importance_seed = s)), error = function(e) e),
      warning = function(x) {
        w <<- c(w, substr(conditionMessage(x), 1, 40))
        invokeRestart("muffleWarning")
      })
    if (inherits(fi, "error")) {
      cat(sprintf("  N %d seed %d STOP %s\n", N, s, substr(conditionMessage(fi), 1, 100)))
      next
    }
    vv <- frmtmb:::varcorr_matrices(fi)[[1L]]
    ex <- -exact(unname(fi$opt$par))
    cat(sprintf("  N %d seed %d ll %.3f exact %.3f err %+.2f (%.1f MCSE) ess %.4f capped %s last move %.3g sd %.2f %.2f cor %.3f slope %.3f warns %d\n",
                N, s, as.numeric(logLik(fi)), ex, as.numeric(logLik(fi)) - ex, (as.numeric(logLik(fi)) - ex) / fi$importance$mcse,
                fi$importance$ess_min, fi$importance$capped, fi$importance$moved, sqrt(vv[1, 1]), sqrt(vv[2, 2]),
                vv[1, 2] / sqrt(vv[1, 1] * vv[2, 2]), fi$opt$par[[2]], length(w)))
  }
}
