# Lane optrap, consolidation fixes: a candidate, measured in-process.
# Newton steps for a fit with random effects, on the outer parameters
# of the Laplace objective, with the Hessian from central differences
# of its (exact) gradient (stats::optimHess(), as sdreport() takes it),
# in natural units; a step needs a Cholesky factor, a predicted gain
# above 1e-6, and an actual gain of at least 0.25 of the prediction.
# Raw polynomials with (1 | g) (dev/optrap-c3-repoly.R's data), against
# the exact profile maximum.
#   Rscript dev/optrap-c3-repolish.R <lib> [degrees]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
degs <- if (length(a) > 1) eval(parse(text = a[2])) else 6:7
ll_at <- function(y, X, Z, sig, tau) {
  V <- sig^2 * diag(length(y)) + tau^2 * tcrossprod(Z)
  R <- chol(V)
  yw <- backsolve(R, y, transpose = TRUE)
  Xw <- backsolve(R, X, transpose = TRUE)
  r <- qr.qty(qr(Xw, LAPACK = TRUE), yw)[-seq_len(ncol(X))]
  -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(R))) + sum(r^2))
}
polish_re <- function(obj, opt, u) {
  p <- opt$par; v <- opt$objective; steps <- 0L
  for (it in 1:5) {
    g <- drop(obj$gr(p))
    H <- stats::optimHess(p, obj$fn, obj$gr)
    gs <- g * u; Hs <- H * outer(u, u)
    R <- tryCatch(chol((Hs + t(Hs)) / 2), error = function(e) NULL)
    if (is.null(R)) break
    st <- -drop(backsolve(R, forwardsolve(t(R), gs)))
    gain <- -sum(gs * st) / 2
    if (!is.finite(gain) || gain <= 1e-6) break
    q <- p + st * u
    vq <- obj$fn(q)
    if (!is.finite(vq) || vq >= v || v - vq < 0.25 * gain) break
    p <- q; v <- vq; steps <- steps + 1L
  }
  list(par = p, objective = v, steps = steps)
}
for (deg in degs) {
  r <- NULL
  for (s in 1:20) {
    set.seed(1000 * deg + s)
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
    for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
    fo <- as.formula(paste("y ~", paste0("x", 1:deg, collapse = " + "),
                           "+ (1 | g)"))
    X <- cbind(1, as.matrix(d[paste0("x", 1:deg)]))
    Z <- model.matrix(~ 0 + g, d)
    f <- suppressWarnings(suppressMessages(frm(fo, data = d)))
    prof <- function(p) -ll_at(d$y, X, Z, exp(p[1]), exp(p[2]))
    o <- optim(log(c(0.2, 0.3)), prof, method = "Nelder-Mead",
               control = list(reltol = 1e-14, maxit = 5000))
    o <- optim(o$par, prof, method = "BFGS", control = list(reltol = 1e-15))
    best <- -o$value
    u <- f$par_units %||% rep(1, length(f$opt$par))
    t0 <- proc.time()[[3]]
    pp <- polish_re(f$obj, f$opt, u)
    r <- rbind(r, data.frame(seed = s, before = -f$opt$objective - best,
                             after = -pp$objective - best, steps = pp$steps,
                             secs = proc.time()[[3]] - t0))
  }
  cat(sprintf("deg %d | short before %2d (min %.3g) | after %2d (min %.3g) | polished %d | secs max %.2f\n",
              deg, sum(r$before < -1e-4), min(r$before), sum(r$after < -1e-4),
              min(r$after), sum(r$steps > 0), max(r$secs)))
}
