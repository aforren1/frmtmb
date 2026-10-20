# Lane illse2, item 2: the 80 sd-at-zero fits of dev/illse-rev-edge4.R
# (raw polynomials with (1 | g), no group effect, n 300, seeds 1-20),
# with each fit's shortfall against the exact maximum: the 256-bit
# profile of the log-likelihood in tau = sd / sigma (beta GLS, sigma^2
# the GLS residual sum of squares over n).
#   ILLSE2_ARM=lane|base Rscript dev/illse2-edge4.R <degrees> [seeds]
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
degs <- eval(parse(text = a[1]))
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:20
prof_max <- function(X, y, g) {
  n <- length(y)
  Xm <- mM(X); ym <- m(y)
  XX <- crossprod(Xm); Xy <- crossprod(Xm, ym); yy <- sum(ym^2)
  Sg <- lapply(split(seq_len(n), g), function(ix) {
    list(s = colSums(Xm[ix, , drop = FALSE]), sy = sum(ym[ix]),
         n = length(ix))
  })
  prof <- function(tau) {
    t2 <- m(tau)^2
    A <- XX; bv <- Xy; c0 <- yy; ld <- m(0)
    for (st in Sg) {
      cg <- t2 / (1 + st$n * t2)
      sm <- col1(st$s)
      A <- A - cg * (sm %*% t(sm))
      bv <- bv - cg * st$sy * sm
      c0 <- c0 - cg * st$sy^2
      ld <- ld + log(1 + st$n * t2)
    }
    bh <- mp_inv(A) %*% bv
    s2 <- (c0 - sum(bv * bh)) / n
    num(-(n / 2) * log(2 * Rmpfr::Const("pi", PREC) * s2) - n / 2 - ld / 2)
  }
  l0 <- prof(0)
  o <- stats::optimize(function(t) -prof(t), c(0, 0.3), tol = 1e-7)
  if (-o$objective > l0) c(tau = o$minimum, ll = -o$objective) else
    c(tau = 0, ll = l0)
}
for (deg in degs) for (seed in seeds) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(1:20, length.out = n)))
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                                "+ (1 | g)"))
  r <- cap(frm(fo, data = d))
  f <- r$v
  if (is_err(f)) { cat(deg, seed, "ERROR\n"); next }
  sd <- exp(f$opt$par[names(f$opt$par) == "theta"])
  l <- ns$sdr_of(f)$se_lost
  V <- suppressWarnings(vcov(f, full = TRUE))
  se_t <- sqrt(V[nrow(V), nrow(V)])
  nb <- sum(!is.finite(diag(suppressWarnings(vcov(f)))))
  X <- stats::model.matrix(stats::as.formula(paste("~",
                                     paste(xn, collapse = " + "))), d)
  pm <- prof_max(X, d$y, d$g)
  cat(sprintf(paste0("deg %d s%2d | code %d | sd %.1e | route %s | ",
                     "theta_1 %s | n lost %d | coef NaN %d | ",
                     "short %+.1e (exact tau %.3g) | escape %s | warn %s\n"),
              deg, seed, f$opt$convergence, sd,
              !is.null(ns$sdr_of(f)$cov_fixed_root),
              if ("theta_1" %in% names(l)) paste0("lost:", l[["theta_1"]])
              else sprintf("SE %.2e", se_t),
              length(l), nb, pm[["ll"]] - as.numeric(logLik(f)),
              pm[["tau"]], if (is.null(f$opt$sd_escape)) "-" else
                sprintf("%.1e", f$opt$sd_escape[["gain"]]),
              wtag(r$w)))
}
