# Re-check of lane illse2, finding 3: which columns reparam_hessian()
# takes as exact (lin_pos), and the route's accuracy, on gaussian models
# with a log-link sigma formula, an offset, and weights. Exact
# reference: 256-bit Hessian of the written-out nll at the fit's point.
#   ILLSE2_ARM=lane Rscript dev/illse2-rev-lincls.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
deg <- 7L
xn <- paste0("x", 1:deg)
rhs <- paste(xn, collapse = " + ")
for (case in c("sigma_z", "offset", "weights")) {
  set.seed(1)
  n <- 300
  d <- data.frame(x = runif(n, 1, 2), z = rnorm(n), off = runif(n),
                  w = runif(n, 0.5, 2))
  for (k in 1:deg) d[[xn[k]]] <- d$x^k
  sz <- if (case == "sigma_z") exp(-1.2 + 0.4 * d$z) else 0.3
  d$y <- sin(3 * d$x) + (if (case == "offset") d$off else 0) +
    rnorm(n, 0, sz)
  fo <- switch(case,
    sigma_z = bf(as.formula(paste("y ~", rhs)), sigma ~ z),
    offset = bf(as.formula(paste("y ~", rhs, "+ offset(off)"))),
    weights = bf(as.formula(paste("y | weights(w) ~", rhs))))
  r <- cap(frm(fo, data = d))
  f <- r$v
  p <- f$opt$par
  X <- model.matrix(as.formula(paste("~", rhs)), d)
  Xm <- mM(X); ym <- m(d$y); zm <- m(d$z); om <- m(d$off); wm <- m(d$w)
  nb <- deg + 1L
  nll <- function(th) {
    mu <- as.vector(Xm %*% col1(th[seq_len(nb)]))
    if (case == "offset") mu <- mu + om
    ls <- if (case == "sigma_z") th[nb + 1L] + th[nb + 2L] * zm else
      th[nb + 1L] + 0 * ym
    e <- ls + (ym - mu)^2 / (2 * exp(2 * ls))
    if (case == "weights") sum(wm * e) else sum(e)
  }
  He <- mp_hessian(nll, p)
  Ve <- mp_inv(He)
  dg <- sqrt(num(diag(Ve)))
  V <- suppressWarnings(vcov(f, full = TRUE))
  E <- if (all(is.finite(V))) max(abs(num(mM(V) - Ve)) / outer(dg, dg)) else
    NaN
  M <- ns$reparam_basis(f$frame, names(p), est = f$estimates)
  g0 <- max(abs(num(mp_grad(nll, p))))
  cat(sprintf("%-8s | code %d | pos %s | lin_pos %s | root %s | cor %.2e | max|grad| %.1e | warn %s\n",
              case, f$opt$convergence,
              paste(attr(M, "pos"), collapse = ","),
              paste(attr(M, "lin_pos"), collapse = ","),
              !is.null(ns$sdr_of(f)$cov_fixed_root), E, g0, wtag(r$w)))
}
