# Review of lane illse2, item 5: the asymmetry guard (se_reparam_asym)
# where its condition should be absent: a statistically well-identified
# nonlinear body (a degree-6 raw polynomial beside A exp(-k z)) whose
# Hessian in the Jacobian's coordinates is mildly asymmetric. Exact
# reference: 256-bit Hessian of the gaussian nll at the fit's point.
#   ILLSE2_ARM=lane|base I2_ASYM=default|off Rscript dev/illse2-rev-asym.R <A> <seeds>
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
a <- commandArgs(TRUE)
A0 <- as.numeric(a[1])
seeds <- eval(parse(text = a[2]))
off <- Sys.getenv("I2_ASYM", "default") == "off"
if (arm == "lane" && off) {
  utils::assignInNamespace("se_reparam_asym", Inf, "frmtmb")
}
deg <- 6L
bn <- paste0("b", 0:deg)
body <- paste(c("b0", paste0("b", 1:deg, " * x", 1:deg),
                "A * exp(-k * z)"), collapse = " + ")
fo <- bf(stats::as.formula(paste("y ~", body)),
         stats::as.formula(paste(paste(c(bn, "A", "k"), collapse = " + "),
                                 "~ 1")), nl = TRUE)
for (seed in seeds) {
  set.seed(seed)
  n <- 2000
  d <- data.frame(x = stats::runif(n, 1, 2), z = stats::runif(n, 0, 10))
  for (k in 1:deg) d[[paste0("x", k)]] <- d$x^k
  d$y <- sin(3 * d$x) + A0 * exp(-0.5 * d$z) + stats::rnorm(n, 0, 0.01)
  l0 <- stats::lm(stats::as.formula(paste("y ~", paste0("x", 1:deg,
                                    collapse = " + "))), data = d)
  st <- list(beta = c(unname(stats::coef(l0)), A0, 0.5))
  r <- cap(frm(fo, data = d, start = st))
  if (is_err(r$v)) { cat(seed, "ERROR", r$v, "\n"); next }
  f <- r$v
  p <- f$opt$par
  X <- stats::model.matrix(stats::as.formula(paste("~", paste0("x", 1:deg,
                                    collapse = " + "))), d)
  Xm <- mM(X); zm <- m(d$z); ym <- m(d$y)
  nb <- deg + 1L
  nll <- function(th) {
    mu <- Xm %*% col1(th[seq_len(nb)])
    mu <- as.vector(mu) + th[nb + 1L] * exp(-th[nb + 2L] * zm)
    s2 <- exp(2 * th[nb + 3L])
    (n * log(s2) + sum((ym - mu)^2) / s2) / 2
  }
  # sanity: the parameter order of the outer vector
  pn <- names(p)
  He <- mp_hessian(nll, p)
  Ve <- mp_inv(He)
  dg <- sqrt(num(diag(Ve)))
  V <- suppressWarnings(vcov(f, full = TRUE))
  E <- if (all(is.finite(V))) max(abs(num(mM(V) - Ve)) / outer(dg, dg)) else
    NaN
  # the asymmetry the guard reads, measured as se_reparam_cov() does
  asy <- NA
  if (arm == "lane") {
    M <- tryCatch(ns$reparam_basis(f$frame, names(p), est = f$estimates),
                  error = function(e) NULL)
    if (!is.null(M)) {
      h <- ns$reparam_hessian_refined(f$obj, p, M)
      if (!is.null(h)) asy <- max(h$E) / max(abs(h$H))
    }
  }
  tk <- p[nb + 2L] / dg[nb + 2L]
  ev <- eigen(cov2cor(num(Ve)), only.values = TRUE)$values
  cat(sprintf(paste0("A %g seed %d | %s asym-guard %s | code %d | kappa %.1e | t(k) %.0f | ",
                     "asym %.2e | lost %s | root %s | cor %.2e | warn %s\n"),
              A0, seed, arm, if (off) "off" else "on", f$opt$convergence,
              max(ev) / min(ev), tk, asy, lost_of(f),
              !is.null(ns$sdr_of(f)$cov_fixed_root), E, wtag(r$w)))
}
