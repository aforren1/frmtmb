# Lane illse2, item 1: is the double-precision reference that the test
# uses (lm()'s design, centered columns) accurate against the 256-bit
# one? Degree-8 gaussian without random effects, seeds 1 to 3.
#   ILLSE2_ARM=lane Rscript dev/illse2-p1ref.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
xn <- paste0("x", 1:8)
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + ")))
mk <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + stats::rnorm(200, 0, 0.3)
  for (k in 1:8) d[[xn[k]]] <- d$x^k
  d
}
for (seed in 1:3) {
  d <- mk(seed)
  fit <- suppressWarnings(frm(fo, data = d))
  p <- fit$opt$par
  X <- stats::model.matrix(fo, d)
  # Cook's: d' H d with H the exact gaussian Hessian in (beta, log s)
  inf <- stats::influence(fit, groups = "g", data = d)
  D <- sweep(inf$fixed, 2, inf$fixed_full)
  b <- p[1:9]; ls <- p[10]; s2 <- exp(2 * ls)
  r <- d$y - drop(X %*% b)
  cd_dbl <- apply(D, 1, function(di) {
    u <- drop(X %*% di[1:9])
    (sum(u^2) + 4 * di[10] * sum(u * r) + 2 * di[10]^2 * sum(r^2)) / s2
  }) / 10
  He <- gauss_hess(X, d$y, p)
  cd_mp <- apply(D, 1, function(di) mp_quad(He, di)) / 10
  cat(sprintf("seed %d cooks: double ref / 256-bit ref %s\n", seed,
              ratio_str(cd_dbl / cd_mp)))
  # Wald form over x1..x8 with the sub-block of V
  Ve <- mp_inv(He)
  q <- b[2:9]
  ex <- mp_quad(mp_inv(Ve[2:9, 2:9]), q)
  Xc <- scale(X[, 2:9], scale = FALSE)
  dbl <- sum(drop(Xc %*% q)^2) / s2
  cat(sprintf("seed %d wald x1..x8: double ref / 256-bit ref %.12f\n",
              seed, dbl / ex))
}
