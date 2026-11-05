# Short check 2, item 3: the cost of one correction at large n, dense v, against the
# fit. Poisson, one series and many groups.
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
options(warn = -1)
for (sh in list(c(1, 1000), c(1, 3000), c(100, 30), c(1, 6000))) {
  G <- sh[1]; K <- sh[2]; n <- G * K
  set.seed(1)
  d <- data.frame(t = rep(seq_len(K), G), g = factor(rep(seq_len(G), each = K)), x = rnorm(n))
  e <- as.vector(replicate(G, arima.sim(list(ar = .5), K, sd = .5)))
  d$y <- rpois(n, exp(.5 + .3 * d$x + e))
  t0 <- proc.time()[[3]]
  f <- frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d)
  tf <- proc.time()[[3]] - t0
  env <- f$obj$env
  bk <- fi("frame_latres")(f$frame)[[1]]
  t0 <- proc.time()[[3]]
  f$obj$fn(f$opt$par)
  H <- env$spHess(env$last.par, random = TRUE)
  pj <- fi("latres_jacobian")(bk, f$estimates$theta[bk$theta_idx])
  A <- Matrix::sparseMatrix(i = pj$i, j = pj$j, x = pj$x, dims = c(n, length(env$random)))
  Ch <- Matrix::Cholesky(H, LDL = FALSE)
  Vv <- as.matrix(A %*% Matrix::solve(Ch, Matrix::t(A)))
  tv <- proc.time()[[3]] - t0
  cat(sprintf("G %d K %d n %d: fit %.2f s; one inner solve and a dense v by sparse factor %.2f s, %.0f MB; 5 such evaluations (theta only) %.1f s\n", G, K, n, tf, tv, 8 * n^2 / 1e6, 5 * tv))
}
