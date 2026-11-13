# Review of lane gpstop, item 5033: the student() step rows. At the default
# optimum and the reference (theta start rho0/sd0, the lane's best), and on
# the segment between them: the Laplace value, the outer gradient, the inner
# Hessian's smallest eigenvalue, and the marginal by importance sampling with
# two proposals (t df 4 floor 1, the lane's; t df 3 floor 0.1, wider in the
# flat direction), two seeds each.
# Args: rows as n:rho0:sd0 (comma list), draws.
# Data: dev/arstep-ng.R's student step rows, set.seed(74000 + 31 seed + n),
# seed 1.
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/reviews/gpstop-review-lib.R")
rv_libs(if (length(a) > 2) a[3] else "lane")
rows <- strsplit(a[1], ",")[[1]]
ndraw <- as.integer(a[2])
inner <- function(obj, par) {
  v <- -obj$fn(par)
  env <- obj$env
  full <- env$last.par
  H <- as.matrix(env$spHess(full, random = TRUE))
  H <- 0.5 * (H + t(H))
  list(lap = v, full = full, H = H, ev = eigen(H, symmetric = TRUE,
                                               only.values = TRUE)$values)
}
is_ll <- function(obj, st, floor, df, seed) {
  env <- obj$env
  r <- env$random
  e <- eigen(st$H, symmetric = TRUE)
  P <- e$vectors %*% diag(pmax(e$values, floor)) %*% t(e$vectors)
  R <- chol(0.5 * (P + t(P)))
  k <- length(r)
  set.seed(seed)
  z <- matrix(rnorm(ndraw * k), ndraw, k)
  s <- sqrt(rchisq(ndraw, df) / df)
  dev <- t(backsolve(R, t(z))) / s
  lq <- lgamma((df + k) / 2) - lgamma(df / 2) - k / 2 * log(df * pi) +
    sum(log(diag(R))) - (df + k) / 2 * log1p(rowSums(z^2) / s^2 / df)
  u0 <- st$full[r]
  x <- st$full
  lj <- vapply(seq_len(ndraw), function(i) {
    x[r] <- u0 + dev[i, ]
    -env$f(x, order = 0)
  }, 0)
  lw <- lj - lq
  m <- max(lw)
  w <- exp(lw - m)
  c(is = m + log(mean(w)), ess = sum(w)^2 / sum(w^2))
}
for (rw in rows) {
  p <- strsplit(rw, ":")[[1]]
  n <- as.integer(p[1]); rho0 <- as.numeric(p[2]); sd0 <- as.numeric(p[3])
  set.seed(74000 + 1 * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  d <- data.frame(x = x, y = 2 * (x > 5) - 1 + 0.02 * rt(n, 5))
  qd <- rv_fit(y ~ gp(x), data = d, family = student())
  qr <- rv_fit(y ~ gp(x), data = d, family = student(),
               start = list(theta = c(log(sd0), log(rho0))))
  fd <- qd$fit; fr <- qr$fit
  cat(sprintf("STU n %d default ll %.4f conv %d nwarn %d | ref ll %.4f conv %d nwarn %d\n",
              n, -fd$opt$objective, fd$opt$convergence, length(qd$warn),
              -fr$opt$objective, fr$opt$convergence, length(qr$warn)))
  cat("STUP n", n, "default", paste(names(fd$opt$par), signif(fd$opt$par, 5),
                                    collapse = " "), "\n")
  cat("STUP n", n, "ref", paste(names(fr$opt$par), signif(fr$opt$par, 5),
                                collapse = " "), "\n")
  obj <- fd$obj
  p0 <- fd$opt$par
  p1 <- fr$opt$par
  for (tt in c(0, 0.5, 0.9, 1)) {
    pt <- p0 + tt * (p1 - p0)
    st <- inner(obj, pt)
    g <- tryCatch(max(abs(obj$gr(pt))), error = function(e) NA)
    st <- inner(obj, pt)
    out <- sprintf("STU n %d t %.2f (default obj) laplace %.4f gmax %.2e ev_min %.3g half-log %.3f",
                   n, tt, st$lap, g, min(st$ev), -log(min(st$ev)) / 2)
    for (cfg in list(c(1, 4, 9001), c(1, 4, 9002), c(0.1, 3, 9001),
                     c(0.1, 3, 9002))) {
      v <- is_ll(obj, st, cfg[1], cfg[2], cfg[3])
      out <- paste0(out, sprintf(" | IS fl%g df%g s%d %.4f ess %.0f",
                                 cfg[1], cfg[2], cfg[3], v[["is"]], v[["ess"]]))
    }
    cat(out, "\n")
    flush.console()
  }
  # the reference on its own tape, as the lane read it
  st <- inner(fr$obj, p1)
  out <- sprintf("STU n %d ref-own laplace %.4f gmax %.2e ev_min %.3g half-log %.3f",
                 n, st$lap, max(abs(fr$obj$gr(p1))), min(st$ev),
                 -log(min(st$ev)) / 2)
  st <- inner(fr$obj, p1)
  for (cfg in list(c(1, 4, 9001), c(0.1, 3, 9001), c(0.1, 3, 9002))) {
    v <- is_ll(fr$obj, st, cfg[1], cfg[2], cfg[3])
    out <- paste0(out, sprintf(" | IS fl%g df%g s%d %.4f ess %.0f",
                               cfg[1], cfg[2], cfg[3], v[["is"]], v[["ess"]]))
  }
  cat(out, "\n")
  flush.console()
}
cat("DONE\n")
