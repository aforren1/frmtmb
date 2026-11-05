# Short check of fix round 2, item 3: Q of a category on a new me() row, against brms's
# quantiles at fixed values, on the lane's build and on mutant band_cat. Seed 7971, the data
# of review-predread-r1-check1.R. Also plain fits, against brms's sampled draws, on which
# the cumulative change also landed.
source("C:/Users/adf44/source/r/frmtmb-wt-predread/dev/review-predread-lib.R")
set.seed(7971)
n <- 400
d <- data.frame(z = rnorm(n), lat = rnorm(n), sx = 0.4)
d$xo <- d$lat + rnorm(n, 0, 0.4)
eta <- 0.4 + 0.7 * d$lat
d$ord <- as.integer(cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf), labels = FALSE))
d$cat <- as.integer(cut(eta + rlogis(n), c(-Inf, 0, 1.5, Inf), labels = FALSE))
nd <- data.frame(xo = c(1, 1), sx = c(0.4, 1.2))
fmt <- function(v) paste(sprintf("%.3f", v), collapse = " ")
for (sp in list(c("bf(ord ~ me(xo, sx))", "cumulative()"), c("bf(cat ~ me(xo, sx))", "categorical()"))) {
  pr <- rev_pair(sp[1], sp[2], d, ndraws = 4000)
  fv <- quiet(fitted(pr$fit, newdata = nd))
  ds <- rev_draws_at(pr$fit, 8000L)
  set.seed(5)
  fe <- quiet(frmtmb.sample::posterior_epred(ds, newdata = nd))
  b <- do.call(abind_rows <- function(...) {
    l <- list(...)
    a <- array(NA_real_, c(sum(sapply(l, function(x) dim(x)[1L])), dim(l[[1L]])[2:3]))
    p <- 0L
    for (x in l) {
      a[p + seq_len(dim(x)[1L]), , ] <- x
      p <- p + dim(x)[1L]
    }
    a
  }, lapply(1:10, function(i) posterior_epred(pr$fx, newdata = nd)))
  cat("\n===", sp, "\n")
  for (r in 1:2) for (k in seq_len(dim(fv)[3L])) {
    cat(sprintf("row %d cat %d  Q %s | draws %s | brms %s\n", r, k, fmt(fv[r, c("Q2.5", "Q97.5"), k]),
                fmt(quantile(fe[, r, k], c(.025, .975))), fmt(quantile(b[, r, k], c(.025, .975)))))
  }
}
cat("\n=== plain cumulative and categorical, brms sampled\n")
d$x <- d$lat
for (sp in list(c("bf(ord ~ x)", "cumulative()"), c("bf(cat ~ x)", "categorical()"))) {
  ff <- quiet(frm(rev_ev(sp[1], "frmtmb"), family = rev_ev(sp[2], "frmtmb"), data = d))
  bb <- quiet(brms::brm(rev_ev(sp[1], "brms"), family = rev_ev(sp[2], "brms"), data = d, chains = 2, iter = 1500, refresh = 0, seed = 3, silent = 2))
  g <- data.frame(x = c(-2, 0, 2))
  fv <- quiet(fitted(ff, newdata = g))
  fb <- fitted(bb, newdata = g)
  cat(sp, "\n")
  for (r in 1:3) for (k in seq_len(dim(fv)[3L])) {
    cat(sprintf("row %d cat %d  frmtmb %.3f %.3f to %.3f | brms %.3f %.3f to %.3f\n", r, k, fv[r, 1, k], fv[r, 3, k], fv[r, 4, k], fb[r, 1, k], fb[r, 3, k], fb[r, 4, k]))
  }
}
