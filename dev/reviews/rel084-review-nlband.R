# rel084 release review, nlband fix round 2: the hole refit is kept where
# its coefficient's SE is below the span of the rows' poles, every run that
# beats the fit is a candidate (highest first), the first refused one is
# said (frmtmb_nl_pole_refused), and a kept refit says its own boundary
# messages. New designs: pole at 1.3 (not 2), x on U(0, 3), noise gaussian,
# student t3, a separated bernoulli, and a small loose design (n 50).
#   REVLIB=r23|r22 Rscript dev/reviews/rel084-review-nlband.R
source("dev/reviews/rel084-review-common.R")
options(width = 170)
form <- bf(y ~ a / (x + b), a + b ~ 1, nl = TRUE)
prof1 <- function(x, y, b) {
  vapply(b, function(bb) {
    z <- 1 / (x + bb)
    a <- sum(y * z) / sum(z * z)
    rss <- sum((y - a * z)^2)
    n <- length(y)
    -n / 2 * (log(2 * pi * rss / n) + 1)
  }, 1)
}
glob1 <- function(x, y) {
  w <- sort(unique(-x))
  r <- diff(range(w))
  u <- seq(0, 1, length.out = 17)[2:16]
  b <- c(as.vector(outer(u, diff(w)) + rep(w[-length(w)], each = 15)),
         w[1] - r * 10^seq(-4, 2, length.out = 60),
         w[length(w)] + r * 10^seq(-4, 2, length.out = 60))
  b <- sort(b)
  p <- prof1(x, y, b)
  p[!is.finite(p)] <- -Inf
  k <- which.max(p)
  o <- optimize(function(bb) prof1(x, y, bb), c(b[max(1, k - 1)],
                b[min(length(b), k + 1)]), maximum = TRUE, tol = 1e-12)
  c(max(o$objective, p[k]), if (o$objective >= p[k]) o$maximum else b[k])
}
pole <- 1.3
mk <- function(seed, n = 300, hole = 0.25, lo = 0, hi = 3) {
  set.seed(seed)
  x <- runif(n, lo, hi)
  x[abs(x - pole) > hole]
}
ll <- function(r) if (is.null(r$v)) NA else as.numeric(logLik(r$v))
mcls <- function(r) {
  m <- unique(unlist(lapply(r$said, function(s) s$cls[[1]])))
  if (length(m)) paste(m, collapse = "+") else "-"
}
bpar <- function(r) if (is.null(r$v)) NA else r$v$opt$par[[2]]
for (s in 1:4) {
  x <- mk(900 + s)
  set.seed(950 + s)
  d <- data.frame(x = x, y = 1 / (x - pole) + rnorm(length(x), 0, 0.2))
  g <- glob1(d$x, d$y)
  r <- conds_of(frm(form, data = d))
  cat(sprintf("GAUSS seed %d | ll %.6f b %.4f code %d [%s] | profile global %.6f at b %.4f\n",
              s, ll(r), bpar(r), r$v$opt$convergence, mcls(r), g[1], g[2]))
}
for (s in 1:4) {
  x <- mk(910 + s)
  set.seed(960 + s)
  d <- data.frame(x = x, y = 1 / (x - pole) + 0.2 * rt(length(x), 3))
  r <- conds_of(frm(form, data = d, family = student()))
  q <- conds_of(frm(form, data = d, family = student(),
                    start = list(beta = c(1, -pole)),
                    control = frmtmb_control(nl_starts = FALSE)))
  cat(sprintf("STUD seed %d | default ll %.6f b %.4f [%s] | from truth ll %.6f b %.4f [%s]\n",
              s, ll(r), bpar(r), mcls(r), ll(q), bpar(q), mcls(q)))
}
for (s in 1:4) {
  x <- mk(920 + s)
  d <- data.frame(x = x, y = as.integer(x > pole))
  r <- conds_of(frm(form, data = d, family = bernoulli()))
  cat(sprintf("SEP seed %d | ll %.6g b %.4f code %s [%s]\n", s, ll(r), bpar(r),
              if (is.null(r$v)) NA else r$v$opt$convergence, mcls(r)))
  rm <- msg_of(r, "frmtmb_nl_pole_refused")
  if (!is.na(rm)) cat("   REFUSED MSG:", substr(gsub("\n", " ", rm), 1, 300), "\n")
}
# loose: few rows, wide noise, narrow hole
for (s in 1:12) {
  set.seed(930 + s)
  x <- runif(400, 0, 3)
  x <- x[abs(x - pole) > 0.5][1:50]
  d <- data.frame(x = x, y = 0.2 / (x - pole) + rnorm(50, 0, 1.5))
  g <- glob1(d$x, d$y)
  r <- conds_of(frm(form, data = d))
  cat(sprintf("LOOSE seed %d | ll %.6f b %.4f code %s [%s] | global %.6f at %.4f\n",
              s, ll(r), bpar(r), if (is.null(r$v)) NA else r$v$opt$convergence,
              mcls(r), g[1], g[2]))
}
cat("DONE\n")
