# Reviewer, punch round 2: log_pbeta_ad(), log_ibeta_pre() and
# log_ibeta_cf(xc =) inside log_ibeta_half().
#  A. value on a grid (shapes 1e-3 to 1e7, plus pairs at a b / (a + b) =
#     15 and 30 and a shape of 12) with x at m, both ties, the clamp and
#     blend edges of log_pbeta_ad(), the log1p/log blend edges of
#     log_ibeta_pre(), and the fraction/pbeta edges; Rmpfr at the worst
#  B. all third partials finite on the grid
#  C. derivatives d/dx to third order (closed forms from the density)
#     and d/da, d/db (central differences) in 256-bit Rmpfr, on a subset
#  D. continuity across every new edge, in x and in the shapes
#  E. the clamp: a fine scan of x across both ties
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(RTMB); library(Rmpfr) })
lib <- frmtmb:::log_ibeta_half
cat("frmtmb from", find.package("frmtmb"), "\n")
err <- function(v, ref) abs(v - ref) / pmax(1, abs(ref))
sdb <- function(a, b) sqrt(a * b / ((a + b)^2 * (a + b + 1)))
prec <- 256
mp_I <- function(x, a, b, maxit = 2e5) {
  x <- mpfr(x, prec); a <- mpfr(a, prec); b <- mpfr(b, prec)
  cf <- function(x, a, b) {
    qab <- a + b; qap <- a + 1; qam <- a - 1
    c <- mpfr(1, prec); d <- 1 / (1 - qab * x / qap); h <- d
    for (m in seq_len(maxit)) {
      m2 <- 2 * m
      aa <- m * (b - m) * x / ((qam + m2) * (a + m2))
      d <- 1 / (1 + aa * d); c <- 1 + aa / c; h <- h * d * c
      aa <- -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
      d <- 1 / (1 + aa * d); c <- 1 + aa / c; del <- d * c; h <- h * del
      if (abs(del - 1) < mpfr(1e-70, prec)) break
    }
    a * log(x) + b * log1p(-x) + lgamma(a + b) - lgamma(a) - lgamma(b) -
      log(a) + log(h)
  }
  if (x < (a + 1) / (a + b + 2)) cf(x, a, b) else log1p(-exp(cf(1 - x, b, a)))
}
# derivatives of g = log I_x(a, b): x to third order in closed form,
# a and b by central differences, all in mpfr
mp_der <- function(x, a, b) {
  X <- mpfr(x, prec); A <- mpfr(a, prec); B <- mpfr(b, prec)
  g <- mp_I(x, a, b)
  lf <- (A - 1) * log(X) + (B - 1) * log1p(-X) + lgamma(A + B) - lgamma(A) - lgamma(B)
  r <- exp(lf - g)
  L <- (A - 1) / X - (B - 1) / (1 - X)
  Lp <- -(A - 1) / X^2 - (B - 1) / (1 - X)^2
  g2 <- r * L - r^2
  g3 <- g2 * L + r * Lp - 2 * r * g2
  h <- mpfr(1e-30, prec)
  ga <- (mp_I(X, A * (1 + h), B) - mp_I(X, A * (1 - h), B)) / (2 * h * A)
  gb <- (mp_I(X, A, B * (1 + h)) - mp_I(X, A, B * (1 - h))) / (2 * h * B)
  as.numeric(c(g, r, g2, g3, ga, gb))
}
Tn <- MakeTape(function(p) lib(p[1], p[2], p[3]), c(0.2, 3, 7))
Jn <- Tn$jacfun(); Hn <- Jn$jacfun(); T3n <- Hn$jacfun()

# points for one shape pair: every edge the new code has
edge_x <- function(a, b) {
  s <- a + b; t1 <- a / s; gap <- b / (s * (s + 1)); m <- (a + 1) / (s + 2)
  sd <- sdb(a, b); mu <- t1; mub <- b / s
  x <- c(m = m, t1 = t1, tie2 = t1 + gap, lo = t1 + 0.25 * gap,
         hi = t1 + 0.75 * gap, mid = t1 + 0.5 * gap,
         t1m = t1 * (1 - 1e-15), t1p = t1 * (1 + 1e-15),
         tie2m = (t1 + gap) * (1 - 1e-15),
         u25m = mu * 0.75, u50m = mu * 0.5, u25p = mu * 1.25, u50p = mu * 1.5,
         v25m = mu - 0.25 * mub, v25p = mu + 0.25 * mub,
         clamp = mu * 0.25,
         e05 = m + 0.5 * sd, e15 = m + 1.5 * sd, e5m = m - 5 * sd,
         e6m = m - 6 * sd, e5p = m + 5 * sd, e6p = m + 6 * sd,
         e65m = m - 6.5 * sd)
  x[x > 0 & x < 0.5]
}
cat("== A. value ==\n")
sh <- 10^seq(-3, 7, by = 0.5)
P <- as.matrix(expand.grid(a = sh, b = sh))
extra <- rbind(c(30, 30), c(60, 60), c(16, 240), c(31.5, 1e6), c(15.02, 1e5),
               c(12, 12), c(12, 1e7), c(11.99, 1e7), c(12.01, 1e7),
               c(150, 150), c(300, 700), c(1e7, 1e7), c(1, 1e7), c(0.3, 8e5))
P <- rbind(P, extra)
rows <- list()
for (i in seq_len(nrow(P))) {
  xs <- edge_x(P[i, 1], P[i, 2])
  if (length(xs)) rows[[i]] <- data.frame(a = P[i, 1], b = P[i, 2], x = xs,
                                         where = names(xs))
}
A <- do.call(rbind, rows)
A$ref <- suppressWarnings(pbeta(A$x, A$a, A$b, log.p = TRUE))
A$ours <- lib(A$x, A$a, A$b)
A <- A[is.finite(A$ref) & A$ref > -700, ]
A$err <- err(A$ours, A$ref)
cat("points", nrow(A), " non-finite:", sum(!is.finite(A$ours)), "\n")
cat("max error by larger shape:\n")
print(signif(tapply(A$err, cut(pmax(A$a, A$b), c(0, 1, 1e2, 1e4, 1e5, 1e6, 1e7 + 1)), max), 3))
cat("max error by location:\n")
print(signif(tapply(A$err, A$where, max), 3))
w <- A[order(-A$err), ][1:8, ]
w$mpfr <- vapply(seq_len(nrow(w)), function(i) as.numeric(mp_I(w$x[i], w$a[i], w$b[i])), 0)
w$ours_vs_mpfr <- err(w$ours, w$mpfr)
w$pbeta_vs_mpfr <- err(w$ref, w$mpfr)
print(w[, c("a", "b", "x", "where", "ref", "err", "ours_vs_mpfr", "pbeta_vs_mpfr")],
      digits = 4, row.names = FALSE)

cat("\n== B. third partials in (x, a, b) finite on the grid ==\n")
fin <- vapply(seq_len(nrow(A)), function(i) {
  p <- c(A$x[i], A$a[i], A$b[i])
  c(all(is.finite(Jn(p))), all(is.finite(Hn(p))), all(is.finite(T3n(p))))
}, logical(3))
cat("points", nrow(A), " gradient", sum(fin[1, ]), " hessian", sum(fin[2, ]),
    " third", sum(fin[3, ]), "\n")
if (!all(fin)) print(head(A[!apply(fin, 2, all), ], 12))

cat("\n== C. derivatives against Rmpfr ==\n")
set_pairs <- list(c(0.3, 7), c(4, 60), c(12, 12), c(16, 240), c(30, 30),
                  c(150, 150), c(300, 700), c(2e3, 8e3), c(3e4, 7e4),
                  c(1e6, 3e6), c(1e7, 1e7), c(1, 1e7))
out <- NULL
for (ab in set_pairs) {
  a <- ab[1]; b <- ab[2]
  xs <- edge_x(a, b)
  xs <- xs[intersect(names(xs), c("m", "t1", "tie2", "lo", "hi", "t1p", "u25m",
                                  "u50m", "v25p", "e05", "e15", "e5m", "e6m"))]
  for (k in seq_along(xs)) {
    x <- xs[[k]]
    if (pbeta(x, a, b, log.p = TRUE) < -600) next
    r <- mp_der(x, a, b)
    p <- c(x, a, b)
    g <- as.vector(Jn(p)); H <- Hn(p); T3 <- T3n(p)
    ours <- c(Tn(p), g[1], H[1, 1], T3[1, 1], g[2], g[3])
    e <- abs(ours - r) / pmax(abs(r), 1e-300)
    out <- rbind(out, data.frame(a = a, b = b, where = names(xs)[k], x = x,
                                 v = e[1], dx = e[2], dxx = e[3], dxxx = e[4],
                                 da = e[5], db = e[6]))
  }
}
op <- options(width = 200)
print(out, digits = 2, row.names = FALSE)
cat("max relative error: value", signif(max(out$v), 3), " d/dx", signif(max(out$dx), 3),
    " d2/dx2", signif(max(out$dxx), 3), " d3/dx3", signif(max(out$dxxx), 3),
    " d/da", signif(max(out$da), 3), " d/db", signif(max(out$db), 3), "\n")

cat("\n== D. continuity across the new edges ==\n")
for (ab in list(c(150, 150), c(300, 700), c(2e3, 8e3), c(3e4, 7e4), c(1e6, 3e6),
                c(1e7, 1e7), c(16, 240), c(60, 60), c(1, 1e7))) {
  a <- ab[1]; b <- ab[2]
  xs <- edge_x(a, b)
  xs <- xs[intersect(names(xs), c("lo", "hi", "u25m", "u50m", "u25p", "u50p",
                                  "v25m", "v25p"))]
  for (k in seq_along(xs)) {
    e <- xs[[k]]
    gap <- b / ((a + b) * (a + b + 1))
    dx <- min(1e-9 * sdb(a, b), 1e-4 * gap, 1e-12 * e)
    p1 <- c(e - dx, a, b); p2 <- c(e + dx, a, b)
    g1 <- as.vector(Jn(p1)); g2 <- as.vector(Jn(p2)); H1 <- Hn(p1)
    jump <- Tn(p2) - Tn(p1)
    pred <- g1[1] * 2 * dx
    dg <- g2 - g1; dgp <- H1[, 1] * 2 * dx
    cat(sprintf("  a %g b %g %-5s: value jump %.3e pred %.3e (diff %.1e) | gradient change vs H dx, worst %.1e of |g|\n",
                a, b, names(xs)[k], jump, pred, jump - pred,
                max(abs(dg - dgp) / pmax(abs(g1), 1e-300))))
  }
}
cat("in the shapes, x at the mean of the pair: across a b / (a + b) = 15, 30 and a = 12\n")
Ta <- MakeTape(function(p) lib(p[3], p[1], p[2]), c(20, 40, 0.3))
Ja <- Ta$jacfun(); Ha <- Ja$jacfun()
for (cfg in list(c(15, 1e3), c(30, 1e3), c(15, 1e7), c(30, 1e7), c(12, 1e7), c(12, 13))) {
  target <- cfg[1]; b <- cfg[2]
  a <- if (target == 12 && b > 1e6) 12 else if (target == 12) 12 else target * b / (b - target)
  x <- a / (a + b)
  if (x >= 0.5) x <- 0.45
  da <- 1e-9 * a
  v1 <- Ta(c(a - da, b, x)); v2 <- Ta(c(a + da, b, x))
  g1 <- Ja(c(a - da, b, x)); g2 <- Ja(c(a + da, b, x)); H1 <- Ha(c(a - da, b, x))
  cat(sprintf("  %s (a %.6g, b %g, x %.4g): value jump %.3e pred %.3e (diff %.1e) | gradient change vs H da %.1e of |g|\n",
              if (target == 12) "a = 12" else paste0("ab/s = ", target), a, b, x,
              v2 - v1, g1[1] * 2 * da, v2 - v1 - g1[1] * 2 * da,
              max(abs(g2 - g1 - H1[, 1] * 2 * da) / pmax(abs(g1), 1e-300))))
}

cat("\n== E. the clamp: x across both ties, 401 points each ==\n")
for (ab in list(c(150, 150), c(300, 700), c(3e4, 7e4), c(1e6, 3e6), c(1e7, 1e7))) {
  a <- ab[1]; b <- ab[2]; s <- a + b; t1 <- a / s; gap <- b / (s * (s + 1))
  xs <- t1 + seq(-1, 2, length.out = 401) * gap
  xs <- xs[xs < 0.5 + 1e-12 & xs > 0]
  ref <- pbeta(xs, a, b, log.p = TRUE)
  v <- lib(xs, a, b)
  # d/dx against dbeta / pbeta in logs, and the tape's d/da, d/db against
  # RTMB::pbeta's own AD at the same points (finite off its tie)
  gref <- exp(dbeta(xs, a, b, log = TRUE) - ref)
  G <- t(vapply(xs, function(x) as.vector(Jn(c(x, a, b))), numeric(3)))
  Tp <- MakeTape(function(p) log(RTMB::pbeta(p[1], p[2], p[3])), c(0.3, 200, 400))
  Jp <- Tp$jacfun()
  Gp <- t(vapply(xs, function(x) as.vector(Jp(c(x, a, b))), numeric(3)))
  okp <- apply(is.finite(Gp), 1, all)
  cat(sprintf("  a %g b %g: value max %.1e | d/dx max rel %.1e | d/da, d/db vs pbeta AD max rel %.1e (%d of %d points where pbeta's AD is finite)\n",
              a, b, max(err(v, ref)), max(abs(G[, 1] - gref) / abs(gref)),
              max(abs(G[okp, 2:3] - Gp[okp, 2:3]) / pmax(abs(Gp[okp, 2:3]), 1e-300)),
              sum(okp), length(xs)))
}
options(op)
