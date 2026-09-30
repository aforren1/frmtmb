# Reviewer, punch round 1, item 2: the new log_ibeta_half().
#  A. value against stats::pbeta(log.p = TRUE) on a grid of shapes 1e-3
#     to 1e7 and x at m, near m, at every blend edge; the worker's
#     bounds; Rmpfr at the worst points
#  B. non-finite derivatives: the whole function on grid A, and
#     RTMB::pbeta() alone over the region the package reads it in
#     (shapes >= 150, x within 6.5 sd of the mean), with exact values
#  C. continuity of value and gradient across each blend edge and the
#     pbeta region's shape boundary
#  D. cost per row of a gradient sweep, against round 1 and pbeta
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({ library(RTMB); library(Rmpfr) })
lib <- frmtmb:::log_ibeta_half
err <- function(v, ref) abs(v - ref) / pmax(1, abs(ref))
sdb <- function(a, b) sqrt(a * b / ((a + b)^2 * (a + b + 1)))
mp_log_ibeta <- function(x, a, b, prec = 300, maxit = 400000) {
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
      if (abs(del - 1) < mpfr(1e-50, prec)) break
    }
    a * log(x) + b * log1p(-x) + lgamma(a + b) - lgamma(a) - lgamma(b) -
      log(a) + log(h)
  }
  if (x < (a + 1) / (a + b + 2)) as.numeric(cf(x, a, b)) else
    as.numeric(log1p(-exp(cf(1 - x, b, a))))
}

cat("== A. accuracy grid ==\n")
sh <- 10^seq(-3, 7, by = 0.5)
extra <- rbind(c(300, 300), c(900, 900), c(151, 1e7), c(450.5, 1e7),
               c(0.268941421 * 1000, 0.731058579 * 1000), c(160, 2400),
               c(3e4, 7e4), c(0.3, 7), c(40, 960), c(60, 8e5), c(0.01, 5e6))
P <- rbind(as.matrix(expand.grid(a = sh, b = sh)), extra)
offs <- c(-6.5, -6, -5.5, -5, -2, -1, -0.5, 0, 0.5, 1, 1.5, 2, 5, 5.5, 6, 6.5)
rows <- list()
for (i in seq_len(nrow(P))) {
  a <- P[i, 1]; b <- P[i, 2]
  m <- (a + 1) / (a + b + 2); s <- sdb(a, b)
  xs <- c(m, m * (1 + c(-1e-3, 1e-3)), m + offs * s,
          m + c(0.5, 1.5, -5, -6, 5, 6) * s * (1 + 1e-12))
  lab <- c("m", "m-1e-3", "m+1e-3", paste0("m", sprintf("%+g", offs), "sd"),
           paste0("edge", c("+0.5", "+1.5", "-5", "-6", "+5", "+6"), "sd+"))
  ok <- xs > 0 & xs < 0.5
  if (any(ok)) rows[[i]] <- data.frame(a = a, b = b, x = xs[ok], where = lab[ok])
}
A <- do.call(rbind, rows)
A$ref <- suppressWarnings(pbeta(A$x, A$a, A$b, log.p = TRUE))
A$ours <- lib(A$x, A$a, A$b)
A <- A[is.finite(A$ref) & A$ref > -700, ]
A$err <- err(A$ours, A$ref)
cat("points", nrow(A), "; non-finite ours:", sum(!is.finite(A$ours)), "\n")
q <- function(v) sprintf("max %.2e q99 %.2e", max(v), quantile(v, 0.99))
cat("all:                       ", q(A$err), "\n")
small <- A$a < 1e4 & A$b < 1e4
cat("both shapes < 1e4:         ", q(A$err[small]), " (worker: < 5e-13)\n")
cat("x == m:                    ", q(A$err[A$where == "m"]), " (worker: 4.97e-13)\n")
cat("x == m, both shapes < 1e4: ", q(A$err[A$where == "m" & small]), "\n")
cat("by location (max):\n")
print(signif(tapply(A$err, A$where, max), 3))
cat("by larger shape (max):\n")
print(signif(tapply(A$err, cut(pmax(A$a, A$b), c(0, 1, 1e2, 1e4, 1e5, 1e6, 1e7 + 1)), max), 3))
cat("worst 8:\n")
w <- A[order(-A$err), ][1:8, ]
w$mpfr <- vapply(seq_len(nrow(w)), function(i) mp_log_ibeta(w$x[i], w$a[i], w$b[i]), 0)
w$err_vs_mpfr <- err(w$ours, w$mpfr)
w$pbeta_vs_mpfr <- err(w$ref, w$mpfr)
print(w, digits = 4, row.names = FALSE)

cat("\n== B1. derivatives of log_ibeta_half in (x, log a, log b), finite on grid A ==\n")
F <- MakeTape(function(p) lib(p[1], exp(p[2]), exp(p[3])), c(0.1, 0, 0))
J <- F$jacfun(); H <- J$jacfun(); T3 <- H$jacfun()
fin <- vapply(seq_len(nrow(A)), function(i) {
  p <- c(A$x[i], log(A$a[i]), log(A$b[i]))
  c(all(is.finite(J(p))), all(is.finite(H(p))), all(is.finite(T3(p))))
}, logical(3))
cat("points", nrow(A), " finite gradient", sum(fin[1, ]), " hessian", sum(fin[2, ]),
    " third", sum(fin[3, ]), "\n")
if (any(!fin)) print(head(A[!apply(fin, 2, all), ], 10))

cat("\n== B2. RTMB::pbeta third derivatives over the region the package reads ==\n")
G3 <- MakeTape(function(p) log(RTMB::pbeta(p[1], p[2], p[3])), c(0.3, 200, 400))
T3p <- G3$jacfun()$jacfun()$jacfun()
Jp <- G3$jacfun()
set.seed(20260930)
nbad <- 0; nn <- 0; badpts <- NULL
for (i in 1:20000) {
  a <- exp(runif(1, log(150), log(1e7))); b <- exp(runif(1, log(150), log(1e7)))
  mm <- a / (a + b); s <- sdb(a, b)
  x <- mm + runif(1, -6.5, 6.5) * s
  if (x <= 0 || x >= 1) next
  nn <- nn + 1
  t3 <- T3p(c(x, a, b))
  if (!all(is.finite(t3)) || !all(is.finite(Jp(c(x, a, b))))) {
    nbad <- nbad + 1
    if (is.null(badpts) || nrow(badpts) < 10) badpts <- rbind(badpts, c(x = x, a = a, b = b))
  }
}
cat("random: points", nn, " non-finite gradient or third:", nbad, "\n")
if (!is.null(badpts)) print(badpts)
cat("exact and edge values:\n")
ex <- list(c(0.268941421, 150, 407.74), c(0.268941421, 268.941421, 731.058579),
           c(0.268941421, 2689.41421, 7310.58579), c(0.268941421, 26894.1421, 73105.8579),
           c(0.5, 150, 150), c(0.5, 1e7, 1e7), c(150 / (150 + 1e7), 150, 1e7),
           c(0.3, 300, 700), c(0.25, 250, 750), c(0.1, 150, 1350))
for (p in ex) {
  a <- p[2]; b <- p[3]; mm <- a / (a + b); s <- sdb(a, b)
  for (k in c(0, -6.5, 6.5, -3, 3)) {
    x <- if (k == 0) p[1] else mm + k * s
    if (x <= 0 || x >= 1) next
    t3 <- T3p(c(x, a, b))
    if (!all(is.finite(t3)))
      cat(sprintf("  NON-FINITE at x %.10g a %g b %g (k %g)\n", x, a, b, k))
  }
}
cat("  (no line above means all finite)\n")
cat("RTMB::pbeta third derivatives at SMALL shapes (outside the region), for contrast:\n")
set.seed(7)
nb <- 0
for (i in 1:2000) {
  a <- exp(runif(1, -4, 5)); b <- exp(runif(1, -4, 5)); x <- runif(1, 0.001, 0.999)
  if (!all(is.finite(T3p(c(x, a, b))))) nb <- nb + 1
}
cat("  non-finite", nb, "of 2000\n")
cat("where the non-finite small-shape points sit, by min shape:\n")
set.seed(7)
mins <- c()
for (i in 1:2000) {
  a <- exp(runif(1, -4, 5)); b <- exp(runif(1, -4, 5)); x <- runif(1, 0.001, 0.999)
  if (!all(is.finite(T3p(c(x, a, b))))) mins <- c(mins, min(a, b))
}
print(summary(mins))

cat("\n== C. continuity across the edges (x) ==\n")
Jx <- F$jacfun()
for (ab in list(c(0.3, 7), c(4, 60), c(40, 960), c(300, 700), c(3e4, 7e4),
                c(160, 2400), c(2e5, 8e5), c(1e6, 3e6))) {
  a <- ab[1]; b <- ab[2]; m <- (a + 1) / (a + b + 2); s <- sdb(a, b)
  mm <- a / (a + b); s2 <- s
  edges <- c(m + 0.5 * s, m + 1.5 * s, m, mm - 6 * s, mm - 5 * s, mm + 5 * s,
             mm + 6 * s)
  for (e in edges) {
    if (e <= 0 || e >= 0.5) next
    dx <- 1e-9 * s
    p1 <- c(e - dx, log(a), log(b)); p2 <- c(e + dx, log(a), log(b))
    g1 <- Jx(p1); g2 <- Jx(p2)
    jump <- F(p2) - F(p1)
    pred <- g1[1] * 2 * dx
    cat(sprintf("  a %g b %g edge %.6g: value jump %.3e pred %.3e (diff %.1e) | grad rel jump %.1e\n",
                a, b, e, jump, pred, jump - pred,
                max(abs(g2 - g1) / pmax(abs(g1), 1e-300))))
  }
}
cat("continuity in a across a b / (a + b) = 150 and 450, x at the mean:\n")
Fa <- MakeTape(function(p) lib(p[3], p[1], p[2]), c(200, 300, 0.4))
Ja <- Fa$jacfun()
for (bb in c(1e3, 1e5)) {
  for (target in c(150, 450)) {
    a <- target * bb / (bb - target)
    x <- a / (a + bb)
    da <- 1e-9 * a
    v1 <- Fa(c(a - da, bb, x)); v2 <- Fa(c(a + da, bb, x))
    g1 <- Ja(c(a - da, bb, x)); g2 <- Ja(c(a + da, bb, x))
    cat(sprintf("  b %g, ab/s = %g (a %.6g): value jump %.2e pred %.2e | grad rel jump %.1e\n",
                bb, target, a, v2 - v1, g1[1] * 2 * da,
                max(abs(g2 - g1) / pmax(abs(g1), 1e-300))))
  }
}

cat("\n== D. cost per row of a gradient sweep ==\n")
# round 1's function, reconstructed from the reviewed diff
sp <- paste0("C:/Users/adf44/AppData/Local/Temp/1/claude/",
             "c--Users-adf44-source-r-frmtmb/",
             "66ed580c-211a-4baa-94bd-45a52ec3082c/scratchpad")
ex1 <- parse(file.path(sp, "r1", "families.R"))
r1 <- new.env(parent = asNamespace("frmtmb"))
for (e in ex1) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c("log_ibeta_half", "log_ibeta_cf")) eval(e, r1)
}
set.seed(5)
n <- 500
xr <- runif(n, 0.02, 0.3); ar <- exp(runif(n, 0, 6)); br <- exp(runif(n, 0, 6))
mk <- function(f) {
  Tp <- MakeTape(function(p) sum(f(xr * 0 + p[1] * xr, p[2] * ar, p[3] * br)), c(1, 1, 1))
  Tp$jacfun()
}
Jn <- mk(lib); J1 <- mk(r1$log_ibeta_half)
Jb <- mk(function(x, a, b) log(RTMB::pbeta(x, a, b)))
tm <- function(Jf) {
  best <- Inf
  for (r in 1:7) best <- min(best, system.time(for (k in 1:20) Jf(c(1, 1, 1)))[["elapsed"]] / 20)
  best
}
tt <- c(new = tm(Jn), round1 = tm(J1), pbeta = tm(Jb), new_again = tm(Jn))
print(signif(tt * 1e3, 3))
cat(sprintf("per row: new %.1f us, round 1 %.1f us, pbeta %.1f us; new / round 1 %.2f, new / pbeta %.1f\n",
            tt["new"] / n * 1e6, tt["round1"] / n * 1e6, tt["pbeta"] / n * 1e6,
            tt["new"] / tt["round1"], tt["new"] / tt["pbeta"]))
