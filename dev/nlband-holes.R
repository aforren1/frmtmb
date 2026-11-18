# Lane nlband: how often a gap between adjacent walls is "wide" under a
# rule that compares each gap with the running median of the gaps around
# it, on designs with no hole (false triggers) and with one (detection).
# Pure R, no frmtmb. Rscript dev/nlband-holes.R > dev/nlband-out/holes.txt
ratios <- function(x, k) {
  w <- sort(unique(x))
  d <- diff(w)
  d <- d[d > 1e-9 * (w[length(w)] - w[1L])]
  m <- length(d)
  win <- 2L * k + 1L
  loc <- if (m <= win) rep(stats::median(d), m) else
    stats::runmed(d, win, endrule = "constant")
  d / loc
}
designs <- list(
  unif = function(n) runif(n),
  norm = function(n) rnorm(n),
  lnorm = function(n) rlnorm(n),
  exp = function(n) rexp(n),
  t3 = function(n) rt(n, 3),
  unif_r3 = function(n) round(runif(n, 0.1, 5), 3),
  norm_r3 = function(n) round(rnorm(n), 3),
  grid = function(n) sample(seq(0, 10, by = 0.5), n, TRUE),
  clust2 = function(n) c(runif(n %/% 2), runif(n - n %/% 2, 3, 4)))
ns <- c(20, 50, 100, 300, 1000, 10000, 100000)
Rs <- c(10, 15, 20, 30, 50)
ks <- c(5L, 10L)
set.seed(20261116)
for (k in ks) for (dn in names(designs)) for (n in ns) {
  reps <- if (n >= 1e5) 20 else if (n >= 1e4) 100 else 400
  cnt <- matrix(0, reps, length(Rs))
  mx <- numeric(reps)
  for (r in seq_len(reps)) {
    q <- ratios(designs[[dn]](n), k)
    cnt[r, ] <- vapply(Rs, function(R) sum(q > R), 0)
    mx[r] <- max(q)
  }
  cat(sprintf("NOHOLE k %d %-7s n %6d reps %3d | any>R %s | mean count %s | max q50 %.1f q99 %.1f\n",
              k, dn, n, reps,
              paste(sprintf("%d:%.3f", Rs, colMeans(cnt > 0)), collapse = " "),
              paste(sprintf("%.2f", colMeans(cnt)), collapse = " "),
              stats::median(mx), stats::quantile(mx, 0.99)))
}
# one hole of half width h in U(0, 4) at 2, n rows before removal
for (k in ks) for (h in c(0.3, 0.1, 0.05, 0.02)) for (n in c(50, 100, 400, 2000)) {
  reps <- 400
  hit <- logical(reps)
  rk <- numeric(reps)
  for (r in seq_len(reps)) {
    x <- runif(n, 0, 4)
    x <- x[abs(x - 2) > h]
    q <- ratios(x, k)
    w <- sort(unique(x))
    d <- diff(w)
    # the gap that holds 2
    j <- which(w[-length(w)] < 2 & w[-1L] > 2)
    rk[r] <- q[j]
  }
  cat(sprintf("HOLE k %d h %.2f n %5d | ratio of the hole q01 %.1f q10 %.1f q50 %.1f | detected at R %s\n",
              k, h, n, stats::quantile(rk, 0.01), stats::quantile(rk, 0.1),
              stats::median(rk),
              paste(sprintf("%d:%.3f", Rs, vapply(Rs, function(R) mean(rk > R), 0)),
                    collapse = " ")))
}
cat("DONE\n")
