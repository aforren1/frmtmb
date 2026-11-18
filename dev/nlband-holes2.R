# Lane nlband: the gap rule, second form. A gap between adjacent walls is
# a hole when it has `side` walls or more on each side and is wider than
# T(m) times the running median of the 2k + 1 gaps around it (m gaps in
# all). T(m) = log2(a * m): under locally uniform rows a gap is wider than
# t times the local mean with probability exp(-t), so m gaps hold one
# wider than log(a m) local means with probability about 1 / a, and the
# median is ln 2 means. Pure R. Rscript dev/nlband-holes2.R
ratios <- function(x, k) {
  w <- sort(unique(x))
  d <- diff(w)
  d <- d[d > 1e-9 * (w[length(w)] - w[1L])]
  m <- length(d)
  win <- 2L * k + 1L
  loc <- if (m <= win) rep(stats::median(d), m) else
    stats::runmed(d, win, endrule = "constant")
  list(q = d / loc, m = m, w = w, d = d)
}
holes <- function(x, k, side, a) {
  r <- ratios(x, k)
  ok <- seq_len(r$m) > side & seq_len(r$m) <= r$m - side
  which(ok & r$q > log2(a * r$m))
}
designs <- list(
  unif = function(n) runif(n),
  norm = function(n) rnorm(n),
  lnorm = function(n) rlnorm(n),
  exp = function(n) rexp(n),
  t3 = function(n) rt(n, 3),
  unif_r3 = function(n) round(runif(n, 0.1, 5), 3),
  norm_r3 = function(n) round(rnorm(n), 3),
  clust2 = function(n) c(runif(n %/% 2), runif(n - n %/% 2, 3, 4)))
ns <- c(20, 50, 100, 300, 1000, 10000, 100000)
vars <- list(c(k = 10, side = 5, a = 1000), c(k = 20, side = 5, a = 1000),
             c(k = 20, side = 10, a = 1000), c(k = 20, side = 5, a = 1e4))
set.seed(20261116)
for (v in vars) for (dn in names(designs)) for (n in ns) {
  reps <- if (n >= 1e5) 20 else if (n >= 1e4) 100 else 400
  cnt <- vapply(seq_len(reps), function(r) {
    length(holes(designs[[dn]](n), v[["k"]], v[["side"]], v[["a"]]))
  }, 0)
  cat(sprintf("NOHOLE k %d side %d a %g %-7s n %6d reps %3d | any %.3f mean %.2f max %d\n",
              v[["k"]], v[["side"]], v[["a"]], dn, n, reps, mean(cnt > 0),
              mean(cnt), max(cnt)))
}
for (v in vars) for (h in c(0.3, 0.1, 0.05)) for (n in c(50, 100, 400, 2000)) {
  reps <- 400
  hit <- vapply(seq_len(reps), function(r) {
    x <- runif(n, 0, 4)
    x <- x[abs(x - 2) > h]
    w <- sort(unique(x))
    j <- which(w[-length(w)] < 2 & w[-1L] > 2)
    j %in% holes(x, v[["k"]], v[["side"]], v[["a"]])
  }, NA)
  cat(sprintf("HOLE k %d side %d a %g h %.2f n %5d | detected %.3f\n",
              v[["k"]], v[["side"]], v[["a"]], h, n, mean(hit)))
}
cat("DONE\n")
