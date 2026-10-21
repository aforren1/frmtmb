# Cost of lgamma_shift_diff() in the dirichlet_multinomial density:
# objective and gradient evaluations of one fit's tape, the final
# density against the plain lgamma difference put back, interleaved in
# one process, minimum of 7 rounds of blocks past 1.2 s, with a control
# arm (the final density against itself) that must report 1.0.
#   Rscript dev/nb2dm-tapecost.R > dev/nb2dm-out/tapecost.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-nb2dm-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(2026)
n <- 400
d <- data.frame(x = rnorm(n), g = gl(20, 20))
P <- exp(cbind(0, 0.4 + 0.6 * d$x, -0.3 - 0.5 * d$x))
P <- P / rowSums(P)
d$tr <- sample(5:40, n, TRUE)
Y <- t(vapply(seq_len(n), function(i) {
  p <- rgamma(3, P[i, ] * 6)
  as.vector(rmultinom(1, d$tr[i], p / sum(p)))
}, numeric(3)))
colnames(Y) <- c("a", "b", "c")
d$Y <- Y
build <- function() {
  frm(Y | trials(tr) ~ x + (1 | g), family = dirichlet_multinomial(),
      data = d)$obj
}
real <- frmtmb:::lgamma_shift_diff
o_shift <- build()
o_ctl <- build()
utils::assignInNamespace("lgamma_shift_diff",
                         function(a, s) lgamma(a + s) - lgamma(a), "frmtmb")
o_plain <- build()
utils::assignInNamespace("lgamma_shift_diff", real, "frmtmb")
p0 <- o_shift$par
blk <- function(o, k) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(k)) {
    o$fn(p0)
    o$gr(p0)
  }
  proc.time()[["elapsed"]] - t0
}
k <- 1L
while (blk(o_shift, k) < 1.2) k <- k * 2L
res <- replicate(7, c(shift = blk(o_shift, k), plain = blk(o_plain, k),
                      ctl = blk(o_ctl, k)))
m <- apply(res, 1, min)
cat("evaluations per block:", k, "(fn + gr, with the Laplace inner",
    "problem)\n")
print(round(m, 3))
cat("ratio shift / plain:", round(m[["shift"]] / m[["plain"]], 3),
    " control shift / shift:", round(m[["ctl"]] / m[["shift"]], 3), "\n")
