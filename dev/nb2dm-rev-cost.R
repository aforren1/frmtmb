# Reviewer: cost of lgamma_shift_diff() in the dirichlet_multinomial
# objective, a second construction beside the lane's: a fixed-effects
# fit (no Laplace inner problem), K = 4, n = 2000, phi ~ z. Arms
# interleaved in one process, blocks grown past 1.2 s, minimum of 9
# rounds, and a control arm (a second shift tape) that must read 1.0.
source("C:/Users/adf44/source/r/frmtmb-wt-nb2dm/dev/nb2dm-rev-hdr.R")
set.seed(31)
n <- 2000
K <- 4
d <- data.frame(x = rnorm(n), z = rnorm(n))
E <- cbind(0, 0.3 + 0.4 * d$x, -0.2 - 0.3 * d$x, 0.1 * d$x)
P <- exp(E) / rowSums(exp(E))
d$N <- sample(5:60, n, TRUE)
d$Y <- t(vapply(seq_len(n), function(i) {
  p <- rgamma(K, P[i, ] * 8)
  as.numeric(rmultinom(1, d$N[i], p / sum(p)))
}, numeric(K)))
colnames(d$Y) <- letters[1:K]
build <- function() {
  frm(bf(Y | trials(N) ~ x, phi ~ z), family = dirichlet_multinomial(),
      data = d)$obj
}
real <- frmtmb:::lgamma_shift_diff
o_shift <- build()
o_ctl <- build()
utils::assignInNamespace("lgamma_shift_diff",
                         function(a, s) lgamma(a + s) - lgamma(a), "frmtmb")
o_plain <- build()
utils::assignInNamespace("lgamma_shift_diff", real, "frmtmb")
o_mn <- frm(Y | trials(N) ~ x, family = multinomial(), data = d)$obj
p0 <- o_shift$par
cat("objective at p0: shift", format(o_shift$fn(p0), digits = 15),
    " plain", format(o_plain$fn(p0), digits = 15), "\n")
blk <- function(o, k, p) {
  t0 <- proc.time()[["elapsed"]]
  for (i in seq_len(k)) {
    o$fn(p)
    o$gr(p)
  }
  proc.time()[["elapsed"]] - t0
}
k <- 1L
while (blk(o_shift, k, p0) < 1.2) k <- k * 2L
pm <- o_mn$par
res <- replicate(9, c(shift = blk(o_shift, k, p0), plain = blk(o_plain, k, p0),
                      ctl = blk(o_ctl, k, p0), mn = blk(o_mn, k, pm)))
m <- apply(res, 1, min)
cat("evaluations per block:", k, "\n")
print(round(m, 3))
cat("ratio shift / plain:", round(m[["shift"]] / m[["plain"]], 3),
    " control / shift:", round(m[["ctl"]] / m[["shift"]], 3),
    " plain / multinomial:", round(m[["plain"]] / m[["mn"]], 3), "\n")
cat("full fit time (s): shift",
    system.time(frm(bf(Y | trials(N) ~ x, phi ~ z),
                    family = dirichlet_multinomial(), data = d))[["elapsed"]],
    "\n")
