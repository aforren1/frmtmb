# Lane famend, fix round 2: gev_parts()'s clamp arithmetic replayed on
# doubles. The old form r = ad_max(z + out, cl) against the new
# zc = r - (ad_pos(r - cl) - ad_pos(r - cl - out)): on
# outside rows (z < cl) how many zc are not exactly cl; on inside rows
# (z > cl) how many zc are not exactly z. (A first form,
# ad_max(r - ad_min(r - cl, out), cl), had the same values but sent an
# inside row's adjoint through one ad_pos() of a positive argument,
# 2 * (0.5 * A), which is not A for an odd subnormal A: inside rows
# with subnormal gradient entries differed from 0.82.0's,
# dev/famend-fix2-xi0diff.R.)
#   Rscript dev/famend-fix2-replay.R
ad_pos <- function(u) 0.5 * (u + abs(u))
ad_min <- function(x, m) x - ad_pos(x - m)
ad_max <- function(x, m) x + ad_pos(m - x)
cl <- -1 + 1e-15
old_zc <- function(z) {
  out <- ad_pos(cl - z)
  ad_max(z + out, cl)
}
new_zc <- function(z) {
  out <- ad_pos(cl - z)
  r <- ad_max(z + out, cl)
  ex <- r - cl
  r - (ad_pos(ex) - ad_pos(ex - out))
}
set.seed(20261114)
eps <- .Machine$double.eps
# outside: a dense sweep of |z| from just below the clamp to 1e20, plus
# the clamp's own neighbors below it
zo <- c(-10^seq(log10(1 + 1e-15), 20, length.out = 2e6),
        -10^stats::runif(1e6, 0, 3),
        cl - (1:200) * eps / 2, -1 - (0:200) * eps)
zo <- zo[zo < cl]
# inside: the clamp's neighbors above it, near 0, and across (-1, 1e20)
zi <- c(cl + (1:400) * eps / 2, stats::runif(1e6, -1, 1),
        10^seq(-300, 20, length.out = 1e6), -10^-seq(0.000001, 300,
                                                       length.out = 1e6))
zi <- zi[zi > cl]
o_old <- old_zc(zo)
o_new <- new_zc(zo)
cat(sprintf("OUTSIDE rows %d: zc != clamp, old %d (above %d, below %d), new %d\n",
            length(zo), sum(o_old != cl), sum(o_old > cl), sum(o_old < cl),
            sum(o_new != cl)))
band <- zo <= cl - 8 & zo > cl - 16
cat(sprintf("  of them cl - z in [8, 16): %d rows, old above the clamp %d, new %d\n",
            sum(band), sum(o_old[band] > cl), sum(o_new[band] != cl)))
big <- abs(zo) > 2^52
cat(sprintf("  |z| > 2^52: %d rows, old zc range [%g, %g], new zc != clamp %d\n",
            sum(big), min(o_old[big]), max(o_old[big]), sum(o_new[big] != cl)))
zf <- -10^seq(15, 40, length.out = 1e6)
cat(sprintf("FAR OUTSIDE |z| 1e15 to 1e40, %d rows: zc != clamp, old %d (zc 0: %d), new %d\n",
            length(zf), sum(old_zc(zf) != cl), sum(old_zc(zf) == 0),
            sum(new_zc(zf) != cl)))
i_old <- old_zc(zi)
i_new <- new_zc(zi)
cat(sprintf("INSIDE rows %d: zc != z, old %d, new %d; new != old %d\n",
            length(zi), sum(i_old != zi), sum(i_new != zi),
            sum(i_new != i_old)))
