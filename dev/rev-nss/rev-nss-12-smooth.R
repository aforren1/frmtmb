# REVIEW of lane nss round 1, attack 2: how smooth is the new shape,
# and where.
#
# The lane measured smoothness through `frm_ode()`, which puts an ODE
# solve between the shape and the number. This measures the shape
# itself. The correction is
#
#     y = yc + d2 * F(r),  F(r) = (r / (um * h)) * step(uc) * gate
#
# so with ya = 0, yb = 1, yc = 1 + r the damping term is negligible,
# the measured ratio IS r, d2 IS r, and F(r) = (y - yc) / r to machine
# precision. Nothing but `ode_ss_extrapolate()` is called.
#
# There are FOUR junctions, not two. The lane's section 4a names the
# stand-down point, r = 1 and "the low gate". `gate` is
# `step(min(max(r / rlow, 0), 1))`, and `step(1) = 1` while
# `step'(1) = -1`, so the gate has a corner where it SATURATES, at
# r = rlow, as well as where it opens, at r = 0. This script asks all
# four.
#
# Script path: dev/rev-nss/rev-nss-12-smooth.R
# No seed: nothing here is random.
source("C:/Users/adf44/source/r/frmtmb-wt-nss/dev/rev-nss/prelude.R")
suppressMessages(library(frmtmb.ode))
ext <- frmtmb.ode:::ode_ss_extrapolate
RCAP <- frmtmb.ode:::ode_ss_rcap
RLOW <- frmtmb.ode:::ode_ss_rlow
cat("\nrcap", RCAP, " rlow", RLOW, " noise",
    frmtmb.ode:::ode_ss_noise, "\n")

# F(r), read off the shipped function and nothing else. The LAST
# difference is held at 1 and the one before it at 1 / r, so the ratio
# the function measures IS r and the correction it returns IS F(r).
# r = 0 needs no probe: the gate is shut at and below it.
mk <- function(r) c(0, 1 / r, 1 / r + 1)
FF <- function(r) {
  if (r == 0) return(0)
  v <- mk(r)
  as.numeric(ext(v[1L], v[2L], v[3L], 1e-14, 1e-14)$y - v[3L])
}
# check the instrument: the ratio the function measures must be r
chk <- vapply(c(-0.3, -1e-4, 0.02, 0.5, 0.9875, 0.995, 1.2),
              function(r) {
                v <- mk(r)
                as.numeric(ext(v[1L], v[2L], v[3L], 1e-14, 1e-14)$r) - r
              }, 0)
cat("max |measured r - requested r| over seven probes:",
    format(max(abs(chk))), "\n")

## ---- A. is F exactly the geometric tail below the cap? -------------
cat("\n== A. F(r) against r / (1 - r) ==\n")
cat(sprintf("%10s %14s %14s %12s\n", "r", "F(r)", "r/(1-r)", "rel"))
for (r in c(0.05, 0.1, 0.3, 0.6, 0.9, 0.95, 0.98, RCAP, 0.99, 0.995,
            0.999, 1, 1.05)) {
  g <- r / (1 - r)
  cat(sprintf("%10.5f %14.6f %14.6f %12.2e\n", r, FF(r), g,
              abs(FF(r) - g) / max(1e-12, abs(g))))
}
cat("\nlargest factor F over r in [0, 1.5]:\n")
rs <- seq(0, 1.5, by = 1e-6)
fs <- vapply(rs, FF, 0)
j <- which.max(fs)
cat(sprintf("  max F = %.4f at r = %.5f\n", fs[[j]], rs[[j]]))
cat(sprintf("  a hard cap at 0.99 would give %.4f\n", 0.99 / 0.01))
cat(sprintf("  max |F(r) - r/(1-r)| over r <= 0.98: %.3e\n",
            max(abs(fs[rs <= 0.98] - rs[rs <= 0.98] /
                      (1 - rs[rs <= 0.98])))))
g <- rs / (1 - rs)
k <- rs > 0 & rs <= RLOW
cat(sprintf("  gate overshoot: max F / (r/(1-r)) on (0, rlow] = %.4f",
            max(fs[k] / g[k])))
cat(sprintf(" at r = %.5f\n", rs[k][which.max(fs[k] / g[k])]))

## ---- B. one-sided derivatives at every junction --------------------
# The gap between the two one-sided k-th difference quotients. If it
# falls like eps the function is C(k); if it settles on a constant, the
# k-th derivative jumps by that constant.
sided <- function(f, x0, eps, k, side) {
  s <- if (side > 0) 1 else -1
  p <- vapply(0:k, function(i) f(x0 + s * i * eps), 0)
  sum(choose(k, 0:k) * (-1)^(0:k) * rev(p)) / (s * eps)^k
}
junc <- list(`r = 0 (gate opens)` = 0,
             `r = rlow (gate saturates)` = RLOW,
             `r = rcap (stand-down begins)` = RCAP,
             `r = 1 (stand-down complete)` = 1)
for (nm in names(junc)) {
  x0 <- junc[[nm]]
  cat(sprintf("\n== B. %s, x0 = %.6f ==\n", nm, x0))
  cat(sprintf("%3s %10s %16s %16s %12s %10s\n", "k", "eps", "left",
              "right", "gap", "gap/eps"))
  for (k in 1:3) {
    prev <- NA_real_
    for (eps in 10^-(2:5)) {
      l <- sided(FF, x0, eps, k, -1)
      r <- sided(FF, x0, eps, k, +1)
      g <- abs(l - r)
      cat(sprintf("%3d %10.0e %16.8f %16.8f %12.3e %10.2e\n", k, eps,
                  l, r, g, g / eps))
    }
    cat("\n")
  }
}

## ---- C. the exactness boundary the cap move cost -------------------
cat("== C. what moving rcap from 0.99 to 0.9875 costs ==\n")
for (ii in c(8, 12, 24)) {
  th <- function(rr) log(2) / (-log(rr) / ii)
  cat(sprintf("  ii %3d h: exact up to a terminal half-life of",  ii))
  cat(sprintf(" %8.1f h at rcap %.4f, %8.1f h at 0.99\n",
              th(RCAP), RCAP, th(0.99)))
}
cat("\ndone\n")

## ---- D. why the gate has a corner, and the one-line repair --------
# `step()` is the degree-7 smootherstep DIVIDED BY x:
#   S(x)  = 35x^4 - 84x^5 + 70x^6 - 20x^7,  S(1)=1, S'(1)=S''(1)=S'''(1)=0
#   step  = S(x)/x = 35x^3 - 84x^4 + 70x^5 - 20x^6,  step'(1) = -1
# That -1 is not a mistake in the STAND-DOWN factor: the branch below
# the cap is r/(1-r), and -1 is exactly the slope that makes the product
# join it smoothly. It IS a mistake in the GATE, which has to saturate
# flat at 1, and a slope of -1 there is a corner.
cat("\n== D. step(x) = S(x)/x, and what that does at each end ==\n")
step <- function(x) x^3 * (35 - x * (84 - x * (70 - 20 * x)))
smoo <- function(x) x^4 * (35 - x * (84 - x * (70 - 20 * x)))
d1 <- function(f, x, e = 1e-6) (f(x + e) - f(x - e)) / (2 * e)
cat(sprintf("  step(1) = %.6f   step'(1) = %+.6f  (a corner in a gate)\n",
            step(1), d1(step, 1)))
cat(sprintf("  S(1)    = %.6f   S'(1)    = %+.6f  (flat, as a gate needs)\n",
            smoo(1), d1(smoo, 1)))
cat(sprintf("  max step(x) on [0,1] = %.4f, so the gate OVERSHOOTS\n",
            max(vapply(seq(0, 1, by = 1e-5), step, 0))))
cat(sprintf("  max S(x)    on [0,1] = %.4f\n",
            max(vapply(seq(0, 1, by = 1e-5), smoo, 0))))
cat("\n  the whole 9.3e-03 disagreement the lane reports for r <= 0.98:\n")
gg <- rs / (1 - rs)
kk <- rs > 0 & rs <= 0.98
w <- which.max(abs(fs[kk] - gg[kk]))
cat(sprintf("    largest at r = %.5f, which is below rlow = %.3f\n",
            rs[kk][w], RLOW))
cat("    below the cap the two stand-down shapes are IDENTICAL, so\n")
cat("    all of it is the new low gate, not the smoothing.\n")
cat("
done
")
