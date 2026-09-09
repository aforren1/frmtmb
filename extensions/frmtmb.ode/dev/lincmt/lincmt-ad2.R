# The GRADIENT identity across the schedule space and the model shapes,
# not just at one shape as lincmt-ad.R does.
#
# For each (ncmt, depot, schedule) the sum of the predictor is taped as
# a function of the log rate constants, and the tape's Jacobian is
# compared with frm_ode()'s adjoint Jacobian at atol = rtol = 1e-12.
# Each shape is also evaluated at a DEGENERATE point where two rate
# constants are exactly equal, which is where the textbook closed form
# has no gradient at all.
source("lincmt-src.R")

mk_dyn <- function(ncmt, depot) {
  if (!depot && ncmt == 1L)
    return(function(t, y, p) list(c(-p[1] * y[1])))
  if (!depot && ncmt == 2L)
    return(function(t, y, p) list(c(-(p[1] + p[2]) * y[1] + p[3] * y[2],
                                    p[2] * y[1] - p[3] * y[2])))
  if (!depot && ncmt == 3L)
    return(function(t, y, p) list(c(
      -(p[1] + p[2] + p[4]) * y[1] + p[3] * y[2] + p[5] * y[3],
      p[2] * y[1] - p[3] * y[2], p[4] * y[1] - p[5] * y[3])))
  if (ncmt == 1L)
    return(function(t, y, p) list(c(-p[6] * y[1],
                                    p[6] * y[1] - p[1] * y[2])))
  if (ncmt == 2L)
    return(function(t, y, p) list(c(
      -p[6] * y[1],
      p[6] * y[1] - (p[1] + p[2]) * y[2] + p[3] * y[3],
      p[2] * y[2] - p[3] * y[3])))
  function(t, y, p) list(c(
    -p[6] * y[1],
    p[6] * y[1] - (p[1] + p[2] + p[4]) * y[2] + p[3] * y[3] + p[5] * y[4],
    p[2] * y[2] - p[3] * y[3], p[4] * y[2] - p[5] * y[4]))
}
st_names <- function(ncmt, depot)
  c(if (depot) "depot", "central",
    if (ncmt >= 2L) "peripheral1", if (ncmt == 3L) "peripheral2")

# th holds log(ke, k12, k21, k13, k31, ka); the unused ones are carried
# so that every shape has the same tape domain
nm6 <- c("ke", "k12", "k21", "k13", "k31", "ka")
pk_list <- function(ncmt, depot, p) {
  z <- list(ke = p[[1]])
  if (ncmt >= 2L) { z$k12 <- p[[2]]; z$k21 <- p[[3]] }
  if (ncmt == 3L) { z$k13 <- p[[4]]; z$k31 <- p[[5]] }
  if (depot) z$ka <- p[[6]]
  z
}

tt <- c(0.5, 1, 2, 4, 8, 12, 18, 24, 36, 48)

f_lin <- function(ncmt, depot, ev) function(th) {
  p <- exp(th)
  sum(frm_lincmt(parms = pk_list(ncmt, depot, p), times = tt,
                 ncmt = ncmt, depot = depot, output = "central",
                 events = ev, n_ss = 20L))
}
f_ode <- function(ncmt, depot, ev) function(th) {
  p <- exp(th)
  n <- ncmt + as.integer(depot)
  sum(frm_ode(mk_dyn(ncmt, depot), init = rep(list(0), n), times = tt,
              parms = as.list(p), states = st_names(ncmt, depot),
              output = "central", events = ev, n_ss = 20L,
              atol = 1e-12, rtol = 1e-12))
}

base <- log(c(ke = 0.2, k12 = 0.4, k21 = 0.1, k13 = 0.05, k31 = 0.01,
              ka = 1.1))
# Every rate constant equal. This coalesces `ka` with a disposition
# eigenvalue, which is lincmt_diff()'s singularity; it does NOT
# coalesce two eigenvalues, because at ke = k12 = k21 = k13 = k31 = r
# the characteristic polynomial factors as (L - r)(L^2 - 4 r L + r^2)
# and the roots r, (2 - sqrt 3) r and (2 + sqrt 3) r are distinct.
degen <- log(rep(0.3, 6))
# The disposition collision, which the arm above does not contain:
# k13 = 0 decouples the third compartment and k31 on the slow root of
# the reduced quadratic makes that root double, so lincmt_disp()'s
# acos() sits where its derivative is infinite. k13 = 0 is written as
# log(1e-300) rather than -Inf so that the tape has a finite input.
# The slow root is taken from the quadratic formula rather than from
# the product of the roots. The two agree to fifteen digits and differ
# in the last two bits, and WHICH of them is used decides whether
# lincmt_disp()'s acos() argument rounds to exactly 1. It does for the
# product form, and there the gradient is NaN; dev/lincmt/lincmt-f6.R
# measures that knife edge.
bq <- 0.2 + 0.4 + 0.1
collide <- log(c(ke = 0.2, k12 = 0.4, k21 = 0.1, k13 = 1e-300,
                 k31 = (bq - sqrt(bq * bq - 4 * 0.2 * 0.1)) / 2,
                 ka = 1.1))

cat(sprintf("%-10s %-26s %-10s %12s %12s\n", "shape", "schedule",
            "point", "rel grad", "rel value"))
worst <- 0
for (ncmt in 1:3) for (dp in c(FALSE, TRUE)) {
  cin <- if (dp) "depot" else "central"
  scheds <- list(
    "single dose" = data.frame(time = 0, state = cin, value = 100),
    "addl/ii" = data.frame(time = 0, state = cin, value = 100, ii = 8,
                           addl = 4L),
    "ss" = data.frame(time = 0, state = cin, value = 100, ii = 8,
                      ss = TRUE),
    "ss + addl" = data.frame(time = c(0, 8), state = cin, value = 100,
                             ii = c(8, 8), addl = c(0L, 4L),
                             ss = c(TRUE, FALSE)),
    "iv infusion, addl" = data.frame(time = 0, state = "central",
                                     value = 100, duration = 2,
                                     ii = 8, addl = 4L),
    "iv infusion, ss" = data.frame(time = 0, state = "central",
                                   value = 100, duration = 2, ii = 8,
                                   ss = TRUE))
  for (snm in names(scheds)) {
    ev <- scheds[[snm]]
    fl <- f_lin(ncmt, dp, ev)
    fo <- f_ode(ncmt, dp, ev)
    pts <- if (ncmt == 3L)
      c("ordinary", "ka on an eigenvalue", "double eigenvalue") else
        c("ordinary", "ka on an eigenvalue")
    for (pt in pts) {
      th <- switch(pt, ordinary = base,
                   "ka on an eigenvalue" = degen, collide)
      gl <- as.numeric(RTMB::MakeTape(fl, th)$jacobian(th))
      go <- suppressWarnings(
        as.numeric(RTMB::MakeTape(fo, th)$jacobian(th)))
      vl <- fl(th); vo <- suppressWarnings(fo(th))
      # only the rate constants this shape uses are compared: the
      # others are dead inputs and both paths report zero for them
      use <- c(TRUE, ncmt >= 2, ncmt >= 2, ncmt == 3, ncmt == 3, dp)
      rg <- max(abs(gl[use] - go[use])) / max(abs(go[use]))
      rv <- abs(vl - vo) / abs(vo)
      worst <- max(worst, rg)
      cat(sprintf("%-10s %-26s %-16s %12.2e %12.2e\n",
                  paste0(ncmt, "cmt", if (dp) "+dep" else ""), snm, pt,
                  rg, rv))
    }
  }
}
cat("\nWORST relative gradient difference:", format(worst), "\n")
