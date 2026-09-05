# Agreement between this package's two new families and EMC2.
#
# Not a testthat file, and deliberately so. Every function in EMC2 that
# computes either of these likelihoods is INTERNAL: `dWald`, `pWald`,
# `dDDM`, `pDDM` and `log_likelihood_ddmgng` are all reached with
# `:::`. A package test that did that would be asserting on another
# package's private surface, which CRAN refuses and which would break
# on an EMC2 release that renamed a helper without changing a single
# number.
#
# So the suite checks the two families against EXPORTED references
# (statmod's inverse Gaussian, WienR's and RWiener's Wiener
# distribution functions), and this script checks them against EMC2
# itself. Run it by hand when EMC2 changes, and copy what it prints
# into dev/rdm-gng-findings.md.
#
# Usage:
#   Rscript dev/rdm-gng-emc2-reference.R
#
# It needs EMC2 and frmtmb.ddm on the library path and prints a table.

if (!requireNamespace("EMC2", quietly = TRUE)) {
  stop("EMC2 is not installed; this script has nothing to compare against.")
}
suppressMessages(library(frmtmb.ddm))

rel <- function(a, b) {
  k <- is.finite(a) & is.finite(b) & b != 0
  if (!any(k)) return(NA_real_)
  max(abs(a[k] - b[k]) / abs(b[k]))
}
say <- function(...) cat(sprintf(...), "\n", sep = "")

# --------------------------------------------------------------- RDM
#
# EMC2's parameter map, read out of `EMC2:::rRDM`, which draws the
# distance to travel as `B + runif(1, 0, A)`:
#
#   EMC2 v -> v_i      EMC2 B -> k      EMC2 A -> A      EMC2 t0 -> ndt
#
# so EMC2's B is the gap ABOVE the start-point range, not the
# threshold. The threshold is B + A.

law <- frmtmb.ddm:::rdm_law
race <- frmtmb.ddm:::lba_race_lpdf
acc <- function(v, A, k) list(v = v, A = A, k = k)

say("== RDM: single-accumulator density vs EMC2:::dWald ==")
gr <- expand.grid(t = c(0.05, 0.2, 0.5, 1, 2, 5, 15), v = c(0.3, 1, 2, 4, 8),
                  k = c(0.3, 1, 2.5), A = c(0.1, 0.5, 1.5))
z <- rep(0, nrow(gr))
d_emc <- EMC2:::dWald(gr$t, v = gr$v, B = gr$k, A = gr$A, t0 = z)
d_me <- exp(law$ldens(gr$t, acc(gr$v, gr$A, gr$k)))

# The independent adjudicator, so that a disagreement can be attributed
# rather than merely noted. statmod knows nothing about either package.
gl <- frmtmb.ddm:::ddm_gauss_legendre(60L)
ref <- vapply(seq_len(nrow(gr)), function(i) {
  d <- gr$k[i] + gr$A[i] * gl$x
  sum(gl$w * statmod::dinvgauss(gr$t[i], mean = d / gr$v[i], shape = d^2))
}, numeric(1))

for (thr in c(1e-3, 1e-8, 1e-12)) {
  k <- d_emc > thr
  say("  dWald > %-7g n=%4d   EMC2 vs mine %-10s | mine vs statmod %-10s | EMC2 vs statmod %s",
      thr, sum(k), format(rel(d_me[k], d_emc[k]), digits = 3),
      format(rel(d_me[k], ref[k]), digits = 3),
      format(rel(d_emc[k], ref[k]), digits = 3))
}

say("")
say("== RDM: survival vs 1 - EMC2:::pWald ==")
s_emc <- 1 - EMC2:::pWald(gr$t, v = gr$v, B = gr$k, A = gr$A, t0 = z)
s_me <- exp(law$lsurv(gr$t, acc(gr$v, gr$A, gr$k)))
for (thr in c(1e-3, 1e-6, 1e-9, 1e-12)) {
  k <- s_emc > thr
  say("  1-pWald > %-7g n=%4d   max rel %s", thr, sum(k),
      format(rel(s_me[k], s_emc[k]), digits = 3))
}
say("  1 - pWald returns EXACTLY ZERO on %d of %d rows; this family returns %d zeros.",
    sum(s_emc == 0), length(s_emc), sum(s_me == 0))
if (any(s_emc == 0)) {
  i <- which(s_emc == 0)[1]
  say("  first such row: t=%g v=%g k=%g A=%g, where this family gives %s",
      gr$t[i], gr$v[i], gr$k[i], gr$A[i], format(s_me[i], digits = 4))
}

say("")
say("== RDM: the whole race, as EMC2 composes it ==")
# EMC2:::log_likelihood_race is
#   log(dfun(winner)) + sum(log(1 - pfun(losers)))
# which is this package's race exactly. Composed here from EMC2's own
# two functions so that the comparison is against EMC2 end to end.
# The times and drifts are kept in the regime where EMC2 can score a
# row at all. Wider ranges are not a harder test of the ALGEBRA, they
# are a test of EMC2's tail, which the survival section above already
# reports: at runif(0.25, 2.5) times and runif(0.5, 3.5) drifts, EMC2
# could score 3 of 40 two-accumulator rows, 1 of 40 with three and 1 of
# 40 with four, because a loser that has almost certainly finished sends
# its `1 - pWald` to exactly zero and the log to -Inf. On the rows it
# COULD score there, it disagreed by up to 12.9 in log units, which is
# the same saturation one step before it reaches zero.
set.seed(7)
A <- 0.4; k <- 0.6
for (nacc in c(2L, 3L, 4L)) {
  m <- 200L
  tt <- stats::runif(m, 0.2, 1.2)
  V <- matrix(stats::runif(m * nacc, 0.5, 3), ncol = nacc)
  win <- rep_len(seq_len(nacc), m)
  pars <- lapply(seq_len(nacc), function(j) acc(V[, j], A, k))
  mine <- race(tt, win, law, pars)
  z1 <- rep(0, m)
  # dWald and pWald are Rcpp and do NOT recycle: every argument has to
  # arrive at the length of `t`. Passing B and A as scalars returns
  # quiet garbage rather than an error, which cost an hour here.
  Bv <- rep(k, m); Av <- rep(A, m)
  ref_emc <- log(EMC2:::dWald(tt, v = V[cbind(seq_len(m), win)], B = Bv,
                              A = Av, t0 = z1))
  # The smallest loser survival on each row decides whether EMC2 can
  # score it at all. `1 - pWald` is a subtraction from one, so it holds
  # nothing below about 1e-16 and is already parts-wrong well above
  # that; the comparison is therefore stratified by it rather than
  # summarized over rows EMC2 was never going to get right.
  worst_s <- rep(Inf, m)
  for (j in seq_len(nacc)) {
    lose <- win != j
    if (any(lose)) {
      sj <- 1 - EMC2:::pWald(tt[lose], v = V[lose, j], B = Bv[lose],
                             A = Av[lose], t0 = z1[lose])
      ref_emc[lose] <- ref_emc[lose] + log(sj)
      worst_s[lose] <- pmin(worst_s[lose], sj)
    }
  }
  say("  %d accumulators, %d rows:", nacc, m)
  for (thr in c(1e-4, 1e-8, 1e-12)) {
    ok <- is.finite(ref_emc) & worst_s > thr
    say("    smallest loser survival > %-7g  n=%3d  max abs log-likelihood difference %s",
        thr, sum(ok),
        if (any(ok)) format(max(abs(mine[ok] - ref_emc[ok])), digits = 3)
        else "n/a")
  }
  # Printed only when there IS such a row. The grid above is chosen to
  # keep EMC2 inside its own regime, so on it the count is zero, and a
  # line reading "cannot score 0 of 200 rows ... e.g. n/a" would claim a
  # comparison it never made. The saturation itself is real and is
  # reported by the survival section above, on a grid built to provoke
  # it.
  bad <- !is.finite(ref_emc)
  if (any(bad)) {
    say("    EMC2 cannot score %d of %d rows at all (a loser survival of exactly zero); this family gives every one of them a finite value, e.g. %s",
        sum(bad), m, format(mine[bad][1], digits = 6))
  } else {
    say("    every row scorable by EMC2 on this grid, which is what it was chosen for")
  }
}

# ----------------------------------------------------------- go/no-go
#
# EMC2:::log_likelihood_ddmgng is
#   go   rows: dDDM(rt, R, pars)
#   nogo rows: 1 - pDDM(TIMEOUT, Rgo, pars)
# and dDDM / pDDM are WienR::dWDM / pWDM. The parameter map is
#   EMC2 v -> mu   a -> bs   t0 -> ndt   Z -> bias   s fixed at 1.

say("")
say("== go/no-go: the no-go probability vs EMC2:::pDDM ==")
gr2 <- expand.grid(t = c(0.05, 0.2, 0.6, 1.5, 4), a = c(0.3, 0.8, 1.4, 2.5, 4),
                   v = c(-3, -1, 0, 1.5, 4), w = c(0.3, 0.5, 0.7))
pars <- cbind(v = gr2$v, a = gr2$a, t0 = 0, Z = gr2$w, sv = 0, SZ = 0,
              st0 = 0, s = 1)
R <- factor(rep("upper", nrow(gr2)), levels = c("lower", "upper"))
emc <- 1 - EMC2:::pDDM(gr2$t, R, pars, precision = 1e-12)
mine <- exp(frmtmb.ddm:::ddm_nogo_lprob(gr2$t, gr2$v, gr2$a, gr2$w))
for (thr in c(1e-2, 1e-6, 1e-10, 0)) {
  k <- emc > thr
  say("  EMC2 nogo > %-7g n=%4d   max rel %s", thr, sum(k),
      format(rel(mine[k], emc[k]), digits = 3))
}

say("")
say("== go/no-go: the whole likelihood, as EMC2 composes it ==")
set.seed(11)
d <- wiener_gng_simulate(400, mu = 1.0, bs = 1.4, ndt = 0.25, bias = 0.45,
                         deadline = 1.5)
dp <- list(mu = 1.0, bs = 1.4, ndt = 0.25, bias = 0.45)
mine_ll <- frmtmb.ddm:::gng_lpdf(
  d$rt, lapply(dp, function(v) rep(v, nrow(d))),
  list(dec = d$responded), 1.5)

n <- nrow(d)
pars2 <- cbind(v = rep(dp$mu, n), a = rep(dp$bs, n), t0 = rep(dp$ndt, n),
               Z = rep(dp$bias, n), sv = 0, SZ = 0, st0 = 0, s = 1)
Rg <- factor(rep("upper", n), levels = c("lower", "upper"))
ref_ll <- numeric(n)
go <- d$responded == 1
ref_ll[go] <- log(EMC2:::dDDM(d$rt[go], Rg[go], pars2[go, , drop = FALSE],
                              precision = 1e-12))
ref_ll[!go] <- log(pmax(0, pmin(1, 1 - EMC2:::pDDM(
  rep(1.5, sum(!go)), Rg[!go], pars2[!go, , drop = FALSE],
  precision = 1e-12))))
say("  go rows   n=%4d   max abs difference %s", sum(go),
    format(max(abs(mine_ll[go] - ref_ll[go])), digits = 3))
say("  nogo rows n=%4d   max abs difference %s", sum(!go),
    format(max(abs(mine_ll[!go] - ref_ll[!go])), digits = 3))
say("  total log likelihood: mine %.10f   EMC2 %.10f   difference %s",
    sum(mine_ll), sum(ref_ll),
    format(abs(sum(mine_ll) - sum(ref_ll)), digits = 3))
