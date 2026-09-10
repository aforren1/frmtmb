# lane gddm, punch round 2: the `poly()` scan, rerun.
#
# The shape most at risk from a tighter tolerance is a computed column
# whose entries are near zero, where the per-pair term nearly vanishes
# and only the floor is left. `poly()` is the only such column anybody
# has found, so it is scanned wider than the sweep's own design: 4 to
# 12 levels by degrees 1 to 4, reporting the worst
# `|row - its condition's first row| / its own tolerance`. A ratio
# below 1 is accepted; the reciprocal is the headroom.
#
# The tolerance is read out of the package rather than retyped, so a
# change to the shipped rule moves this measurement with it.
#
# Seed 202609, the sweep's own. Arm chosen by GDDM_LIB.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, "  frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")

set.seed(202609)
n <- 480L
worst <- 0
worst_at <- ""
flagged <- 0L
cells <- 0L
cat(sprintf("\n%-8s %-8s %-13s %-13s %s\n", "levels", "degree",
            "worst dev", "dev / its tol", "flagged"))
for (L in c(4L, 6L, 8L, 12L)) {
  z <- rep(seq_len(L), length.out = n)
  gj <- match(z, sort(unique(z)))
  fj <- match(seq_len(L), gj)
  ref <- fj[gj]
  for (dg in seq_len(min(4L, L - 1L))) {
    P <- stats::poly(z, dg)
    dv <- abs(P - P[ref, , drop = FALSE])
    # the shipped rule, per pair plus the column floor
    colf <- apply(P, 2L, function(v) frmtmb.eam:::gd_col_floor(v))
    tol <- 1e-8 * pmax(abs(P), abs(P[ref, , drop = FALSE])) +
      rep(colf, each = nrow(P))
    r <- dv / tol
    r[dv == 0] <- 0
    fl <- length(frmtmb.eam:::gd_varying_groups(P, gj, fj))
    cells <- cells + 1L
    flagged <- flagged + (fl > 0L)
    if (max(r) > worst) {
      worst <- max(r)
      worst_at <- paste0(L, " levels, degree ", dg)
    }
    cat(sprintf("%-8d %-8d %-13.3g %-13.3g %d\n", L, dg, max(dv),
                max(r), fl))
  }
}
cat(sprintf("\nworst ratio anywhere: %.3g, at %s\n", worst, worst_at))
cat(sprintf("headroom, the reciprocal: %.3g x\n", 1 / worst))
cat(sprintf("cells flagged: %d of %d (0 is right; every one of these",
            flagged, cells))
cat(" is a correct design)\n")

# The review's own design of the same family, reproduced, because two
# scans of one family that disagree mean one of them has not been
# understood. It differs in three ways: 120 rows rather than 480,
# levels equally spaced on (0, 1) rather than the integers, and blocked
# rather than interleaved. All three make it gentler, and the harsher
# arm above is the one to read.
cat("\n== the review's design, reproduced\n")
n2 <- 120L
w2 <- 0
fl2 <- 0L
c2 <- 0L
for (nl in c(4L, 6L, 8L, 12L)) {
  for (deg in 1:4) {
    if (deg >= nl) next
    x <- rep(seq(0, 1, length.out = nl), each = n2 / nl)
    P <- stats::poly(x, deg)
    gj <- match(x, sort(unique(x)))
    fj <- match(seq_len(nl), gj)
    ref <- fj[gj]
    colf <- apply(P, 2L, function(v) frmtmb.eam:::gd_col_floor(v))
    tol <- 1e-8 * pmax(abs(P), abs(P[ref, , drop = FALSE])) +
      rep(colf, each = nrow(P))
    dv <- abs(P - P[ref, , drop = FALSE])
    r <- dv / tol
    r[dv == 0] <- 0
    c2 <- c2 + 1L
    fl2 <- fl2 + (length(frmtmb.eam:::gd_varying_groups(P, gj, fj)) > 0L)
    w2 <- max(w2, max(r))
    cat(sprintf("   %2d levels, poly degree %d: worst/tol = %.3e\n",
                nl, deg, max(r)))
  }
}
cat(sprintf("worst ratio: %.3g, headroom %.3g x, cells flagged %d of %d\n",
            w2, 1 / w2, fl2, c2))
