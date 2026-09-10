# lane gddm, punch round 2: does the computed-column predicate hold?
#
# The floor `1e-11 * max|column|` is now charged only where the model
# frame RAN arithmetic, which the check reads off the column name:
# `is.call(str2lang(cn))`. Two things have to be true before that is
# safe, and neither is safe to assume.
#
#  1. every term whose values are computed must still get the floor, in
#     particular scale(), bs(), ns(), log1p() and cut(), which an
#     earlier pass showed are row-deterministic. Row-deterministic is
#     not the same as needing no floor, so the predicate is checked on
#     the NAME rather than on the earlier measurement.
#  2. a bare symbol whose column is nevertheless computed upstream
#     would lose its floor. This looks for one.
#
# Section 3 measures what the predicate actually buys: the band on a
# wide-range bare column, before and after.
#
# Seed 4242. Arm chosen by GDDM_LIB.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, "  frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")

# The two rules. `with_floor` is what the package ships; `no_floor` is
# the PROPOSED form, written out here because the package does not
# carry it, so both arms are measured in one process against the same
# columns.
with_floor <- function(v, gj, fj)
  frmtmb.eam:::gd_varying_groups(v, gj, fj)
no_floor <- function(v, gj, fj) {
  ref <- fj[gj]
  tol <- 1e-8 * pmax(abs(v), abs(v[ref]))
  sort(unique(gj[which(abs(v - v[ref]) > tol)]))
}

set.seed(4242)
n <- 120L
ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$z <- rep(c(0, 0.128, 0.256, 0.512), length.out = n)
d$g <- factor(rep(1:4, length.out = n))
d$ord <- factor(d$z, ordered = TRUE)
d$off <- d$z / 10
d$cond <- gddm_conditions(d, z)
# a bare symbol whose values were computed BEFORE frm() saw them: the
# case the predicate could get wrong
d$pre <- as.numeric(stats::poly(d$z, 2)[, 2])
d$prescale <- as.numeric(scale(d$z))
d$preave <- stats::ave(as.numeric(d$rt), d$cond, FUN = mean)

# 1. which model-frame columns the predicate calls computed, and
#    whether the design is accepted.
probe <- function(label, mu) {
  f <- stats::as.formula(call("~", quote(rt | vint(upper, cond)), mu))
  sp <- bf(f, bs ~ 1, ndt ~ 1, bias = 0.5)
  fr <- tryCatch(frm(sp, family = gddm(control = ctl), data = d,
                     dry_run = "frame"),
                 error = function(e) conditionMessage(e))
  if (is.character(fr)) {
    cat(sprintf("%-26s REFUSED: %s\n", label, substr(fr, 8, 60)))
    return(invisible(NULL))
  }
  cn <- names(fr[["data_frame"]])
  keep <- setdiff(cn, c("rt", "upper", "cond"))
  isc <- vapply(keep, function(z)
    is.call(tryCatch(str2lang(z), error = function(e) NULL)), TRUE)
  cat(sprintf("%-26s accepted; columns: %s\n", label,
              paste0(keep, if (length(keep)) "=" else "",
                     ifelse(isc, "computed", "as supplied"),
                     collapse = "  ")))
}

cat("\n== 1. what the predicate calls computed\n")
probe("bare symbol", quote(z))
probe("poly()", quote(poly(z, 2)))
probe("scale()", quote(scale(z)))
probe("bs()", quote(splines::bs(z, df = 3)))
probe("ns()", quote(splines::ns(z, df = 3)))
probe("log1p()", quote(log1p(z)))
probe("cut()", quote(cut(z, 3)))
probe("I()", quote(I(z^2)))
probe("offset()", quote(offset(off)))
probe("mo()", quote(mo(ord)))
probe("s()", quote(s(z, k = 4)))
probe("factor grouping", quote(1 + (1 | g)))
probe("precomputed poly column", quote(pre))
probe("precomputed scale column", quote(prescale))
probe("precomputed ave column", quote(preave))


# 2. the five calls that must KEEP the floor, checked on the tolerance
#    itself rather than on a verdict that could come out right for
#    another reason. A ratio of 0/0 is reported as 0: where both the
#    deviation and the tolerance are zero the entries are bitwise equal
#    and there is nothing to carry.
cat("\n== 2. the floor is still charged to a call, and the slack each",
    "one has\n")
gi <- match(d$cond, sort(unique(d$cond)))
first <- match(seq_along(sort(unique(d$cond))), gi)
ref <- first[gi]
ratio <- function(dv, tol) {
  r <- dv / tol
  r[dv == 0] <- 0
  max(r)
}
slack <- function(label, M) {
  M <- as.matrix(M)
  r <- M[ref, , drop = FALSE]
  dv <- abs(M - r)
  colmax <- apply(abs(M), 2L, function(z) max(z[is.finite(z)]))
  tol_no <- 1e-8 * pmax(abs(M), abs(r))
  tol_yes <- tol_no + rep(1e-11 * colmax, each = nrow(M))
  cat(sprintf("  %-24s worst dev %-10.3g dev/tol %.3g", label,
              max(dv), ratio(dv, tol_yes)))
  cat(sprintf("   without the floor %.3g\n", ratio(dv, tol_no)))
}
slack("poly(z, 2)", stats::poly(d$z, 2))
slack("scale(z)", scale(d$z))
slack("splines::bs(z, df=3)", splines::bs(d$z, df = 3))
slack("splines::ns(z, df=3)", splines::ns(d$z, df = 3))
slack("log1p(z)", log1p(d$z))
slack("z, a bare symbol", d$z)

# 2b. The case the predicate could get wrong: a column computed
#     UPSTREAM and stored under a bare name, so it carries arithmetic
#     the model frame did not do and loses the floor. The dangerous
#     shape is a CENTERED one, where an entry sits at zero and the
#     per-pair term goes to zero with it.
cat("\n== 2b. a bare symbol whose values were computed before frm()\n")
zc <- rep(c(0, 1, 2), length.out = n)     # mean 1, so poly hits zero
dc <- d
dc$cond <- gddm_conditions(dc, zc = zc)
dc$zc <- zc
dc$p1 <- as.numeric(stats::poly(zc, 1))
dc$sc <- as.numeric(scale(zc))
cat(sprintf("  poly(zc, 1) at the middle level: %.3g  (entries near",
            max(abs(dc$p1[zc == 1]))))
cat(" zero are where a per-pair tolerance vanishes)\n")
for (nm in c("zc", "p1", "sc")) {
  gj <- match(dc$cond, sort(unique(dc$cond)))
  fj <- match(seq_along(sort(unique(dc$cond))), gj)
  cat(sprintf("  %-4s stored, treated as supplied: %d conditions",
              nm, length(no_floor(dc[[nm]], gj, fj))))
  cat(sprintf("   with a floor: %d\n",
              length(with_floor(dc[[nm]], gj, fj))))
}
r <- tryCatch({
  frm(bf(rt | vint(upper, cond) ~ p1 + sc, bs ~ 1, ndt ~ 1,
         bias = 0.5),
      family = gddm(control = ctl), data = dc, dry_run = "frame")
  "accepted"
}, error = function(e) paste("REFUSED:", substr(conditionMessage(e),
                                                8, 70)))
cat("  a model on both precomputed columns:", r, "\n")

# 3. what the predicate buys: the review's sentinel column. The
#    sentinel has to sit in ANOTHER condition, so that the flagged
#    comparison is the 1-against-3 pair and the column maximum is the
#    sentinel. A first version put the sentinel inside the same
#    condition, where its own rows are refused and the band never
#    shows.
cat("\n== 3. the band the predicate closes\n")
gi2 <- rep(1:2, each = 20L)
f2 <- c(1L, 21L)
cat(sprintf("  %-10s %-10s %-22s %s\n", "sentinel", "floor",
            "as a call (floor)", "as a bare symbol (none)"))
for (s in c(1e6, 1e9, 1e11, 1e12, 1e15, 1e18)) {
  v <- c(rep(c(1, 3), each = 10L), rep(s, 20L))
  a <- length(with_floor(v, gi2, f2))
  b <- length(no_floor(v, gi2, f2))
  cat(sprintf("  %-10.0e %-10.3g %-22s %s\n", s, 1e-11 * s,
              if (a) "refused" else "ACCEPTED",
              if (b) "refused" else "ACCEPTED"))
}
cat("  values 1 and 3 inside condition 1, sentinel alone in",
    "condition 2\n")

# 2c. How wide is that false alarm? A precomputed orthogonal
#     polynomial stored under a bare name, over symmetric level counts
#     and degrees, reporting the relative size of the deviation the
#     floor would have to carry.
cat("\n== 2c. precomputed poly() stored bare, the shape 2b found\n")
cat(sprintf("  %-8s %-7s %-12s %-12s %-9s %s\n", "levels", "degree",
            "worst dev", "rel to colmax", "no floor", "with floor"))
worst_rel <- 0
nfa <- 0L; ntot <- 0L
for (L in c(3L, 5L, 7L, 9L)) {
  zz <- rep(seq_len(L), length.out = 120L)
  gj <- match(zz, sort(unique(zz)))
  fj <- match(seq_len(L), gj)
  rj <- fj[gj]
  for (dg in seq_len(min(3L, L - 1L))) {
    P <- stats::poly(zz, dg)
    for (k in seq_len(dg)) {
      v <- as.numeric(P[, k])
      dv <- max(abs(v - v[rj]))
      rel <- dv / max(abs(v))
      a <- length(no_floor(v, gj, fj))
      b <- length(with_floor(v, gj, fj))
      ntot <- ntot + 1L
      if (a > 0L) nfa <- nfa + 1L
      if (dv > 0) worst_rel <- max(worst_rel, rel)
      cat(sprintf("  %-8d %-7s %-12.3g %-12.3g %-9s %s\n", L,
                  paste0(dg, ".", k), dv, rel,
                  if (a) "REFUSED" else "accepted",
                  if (b) "REFUSED" else "accepted"))
    }
  }
}
cat(sprintf("\n  false alarms without the floor: %d of %d columns\n",
            nfa, ntot))
cat(sprintf("  worst relative deviation to carry: %.3g, that is %.2f",
            worst_rel, worst_rel / .Machine$double.eps))
cat(" ulp of the column maximum\n")
