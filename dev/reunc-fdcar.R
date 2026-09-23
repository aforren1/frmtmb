# Lane wt-reunc, punch round 2 (B1): the finite-difference route on
# car() blocks, where the batch's premise fails for one type only.
#
# `car(type = "esicar")` is the hard-constrained intrinsic CAR, whose
# coefficient vector is `car_center(b, a) = b - m[comp]`, the parameter
# with its connected component's mean removed. So perturbing every
# level of the block at once moves each row by its own step MINUS the
# mean of all the steps, and the difference read off that row is not
# the one a single-coefficient perturbation gives. `length(c_idx)` and
# `length(b_idx)` are equal there, so the batch's old proxy for "not
# positionwise" let it through. `re_used_b()`'s bound fails the same
# way: under centering every `b` of a component reaches every
# coefficient of that component, so a bound taken from the loaded
# columns drops Jacobian columns that are not zero.
#
# icar and escar do not center and are unaffected, which is what makes
# this a property of the TYPE and not of car().
#
# The reference is fit_fd_se() driven with the full governed set and no
# batches, which is one coefficient at a time. It is VALIDATED first,
# against the pre-batch route, rather than assumed to be right.
#
#   Rscript dev/reunc-fdcar.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc2-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

relmax <- function(a, b) {
  a <- unname(as.vector(a)); b <- unname(as.vector(b))
  max(abs(a - b) / abs(a))
}
ulpmax <- function(a, b) {
  a <- unname(as.vector(a)); b <- unname(as.vector(b))
  max(abs(a - b) / (.Machine$double.eps * abs(a)))
}

# a 3 x 3 rook lattice, 20 rows per location
lat <- expand.grid(r = 1:3, c = 1:3)
loc <- factor(seq_len(nrow(lat)))
W <- matrix(0L, nrow(lat), nrow(lat))
for (i in seq_len(nrow(lat))) {
  for (j in seq_len(nrow(lat))) {
    if (abs(lat$r[i] - lat$r[j]) + abs(lat$c[i] - lat$c[j]) == 1L) {
      W[i, j] <- 1L
    }
  }
}
dimnames(W) <- list(levels(loc), levels(loc))

make_car <- function(type, extra = FALSE, seed = 31) {
  set.seed(seed)
  m <- 20
  d <- data.frame(g = factor(rep(levels(loc), each = m), levels = levels(loc)),
                  x = stats::rnorm(nrow(lat) * m))
  d$h <- factor(rep(1:4, length.out = nrow(d)))
  u <- stats::rnorm(nrow(lat), 0, 0.8)
  u <- u - mean(u)
  eta <- 0.6 * d$x + u[as.integer(d$g)] +
    if (extra) stats::rnorm(4, 0, 0.5)[as.integer(d$h)] else 0
  y <- factor(cut(eta + stats::rlogis(nrow(d)), c(-Inf, -0.8, 0.8, Inf),
                  labels = FALSE), ordered = TRUE)
  d$y <- y
  f <- if (extra) {
    bf(y ~ x + car(W, gr = g, type = type) + (1 | h))
  } else {
    bf(y ~ x + car(W, gr = g, type = type))
  }
  list(fit = suppressWarnings(frm(f + cumulative(), data = d)), d = d)
}

# the reference, one coefficient at a time, and its validation
ref <- function(fit, nd = NULL) {
  gov <- frmtmb:::re_governed_b(fit)
  frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x, nd),
                     b_idx = gov, b_batch = NULL)
}
shipped <- function(fit, nd = NULL) {
  a <- if (is.null(nd)) fitted(fit) else fitted(fit, newdata = nd)
  a[, "Est.Error", ]
}

cat("\n0. the reference is one coefficient at a time, checked\n")
z <- make_car("esicar")
gov <- frmtmb:::re_governed_b(z$fit)
bt <- frmtmb:::re_b_batches(z$fit, NULL, "y", FALSE, gov)
cat("   esicar: kept ", length(gov), " effects; re_b_batches() says ",
    if (is.null(bt)) "REFUSED" else paste(length(bt), "batches"), "\n",
    sep = "")
# driving fit_fd_se() with an explicit empty batch list is the same
# code path as a refusal, so the two must agree bit for bit
one <- frmtmb:::fit_fd_se(z$fit, function(x) frmtmb:::fitted_point(x),
                          b_idx = gov, b_batch = list())
cat("   b_batch = NULL equals b_batch = list(): ",
    identical(unname(as.vector(ref(z$fit))), unname(as.vector(one))),
    "\n", sep = "")

cat("\n1. every car type, in sample and on 3 newdata rows\n")
cat(sprintf("%-22s %10s %10s %8s\n", "fit", "in (rel)", "nd (rel)",
            "batches"))
for (ty in c("icar", "escar", "esicar")) {
  for (extra in c(FALSE, TRUE)) {
    z <- tryCatch(make_car(ty, extra), error = function(e) NULL)
    if (is.null(z)) {
      cat(sprintf("%-22s %10s\n", paste0(ty, if (extra) " + (1|h)"),
                  "no fit"))
      next
    }
    nd <- z$d[c(1, 25, 45), c("x", "g", "h")]
    g <- frmtmb:::re_governed_b(z$fit)
    b <- frmtmb:::re_b_batches(z$fit, NULL, "y", FALSE, g)
    tag <- paste0(ty, if (extra) " + (1|h)")
    cat(sprintf("%-22s %10.3g %10.3g %8s\n", tag,
                relmax(ref(z$fit), shipped(z$fit)),
                relmax(ref(z$fit, nd), shipped(z$fit, nd)),
                if (is.null(b)) "refused" else length(b)))
  }
}

cat("\n2. the same four numbers as ulps, where they should be exact\n")
for (ty in c("icar", "escar", "esicar")) {
  z <- make_car(ty)
  nd <- z$d[c(1, 25, 45), c("x", "g", "h")]
  cat(sprintf("   %-8s in %6.1f ulp, nd %6.1f ulp\n", ty,
              ulpmax(ref(z$fit), shipped(z$fit)),
              ulpmax(ref(z$fit, nd), shipped(z$fit, nd))))
}

cat("\n3. the bound alone, which fails the same way under centering\n")
z <- make_car("esicar")
nd <- z$d[c(1, 25, 45), c("x", "g", "h")]
g <- frmtmb:::re_governed_b(z$fit)
u <- frmtmb:::re_used_b(z$fit, nd, "y", FALSE)
cat("   kept ", length(g), ", bound to ", length(u),
    " (equal means the whole block joins, as it must)\n", sep = "")

cat("\n4. the defect, reconstructed: the OLD proxy for 'not positionwise'\n")
# block_b_positionwise() replaced by the length test the batch used
# before, which is TRUE for esicar. The swap must change the BATCH and
# not the MODEL, so the point estimate is checked across it: expand_b()
# reaches the same predicate through frame_needs_expand(), and the
# frame's own has_expand flag is what keeps the centering on.
old_pred <- function(bk) {
  length(bk[["c_idx"]]) == length(bk[["b_idx"]])
}
for (extra in c(FALSE, TRUE)) {
  z <- make_car("esicar", extra)
  nd <- z$d[c(1, 25, 45), c("x", "g", "h")]
  pt_in <- frmtmb:::fitted_point(z$fit)
  pt_nd <- frmtmb:::fitted_point(z$fit, nd)
  r_in <- ref(z$fit)
  r_nd <- ref(z$fit, nd)
  assignInNamespace("block_b_positionwise", old_pred, ns = "frmtmb")
  same_model <- identical(pt_in, frmtmb:::fitted_point(z$fit)) &&
    identical(pt_nd, frmtmb:::fitted_point(z$fit, nd))
  g <- frmtmb:::re_governed_b(z$fit)
  b <- frmtmb:::re_b_batches(z$fit, NULL, "y", FALSE, g)
  u <- frmtmb:::re_used_b(z$fit, nd, "y", FALSE)
  old_in <- shipped(z$fit)
  old_nd <- shipped(z$fit, nd)
  assignInNamespace("block_b_positionwise",
                    function(bk) {
                      !identical(bk[["covstruct"]], "rr") &&
                        !frmtmb:::block_is_esicar(bk)
                    }, ns = "frmtmb")
  cat("   esicar", if (extra) " + (1|h)" else "", ": model unchanged by ",
      "the swap ", same_model, "\n", sep = "")
  cat("     old batches ", if (is.null(b)) "refused" else length(b),
      ", old bound ", length(u), " of ", length(g), "\n", sep = "")
  cat(sprintf("     old in %.3g rel, old nd %.3g rel; fixed %.1f, %.1f ulp\n",
              relmax(r_in, old_in), relmax(r_nd, old_nd),
              ulpmax(r_in, shipped(z$fit)),
              ulpmax(r_nd, shipped(z$fit, nd))))
  # which way the error goes is not fixed: the variance carries the
  # cross term 2 Jd' Vdb Jb, so a mis-attributed b column can move a
  # cell either way. Report both extremes of the signed ratio.
  rr_in <- unname(as.vector(old_in)) / unname(as.vector(r_in))
  cat(sprintf("     old/reference in sample: min %.4f, max %.4f\n",
              min(rr_in), max(rr_in)))
}

cat("\n5. rr, which the old proxy caught only below full rank\n")
set.seed(41)
nd2 <- data.frame(g = factor(rep(1:8, each = 12)),
                  x = stats::rnorm(96))
nd2$y <- factor(cut(0.5 * nd2$x + stats::rnorm(8, 0, 0.8)[nd2$g] +
                      stats::rlogis(96), c(-Inf, -0.8, 0.8, Inf),
                    labels = FALSE), ordered = TRUE)
for (rk in c(2, 1)) {
  f <- tryCatch(suppressWarnings(
    frm(bf(y ~ x + rr(1 + x | g, rank = rk)) + cumulative(), data = nd2)),
    error = function(e) NULL)
  if (is.null(f)) { cat("   rank ", rk, ": no fit\n", sep = ""); next }
  g2 <- frmtmb:::re_governed_b(f)
  b2 <- frmtmb:::re_b_batches(f, NULL, "y", FALSE, g2)
  cat("   rank ", rk, ": kept ", length(g2), ", batches ",
      if (is.null(b2)) "refused" else length(b2),
      ", shipped vs reference ",
      sprintf("%.3g", relmax(ref(f), shipped(f))), " rel\n", sep = "")
}
