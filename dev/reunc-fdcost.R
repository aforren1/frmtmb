# Lane wt-reunc, punch round 1 (M1): what the ordinal and categorical
# fitted() route costs, and the proof that bounding it changes no
# answer.
#
# The finite-difference route differences the group effects too. It used
# to difference EVERY kept level, which is one model evaluation per
# level whether or not the rows being predicted load it;
# `re_used_b()` bounds the set to the levels the rows actually load.
#
# Two things are measured here:
#  1. THE ANSWER IS THE SAME. On a fit where the two sets differ (a few
#     newdata rows out of many levels), fitted() is compared against
#     fit_fd_se() driven with the FULL set, one coefficient at a time.
#     Two things are being checked at once and they differ in kind. The
#     BOUND is exact: an unloaded level's Jacobian column is exactly
#     zero and dropping it changes no sum. The BATCH (the levels of one
#     block differenced together and attributed by row) is exact too,
#     but it lets the quadratic form be taken over the block structure
#     rather than over a dense Jacobian, and that adds the same terms in
#     a different ORDER. So the comparison is reported as a relative
#     difference and a count of ulps rather than as identity.
#  2. THE COST. Timed at 20, 100, 200 and 1000 grouping levels, in
#     sample and on 5 newdata rows, against the same fit's
#     re_formula = NA call, which does not difference b at all.
#
#   Rscript dev/reunc-fdcost.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")
# the base build is the ONLY one under rellib-r3, so this does not
# depend on what a later worker names its own private library
lane <- !grepl("rellib-r3", find.package("frmtmb"), fixed = TRUE)

make_fit <- function(ng, m = 6, seed = 5) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(seq_len(ng), each = m)),
                  x = stats::rnorm(ng * m))
  lat <- 0.8 * d$x + stats::rnorm(ng, 0, 1)[d$g] + stats::rlogis(ng * m)
  d$y <- factor(cut(lat, c(-Inf, -1, 0.5, Inf), labels = FALSE),
                ordered = TRUE)
  list(fit = suppressWarnings(frm(bf(y ~ x + (1 | g)) + cumulative(),
                                  data = d)), d = d)
}
# a FUNCTION, not an expression: a promise evaluates once, so the first
# spelling of this timer ran the call in round 1 and nothing in rounds 2
# and 3, and reported 0.00 s for every cell
tm <- function(fun, rounds = 3L) {
  best <- Inf
  for (i in seq_len(rounds)) {
    t0 <- proc.time()[["elapsed"]]
    invisible(fun())
    best <- min(best, proc.time()[["elapsed"]] - t0)
  }
  best
}

cat("\n1. the bound changes no answer (lane build only)\n")
if (lane) {
  z <- make_fit(20)
  fit <- z$fit
  nd <- z$d[c(1, 61), c("x", "g")]
  gov <- frmtmb:::re_governed_b(fit)
  used <- frmtmb:::re_used_b(fit, nd, "y", FALSE)
  cat(sprintf("  20 levels: kept %d, loaded by the 2 rows %d (%s)\n",
              length(gov), length(used), paste(used, collapse = ",")))
  f <- function(x) frmtmb:::fitted_point(x, nd)
  full <- frmtmb:::fit_fd_se(fit, f, b_idx = gov)
  bounded <- fitted(fit, newdata = nd)[, "Est.Error", ]
  # The BOUND is exact arithmetic: an unloaded level's column is zero,
  # and dropping a zero column changes no sum. The BATCH is exact too,
  # but the sparse quadratic form it allows adds the terms in a
  # different order, so the comparison is reported as a relative
  # difference and a count of ulps, not as identity.
  cmp <- function(a, b) {
    a <- unname(as.vector(a))
    b <- unname(as.vector(b))
    cat(sprintf(paste0("    identical %s; max abs %.3g, max relative %.3g",
                       ", max ulp %.1f\n"),
                identical(a, b), max(abs(a - b)),
                max(abs(a - b) / abs(a)),
                max(abs(a - b) / (.Machine$double.eps * abs(a)))))
  }
  cat("  against differencing every kept level, one at a time:\n")
  cmp(full, bounded)
  # and in sample, where the two sets are the SAME, the saving comes
  # from batching the block instead, which must also change no answer
  usedi <- frmtmb:::re_used_b(fit, NULL, "y", FALSE)
  cat(sprintf("  in sample: loaded %d of %d kept\n", length(usedi),
              length(gov)))
  fin <- function(x) frmtmb:::fitted_point(x)
  one <- frmtmb:::fit_fd_se(fit, fin, b_idx = usedi)
  bat <- fitted(fit)[, "Est.Error", ]
  cat("  in sample, batched against one at a time:\n")
  cmp(one, bat)
  # two blocks, so two batches, and a multi-membership block, which
  # cannot be batched and must fall back
  set.seed(11)
  d2 <- z$d
  d2$h <- factor(rep(1:5, length.out = nrow(d2)))
  f2 <- suppressWarnings(frm(bf(y ~ x + (1 | g) + (1 | h)) + cumulative(),
                             data = d2))
  b2 <- frmtmb:::re_governed_b(f2)
  o2 <- frmtmb:::fit_fd_se(f2, function(x) frmtmb:::fitted_point(x),
                           b_idx = b2)
  cat("  two blocks, batched against one at a time:\n")
  cmp(o2, fitted(f2)[, "Est.Error", ])
  # a random SLOPE block: two columns per row, so the batch has to be
  # per column position rather than per block
  fs2 <- suppressWarnings(frm(bf(y ~ x + (1 + x | g)) + cumulative(),
                              data = d2))
  bs2 <- frmtmb:::re_governed_b(fs2)
  bts <- frmtmb:::re_b_batches(fs2, NULL, "y", FALSE, bs2)
  os2 <- frmtmb:::fit_fd_se(fs2, function(x) frmtmb:::fitted_point(x),
                            b_idx = bs2)
  cat(sprintf("  random slope: %d coefficients in %s batches\n",
              length(bs2),
              if (is.null(bts)) "no" else paste(length(bts))))
  cmp(os2, fitted(fs2)[, "Est.Error", ])
  d3 <- d2
  d3$h2 <- factor(rep(c(2:5, 1), length.out = nrow(d3)))
  f3 <- tryCatch(suppressWarnings(
    frm(bf(y ~ x + (1 | mm(h, h2))) + cumulative(), data = d3)),
    error = function(e) NULL)
  if (!is.null(f3)) {
    b3 <- frmtmb:::re_governed_b(f3)
    bt3 <- frmtmb:::re_b_batches(f3, NULL, "y", FALSE, b3)
    o3 <- frmtmb:::fit_fd_se(f3, function(x) frmtmb:::fitted_point(x),
                             b_idx = b3)
    cat(sprintf("  multi-membership: batches %s, identical: %s\n",
                if (is.null(bt3)) "refused (as it must)" else
                  paste(length(bt3)),
                identical(unname(as.vector(o3)),
                          unname(as.vector(fitted(f3)[, "Est.Error", ])))))
  }
} else {
  cat("  (base build has no bound to prove)\n")
}

# The clock on this machine measures LOAD as much as work: a lane
# running its own suite beside this one moved every cell of the table
# below by three to five times, including the nd_NA control, which is
# 0.61.0's own route and is not touched here. So the primary figure is
# the one thing the fix actually changes and the machine cannot move:
# the NUMBER OF MODEL EVALUATIONS the finite-difference route makes.
# fit_fd_se() calls its function once for the point, twice per outer
# parameter, and then twice per BATCH, or twice per coefficient when
# there is no batch. Counting them is exact and reproducible.
cat("\n2a. model evaluations, which do not depend on the machine\n")
if (lane) {
  cat(sprintf("%6s %8s %8s %8s %8s\n", "levels", "in", "nd", "in_full",
              "nd_full"))
  for (ng in c(20, 100, 200, 1000)) {
    z <- make_fit(ng)
    fit <- z$fit
    nd <- z$d[c(1, 7, 13, 19, 25), c("x", "g")]
    count <- function(newdata, bounded) {
      n <- 0L
      f <- function(x) { n <<- n + 1L; frmtmb:::fitted_point(x, newdata) }
      gov <- frmtmb:::re_governed_b(fit)
      bi <- gov
      bt <- NULL
      if (bounded) {
        used <- frmtmb:::re_used_b(fit, newdata, "y", FALSE)
        if (!is.null(used)) bi <- intersect(gov, used)
        bt <- frmtmb:::re_b_batches(fit, newdata, "y", FALSE, bi)
      }
      frmtmb:::fit_fd_se(fit, f, b_idx = bi, b_batch = bt)
      n
    }
    cat(sprintf("%6d %8d %8d %8d %8d\n", ng,
                count(NULL, TRUE), count(nd, TRUE),
                count(NULL, FALSE), count(nd, FALSE)))
  }
}

cat("\n2b. cost, best of 3, seconds. READ THE CONTROL FIRST: nd_NA is\n")
cat("   0.61.0's own route and this lane does not touch it, so a\n")
cat("   nd_NA that moves between runs says the machine moved.\n")
cat("   in/nd: fitted() in sample and on 5 newdata rows, as shipped\n")
cat("   in_full/nd_full: the same two calls differencing EVERY kept\n")
cat("   level, one coefficient at a time, which is what this lane did\n")
cat("   before the bound and the batch. in_full at 1000 levels is one\n")
cat("   round rather than the best of three, because it is minutes.\n")
cat("   nd_NA: re_formula = NA, which differences no group effect at\n")
cat("   all and is what 0.61.0 cost\n")
cat(sprintf("%6s %6s %8s %8s %8s %8s %8s\n", "levels", "rows", "in",
            "nd", "in_full", "nd_full", "nd_NA"))
for (ng in c(20, 100, 200, 1000)) {
  z <- make_fit(ng)
  fit <- z$fit
  nd <- z$d[c(1, 7, 13, 19, 25), c("x", "g")]
  t_in <- tm(function() fitted(fit))
  t_nd <- tm(function() fitted(fit, newdata = nd))
  t_na <- tm(function() fitted(fit, newdata = nd, re_formula = NA))
  t_full <- NA_real_
  t_infull <- NA_real_
  if (lane) {
    gov <- frmtmb:::re_governed_b(fit)
    t_full <- tm(function() {
      frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x, nd),
                         b_idx = gov)
    })
    t_infull <- tm(function() {
      frmtmb:::fit_fd_se(fit, function(x) frmtmb:::fitted_point(x),
                         b_idx = gov)
    }, rounds = if (ng >= 1000) 1L else 3L)
  }
  cat(sprintf("%6d %6d %8.2f %8.2f %8.2f %8.2f %8.2f\n", ng, nrow(z$d),
              t_in, t_nd, t_infull, t_full, t_na))
}
