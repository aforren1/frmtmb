# Lane wt-arcovsample, punch round 1: which rows of a group are shifted,
# and by how many lags.
#
# The review (dev/reviews/2026-09-29-arcovsample.md, section 6) measured
# that only the FIRST row of a group is unshifted, at every order, where
# ?sample-log_lik said the first max(p, q) rows were. This reproduces
# that on its own seed and then measures the RULE the replacement
# sentence has to state: which lag first reaches which within-group
# position. Lag i is detected by perturbing coefficient i alone, which
# is a statement about the shift and not about the code that makes it.
#
#   Rscript dev/arcovsample-firstrows.R > dev/arcovsample-log/firstrows.txt

LIB <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
stopifnot("arma_cond_resp" %in% getNamespaceExports("frmtmb"))
cat("frmtmb", format(packageVersion("frmtmb")), "\n")

SEED <- 6162L
set.seed(SEED)
# groups of UNEQUAL length, including one of length 1 and one shorter
# than max(p, q) + 1, so the rule is measured where it is awkward
lens <- c(9L, 7L, 3L, 1L, 6L)
dd <- do.call(rbind, lapply(seq_along(lens), function(i) {
  data.frame(g = factor(i, levels = seq_along(lens)),
             t = seq_len(lens[i]))
}))
dd$x <- stats::rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + stats::rnorm(nrow(dd))
dd <- dd[sample(nrow(dd)), ]
rownames(dd) <- NULL
cat("N =", nrow(dd), " group lengths", paste(lens, collapse = ","),
    " seed =", SEED, "\n")

# within-group position of every row, in the (gr, time) order the block
# is built on
ord <- order(dd$g, dd$t)
pos <- integer(nrow(dd))
pos[ord] <- unlist(lapply(split(seq_along(ord), dd$g[ord]), seq_along))

shift_at <- function(f, th) {
  ac <- f$frame[["autocor"]][[1L]]
  d0 <- eval_dpars(f)
  mu0 <- as.numeric(d0$y$mu)
  mu0 + frmtmb:::autocor_cond_shift(as.numeric(dd$y) - mu0, th, ac) - mu0
}

report <- function(lab, form) {
  f <- frm(form, family = gaussian(), data = dd)
  ac <- f$frame[["autocor"]][[1L]]
  p <- ac[["p"]]
  q <- ac[["q"]]
  k <- p + q
  th <- seq(0.31, by = 0.07, length.out = k)
  s0 <- shift_at(f, th)
  cat("\n---- ", lab, ": p = ", p, ", q = ", q, ", max(p, q) = ",
      max(p, q), "\n", sep = "")
  # (i) which positions have NO shift at all
  nz <- tapply(abs(s0) > 1e-12, pos, sum)
  n <- tapply(abs(s0) > 1e-12, pos, length)
  none <- as.integer(names(nz))[vapply(nz, function(v) v == 0L, TRUE)]
  cat("  positions with NO shift in ANY group : ",
      paste(none, collapse = ","), "\n", sep = "")
  cat("  shifted groups by position           : ",
      paste(sprintf("%d:%d/%d", as.integer(names(nz)), unlist(nz),
                    unlist(n)), collapse = " "), "\n", sep = "")
  # (ii) the first position each COEFFICIENT reaches. Perturb one
  # coefficient and see which positions move: that names the lag
  # without reading the recursion.
  for (j in seq_len(k)) {
    th2 <- th
    th2[j] <- th2[j] + 0.13
    moved <- abs(shift_at(f, th2) - s0) > 1e-12
    first <- if (any(moved)) min(pos[moved]) else NA_integer_
    nm <- if (j <= p) paste0("ar[", j, "]") else paste0("ma[", j - p, "]")
    lag <- if (j <= p) j else j - p
    cat("    ", nm, " (lag ", lag, ") first moves position ", first,
        "\n", sep = "")
  }
}

report("ar(t, g, p = 1)", bf(y ~ x + ar(t, g, p = 1)))
report("ar(t, g, p = 2)", bf(y ~ x + ar(t, g, p = 2)))
report("ar(t, g, p = 3)", bf(y ~ x + ar(t, g, p = 3)))
report("ma(t, g, q = 2)", bf(y ~ x + ma(t, g, q = 2)))
report("arma(t, g, p = 2, q = 2)", bf(y ~ x + arma(t, g, p = 2, q = 2)))
report("arma(t, g, p = 3, q = 1)", bf(y ~ x + arma(t, g, p = 3, q = 1)))
cat("\nDONE\n")
