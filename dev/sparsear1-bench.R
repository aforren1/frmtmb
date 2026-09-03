# dev/sparsear1-bench.R
#
# The scaling of the ar1() block density in the series length d: the
# dense d x d covariance the registry carried through v0.42.0 against
# the O(d) innovation form that replaced it.
#
# The dense route is reconstructed INLINE here, exactly as it stood, so
# the comparison needs nothing but the installed package. Each cell is
# tape construction plus one fn() and one gr() call on the block
# density alone, which is the unit the audit's roadmap item 1 quotes.
#
# Usage:
#   Rscript dev/sparsear1-bench.R [reps]

suppressMessages(library(RTMB))
pkgload::load_all(".", quiet = TRUE)

reps <- {
  a <- commandArgs(trailingOnly = TRUE)
  if (length(a)) as.integer(a[[1]]) else 3L
}

# ---- the pre-v0.43 dense density, verbatim -------------------------

dense_ar1_nll <- function(b, theta, blk) {
  "[<-" <- RTMB::ADoverload("[<-")
  d <- blk$dim
  sd1 <- exp(theta[1])
  rho <- theta[2] / sqrt(1 + theta[2]^2)
  pows <- rep(rho, d)
  pows[1] <- 1
  for (k in seq_len(d - 1L) + 1L) pows[k] <- pows[k - 1L] * rho
  M <- abs(outer(seq_len(d), seq_len(d), "-")) + 1L
  C <- RTMB::matrix(pows[as.vector(M)], d, d)
  Sigma <- sd1^2 * C
  dim(b) <- c(d, length(b) %/% d)
  sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
}

dense_hetar1_nll <- function(b, theta, blk) {
  "[<-" <- RTMB::ADoverload("[<-")
  d <- blk$dim
  sdv <- exp(theta[seq_len(d)])
  rho <- theta[d + 1L] / sqrt(1 + theta[d + 1L]^2)
  pows <- rep(rho, d)
  pows[1] <- 1
  for (k in seq_len(d - 1L) + 1L) pows[k] <- pows[k - 1L] * rho
  M <- abs(outer(seq_len(d), seq_len(d), "-")) + 1L
  C <- RTMB::matrix(pows[as.vector(M)], d, d)
  Sigma <- C * (RTMB::matrix(sdv, ncol = 1) %*% RTMB::matrix(sdv, nrow = 1))
  dim(b) <- c(d, length(b) %/% d)
  sum(RTMB::dmvnorm(t(b), 0, Sigma, log = TRUE))
}

# ---- the measurement -----------------------------------------------

# One cell: build the tape, call fn once and gr once. `proc.time()` on
# Windows quantizes to 1/100 s, which is coarser than the whole O(d)
# cell, so the batch is repeated until at least `min_secs` of work has
# accumulated and the reported number is the mean. The best of `reps`
# batches is taken, so one garbage collection does not become the
# number.
time_cell <- function(fn, theta, reps, min_secs = 0.5) {
  best <- Inf
  val <- NA_real_
  for (r in seq_len(reps)) {
    k <- 0L
    t0 <- Sys.time()
    repeat {
      tp <- RTMB::MakeTape(fn, theta)
      val <- tp(theta)
      tp$jacobian(theta)
      k <- k + 1L
      el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      if (el >= min_secs || k >= 200L) break
    }
    best <- min(best, el / k)
  }
  list(t = best, v = val)
}

fmt <- function(x) {
  if (is.na(x)) return("not run")
  if (x < 0.01) sprintf("%.1f ms", x * 1000) else sprintf("%.3f s", x)
}
fmt_speedup <- function(old, new) {
  if (is.na(old)) return("-")
  r <- old / new
  if (r < 10) sprintf("%.1fx", r) else sprintf("%.0fx", r)
}

run <- function(dims, n_levels = 1L, dense_cap = Inf, reps = 3L) {
  set.seed(9001)
  rows <- lapply(dims, function(d) {
    b <- stats::rnorm(d * n_levels)
    blk <- list(dim = d, n_levels = n_levels)
    th_a <- c(0.2, 0.9)
    th_h <- c(stats::rnorm(d, 0, 0.3), 0.9)

    new_a <- time_cell(function(th) -covstruct_registry$ar1$nll(b, th, blk),
                       th_a, reps)
    new_h <- time_cell(function(th) -covstruct_registry$hetar1$nll(b, th, blk),
                       th_h, reps)
    if (d <= dense_cap) {
      old_a <- time_cell(function(th) -dense_ar1_nll(b, th, blk), th_a, reps)
      old_h <- time_cell(function(th) -dense_hetar1_nll(b, th, blk), th_h,
                         reps)
    } else {
      old_a <- list(t = NA_real_, v = NA_real_)
      old_h <- list(t = NA_real_, v = NA_real_)
    }
    data.frame(
      d = d, n_levels = n_levels,
      ar1_old = old_a$t, ar1_new = new_a$t,
      het_old = old_h$t, het_new = new_h$t,
      ar1_reldiff = abs(new_a$v - old_a$v) / max(1, abs(old_a$v)),
      het_reldiff = abs(new_h$v - old_h$v) / max(1, abs(old_h$v))
    )
  })
  do.call(rbind, rows)
}

cat("frmtmb ar1(): dense (pre-v0.43) against O(d), tape + fn + gr\n")
cat("R", format(getRversion()), " RTMB", format(packageVersion("RTMB")),
    " best of", reps, "\n\n")

tab <- run(c(50L, 200L, 800L, 2000L), n_levels = 1L, reps = reps)

cat("| block dimension d | ar1 dense | ar1 O(d) | speedup |",
    "hetar1 dense | hetar1 O(d) | speedup |\n")
cat("| --- | --- | --- | --- | --- | --- | --- |\n")
for (i in seq_len(nrow(tab))) {
  r <- tab[i, ]
  cat(sprintf("| %d | %s | %s | %s | %s | %s | %s |\n", r$d,
              fmt(r$ar1_old), fmt(r$ar1_new),
              fmt_speedup(r$ar1_old, r$ar1_new),
              fmt(r$het_old), fmt(r$het_new),
              fmt_speedup(r$het_old, r$het_new)))
}

cat("\nvalue agreement, dense against O(d) (relative):\n")
for (i in seq_len(nrow(tab))) {
  cat(sprintf("  d = %-5d ar1 %.3e   hetar1 %.3e\n", tab$d[i],
              tab$ar1_reldiff[i], tab$het_reldiff[i]))
}

# The other axis: many short blocks, which is the repeated-measures
# shape the dense route was written for. The O(d) form must not lose
# there, because that is where most fits live.
cat("\nmany short blocks (the repeated-measures shape):\n")
tab2 <- run(c(4L, 10L), n_levels = 500L, reps = reps)
cat("| d | levels | ar1 dense | ar1 O(d) |\n| --- | --- | --- | --- |\n")
for (i in seq_len(nrow(tab2))) {
  cat(sprintf("| %d | %d | %s | %s |\n", tab2$d[i], tab2$n_levels[i],
              fmt(tab2$ar1_old[i]), fmt(tab2$ar1_new[i])))
}
