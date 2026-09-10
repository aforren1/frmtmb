# rev-gddm, second pass: the residual band under the new tolerance
#
#     tol = 1e-8 * max(|a|, |b|)  +  1e-11 * max|finite column|
#
# The second term still scales with the column maximum, so a column of
# dynamic range R hides a within-condition difference of 1e-11 * R.
# Two questions, both measured:
#
#   1. how wide is that band in the units of the column, on the same
#      delay construction that broke the first form;
#   2. can a design carry a MEANINGFUL difference that small? A column
#      is only estimable at all if its total range clears frmtmb's rank
#      test, so the band is bounded from the other side too. That
#      threshold is measured here rather than assumed.
#
# Seed 909, 120 rows, grid dt = 0.05, ny = 51, t_max = 2.
# Arm from GDDM_LIB; default is this review's own install.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
`%||%` <- function(x, y) if (is.null(x)) y else x
cat("arm:", lib, " frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n\n")

ctl <- gddm_control(t_max = 2, dt = 0.05, ny = 51L)
set.seed(909)
n <- 120L
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$blk <- rep(1:4, each = 30L)
d$cond <- d$blk

form <- bf(rt | vint(upper, cond) ~ delay, bs ~ 1, ndt ~ 1, bias = 0.5)
run <- function(dat) {
  w <- NULL
  r <- withCallingHandlers(
    tryCatch({
      frm(form, family = gddm(control = ctl), data = dat,
          dry_run = "frame")
      ""
    }, error = function(e) conditionMessage(e)),
    warning = function(z) { w <<- c(w, conditionMessage(z))
                            invokeRestart("muffleWarning") })
  list(msg = r, warn = w)
}
report <- function(label, delay, want) {
  dat <- d; dat$delay <- delay
  mx <- max(abs(delay))
  sp <- max(tapply(delay, dat$cond, function(z) max(z) - min(z)))
  floor_t <- 1e-11 * mx
  r <- run(dat)
  got <- if (nzchar(r$msg)) "REFUSED" else "accepted"
  cat(sprintf("%s %-38s max %.3g  floor %.3g  spread %.3g  %s\n",
              if (got == want) "   " else ">>>", label, mx, floor_t, sp,
              got))
  if (!is.null(r$warn))
    cat("       warning: ", substr(r$warn[[1L]], 1, 78), "\n", sep = "")
  invisible(got)
}

cat("== 1. the review's delay column under the new form\n")
short <- c(1, 2, 3)
report("1/2/3 s beside ten years",
       ifelse(d$blk <= 2L, sample(short, n, replace = TRUE), 3.15e8),
       "REFUSED")
report("the same column, centered",
       {v <- ifelse(d$blk <= 2L, sample(short, n, replace = TRUE),
                    3.15e8); v - mean(v)}, "REFUSED")

cat("\n== 2. the residual band, priced in the column's own units\n")
cat("   the floor on a ten-year column is 1e-11 * 3.15e8 = 3.15e-3 s\n")
for (step in c(1e-1, 1e-2, 3.2e-3, 3.1e-3, 1e-3, 1e-6)) {
  v <- ifelse(d$blk <= 2L,
              1 + step * rep(c(0, 1), length.out = n), 3.15e8)
  report(sprintf("within-condition step of %.1e s", step),
         v, if (step > 3.15e-3) "REFUSED" else "accepted")
}

cat("\n== 3. the other side: what range is estimable at all?\n")
cat("   a column is dropped as rank deficient below some relative\n")
cat("   total range. Measured on the same model, base 1e9:\n")
for (rel in 10^-(4:12)) {
  v <- 1e9 * (1 + rel * (d$blk - 1L) / 3)     # constant within blk
  dat <- d; dat$delay <- v
  r <- run(dat)
  dropped <- any(grepl("rank deficient", r$warn %||% "", fixed = TRUE))
  cat(sprintf("   relative total range %.0e: %s%s\n", rel,
              if (nzchar(r$msg)) "REFUSED" else "accepted",
              if (dropped) ", COLUMN DROPPED as rank deficient" else
                ", column kept"))
}

cat("\n== 4. a column with a sentinel, which is how 1e11 of dynamic\n")
cat("      range can appear without 1e11 of meaningful range\n")
for (sent in c(1e6, 1e9, 1e11, 1e12)) {
  v <- ifelse(d$blk <= 2L, rep(c(1, 3), length.out = n), sent)
  report(sprintf("real values 1 and 3, sentinel %.0e", sent),
         v, if (sent < 1e11 * 2 / 1e-11 / 1e11) "REFUSED" else "?")
}

cat("\n== 5. false-alarm risk under the TIGHTER form: poly() near zero\n")
cat("   the per-pair term nearly vanishes at an entry near zero, so\n")
cat("   only the floor carries it. Worst |row - first| against its own\n")
cat("   tolerance, over degree and number of levels:\n")
for (nl in c(4L, 6L, 8L, 12L)) {
  for (deg in 1:4) {
    if (deg >= nl) next
    x <- rep(seq(0, 1, length.out = nl), each = n / nl)
    P <- poly(x, deg)
    gi <- match(x, sort(unique(x)))
    f1 <- match(seq_len(nl), gi)
    worst <- 0
    for (k in seq_len(ncol(P))) {
      a <- P[, k]; b <- a[f1[gi]]
      tol <- 1e-8 * pmax(abs(a), abs(b)) + 1e-11 * max(abs(a))
      worst <- max(worst, max(abs(a - b) / tol))
    }
    cat(sprintf("   %2d levels, poly degree %d: worst/tol = %.3e\n",
                nl, deg, worst))
  }
}
