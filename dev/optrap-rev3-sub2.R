# Reviewer of lane optrap, final check: the one design where the
# subsample's unit moved (dev/optrap-rev3-sub.R, rows sorted by x with
# sd exp(1.5 x)). Is it the sorting or the subsample's size? Arms: rows
# sorted and unsorted, subsample of 100 and of 1000 rows (n = 1e3 and
# 1e4, stride 10), against the unit from every row; then a fit of the
# worst seed with the subsample's unit forced, against the full unit.
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
pu <- frmtmb:::predictor_units
unit_at <- function(d, cells) {
  old <- options(frmtmb.unit_ls_cells = cells)
  on.exit(options(old))
  u <- pu(frm(y ~ x, data = d, dry_run = "frame"))
  if (length(u)) u[[1]] else 1
}
X <- list()
for (n in c(1e3, 1e4)) for (sorted in c(TRUE, FALSE)) for (s in 1:20) {
  set.seed(s)
  x <- rnorm(n)
  if (sorted) x <- sort(x)
  d <- data.frame(x = x)
  d$y <- 1 + 0.5 * d$x + rnorm(n, 0, exp(1.5 * d$x))
  X[[length(X) + 1L]] <- data.frame(n = n, sorted = sorted, seed = s,
    full = unit_at(d, Inf), sub = unit_at(d, 2 * n / 10))
}
X <- do.call(rbind, X)
X$l2 <- abs(log2(X$sub / X$full))
for (n in unique(X$n)) for (so in c(TRUE, FALSE)) {
  x <- X[X$n == n & X$sorted == so, ]
  cat(sprintf("n %g sorted %s | subsample rows %g | |log2 sub/full| > 1: %d of 20, max %g\n",
              n, so, n / 10, sum(x$l2 > 1), max(x$l2)))
}
# the fit at the worst sorted seed of n = 1e3, with each unit
w <- X[X$n == 1e3 & X$sorted, ]
s <- w$seed[which.max(w$l2)]
set.seed(s)
d <- data.frame(x = sort(rnorm(1e3)))
d$y <- 1 + 0.5 * d$x + rnorm(1e3, 0, exp(1.5 * d$x))
ll <- as.numeric(logLik(lm(y ~ x, data = d)))
for (cells in c(Inf, 200)) {
  options(frmtmb.unit_ls_cells = cells)
  f <- suppressWarnings(frm(y ~ x, data = d))
  cat(sprintf("seed %d cells %g unit %s | code %d | logLik - lm %.3g\n", s,
              cells, paste(f$par_units %||% 1, collapse = ","),
              f$opt$convergence, as.numeric(logLik(f)) - ll))
}
