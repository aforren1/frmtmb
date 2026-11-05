# The warm start's own cost, alone: setup, grid, descent, against the whole frm(), with
# an Rprof of the marginal at 500. 60, 200, 500 positions, and 1000 rows on 500 positions.
# Same design as dev/ncgp-perf.R, seed 20261103. Run with nothing else of mine on the box.
#   Rscript dev/ncgp-warmcost.R
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
tick <- function(e) { t0 <- proc.time()[[3]]; force(e); proc.time()[[3]] - t0 }
for (n in c(60L, 200L, 500L, 1000L)) {
  set.seed(20261103)
  x <- if (n == 1000L) sample(sort(stats::runif(500, 0, 10)), n, TRUE) else sort(stats::runif(n, 0, 10))
  d <- data.frame(x = x, y = sin(x) + 0.3 * stats::rnorm(n))
  fr <- suppressWarnings(frm(y ~ gp(x), d))$frame
  for (rep in 1:3) {
    ts <- tick(s <- ns$gp_marginal_setup(fr))
    tm <- tick(m <- ns$gp_marginal_opt(s, 3L))
    t1 <- tick(ns$gp_marginal_opt(s, 1L))
    tf <- tick(f <- suppressWarnings(frm(y ~ gp(x), d)))
    cat(sprintf("WC n %d u %d rep %d setup %.3fs marginal keep 3 %.2fs keep 1 %.2fs | frm %.2fs conv %d iter %d\n",
                n, s$u, rep, ts, tm, t1, tf, f$opt$convergence, f$opt$iterations))
  }
}
set.seed(20261103)
x <- sort(stats::runif(500, 0, 10))
d <- data.frame(x = x, y = sin(x) + 0.3 * stats::rnorm(500))
fr <- suppressWarnings(frm(y ~ gp(x), d))$frame
s <- ns$gp_marginal_setup(fr)
pf <- tempfile()
Rprof(pf, interval = 0.02)
invisible(ns$gp_marginal_opt(s, 3L))
Rprof(NULL)
sm <- summaryRprof(pf)$by.total
print(utils::head(sm[, c("total.time", "self.time")], 16))
