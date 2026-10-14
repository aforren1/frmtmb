# Reviewer of lane optstart: cost of the least-squares start at n = 1e6.
#   Rscript dev/optstart-rev-cost.R <lib> <label>
# A gaussian GLM (no random effects, so the fit itself is cheap) and the
# same with (1 | g). Elapsed seconds, minimum of 3 rounds: the whole
# frm() call, and make_start() alone (frame from dry_run = "frame").
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
set.seed(1)
n <- 1e6
d <- data.frame(matrix(rnorm(n * 5), n, 5), g = factor(sample(1:200, n, TRUE)))
d$y <- 1 + 50 * d$X1 - 3 * d$X2 + rnorm(200)[d$g] + rnorm(n)
tm <- function(fn) min(vapply(1:3, function(i) system.time(fn())[["elapsed"]], 0))
for (fo in list(y ~ X1 + X2 + X3 + X4 + X5,
                y ~ X1 + X2 + X3 + X4 + X5 + (1 | g))) {
  fr <- frm(fo, data = d, dry_run = "frame")
  t_ms <- tm(function() ns$make_start(fr, NULL))
  t_fit <- tm(function() f <<- suppressWarnings(frm(fo, data = d)))
  cat(sprintf("%s: make_start %.2f s, frm %.2f s, logLik %.6f code %d\n",
              deparse(fo), t_ms, t_fit, as.numeric(logLik(f)),
              f$opt$convergence))
}
