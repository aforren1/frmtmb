# Short review of lane selimit2's fix round, F3: the cost of
# jc_limit_extrap() moving every predictor along each limit direction, on
# a large fit with many limits. Poisson y ~ f + x, 60 levels of 1000 rows,
# 45 of them all zero (45 ends), and a ZIP with zi ~ x + f whose zi runs
# off; fitted() on every row, timed with peak memory.
#   Rscript selimit2-short-cost.R <base|lane>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "lib", dirname(find.package("frmtmb")), "\n")
quiet <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    cat("   W:", substr(conditionMessage(w), 1, 140), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
}
timed <- function(label, fit) {
  sl <- frmtmb:::sdr_of(fit)$se_lost
  cat(label, "| lost", sum(sl == "limit"), "limit of", length(sl), "\n")
  invisible(gc(reset = TRUE))
  t0 <- proc.time()[["elapsed"]]
  fp <- quiet(stats::fitted(fit))
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("   fitted() on %d rows: %.2f s, max used %.0f MB, NaN SE %d\n",
              nrow(fp), el, sum(gc()[, 6]), sum(is.nan(fp[, "Est.Error"]))))
}
set.seed(12)
nl <- 60
d <- data.frame(f = factor(rep(sprintf("L%02d", 1:nl), each = 1000)),
                x = stats::rnorm(nl * 1000))
lev <- as.integer(d$f)
d$y <- ifelse(lev > 15, 0L, stats::rpois(nrow(d), exp(0.5 + 0.3 * d$x)))
fit <- quiet(frm(y ~ f + x, data = d, family = poisson()))
timed("poisson y ~ f + x, 45 levels all zero", fit)
set.seed(13)
d2 <- data.frame(x = stats::rnorm(50000),
                 f = factor(sample(letters[1:10], 50000, TRUE)))
d2$y <- stats::rpois(50000, exp(0.8 + 0.3 * d2$x))
fit2 <- quiet(frm(bf(y ~ x, zi ~ x + f), data = d2,
                  family = zero_inflated_poisson()))
timed("zip zi ~ x + f, no excess zeros", fit2)
