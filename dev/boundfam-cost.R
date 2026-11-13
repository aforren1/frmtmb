# Lane boundfam: what se_edge_limit() costs, against the fit it runs in.
#   Rscript dev/boundfam-cost.R <lane lib>
# Each fit: the whole frm() call's elapsed time, the time spent inside
# se_edge_limit() (traced), and how many objective and gradient
# evaluations it made (counted by wrapping fit$obj$fn and $gr inside
# the trace). Three replicates each, in one process.
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib", dirname(find.package("frmtmb")), "\n")
acc <- new.env()
acc$t <- 0
acc$fn <- 0L
acc$gr <- 0L
acc$calls <- 0L
real <- ns$se_edge_limit
counted <- function(fit, H, skip = logical(NROW(H))) {
  obj <- fit$obj
  f0 <- obj$fn
  g0 <- obj$gr
  obj$fn <- function(...) { acc$fn <- acc$fn + 1L; f0(...) }
  obj$gr <- function(...) { acc$gr <- acc$gr + 1L; g0(...) }
  on.exit({ obj$fn <- f0; obj$gr <- g0 }, add = TRUE)
  fit$obj <- obj
  t0 <- proc.time()[[3]]
  out <- real(fit, H, skip)
  acc$t <- acc$t + proc.time()[[3]] - t0
  acc$calls <- acc$calls + 1L
  out
}
utils::assignInNamespace("se_edge_limit", counted, "frmtmb")
cases <- list(
  nb_limit_re = function() {
    set.seed(5)
    n <- 2e4
    d <- data.frame(x = rnorm(n), g = factor(sample(200, n, TRUE)))
    d$y <- rpois(n, exp(0.5 + 0.3 * d$x + rnorm(200, 0, 0.3)[d$g]))
    frm(y ~ x + (1 | g), family = negbinomial(), data = d)
  },
  nb_max_re = function() {
    set.seed(5)
    n <- 2e4
    d <- data.frame(x = rnorm(n), g = factor(sample(200, n, TRUE)))
    mu <- exp(0.5 + 0.3 * d$x + rnorm(200, 0, 0.3)[d$g])
    d$y <- MASS::rnegbin(n, mu, theta = 5)
    frm(y ~ x + (1 | g), family = negbinomial(), data = d)
  },
  gauss_re = function() {
    set.seed(5)
    n <- 2e4
    d <- data.frame(x = rnorm(n), g = factor(sample(200, n, TRUE)))
    d$y <- 1 + 0.3 * d$x + rnorm(200, 0, 0.5)[d$g] + rnorm(n)
    frm(y ~ x + (1 | g), data = d)
  })
for (nm in names(cases)) {
  for (r in 1:3) {
    acc$t <- 0; acc$fn <- 0L; acc$gr <- 0L; acc$calls <- 0L
    t0 <- proc.time()[[3]]
    f <- suppressMessages(suppressWarnings(cases[[nm]]()))
    se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
    tt <- proc.time()[[3]] - t0
    cat(sprintf("%-12s rep %d: fit+vcov %.2f s, se_edge_limit %.2f s in %d calls, fn %d, gr %d, lost %s\n",
                nm, r, tt, acc$t, acc$calls, acc$fn, acc$gr,
                paste(names(ns$sdr_of(f)$se_lost), collapse = ",")))
  }
}
