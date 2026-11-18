# Lane drawdir: the three fits of the new tests in
# tests/testthat/test-predict-draw-check.R, with the data perturbed by a
# few ulps (x for the poisson mixture, whose counts must stay whole; y
# for the gaussian GLMMs) times 1 + 2e-15 k, k = 0 to 20. Per rounding:
# the summed check's excess and ratio, the per-direction reading, its
# evaluation count and the verdict, so that no assertion is seen to sit
# on a rounding-decided branch.
# Usage: Rscript dev/drawdir-ulp.R <lib or "base"> <tag>
args <- commandArgs(TRUE)
lib <- if (args[1] != "base") args[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), " tag", args[2],
    " tol", ns$draw_dir_tol, "\n")
q <- function(expr) suppressWarnings(suppressMessages(expr))
mixk <- function(k) {
  set.seed(25)
  x <- stats::rnorm(300)
  invisible(stats::runif(300))
  d <- data.frame(x = x, y = stats::rpois(300, exp(1 + 0.3 * x)))
  d$x <- d$x * (1 + 2e-15 * k)
  set.seed(5025)
  for (j in 1:20) d[[paste0("z", j)]] <- stats::rnorm(300)
  fo <- stats::as.formula(paste("y ~", paste(c("x", paste0("z", 1:20)),
                                             collapse = " + ")))
  q(frm(bf(fo), family = mixture(poisson(), poisson()), data = d))
}
slope <- function(k) {
  set.seed(2014)
  d <- data.frame(g = factor(rep(1:15, each = 8)), x = stats::rnorm(120))
  ug <- stats::rnorm(15, 0, 0.8)
  sg <- stats::rnorm(15, 0, 0.1)
  d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + stats::rnorm(120)
  d$y <- d$y * (1 + 2e-15 * k)
  q(frm(bf(y ~ x + (1 + x | g)), data = d))
}
grp <- function(k) {
  set.seed(1)
  d <- data.frame(g = factor(rep(1:30, each = 4)), x = stats::rnorm(120))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(30, 0, 0.3)[d$g] + stats::rnorm(120)
  d$y <- d$y * (1 + 2e-15 * k)
  q(frm(bf(y ~ x + (1 | g)), data = d))
}
for (nm in c("mixk", "slope", "grp")) {
  f <- get(nm)
  for (k in 0:20) {
    fit <- f(k)
    set.seed(1)
    q(predict(fit, ndraws = 20, summary = FALSE))
    chk <- fit$cache[["draw_check"]]
    cat(sprintf("%-6s k %2d np %2d fire %-5s excess %9.5f rel %8.5f | dir %10.5f (%s) n %3d | logLik %.10f\n",
                nm, k, length(fit$opt$par), chk$fire, chk$excess, chk$rel,
                chk$dir %||% NA, chk$dir_par %||% "-", chk$dir_n %||% 0L,
                -fit$opt$objective))
  }
}
