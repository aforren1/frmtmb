# Lane ci701: the body of test-ordinal-mixture.R "simulate() draws from
# the fitted category distribution", as revised, run on data moved by a
# few ulps (seed 2000 + k, as dev/ci701-omx-probe.R), so that both sides
# of the "singular convergence (7)" branch are seen to pass. Reports per
# fit the optimizer code, every warning that would escape the test's
# allow_warnings() lists, and the frequency assertion.
#
#   Rscript dev/ci701-omx-body.R <lib|base> <k0> <k1>
args <- commandArgs(TRUE)
lib <- args[1]; k0 <- as.integer(args[2]); k1 <- as.integer(args[3])
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
omx_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) +
    stats::rlogis(n) / exp(0.3 * d$z * cls)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < stats::plogis(-1.2 + 0.5 * d$z), 0L, d$y)
  d
}
escaped <- character(0)
seen <- character(0)
allow <- function(expr, allowed) {
  withCallingHandlers(expr, warning = function(w) {
    m <- conditionMessage(w)
    if (any(vapply(allowed, function(a) grepl(a, m, fixed = TRUE), NA))) {
      seen <<- c(seen, m)
    } else {
      escaped <<- c(escaped, m)
    }
    invokeRestart("muffleWarning")
  })
}
d0 <- omx_data(20261011, n = 300)
for (k in k0:k1) {
  d <- d0
  if (k > 0) {
    set.seed(2000 + k)
    d$x <- d$x * (1 + sample(0:4, nrow(d), TRUE) * 2^-52)
  }
  escaped <- character(0)
  seen <- character(0)
  fit <- allow(
    frm(bf(y | thres(gr = g) ~ x), family = mixture(cumulative(), sratio()),
        data = d),
    c("degenerate boundary", "did not report convergence"))
  req <- any(grepl("degenerate boundary", seen, fixed = TRUE))
  P <- allow(fitted(fit), "move along a direction")[, "Estimate", ]
  sims <- simulate(fit, nsim = 400, seed = 7)
  codes <- vapply(sims, as.integer, integer(nrow(d)))
  freq <- vapply(1:4, function(j) mean(codes == j), 0)
  pbar <- colMeans(P)
  se <- sqrt(colSums(P * (1 - P))) / (nrow(d) * sqrt(400))
  ok <- all(abs(freq - pbar) < 5 * se)
  cat(sprintf("RES k=%d code=%d required=%s freq_ok=%s max_z=%.3f escaped=%d\n",
              k, fit$opt$convergence, req, ok, max(abs(freq - pbar) / se),
              length(escaped)))
  for (m in unique(escaped)) cat("   ESC:", substr(m, 1, 200), "\n")
}
