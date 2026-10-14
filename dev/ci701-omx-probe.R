# Lane ci701: the fit of test-ordinal-mixture.R:359 ("simulate() draws
# from the fitted category distribution"), omx_data(20261011, n = 300),
# refitted with x and z moved by a few ulps (seed 2000 + k), to see how
# often nlminb ends with "singular convergence (7)" and what the fit
# looks like on each side.
#
#   Rscript dev/ci701-omx-probe.R <lib|base> <k0> <k1>
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
d0 <- omx_data(20261011, n = 300)
for (k in k0:k1) {
  d <- d0
  if (k > 0) {
    set.seed(2000 + k)
    d$x <- d$x * (1 + sample(0:4, nrow(d), TRUE) * 2^-52)
  }
  w <- character(0)
  fit <- withCallingHandlers(
    frm(bf(y | thres(gr = g) ~ x), family = mixture(cumulative(), sratio()),
        data = d),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
    })
  gr <- drop(fit$obj$gr(fit$opt$par))
  cat(sprintf("RES k=%d code=%d msg=%s ll=%.10f gmax=%.3g nwarn=%d\n", k,
              fit$opt$convergence, gsub(" ", "_", fit$opt$message),
              as.numeric(logLik(fit)), max(abs(gr * (fit$par_units %||% 1))),
              length(w)))
  for (m in unique(w)) cat("   W:", substr(m, 1, 160), "\n")
}
