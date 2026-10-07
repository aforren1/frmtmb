# Reviewer: does test-perf.R's own perf_census() catch an
# observation-length elementwise sub-assignment loop placed INSIDE the
# objective that build_objective() returns (dev/rtmb-pitfalls.md item
# 14)? The lane measured the loop on a separate RTMB objective.
.libPaths(c("C:/Users/adf44/source/r/wt-ciharden-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
orig <- get("build_objective", ns)
nobs_now <- 0L
mutant <- function(fr, ...) {
  nll <- orig(fr, ...)
  n <- nobs_now
  function(pars) {
    val <- nll(pars)
    "[<-" <- RTMB::ADoverload("[<-")
    b <- pars[[1L]][1L]
    m <- b * rep(1, n)
    for (i in seq_len(n)) m[i] <- b * 2
    val + 1e-12 * sum(m)
  }
}
env <- new.env(parent = ns)
env$build_objective <- mutant
ex <- parse("C:/Users/adf44/source/r/frmtmb-wt-ciharden/tests/testthat/test-perf.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, env)
}
make <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rpois(n, exp(0.3 + 0.4 * d$x))
  d
}
f <- bf(y ~ x) + poisson()
nobs_now <- 1000L
invisible(env$perf_census(f, make(1000L, 70)))
s <- env$perf_census(f, make(1000L, 71))
for (nl in c(10000L, 100000L)) {
  nobs_now <- nl
  t0 <- proc.time()[["elapsed"]]
  l <- env$perf_census(f, make(nl, 72))
  cat(sprintf(paste0("mutant n %d/%d: nodes %d -> %d ratio %.2f (bound %d);",
                     " bytes %.4g -> %.4g ratio %.1f (bound %d); %.0f s\n"),
              1000L, nl, s[["nodes"]], l[["nodes"]],
              l[["nodes"]] / s[["nodes"]], nl / 1000L, s[["bytes"]],
              l[["bytes"]], l[["bytes"]] / s[["bytes"]], 2L * nl / 1000L,
              proc.time()[["elapsed"]] - t0))
}
