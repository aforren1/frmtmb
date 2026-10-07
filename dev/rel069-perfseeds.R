# The ciharden review's script on the 0.69.0 release library (rellib-r7)
# and the release tree's test-perf.R (the merged-tree rerun it asks).
# Reviewer: how much headroom does test-perf.R's node bound of 100 have
# when the data are drawn with other seeds? The test pins seeds 71 and
# 72 (GLM) and 7 and 8 (GLMM); this draws 30 other pairs of each.
#   Rscript dev/rel069-perfseeds.R <glm|glmm> <first seed> <pairs>
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c("C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
env <- new.env(parent = asNamespace("frmtmb"))
ex <- parse("C:/Users/adf44/source/r/frmtmb-wt-release/tests/testthat/test-perf.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, env)
}
perf_census <- env$perf_census
kind <- a[1]
s0 <- as.integer(a[2])
k <- as.integer(a[3])
make_glm <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rpois(n, exp(0.3 + 0.4 * d$x))
  d
}
make_glmm <- function(n, seed) {
  set.seed(seed)
  n_g <- 500L
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(sample.int(n_g, n, replace = TRUE),
                             levels = seq_len(n_g)))
  d$y <- stats::rpois(n, exp(0.2 + 0.3 * d$x +
                               stats::rnorm(n_g, 0, 0.4)[d$g]))
  d
}
if (kind == "glm") {
  f <- bf(y ~ x) + poisson()
  mk <- make_glm
  rnd <- NULL
} else {
  f <- bf(y ~ x + (1 | g)) + poisson()
  mk <- make_glmm
  rnd <- "b"
}
invisible(perf_census(f, mk(1000L, 1), random = rnd))
for (i in seq_len(k)) {
  ss <- s0 + 2L * (i - 1L)
  ds <- mk(1000L, ss)
  dl <- mk(100000L, ss + 1L)
  sm <- perf_census(f, ds, random = rnd)
  lg <- perf_census(f, dl, random = rnd)
  cat(sprintf(paste0("%s seeds %d,%d nodes %d %d ratio %.3f bytes ratio %.2f",
                     " zeros_small %.4f zeros_large %.4f ones_small %.4f",
                     " ones_large %.4f\n"),
              kind, ss, ss + 1L, sm[["nodes"]], lg[["nodes"]],
              lg[["nodes"]] / sm[["nodes"]], lg[["bytes"]] / sm[["bytes"]],
              mean(ds$y == 0), mean(dl$y == 0), mean(ds$y == 1),
              mean(dl$y == 1)))
}
