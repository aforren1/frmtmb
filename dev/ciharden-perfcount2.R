# The counts test-perf.R asserts on, printed: tape operations and bytes
# allocated at n = 1000 and 100000, for the GLM and the GLMM, with the
# test file's own perf_census() (parsed out of it, so this measures what
# the test measures). Run under each BLAS to see the counts do not move.
# Usage: Rscript dev/ciharden-perfcount2.R <lib or base>
a <- commandArgs(TRUE)
.libPaths(c(if (!identical(a[1], "base")) a[1],
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
root <- "C:/Users/adf44/source/r/frmtmb-wt-ciharden"
env <- new.env(parent = asNamespace("frmtmb"))
for (e in parse(file.path(root, "tests/testthat/test-perf.R"))) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      identical(e[[2]], as.name("perf_census"))) eval(e, env)
}
cat("R:", R.home(), " threads:", Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
glm_d <- function(n, seed) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(n))
  d$y <- stats::rpois(n, exp(0.3 + 0.4 * d$x))
  d
}
glmm_d <- function(n, seed) {
  set.seed(seed)
  n_g <- 500L
  d <- data.frame(x = stats::rnorm(n),
                  g = factor(sample.int(n_g, n, replace = TRUE),
                             levels = seq_len(n_g)))
  d$y <- stats::rpois(n, exp(0.2 + 0.3 * d$x +
                               stats::rnorm(n_g, 0, 0.4)[d$g]))
  d
}
pc <- env$perf_census
f1 <- bf(y ~ x) + poisson()
invisible(pc(f1, glm_d(1000L, 70)))
s1 <- pc(f1, glm_d(1000L, 71))
l1 <- pc(f1, glm_d(100000L, 72))
f2 <- bf(y ~ x + (1 | g)) + poisson()
invisible(pc(f2, glmm_d(1000L, 6), random = "b"))
s2 <- pc(f2, glmm_d(1000L, 7), random = "b")
l2 <- pc(f2, glmm_d(100000L, 8), random = "b")
out <- rbind(glm_small = s1, glm_large = l1, glmm_small = s2,
             glmm_large = l2)
print(out)
cat(sprintf("GLM  ratio: nodes %.4f (bound 100), bytes %.4f (bound 200)\n",
            l1[["nodes"]] / s1[["nodes"]], l1[["bytes"]] / s1[["bytes"]]))
cat(sprintf("GLMM ratio: nodes %.4f (bound 100), bytes %.4f (bound 200)\n",
            l2[["nodes"]] / s2[["nodes"]], l2[["bytes"]] / s2[["bytes"]]))
