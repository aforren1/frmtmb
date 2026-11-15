# Backlog 4483: what default_prior() and frm_sample()'s prior stack give a
# car() block's rho (escar, bym2), on the route of the sampler.
lib <- Sys.getenv("ARSPAT2_LIB", "C:/Users/adf44/source/r/wt-arspat2-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(frmtmb.sample)
})
k <- 4
g <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(g)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
rownames(A) <- colnames(A) <- seq_len(n)
set.seed(1)
d <- data.frame(loc = factor(rep(seq_len(n), 3)), x = rnorm(3 * n))
d$y <- 1 + 0.5 * d$x + rnorm(3 * n)
for (ty in c("escar", "bym2")) {
  fit <- suppressWarnings(frm(y ~ x + car(A, gr = loc, type = ty), data = d,
                              data2 = list(A = A)))
  cat("\n=== ", ty, "\n")
  print(frmtmb:::outer_par_labels(fit))
  dp <- default_prior(fit, route = "sample")
  print(dp)
}
