# Lane sampleprior: where the per-draw cost of the derived sd_/cor_
# columns goes (Rprof over draws_cov_values() on 2000 synthetic draws of
# the dev/sampleprior-names.R model; data seed 11).
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
set.seed(11)
n <- 160
dd <- data.frame(x = stats::rnorm(n), g = factor(rep(1:16, 10)),
                 h = factor(rep(1:8, each = 20)))
u <- cbind(stats::rnorm(16, 0, 0.7), stats::rnorm(16, 0, 0.4))
dd$y <- stats::rnorm(n, 1 + 0.5 * dd$x + u[dd$g, 1] + u[dd$g, 2] * dd$x +
                       stats::rnorm(8, 0, 0.5)[dd$h], 1)
fit <- frm(bf(y ~ x + (1 + x | g) + (1 | h)), family = gaussian(),
           data = dd)
nc <- frmtmb.sample:::draws_natural_cols(fit)
dv <- nc$derived
lab <- unlist(dv$cols)
m <- matrix(rnorm(2000 * length(lab), 0, 0.3), 2000,
            dimnames = list(NULL, lab))
# the slow path, every draw through hyp_env_vals(), as the first lane
# build ran it; profiled
slow <- dv
slow$fast <- NULL
f <- tempfile()
Rprof(f, interval = 0.002)
t0 <- proc.time()[["elapsed"]]
v <- frmtmb.sample:::draws_cov_values(m, fit, slow)
el <- proc.time()[["elapsed"]] - t0
Rprof(NULL)
cat(sprintf("slow path, %d draws, %d columns: %.2f s\n", nrow(v), ncol(v),
            el))
s <- summaryRprof(f)$by.total
print(utils::head(s[order(-s$total.time), ], 12))
# the two paths interleaved in one process, three rounds each, minimum
# reported, and their largest relative difference
tm <- function(d) {
  t0 <- proc.time()[["elapsed"]]
  frmtmb.sample:::draws_cov_values(m, fit, d)
  proc.time()[["elapsed"]] - t0
}
el <- sapply(1:3, function(k) c(fast = tm(dv), slow = tm(slow)))
cat(sprintf("2000 draws: fast path min %.2f s, slow path min %.2f s\n",
            min(el["fast", ]), min(el["slow", ])))
vf <- frmtmb.sample:::draws_cov_values(m, fit, dv)
cat("max |fast - slow| / max |slow|:", max(abs(vf - v)) / max(abs(v)), "\n")
