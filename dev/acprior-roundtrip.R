# Does a default_prior() table passed back as prior = fit? A sar() row
# carries lb/ub at the ends of its interval. Arg 1: "before"/"after".
# Log: dev/acprior-out/roundtrip-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- c(if (identical(arm, "after")) "C:/Users/adf44/source/r/wt-acprior-lib",
         "C:/Users/adf44/source/r/rellib-r13",
         "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(lib)
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
set.seed(2)
k <- 5
g <- expand.grid(r = seq_len(k), c = seq_len(k))
A <- outer(seq_len(k^2), seq_len(k^2), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- A / rowSums(A)
d <- data.frame(x = rnorm(k^2))
d$y <- as.vector(solve(diag(k^2) - 0.4 * W, 1 + 0.5 * d$x + rnorm(k^2)))
d$t <- seq_len(k^2)
d$gg <- rep(1:5, each = 5)
rt <- function(lab, f, d2 = list()) {
  cat("\n---", lab, "\n")
  tab <- default_prior(f, data = d, data2 = d2)
  print(as.data.frame(tab)[, c("prior", "class", "coef", "lb", "ub")])
  r <- tryCatch({
    fit <- frm(f, data = d, data2 = d2, prior = tab)
    fit0 <- frm(f, data = d, data2 = d2)
    cat("  round trip fits; logLik", format(as.numeric(logLik(fit)),
                                          digits = 12),
        "plain", format(as.numeric(logLik(fit0)), digits = 12), "\n")
  }, error = function(e) cat("  ERROR:", conditionMessage(e), "\n"))
}
rt("sar lag", y ~ x + sar(W), list(W = W))
rt("ar cov", y ~ x + ar(time = t, gr = gg, cov = TRUE))
rt("cosy", y ~ x + cosy(time = t, gr = gg))
