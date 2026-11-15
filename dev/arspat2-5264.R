# Backlog 5264 (= 3992): confint_varcorr() of sar()'s rho. Over 20 seeds of
# a 5 x 5 rook lattice (row-standardized, so rho's upper end is 1) with
# data drawn at rho = 0.9, lag type: how often the interval passes an end
# of rho's interval, and the summary()'s cor_pars row. Base build when
# ARSPAT2_LIB="".
lib <- Sys.getenv("ARSPAT2_LIB", "C:/Users/adf44/source/r/wt-arspat2-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
cat("frmtmb from", dirname(dirname(system.file(package = "frmtmb"))), "\n")
k <- 5
g <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(g)
Aj <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- Aj / rowSums(Aj)
past <- 0
for (seed in 1:20) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  d$y <- as.vector(solve(diag(n) - 0.9 * W, 1 + 0.5 * d$x + rnorm(n)))
  fit <- suppressWarnings(suppressMessages(
    frm(y ~ x + sar(W), data = d, data2 = list(W = W))))
  ac <- fit$frame$autocor[[1L]]
  cv <- suppressWarnings(confint_varcorr(fit))
  s <- summary(fit)$cor_pars
  out <- cv$upr >= ac$aux$hi || cv$lwr <= ac$aux$lo
  past <- past + out
  cat(sprintf("seed %2d: rho %.4f  interval (%.4f, %.4f)  ends (%.4f, %.4f)%s  summary (%.4f, %.4f)\n",
              seed, cv$estimate, cv$lwr, cv$upr, ac$aux$lo, ac$aux$hi,
              if (out) "  PAST AN END" else "", s[1, 3], s[1, 4]))
}
cat(sprintf("intervals past an end: %d of 20\n", past))
