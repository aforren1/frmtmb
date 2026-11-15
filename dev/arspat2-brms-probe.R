# Which second residual term brms 2.23.0 takes beside ar(cov = FALSE)
# (backlog 5256 = 3994): make_stancode() only, no compile.
.libPaths(c("C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
set.seed(1)
k <- 4
g <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(g)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- A / rowSums(A)
M <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
d <- data.frame(y = rnorm(n), x = rnorm(n), t = rep(1:4, 4),
                gr = factor(rep(1:4, each = 4)))
cases <- list(
  sar_lag = y ~ x + ar(t, gr, cov = FALSE) + sar(W, type = "lag"),
  sar_err = y ~ x + ar(t, gr, cov = FALSE) + sar(W, type = "error"),
  fcor = y ~ x + ar(t, gr, cov = FALSE) + fcor(M),
  cosy = y ~ x + ar(t, gr, cov = FALSE) + cosy(t, gr),
  unstr = y ~ x + ar(t, gr, cov = FALSE) + unstr(t, gr),
  ma_cov = y ~ x + ar(t, gr, cov = FALSE) + ma(t, gr, cov = TRUE),
  sar_ma = y ~ x + ma(t, gr, cov = FALSE) + sar(W)
)
for (nm in names(cases)) {
  for (fam in c("gaussian", "student")) {
    cat("\n=== ", nm, " ", fam, "\n", sep = "")
    r <- tryCatch({
      sc <- make_stancode(cases[[nm]], data = d, family = fam,
                          data2 = list(W = W, M = M))
      lines <- strsplit(sc, "\n")[[1]]
      keep <- grep("Err|ar\\[|_lpdf\\(Y|target \\+=.*Y|mu\\[n\\] \\+=",
                   lines, value = TRUE)
      cat(paste(trimws(keep), collapse = "\n"), "\n")
      "fits"
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(r, "\n")
  }
}
