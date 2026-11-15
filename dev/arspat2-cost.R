# Cost of the student() margin route (backlog 5289): one whiten() call per
# missing value on the tape. Seconds per frm() fit, minimum of 3 rounds,
# for fcor() and sar(type = "error") with (1 | h) and m missing values,
# beside the same model fit on the observed rows (fcor(M_oo), the cheap
# exact form for fcor). Base build when ARSPAT2_LIB="" (its fits are the
# wrong ones; their time is the old route's). Data seed 11.
lib <- Sys.getenv("ARSPAT2_LIB", "C:/Users/adf44/source/r/wt-arspat2-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
cat("frmtmb from", dirname(dirname(system.file(package = "frmtmb"))), "\n")
k <- 20
g <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(g)
Aj <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- Aj / rowSums(Aj)
M <- solve(crossprod(diag(n) - 0.4 * W))
set.seed(11)
d <- data.frame(x = rnorm(n), h = factor(rep(1:8, length.out = n)))
d$y <- 1 + 0.5 * d$x + rnorm(8)[d$h] +
  as.vector(t(chol(M)) %*% rnorm(n)) * sqrt(4 / rchisq(1, 4))
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
tm <- function(expr) {
  best <- Inf
  for (i in 1:3) {
    t0 <- proc.time()[["elapsed"]]
    force(expr())
    best <- min(best, proc.time()[["elapsed"]] - t0)
  }
  best
}
for (m in c(10, 40)) {
  dm <- d
  set.seed(m)
  miss <- sort(sample(n, m))
  dm$y[miss] <- NA
  ok <- !is.na(dm$y)
  pr <- set_prior("gamma(2, 0.1)", class = "nu")
  t_f <- tm(function() quiet(frm(y | mi() ~ x + (1 | h) + fcor(M),
                                 data = dm, data2 = list(M = M),
                                 family = student(), prior = pr)))
  t_s <- tm(function() quiet(frm(y | mi() ~ x + (1 | h) +
                                   sar(W, type = "error"),
                                 data = dm, data2 = list(W = W),
                                 family = student(), prior = pr)))
  t_o <- tm(function() quiet(frm(y ~ x + (1 | h) + fcor(Mo),
                                 data = dm[ok, ],
                                 data2 = list(Mo = M[ok, ok]),
                                 family = student(), prior = pr)))
  cat(sprintf("n %d m %d: fcor mi %.2f s, sar mi %.2f s, fcor observed rows %.2f s\n",
              n, m, t_f, t_s, t_o))
}
