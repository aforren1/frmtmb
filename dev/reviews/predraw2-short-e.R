# Short check of predraw2's fix round, item 2: the review's 1-D exact test
# (dev/reviews/predraw2-review-a.R jac) over 40 seeds, to tell Monte
# Carlo error from bias. Per seed: the self-normalized importance
# estimate of P(t < 0) from r$lw, and the share after draw_resample().
# Usage (from the lane tree): Rscript dev/reviews/predraw2-short-e.R
local({
.libPaths(c("C:/Users/adf44/source/r/wt-predraw2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
for (a in c(0.3, 0.1)) {
  obj <- list(fn = function(p) {
    t <- p[1]
    if (t >= 0) t^2 / 2 else (t / a)^2 / 2
  }, env = new.env())
  object <- list(obj = obj, opt = list(par = 0))
  law <- list(lo = list(a * (1:4)), hi = list(1), s = 1, Lc = NULL,
              edge = integer(0), narrowed = 1L)
  L <- matrix(1, 1, 1)
  N <- 20000
  is <- rs_sh <- ks <- numeric(0)
  for (sd in 1:40) {
    set.seed(sd)
    Z <- matrix(stats::rnorm(N), 1, N)
    D <- ns$draw_edge_apply(law, L, Z)
    dw <- list(D = D, Z = Z, L = L, law = law, v0 = 0, held = integer(0),
               off = 0L)
    r <- ns$draw_logratio(object, dw, D, Z)
    w <- exp(r$lw - max(r$lw))
    is <- c(is, sum(w * (D[1, ] < 0)) / sum(w))
    set.seed(1000 + sd)
    rs <- ns$draw_resample(r$lw, r$stratum)
    rs_sh <- c(rs_sh, mean(D[1, rs$idx] < 0))
    ks <- c(ks, rs$k)
  }
  ex <- a / (1 + a)
  cat(sprintf("a %.2f exact %.4f | IS mean %.5f (se %.5f, z %.2f) | resampled mean %.5f (se %.5f, z %.2f) | per-seed sd IS %.5f resampled %.5f | k %s\n",
              a, ex, mean(is), stats::sd(is) / sqrt(40),
              (mean(is) - ex) / (stats::sd(is) / sqrt(40)), mean(rs_sh),
              stats::sd(rs_sh) / sqrt(40),
              (mean(rs_sh) - ex) / (stats::sd(rs_sh) / sqrt(40)),
              stats::sd(is), stats::sd(rs_sh),
              paste(unique(ks), collapse = ",")))
}
})
