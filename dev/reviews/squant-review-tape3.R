# Review of lane squant: the pieces of elf_cdf_log()'s lF on the tape,
# to find which one zeroes the gradient for -690 < t < about -260.
# Run (worktree root):
#   Rscript dev/reviews/squant-review-tape3.R > dev/reviews/squant-review-tape3.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-squant-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
for (nm in c("ad_min", "ad_max", "ad_pos", "log1m_exp_ad", "log_ibeta_two",
             "cdf_rows")) assign(nm, get(nm, ns))
a <- 0.125
b <- 0.125
piece <- function(t, k) {
  tc <- ad_max(ad_min(t, 690), -690)
  dl <- t - tc
  x <- 1 / (1 + exp(-tc))
  xc <- 1 / (1 + exp(tc))
  lp <- log_ibeta_two(x, xc, a, b)
  lf0 <- lp[["lF"]]
  ls0 <- lp[["lS"]]
  lf_low <- lf0 + a * ad_min(dl, 0)
  ls_high <- ls0 - b * ad_pos(dl)
  cap <- -1e-300
  t1 <- log1m_exp_ad(ad_min(ls_high, cap))
  t2 <- log1m_exp_ad(ad_min(ls0, cap))
  switch(k, lf0 = lf0, ls0 = ls0, lf_low = lf_low, t1 = t1, t2 = t2,
         t1mt2 = t1 - t2, lf = lf_low + t1 - t2)
}
for (t in c(-250, -270, -280, -290, -300)) {
  out <- vapply(c("lf0", "ls0", "lf_low", "t1", "t2", "t1mt2", "lf"),
                function(k) {
                  tp <- RTMB::MakeTape(function(tt) piece(tt, k), t)
                  c(tp(t), as.numeric(tp$jacobian(t)))
                }, c(0, 0))
  cat(sprintf("t %g\n", t))
  print(signif(out, 6))
}
