# Reviewer round 2: high-rank kriging factor cost against the dense
# route it replaced (gp_krig_cov() then chol()), at 100, 400 and 1000
# unseen positions with the length scale cut to 0.05.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(5)
d <- data.frame(x = round(runif(60, 0, 6), 1))
d$y <- sin(d$x) + rnorm(60, 0, 0.3)
f1 <- frm(bf(y ~ gp(x)), data = d)
for (ell in c(1, 0.05)) for (n in c(100, 400, 1000)) {
  fs <- f1
  bk <- Filter(function(b) b$covstruct == "gp", fs$frame$re_blocks)[[1]]
  if (ell != 1) fs$estimates$theta[bk$theta_idx[2]] <- log(ell)
  nd <- data.frame(x = seq(-1, 9, length.out = n) + 1e-3)
  ed <- ns$lp_eta_design(fs, fs$frame$linpreds[["y.mu"]], nd, FALSE, FALSE)
  kg <- Filter(Negate(is.null), lapply(ed$sm_parts, `[[`, "krig"))[[1]]
  k <- 3
  tf <- system.time(for (i in 1:k) fc <- ns$gp_krig_factor(kg))[[3]] / k
  td <- system.time(for (i in 1:k) {
    S <- ns$gp_krig_cov(kg); ch <- chol(S, pivot = TRUE)
  })[[3]] / k
  cat(sprintf("ell %-5s unseen %4d rank %4d | pivoted factor %.3f s | dense cov + chol %.3f s\n",
              format(if (ell == 1) exp(f1$estimates$theta[bk$theta_idx[2]]) else ell, digits = 3),
              length(kg$rows), ncol(fc$L), tf, td))
}
