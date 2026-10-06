# Punch round 1, m4: how closely the pivoted factor of the kriging draw
# reproduces gp_krig_cov(), and its rank, on a grid past the data, a
# grid inside it, and near-duplicate positions.
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d$y <- 0.5 + sin(d$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d)
lp <- fit$frame$linpreds[["y.mu"]]
grids <- list(past = seq(6.05, 9, length.out = 300),
              inside = seq(0.03, 5.97, length.out = 300),
              neardup = c(6.7, 6.7, 7.4, 7.4 + 1e-9, 7.4 + 1e-6, 2.55))
for (g in names(grids)) {
  nd <- data.frame(x = grids[[g]])
  ed <- ns$lp_eta_design(fit, lp, nd, TRUE, FALSE)
  krig <- Filter(Negate(is.null), lapply(ed$sm_parts, `[[`, "krig"))[[1]]
  S <- ns$gp_krig_cov(krig)
  tt <- system.time(for (i in 1:50) fc <- ns$gp_krig_factor(krig))
  C <- tcrossprod(fc$L) + diag(fc$white, length(fc$white))
  C <- (krig$sd2 * C[fc$idx, fc$idx]) * outer(krig$w, krig$w)
  cat(sprintf(paste0("%-8s rows %3d | rank %3d | max |C - S| / sd^2 ",
                     "%.3e | diag max rel %.3e | %.2f ms per factor\n"),
              g, nrow(S), ncol(fc$L), max(abs(C - S)) / krig$sd2,
              max(abs(diag(C) / diag(S) - 1)), 1000 * tt[["elapsed"]] / 50))
}
