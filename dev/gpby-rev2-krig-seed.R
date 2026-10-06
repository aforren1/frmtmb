# Reviewer round 2: dev/gpby-rev2-krig.R with the draw seed read from REVSEED.
# Reviewer round 2: the kriging draw's pivoted factor (gp_krig_factor()).
# (1) Analytic: sd^2 w w' (L L' + diag(white)) at idx against
#     gp_krig_cov(), which carries the nugget on its diagonal and on rows
#     at one position.
# (2) Empirical: 40000 draws of gp_krig_draw(), sample covariance against
#     gp_krig_cov(), in units of each entry's Monte Carlo sd.
# (3) A high-rank case: the fitted length scale cut to 0.05 and 400
#     unseen positions spread over [-1, 9]: rank, accuracy and time per
#     factor against a dense Cholesky of gp_krig_cov().
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "\n")
set.seed(5)
d <- data.frame(x = round(runif(60, 0, 6), 1),
                f = factor(rep(c("a", "b"), 30)), w = runif(60, 0.5, 2))
d$y <- sin(d$x) + rnorm(60, 0, 0.3)
krig_of <- function(fit, nd, j = 1L) {
  lp <- fit$frame$linpreds[["y.mu"]]
  ed <- ns$lp_eta_design(fit, lp, nd, FALSE, FALSE)
  Filter(Negate(is.null), lapply(ed$sm_parts, `[[`, "krig"))[[j]]
}
draw_cov <- function(kg) {
  fc <- ns$gp_krig_factor(kg)
  C <- fc$L %*% t(fc$L) + diag(fc$white, length(fc$white))
  list(C = kg$sd2 * C[fc$idx, fc$idx] * outer(kg$w, kg$w),
       rank = ncol(fc$L))
}
check <- function(lab, fit, nd, j = 1L, nsim = 40000) {
  kg <- krig_of(fit, nd, j)
  S <- ns$gp_krig_cov(kg)
  dc <- draw_cov(kg)
  set.seed(as.integer(Sys.getenv("REVSEED", "1")))
  Z <- vapply(seq_len(nsim), function(i) ns$gp_krig_draw(kg),
              numeric(length(kg$rows)))
  Z <- if (is.null(dim(Z))) matrix(Z, 1) else Z
  Se <- stats::cov(t(Z))
  # Monte Carlo sd of a sample covariance entry: sqrt((S_ii S_jj +
  # S_ij^2) / nsim)
  mcsd <- sqrt((outer(diag(S), diag(S)) + S^2) / nsim)
  same <- outer(ns$pos_rowkey(kg$X), ns$pos_rowkey(kg$X), `==`)
  cat(sprintf(paste0("%-34s rows %3d rank %3d | analytic max |C - S| / ",
                     "sd^2 %.2e | empirical max |Se - S| / mcsd %.2f ",
                     "(%d entries) | same-position rows identical in ",
                     "draws: %s\n"),
              lab, length(kg$rows), dc$rank,
              max(abs(dc$C - S)) / kg$sd2, max(abs(Se - S) / mcsd),
              length(S),
              all(vapply(which(rowSums(same) > 1), function(i) {
                k <- which(same[i, ])
                all(Z[k[1], ] == Z[k[2], ] * kg$w[k[1]] / kg$w[k[2]])
              }, NA))))
}
f1 <- frm(bf(y ~ gp(x)), data = d)
check("gp(x) past the data", f1, data.frame(x = seq(6.05, 8, length.out = 12)))
check("gp(x) inside, repeated, near-dup", f1,
      data.frame(x = c(0.55, 0.55, 2.05, 2.05 + 1e-9, 3.33, 4.44, 7)))
fb <- frm(bf(y ~ gp(x, by = f)), data = d)
check("gp(x, by = f) sub-GP 2", fb,
      data.frame(x = c(0.55, 6.5, 7, 0.55, 6.5, 7),
                 f = factor(c("a", "a", "a", "b", "b", "b"))), j = 2L)
fw <- frm(bf(y ~ gp(x, by = w)), data = d)
check("gp(x, by = w), repeated x, w 1 and 2", fw,
      data.frame(x = c(6.5, 6.5, 7.2), w = c(1, 2, 0.5)))
# (3) high rank
fs <- f1
bk <- Filter(function(b) b$covstruct == "gp", fs$frame$re_blocks)[[1]]
fs$estimates$theta[bk$theta_idx[2]] <- log(0.05)
fs$cache <- new.env(parent = emptyenv())
nd <- data.frame(x = seq(-1, 9, length.out = 400))
kg <- krig_of(fs, nd)
S <- ns$gp_krig_cov(kg)
t_f <- system.time(for (i in 1:5) fc <- ns$gp_krig_factor(kg))[["elapsed"]] / 5
t_c <- system.time(for (i in 1:5) ch <- chol(S + diag(1e-14, nrow(S)),
                                            pivot = TRUE))[["elapsed"]] / 5
dc <- draw_cov(kg)
g0 <- gc(reset = TRUE); fc <- ns$gp_krig_factor(kg); g1 <- gc()
cat(sprintf(paste0("high rank: rows %d (unseen %d) rank %d | max |C - S| / ",
                   "sd^2 %.2e | factor %.3f s, dense pivoted chol %.3f s | ",
                   "peak Vcells during factor %.1f MB, L %.1f MB\n"),
            nrow(nd), length(kg$rows), dc$rank, max(abs(dc$C - S)) / kg$sd2,
            t_f, t_c, g1[2, 6] - g0[2, 2],
            as.numeric(object.size(fc$L)) / 2^20))
cat("DONE\n")
