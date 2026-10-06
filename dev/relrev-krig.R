# Reviewer, item 4: gp_krig_factor()'s dense fallback at its boundary.
# On dev/rel068-krig-rank.R's fit (seed 5), a 400-position grid whose
# lengthscale is tuned so that the lane's pivoted-only factor has rank
# 99, 100, 101 and 102 (n / 4 = 100): which route the release takes,
# whether the pivoted route is bitwise the lane's, the law against
# gp_krig_cov(), and two calls (and two seeded draws) for determinism.
#   Rscript dev/relrev-krig.R > dev/relrev-log/krig.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
lane_src <- "C:/Users/adf44/source/r/frmtmb-wt-gpby/R/predict.R"
ex <- parse(lane_src, keep.source = FALSE)
def <- Filter(function(e) {
  is.call(e) && identical(as.character(e[[2]]), "gp_krig_factor")
}, as.list(ex))[[1]]
old <- eval(def[[3]], ns); environment(old) <- ns
new <- ns$gp_krig_factor
set.seed(5)
d <- data.frame(x = round(runif(60, 0, 6), 1))
d$y <- sin(d$x) + rnorm(60, 0, 0.3)
f1 <- frm(bf(y ~ gp(x)), data = d)
bk <- Filter(function(b) b$covstruct == "gp", f1$frame$re_blocks)[[1]]
krig_at <- function(ell, n) {
  fs <- f1
  fs$estimates$theta[bk$theta_idx[2]] <- log(ell)
  nd <- data.frame(x = seq(-1, 9, length.out = n) + 1e-3)
  ed <- ns$lp_eta_design(fs, fs$frame$linpreds[["y.mu"]], nd, FALSE, FALSE)
  Filter(Negate(is.null), lapply(ed$sm_parts, `[[`, "krig"))[[1]]
}
law_err <- function(fc, kg) {
  S <- ns$gp_krig_cov(kg)
  M <- tcrossprod(fc$L) + diag(fc$white, nrow(fc$L))
  M <- kg$sd2 * M[fc$idx, fc$idx] * outer(kg$w, kg$w)
  max(abs(M - S)) / kg$sd2
}
n <- 400
rank_old <- function(ell) ncol(old(krig_at(ell, n))$L)
# bisection on log(ell) for each target rank (rank falls as ell grows)
find_ell <- function(target) {
  lo <- log(0.01); hi <- log(5)
  for (i in 1:60) {
    mid <- (lo + hi) / 2
    r <- rank_old(exp(mid))
    if (r == target) return(exp(mid))
    if (r > target) lo <- mid else hi <- mid
  }
  NA
}
cat(sprintf("%-10s %-8s %-8s %-7s %-10s %-10s %-10s %-6s %-6s\n", "ell", "rank_ln",
            "rank_new", "route", "pivot==ln", "law new", "law lane", "det_f", "det_draw"))
for (tg in c(99, 100, 101, 102)) {
  ell <- find_ell(tg)
  if (is.na(ell)) { cat("rank", tg, "not hit\n"); next }
  kg <- krig_at(ell, n)
  fo <- old(kg); fn <- new(kg); fn2 <- new(kg)
  same <- identical(fo, fn)
  set.seed(1); a <- ns$gp_krig_draw(kg); set.seed(1); b <- ns$gp_krig_draw(kg)
  cat(sprintf("%-10.6g %-8d %-8d %-7s %-10s %-10.2g %-10.2g %-6s %-6s\n", ell,
              ncol(fo$L), ncol(fn$L), if (same) "pivot" else "dense", same,
              law_err(fn, kg), law_err(fo, kg), identical(fn, fn2), identical(a, b)))
}
# the smallest grids, where n / 4 < 1
for (m in 1:5) {
  kg <- krig_at(0.05, m)
  fo <- old(kg); fn <- new(kg)
  cat(sprintf("n=%d positions %d: rank lane %d new %d, route %s, law new %.2g lane %.2g\n",
              m, nrow(fn$L), ncol(fo$L), ncol(fn$L),
              if (identical(fo, fn)) "pivot" else "dense", law_err(fn, kg), law_err(fo, kg)))
}
