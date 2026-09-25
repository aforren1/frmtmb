# Reviewer 2, item 2/3: the interval mass ddm_rt_linterval_b() against
# an Rmpfr reference, P(t1 < T <= t2, lower) = tail(t1) - tail(t2),
# with the eigenfunction tail at 1200 bits (r2-mpfr-lib.R). Focus: the
# lower edge inside the blend's region (u0B = 0.2, scale 0.12), where
# the early route is F(t2) (1 - F(t1) / F(t2)) in double and loses every
# digit once the mass is below 1e-16 of F.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(Rmpfr); library(frmtmb.eam)})
source("dev/phase3b-review2/r2-mpfr-lib.R")
ns <- asNamespace("frmtmb.eam")
ref_int <- function(t1, t2, v, a, w) {
  as.numeric(log(tail_eigen(t1, v, a, w) - tail_eigen(t2, v, a, w)))
}
ref_F <- function(t, v, a, w) as.numeric(log(p_lower(v, a, w) - tail_eigen(t, v, a, w)))
set.seed(99)
pts <- rbind(
  data.frame(v = -5.847174087, a = 2.647920, w = 0.5609502, u1 = 2.133821 / 2.647920^2,
             r = 1.3),
  expand.grid(v = c(-6, -3, -1, 1, 3), a = c(1, 2.5), w = c(0.3, 0.6),
              u1 = c(0.1, 0.15, 0.2, 0.3, 0.45, 0.7), r = c(1.05, 1.3, 2)),
  data.frame(v = runif(150, -8, 8), a = exp(runif(150, log(0.5), log(3))),
             w = runif(150, 0.1, 0.9), u1 = exp(runif(150, log(0.05), log(3))),
             r = 1 + exp(runif(150, log(0.02), log(1)))))
pts$t1 <- pts$u1 * pts$a^2
pts$t2 <- pts$t1 * pts$r
pts$ref <- vapply(seq_len(nrow(pts)), function(i) {
  ref_int(pts$t1[i], pts$t2[i], pts$v[i], pts$a[i], pts$w[i])
}, 0)
pts$refF1 <- vapply(seq_len(nrow(pts)), function(i) {
  ref_F(pts$t1[i], pts$v[i], pts$a[i], pts$w[i])
}, 0)
pts$got <- ns$ddm_rt_linterval_b(pts$t1, pts$t2, pts$v, pts$a, pts$w, 0)
pts$abs_err <- abs(pts$got - pts$ref)
pts$mass_over_F1 <- exp(pts$ref - pts$refF1)
saveRDS(pts, "dev/phase3b-review2/interval.rds")
options(width = 180, digits = 5)
cat("points:", nrow(pts), "\n")
cat("abs error of the log interval mass: max", max(pts$abs_err), " median",
    median(pts$abs_err), "\n")
cat("by mass / F(t1) band:\n")
b <- cut(log10(pts$mass_over_F1), c(-Inf, -16, -12, -8, -4, 0.5))
print(tapply(pts$abs_err, b, max))
print(table(b))
cat("\nworst 12:\n")
print(head(pts[order(-pts$abs_err), c("v", "a", "w", "u1", "r", "ref", "got",
                                      "abs_err", "mass_over_F1")], 12))
