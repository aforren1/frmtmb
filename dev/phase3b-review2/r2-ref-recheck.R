# Reviewer 2: ref-points.rds was written at 18:35:17, after the 18:30
# cutoff of the battery loss. Read it back, then recompute on the
# restored library, from scratch, the Rmpfr reference and the lane's
# values at the 40 worst survival rows, the 20 worst defective rows and
# 40 random rows, and compare with the stored ones.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(Rmpfr); library(frmtmb.eam)})
cat("Rmpfr", as.character(packageVersion("Rmpfr")), "from", find.package("Rmpfr"), "\n")
cat("frmtmb.eam from", find.package("frmtmb.eam"), "\n")
source("dev/phase3b-review2/r2-mpfr-lib.R")
ns <- asNamespace("frmtmb.eam")
x <- readRDS("dev/phase3b-review2/ref-points.rds")
cat("read back:", nrow(x), "rows; all reference columns finite:",
    all(is.finite(as.matrix(x[, c("lS", "lF", "lFl", "lFu")]))), "\n")
eS <- abs(x$g_lS - x$lS)
eB <- pmax(ifelse(x$lFl > log(1e-300), abs(x$g_lFl - x$lFl), 0),
           ifelse(x$lFu > log(1e-300), abs(x$g_lFu - x$lFu), 0))
set.seed(5)
pick <- unique(c(order(-eS)[1:40], order(-eB)[1:20], sample(nrow(x), 40)))
cat("rows rechecked:", length(pick), "\n")
out <- t(vapply(pick, function(i) {
  t <- x$t[i]; v <- x$v[i]; a <- x$a[i]; w <- x$w[i]
  tl <- tail_eigen(t, v, a, w); tu <- tail_eigen(t, -v, a, 1 - w)
  lS <- as.numeric(log(tl + tu))
  u <- t / a^2
  Fl <- if (u <= 1) lower_images(t, v, a, w) else p_lower(v, a, w) - tl
  Fu <- if (u <= 1) lower_images(t, -v, a, 1 - w) else p_lower(-v, a, 1 - w) - tu
  c(lS = lS, lFl = as.numeric(log(Fl)), lFu = as.numeric(log(Fu)))
}, numeric(3)))
g <- ns$ddm_rt_lcdf2(x$t[pick], x$v[pick], x$a[pick], x$w[pick])$lS
gl <- ns$ddm_rt_lcdf_b(x$t[pick], x$v[pick], x$a[pick], x$w[pick], 0)
gu <- ns$ddm_rt_lcdf_b(x$t[pick], x$v[pick], x$a[pick], x$w[pick], 1)
cat("stored vs recomputed reference, max |diff| of the log: S",
    max(abs(out[, "lS"] - x$lS[pick])), " Fl", max(abs(out[, "lFl"] - x$lFl[pick])),
    " Fu", max(abs(out[, "lFu"] - x$lFu[pick])), "\n")
cat("stored vs recomputed lane values: identical S", identical(g, x$g_lS[pick]),
    " Fl", identical(gl, x$g_lFl[pick]), " Fu", identical(gu, x$g_lFu[pick]), "\n")
va <- abs(x$v[pick] * x$a[pick])
e <- abs(g - out[, "lS"])
cat("lane S error on the recheck rows by |v| a band:\n")
print(tapply(e, cut(va, c(-1, 1, 24, 72, 120, 400)), max))
