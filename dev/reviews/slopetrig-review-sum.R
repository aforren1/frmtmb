# Review of lane slopetrig: an independent summary of the lane's
# calibration tsv files (dev/slopetrig-log/calib*-{gau,pois}-{base,lane}.tsv).
# Recomputes the ratio over every row from the design itself (regenerates
# x with the lane's own generator, seeds and cells), not from the lane's
# closed-form approximation, and scores short fits against two references:
# (a) the best of the four arms (the lane's), (b) the ref arm alone (the
# centered, standardized fit), each at tolerances 1e-3 and 1e-6 relative.
# Usage: Rscript slopetrig-review-sum.R <logdir>
a <- commandArgs(TRUE)
ld <- if (length(a)) a[1] else "dev/slopetrig-log"
rd <- function(p) utils::read.delim(file.path(ld, p))
B <- rbind(rd("calib-gau-base.tsv"), rd("calib-pois-base.tsv"),
           rd("calib2-gau-base.tsv"), rd("calib2-pois-base.tsv"))
L <- rbind(rd("calib-gau-lane.tsv"), rd("calib-pois-lane.tsv"),
           rd("calib2-gau-lane.tsv"), rd("calib2-pois-lane.tsv"))
cat("rows base", nrow(B), "lane", nrow(L), "\n")
key <- function(d) paste(d$design, d$c0, signif(d$r, 6), d$fam, d$form, d$seed)
stopifnot(identical(key(B), key(L)))
# the F, T and ref arms run the same code on both libraries
cat("F arm identical base/lane:", sum(B$ll_F == L$ll_F), "of", nrow(B),
    "; max |diff|", signif(max(abs(B$ll_F - L$ll_F)), 3), "\n")
cat("T arm identical base/lane:", sum(B$ll_T == L$ll_T), "of", nrow(B), "\n")
cat("ref identical base/lane:", sum(B$ll_ref == L$ll_ref), "of", nrow(B), "\n")
# ratio over every row, regenerated
gen <- function(design, c0, r, seed) {
  set.seed(seed)
  n <- 300
  trt <- rep(c(TRUE, FALSE), length.out = n)
  z <- rnorm(n)
  rnorm(20); rnorm(20)
  switch(design,
    unc = c0 + r * c0 * z, year = c0 + r * c0 * z,
    dose = ifelse(trt, 1 + r * z, 0),
    twonz = ifelse(trt, sample(c(1, 1.01), n, TRUE), 0),
    ulpind = ifelse(trt, sample(c(0.3, 0.1 * 3), n, TRUE), 0),
    binlev = c0 + trt + 0.01 * z)
}
ra <- t(vapply(seq_len(nrow(B)), function(i) {
  x <- gen(B$design[i], B$c0[i], B$r[i], B$seed[i])
  nz <- x[x != 0]
  c(all = sqrt(mean(x^2)) / sd(x), nz = sqrt(mean(nz^2)) / sd(nz),
    sp = sd(nz), rms = sqrt(mean(x^2)))
}, numeric(4)))
cat("regenerated nz ratio matches tsv (rel 1e-9):",
    sum(abs(ra[, "nz"] / B$ratio - 1) < 1e-9), "of", nrow(B), "\n")
B$ratio_all <- ra[, "all"]
B$rms_all <- ra[, "rms"]
best <- pmax(B$ll_ref, B$ll_F, B$ll_T, B$ll_D, L$ll_D, na.rm = TRUE)
sp <- B$spread < 0.05
for (ref in c("best", "ref")) for (tol in c(1e-3, 1e-6)) {
  top <- if (ref == "best") best else B$ll_ref
  shF <- top - B$ll_F > tol
  shDb <- top - B$ll_D > tol
  shDl <- top - L$ll_D > tol
  shT <- top - B$ll_T > tol
  cat(sprintf("\n== reference %s, tol %g: short F %d, T %d, D base %d (silent %d), D lane %d (silent %d)\n",
              ref, tol, sum(shF), sum(shT), sum(shDb), sum(shDb & B$w_D == 0),
              sum(shDl), sum(shDl & L$w_D == 0)))
  if (ref == "best") cat("   ref arm itself short:", sum(best - B$ll_ref > tol),
                         "max", signif(max(best - B$ll_ref), 3), "\n")
  bands <- c(1, 1.8, 2, 2.5, 3.5, 5, 10, 40, Inf)
  for (nm in c("all", "nz")) {
    rr <- if (nm == "all") B$ratio_all else B$ratio
    bb <- cut(rr[!sp], bands, right = FALSE)
    cat("  spread>=0.05, short F by ratio_", nm, " band:\n", sep = "")
    print(rbind(short = tapply(shF[!sp], bb, sum), n = table(bb)))
  }
  lft <- shF & !sp
  cat("  short F left by spread rule alone:", sum(lft), "silent",
      sum(lft & B$w_F == 0), "; max", signif(max((top - B$ll_F)[lft]), 4), "\n")
  eng <- L$as_D
  cat("  lane engages", sum(eng), "; by ratio rule only", sum(eng & !sp),
      "; of those not needed (F not short)", sum(eng & !sp & !shF), "\n")
  cat("  engaged-by-rule check: ratio_all>2|sp vs as_D mismatches",
      sum((B$ratio_all > 2 | sp) != eng), "\n")
}
# the 36 warned-at-max fits
best <- pmax(B$ll_ref, B$ll_F, B$ll_T, B$ll_D, L$ll_D, na.rm = TRUE)
nr <- L$as_D & !sp & best - L$ll_D <= 1e-3 & L$w_D > 0
cat("\nnewly engaged by ratio, not short, warned:", sum(nr), "\n")
print(table(design = B$design[nr], c0 = B$c0[nr]))
cat("ratio_all range of those:", signif(range(B$ratio_all[nr]), 4), "\n")
cat("of them T arm warns too:", sum(B$w_T[nr] > 0), "; ref warns:",
    sum(B$w_ref[nr] > 0), "\n")
cat("base D on those: short", sum(best[nr] - B$ll_D[nr] > 1e-3), "warned",
    sum(B$w_D[nr] > 0), "\n")
# base D: warned and not short
cat("base D not short but warned:", sum(best - B$ll_D <= 1e-3 & B$w_D > 0),
    "; lane D:", sum(best - L$ll_D <= 1e-3 & L$w_D > 0), "\n")
# cost: time ratio D/F where newly engaged and not needed
nn <- L$as_D & !sp & best - B$ll_F <= 1e-3
cat("time D_lane/F where newly engaged and not needed: median",
    signif(median(L$t_D[nn] / pmax(L$t_F[nn], 0.01)), 3), "q90",
    signif(quantile(L$t_D[nn] / pmax(L$t_F[nn], 0.01), 0.9), 3),
    "; total seconds D", sum(L$t_D[nn]), "F", sum(L$t_F[nn]), "\n")
# rms<2 table
cat("\nshort F (best, 1e-3) among spread>=0.05 by rms_all band and ratio_all>2:\n")
shF <- best - B$ll_F > 1e-3
print(table(rms = cut(B$rms_all[!sp], c(0, 2, 5, 20, Inf)),
            ratio_gt2 = B$ratio_all[!sp] > 2, short = shF[!sp]))
