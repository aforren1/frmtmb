suppressMessages(library(frmtmb))
sq <- frmtmb:::frmtmb_links$sqrt
lo <- c(-3, 1); hi <- c(-1, 3)
dm <- sq$mu_eta((lo + hi) / 2)
a <- sq$linkinv(lo); b <- sq$linkinv(hi)
cat("dm:", dm, " (mixed sign:", any(dm < 0) && any(dm > 0), ")\n")
fin <- is.finite(dm)
old_swaps <- any(fin) && all(dm[fin] < 0)
cat("OLD global rule swaps?", old_swaps, "\n")
oa <- a; ob <- b
if (old_swaps) { t <- oa; oa <- ob; ob <- t }
cat("OLD lower:", oa, " upper:", ob, " inverted rows:", sum(oa > ob), "\n")
n <- frmtmb:::ce_band_ends(sq$linkinv, lo, hi, dm)
cat("NEW lower:", n$lower, " upper:", n$upper, " inverted rows:",
    sum(n$lower > n$upper), "\n")
