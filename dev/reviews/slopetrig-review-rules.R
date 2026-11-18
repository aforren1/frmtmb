# Review of lane slopetrig, check 7: the three candidate rules on the
# lane's grids, counted against fits newly engaged (beyond the spread
# rule): short left, silent left, newly engaged, of them not needed
# (autoscale = FALSE within 1e-3 of the best arm), and the time ratio
# D/F there. Usage: Rscript slopetrig-review-rules.R
rd <- function(p) utils::read.delim(file.path("dev/slopetrig-log", p))
B <- rbind(rd("calib-gau-base.tsv"), rd("calib-pois-base.tsv"),
           rd("calib2-gau-base.tsv"), rd("calib2-pois-base.tsv"))
L <- rbind(rd("calib-gau-lane.tsv"), rd("calib-pois-lane.tsv"),
           rd("calib2-gau-lane.tsv"), rd("calib2-pois-lane.tsv"))
best <- pmax(B$ll_ref, B$ll_F, B$ll_T, B$ll_D, L$ll_D)
shF <- best - B$ll_F > 1e-3
sp <- B$spread < 0.05
zr <- B$design %in% c("dose", "twonz", "ulpind")
m_nz <- sqrt(pmax(B$rms^2 - B$spread^2, 0))
ra <- ifelse(zr, sqrt(B$rms^2 / 2) / sqrt(B$rms^2 / 2 - (m_nz / 2)^2), B$ratio)
rms_all <- ifelse(zr, sqrt(B$rms^2 / 2), B$rms)
tr <- L$t_T / pmax(B$t_F, 0.01)
rules <- list("ratio > 2" = ra > 2, "ratio > 2 & rms > 2 (nonzero rows)" = ra > 2 & B$rms > 2,
              "ratio > 3" = ra > 3, "ratio > 5" = ra > 5, "ratio > 10" = ra > 10)
for (nm in names(rules)) {
  e <- !sp & rules[[nm]]
  left <- shF & !sp & !rules[[nm]]
  cat(sprintf("%-36s newly engaged %4d, not needed %4d (%.0f%%), median T/F time there %.2f, q90 %.2f ; short left %3d (silent %3d, max %.3g)\n",
              nm, sum(e), sum(e & !shF), 100 * mean(!shF[e]),
              median(tr[e & !shF]), quantile(tr[e & !shF], 0.9), sum(left),
              sum(left & B$w_F == 0), if (any(left)) max((best - B$ll_F)[left]) else 0))
}
cat("rms in x's units: scaling x by 1/1000 (years to millennia) moves rms > 2 but not the ratio\n")
