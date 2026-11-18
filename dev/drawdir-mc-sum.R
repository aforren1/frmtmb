# Lane drawdir: dev/drawdir-mixcov.R's lane arm against its base arm
# (rellib-r22), per design, split by verdict: fired on both (the summed
# check), fired on the lane only (the per-direction reading), fired on
# neither. Means of the 90% and 98% coverage of the true predictive law,
# the 98% width, the predictive means' error against the truth, and
# |z| (predict()'s mean against fitted() in Monte Carlo SEs).
# Usage: Rscript dev/drawdir-mc-sum.R <design> [<design> ...]
for (des in commandArgs(TRUE)) {
  b <- utils::read.delim(sprintf("dev/drawdir-mc-%s-base.tsv", des))
  l <- utils::read.delim(sprintf("dev/drawdir-mc-%s-lane.tsv", des))
  m <- merge(b, l, by = "seed", suffixes = c(".b", ".l"))
  cat(sprintf("\n== %s: %d seeds with both arms\n", des, nrow(m)))
  grp <- ifelse(m$fire.b, "fired both ", ifelse(m$fire.l, "lane only  ",
                                                "neither    "))
  if (any(m$fire.b & !m$fire.l)) cat("  base fired, lane not:",
                                     m$seed[m$fire.b & !m$fire.l], "\n")
  f <- function(x) if (all(is.na(x))) "    -" else sprintf("%.4f", mean(x))
  for (g in unique(grp)) {
    j <- grp == g
    cat(sprintf(paste("  %s n %2d | cov90 base %s lane %s | cov98 base %s",
                      "lane %s | w98 base %s lane %s | mae base %s lane %s",
                      "plug %s | |z| base %s lane %s\n"),
                g, sum(j), f(m$cov90.b[j]), f(m$cov90.l[j]), f(m$cov98.b[j]),
                f(m$cov98.l[j]), f(m$w98.b[j]), f(m$w98.l[j]),
                f(m$mae_pred.b[j]), f(m$mae_pred.l[j]), f(m$mae_plug.l[j]),
                f(abs(m$z.b[j])), f(abs(m$z.l[j]))))
  }
  j <- grp == "lane only  "
  if (any(j)) {
    cat("  lane-only fits: seed dir k | mae base lane | cov90 base lane\n")
    for (i in which(j)) {
      cat(sprintf("    %2d %8.2f %5s | %.4f %.4f | %s %s\n", m$seed[i],
                  m$dir.l[i], m$k.l[i], m$mae_pred.b[i], m$mae_pred.l[i],
                  format(m$cov90.b[i], digits = 3),
                  format(m$cov90.l[i], digits = 3)))
    }
  }
  j <- grp == "neither    "
  same <- sum(abs(m$mae_pred.b[j] - m$mae_pred.l[j]) <= 1e-9 *
                pmax(1, abs(m$mae_pred.b[j])))
  cat(sprintf("  neither: %d of %d with mae equal to 1e-9 base/lane\n",
              same, sum(j)))
}
