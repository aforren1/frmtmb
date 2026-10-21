.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-lib", "C:/Users/adf44/source/r/rellib-r10", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(RTMB))
ns <- asNamespace("frmtmb")
for (z in c(-38, -30, -12, -6, -3, -1, -0.2)) for (a in c(-50, -5, -0.5, 0.5, 3, 20, 300)) {
  v1 <- ns$owen_g_log(z, a); v2 <- ns$owen_g_log(-z, -a)
  if (!is.finite(v1) || !is.finite(v2)) cat("z", z, "a", a, v1, v2, "\n")
}
