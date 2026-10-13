# dev/moparity-gp-ncp.R: each frmtmb run against brms 2.23.0's, per
# quantity: means and sds over their joint Monte Carlo standard error.
a <- c(readRDS("dev/moparity-log/gp-ncp-base.rds"),
       readRDS("dev/moparity-log/gp-ncp-lane.rds"))
names(a) <- paste(rep(c("base", "lane"), c(length(readRDS(
  "dev/moparity-log/gp-ncp-base.rds")), length(readRDS(
  "dev/moparity-log/gp-ncp-lane.rds")))), names(a))
b <- readRDS("dev/moparity-log/gp-ncp-brms.rds")
for (nm in names(a)) {
  case <- strsplit(nm, " ")[[1]][2]
  B <- b[[paste(case, "brms")]]$S
  A <- a[[nm]]$S
  v <- intersect(A$variable, B$variable)
  A <- A[match(v, A$variable), ]
  B <- B[match(v, B$variable), ]
  zm <- (A$mean - B$mean) / sqrt(A$mcse_mean^2 + B$mcse_mean^2)
  zs <- (A$sd - B$sd) / sqrt(A$mcse_sd^2 + B$mcse_sd^2)
  cat(sprintf("== %s: divergences %d, at max tree depth %d, seconds %.0f | max R-hat %.3f, min bulk ESS %.0f | max |z| mean %.2f sd %.2f\n",
              nm, a[[nm]]$div, a[[nm]]$tdh, a[[nm]]$secs, max(A$rhat),
              min(A$ess_bulk), max(abs(zm)), max(abs(zs))))
}
for (nm in names(b)) {
  B <- b[[nm]]$S
  cat(sprintf("== %s: divergences %d, seconds %.0f | max R-hat %.3f, min bulk ESS %.0f\n",
              nm, b[[nm]]$div, b[[nm]]$secs, max(B$rhat), min(B$ess_bulk)))
}
