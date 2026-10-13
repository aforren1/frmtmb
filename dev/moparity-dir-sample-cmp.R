# frm_sample() against brms 2.23.0 (dev/moparity-dir-sample.R output):
# per quantity the difference of the means over the joint Monte Carlo
# standard error, and the same for the sds, with divergences and R-hat.
a <- readRDS("dev/moparity-log/dir-sample-lane.rds")
b <- readRDS("dev/moparity-log/dir-sample-brms.rds")
worst <- 0
for (nm in names(a)) {
  A <- a[[nm]]$S
  B <- b[[nm]]$S
  v <- intersect(A$variable, B$variable)
  A <- A[match(v, A$variable), ]
  B <- B[match(v, B$variable), ]
  zm <- (A$mean - B$mean) / sqrt(A$mcse_mean^2 + B$mcse_mean^2)
  zs <- (A$sd - B$sd) / sqrt(A$mcse_sd^2 + B$mcse_sd^2)
  worst <- max(worst, abs(zm), abs(zs))
  cat(sprintf("== %s: divergences frmtmb %d brms %d; only in brms: %s\n",
              nm, a[[nm]]$div, b[[nm]]$div,
              paste(setdiff(b[[nm]]$S$variable, v), collapse = " ")))
  print(data.frame(variable = v, mean_frmtmb = A$mean, mean_brms = B$mean,
                   z_mean = zm, sd_frmtmb = A$sd, sd_brms = B$sd, z_sd = zs,
                   rhat_frmtmb = A$rhat, rhat_brms = B$rhat),
        digits = 4, row.names = FALSE)
}
cat("largest |z| over every mean and sd:", format(worst, digits = 3), "\n")
