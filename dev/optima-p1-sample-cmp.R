# Lane optima, punch round 1, B1: frm_sample() against brms 2.23.0 per
# simplex weight (dev/optima-p1-sample.R output): the difference of the
# means over the joint Monte Carlo standard error, and the same for the
# sds, with divergences and R-hat on both sides.
a <- readRDS("dev/optima-log/p1-sample-lane.rds")
b <- readRDS("dev/optima-log/p1-sample-brms.rds")
worst <- 0
for (nm in names(a)) {
  A <- a[[nm]]$S
  B <- b[[nm]]$S
  zm <- (A$mean - B$mean) / sqrt(A$mcse_mean^2 + B$mcse_mean^2)
  zs <- (A$sd - B$sd) / sqrt(A$mcse_sd^2 + B$mcse_sd^2)
  worst <- max(worst, abs(zm), abs(zs))
  cat(sprintf("== %s: divergences frmtmb %d brms %d\n", nm, a[[nm]]$div,
              b[[nm]]$div))
  print(data.frame(weight = A$variable, mean_frmtmb = A$mean,
                   mean_brms = B$mean, z_mean = zm, sd_frmtmb = A$sd,
                   sd_brms = B$sd, z_sd = zs, rhat_frmtmb = A$rhat,
                   rhat_brms = B$rhat), digits = 4, row.names = FALSE)
}
cat("largest |z| over every weight's mean and sd:", format(worst, digits = 3),
    "\n")
