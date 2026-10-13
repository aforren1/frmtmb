# Reviewer: dev/moparity-rev-cox-long.R's two arms, z of each mean and sd
# difference over the joint Monte Carlo error.
a <- readRDS("dev/moparity-rev-log/cox-long-lane.rds")[["5"]]$S
b <- readRDS("dev/moparity-rev-log/cox-long-brms.rds")[["5"]]$S
m <- merge(a, b, by = "variable", suffixes = c("_frmtmb", "_brms"))
m$z_mean <- (m$mean_frmtmb - m$mean_brms) /
  sqrt(m$mcse_mean_frmtmb^2 + m$mcse_mean_brms^2)
m$z_sd <- (m$sd_frmtmb - m$sd_brms) / sqrt(m$mcse_sd_frmtmb^2 + m$mcse_sd_brms^2)
m$sd_ratio <- m$sd_frmtmb / m$sd_brms
print(m[, c("variable", "z_mean", "z_sd", "sd_ratio")], digits = 4, row.names = FALSE)
