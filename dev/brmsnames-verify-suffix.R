## Punch round 1, MAJOR 2: which coefficient brms gives the __1 suffix
## when mu's covariate sigma_z meets sigma's z. Reads the reviewer's C6
## brmsfit (dev/brmsnames-rev-brms.R, data dev/brmsnames-rev-data.R) and
## fits the same data by maximum likelihood on the lane build, so the
## content under each name is compared, not only the names.
##   Rscript dev/brmsnames-verify-suffix.R
source("dev/brmsnames-libs.R")
brmsnames_libs("lane")
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
b <- readRDS("dev/stan-cache/brmsnames-rev-brms-C6.rds")
m <- as.matrix(b$fit)
cat("brms means:\n")
print(round(colMeans(m[, c("b_sigma_z", "b_sigma_z__1",
                           "b_sigma_Intercept")]), 4))
cat("data: b$data, the brmsfit's own\n")
d <- b$data
fit <- q(frm(bf(y ~ sigma_z + (1 | g), sigma ~ z) + gaussian(), data = d))
fe <- fixef(fit)
cat("frmtmb ML: mu sigma_z", round(fe$mu[["sigma_z"]], 4),
    " sigma z", round(fe$sigma[["z"]], 4), "\n")
h <- function(s) q(hypothesis(fit, paste(s, "= 0")))$hypothesis$Estimate
cat("frmtmb hypothesis: sigma_z", round(h("sigma_z"), 4),
    " sigma_z__1", round(h("sigma_z__1"), 4), "\n")
same <- sign(mean(m[, "b_sigma_z"]) - mean(m[, "b_sigma_z__1"])) ==
  sign(h("sigma_z") - h("sigma_z__1"))
cat("suffix on the same coefficient as brms:", same, "\n")
