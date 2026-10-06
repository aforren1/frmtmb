# Punch round 1, B2: how much of the gap between frmtmb's ML residual SD
# of bmi (the two mi() models of brms_missings) and brms's posterior
# mean the degrees-of-freedom factor of ML explains.
#
#   Rscript dev/vigport-pr1-nhanes.R
#
# The missing bmi rows are latent in an mi() model, so the residual SD
# is estimated from the observed rows; the bmi formula has 4
# coefficients (Intercept, age, chl, age:chl).
here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  normalizePath(dirname(sub("^--file=", "", a[1])), winslash = "/")
})
P <- readRDS(file.path(here, "vigport-port-out", "r5", "plausibility.rds"))
data("nhanes", package = "mice")
n <- sum(!is.na(nhanes$bmi))
p <- 4
k <- sqrt((n - p) / n)
cat(sprintf("observed bmi rows %d of %d; sqrt((n - p) / n) = %.4f\n",
            n, nrow(nhanes), k))
for (id in c("brms_missings.8.2", "brms_missings.11.2")) {
  s <- P$detail[[id]]$sd
  s <- s[s$name == "sd_residual____bmi", ]
  adj <- s$frmtmb / k
  cat(sprintf(paste0("%s: frmtmb ML %.4f, divided by the factor %.4f; ",
                     "brms %.4f (sd %.4f); z %.3f before, %.3f after\n"),
              id, s$frmtmb, adj, s$brms, s$brms_sd, s$z,
              (adj - s$brms) / s$brms_sd))
}
