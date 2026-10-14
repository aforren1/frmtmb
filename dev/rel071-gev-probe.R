# bin-2P row log_lik:349 (brms's tests.log_lik.R, gen_extreme_value):
# what frmtmb's density and brms's dgen_extreme_value() give at a row
# inside and at a row outside the support (y = 0 below mu - sigma / xi).
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
fam <- gen_extreme_value()
dp <- list(mu = c(1, -1), sigma = c(0.1, 1), xi = c(0.3, 0.3))
cat("rows: y = 0 with mu = 1, sigma = 0.1, xi = 0.3 (outside: support y > 0.667);",
    "y = 0 with mu = -1, sigma = 1, xi = 0.3 (inside)\n")
cat("frmtmb lpdf:", format(fam$lpdf(c(0, 0), dp, list()), digits = 6), "\n")
cat("brms dgen_extreme_value:", format(suppressWarnings(brms::dgen_extreme_value(
  c(0, 0), mu = dp$mu, sigma = dp$sigma, xi = dp$xi, log = TRUE)), digits = 6), "\n")
