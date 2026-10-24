# Probe: how a fit registered through the seam reports itself.
.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
library(testthat)
ns <- asNamespace("frmtmb")
env <- new.env(parent = ns)
sys.source("tests/testthat/helper-registry.R", envir = env)
src <- readLines("tests/testthat/test-covstruct-seam.R")
# the definitions above the first test only
first <- grep("^test_that", src)[1]
eval(parse(text = src[seq_len(first - 1)]), envir = env)
local({
  register_xcs()
  d <- xcs_data()
  fit <- frm(bf(y ~ x + xcs(0 + t | g)), family = gaussian(), data = d)
  print(str(fit$estimates))
  print(VarCorr(fit))
  print(confint_varcorr(fit))
  print(log_sd_theta_index(fit))
  print(try(sigma(fit)))
}, envir = env)
