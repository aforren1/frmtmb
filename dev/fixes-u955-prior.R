# Lane fixes, punch round, m1: the update of brmsfit-methods:955 under
# brms's fit2 priors (normal(2, 2) on a, normal(0, 3) on b), as brms's
# update() keeps them; the default start, and the recovery start.
#   Rscript dev/fixes-u955-prior.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
d <- brms_fixture_data(2)
d$Trt <- as.numeric(as.character(d$Trt))
fo <- bf(count | weights(AgeSD) ~ a + b, a ~ Age + (1 | ID1 | patient),
         b ~ Age + (1 | ID1 | patient), nl = TRUE)
pr <- c(set_prior("normal(2, 2)", nlpar = "a"),
        set_prior("normal(0, 3)", nlpar = "b"))
fit <- tryCatch(withCallingHandlers(
  frm(fo, data = d, family = Gamma("identity"), prior = pr,
      control = frmtmb_control(verbose = TRUE)),
  message = function(m) {
    if (grepl("restarting|placed", conditionMessage(m))) {
      cat("  [message]", conditionMessage(m))
    }
    invokeRestart("muffleMessage")
  }), error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n"); NULL
  })
if (!is.null(fit)) {
  cat("conv:", fit$opt$convergence, fit$opt$message, "\n")
  print(round(fixef(fit)[, 1:2], 3))
  cat("brms posterior means (dev/fixes-rev-log/u955-brms.txt): a_Intercept",
      "6.883, a_Age 1.603, b_Intercept 10.963, b_Age -0.762\n")
}
# and through update() of the ported fixture, given the same priors
fit2 <- brms_fixture(2)
up <- tryCatch(suppressMessages(
  update(fit2, formula. = bf(count ~ a + b, nl = TRUE), prior = pr)),
  error = function(e) conditionMessage(e))
cat("update(fit2, ..., prior = fit2's priors):",
    if (is.character(up)) up else paste("class", class(up)[1], "conv",
                                        up$opt$convergence), "\n")
