# Lane fixes, item 2: the frame of the 955 model, and mu at the starts.
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(testthat); library(frmtmb)})
sys.source("tests/testthat/helper-brms-suite.R", envir = environment())
d <- brms_fixture_data(2)
d$Trt <- as.numeric(as.character(d$Trt))
cat("Age range:", range(d$Age), "\n")
fo <- bf(count | weights(AgeSD) ~ a + b, a ~ Age + (1 | ID1 | patient),
         b ~ Age + (1 | ID1 | patient), nl = TRUE)
fr <- frm(fo, data = d, family = Gamma("identity"), dry_run = "frame")
for (k in names(fr$linpreds)) {
  lp <- fr$linpreds[[k]]
  cat("==", k, "\n"); print(names(lp))
  cat(" nl_body:", deparse(lp$nl_body), " nl_pars:", lp$nl_pars,
      " par:", lp$par, " idx:", lp$idx, "\n")
}
cat("mu at the prior-located start a = 2 + 2 Age:",
    range(2 + 2 * d$Age), "\n")
