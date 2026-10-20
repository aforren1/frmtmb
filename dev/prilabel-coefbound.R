# Lane prilabel: does frmtmb accept a bound on one coefficient's prior
# (brms refuses "coef may not be specified when using boundaries")?
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
source("C:/Users/adf44/source/r/frmtmb-wt-prilabel/dev/dblbar-data.R")
suppressPackageStartupMessages(library(frmtmb))
dd <- dblbar_data()
r <- tryCatch({
  fit <- frm(bf(y ~ x + z), family = gaussian(), data = dd,
             prior = set_prior("normal(0, 1)", class = "b", coef = "x",
                               lb = 0.6))
  paste("accepted, b_x =", signif(fixef(fit)["x", 1], 4))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
cat(r, "\n")
