# Review of lane mfxse, attack 5: the kriging warning on gp(x) + (1 | g).
#   Rscript dev/review-mfxse-krig.R
# Data seed 8150 (follow.R's design). Rows x = 7, 8, 9, 10, 12 at seen level 2.
lane <- "C:/Users/adf44/source/r/wt-mfxse-lib"
.libPaths(c(lane, "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(marginaleffects)})
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
quiet <- function(e) suppressWarnings(suppressMessages(e))
set.seed(8150)
mg <- 90
dg <- data.frame(x = runif(mg, 0, 6), g = factor(sample(1:6, mg, TRUE)))
dg$y <- sin(dg$x) + rnorm(6, 0, 0.5)[dg$g] + rnorm(mg, sd = 0.3)
for (form in c("y ~ gp(x) + (1 | g)", "y ~ gp(x)")) {
  fg <- quiet(frm(as.formula(form), data = dg))
  nd <- data.frame(x = c(7, 8, 9, 10, 12), g = factor("2", levels = levels(dg$g)))
  ws <- character(0)
  p <- withCallingHandlers(predictions(fg, newdata = nd), warning = function(w) {
    ws <<- c(ws, paste(class(w)[1], substr(conditionMessage(w), 1, 90))); invokeRestart("muffleWarning")})
  cat("\n", form, ": warnings:\n  ", paste(ws, collapse = "\n   "), "\n")
  for (bk in fg$frame$re_blocks) cat("  block", bk$covstruct, "theta_idx", bk$theta_idx, "\n")
  lp <- fg$frame$linpreds[[1]]
  ed <- ns$lp_eta_design(fg, lp, nd, TRUE, FALSE)
  for (sp in ed$sm_parts) cat("  sm_part extra_var:", if (is.null(sp$extra_var)) "NULL" else sprintf("%.3g", sp$extra_var),
                              " block", sp$bk$covstruct, "\n")
}
