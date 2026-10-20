# Reviewer: what a refcat = NA ML fit reports for the quantities it does
# determine (category probabilities) and for the df it counts.
# Log dev/catmulti-rev-out/na.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
suppressPackageStartupMessages(library(frmtmb))
source("dev/catmulti-rev-data.R")
dc <- rev_cat_data(23, 3, n = 300)
w <- character(0)
fn <- withCallingHandlers(frm(bf(y ~ x), family = categorical(refcat = NA), data = dc),
  warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
cat("fit warnings:\n"); cat(w, sep = "\n")
fa <- frm(bf(y ~ x), family = categorical(), data = dc)
w <- character(0)
Pn <- withCallingHandlers(fitted(fn),
  warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
Pa <- fitted(fa)
cat("fitted() warnings:", w, sep = "\n")
cat("refcat NA fitted SE finite:", sum(is.finite(Pn[, "Est.Error", ])), "of", length(Pn[, "Est.Error", ]), "\n")
cat("refcat a  fitted SE finite:", sum(is.finite(Pa[, "Est.Error", ])), "of", length(Pa[, "Est.Error", ]), "\n")
cat("head SE NA:", head(Pn[, "Est.Error", 1]), "\nhead SE a:", head(Pa[, "Est.Error", 1]), "\n")
cat("logLik df: NA", attr(logLik(fn), "df"), " ref a", attr(logLik(fa), "df"), "\n")
cat("AIC: NA", AIC(fn), " ref a", AIC(fa), "\n")
cat("anova(fa, fn):\n"); print(tryCatch(anova(fa, fn), error = function(e) conditionMessage(e)))
cat("DONE\n")
