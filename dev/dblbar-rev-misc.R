# Reviewer, lane dblbar: an lkj prior on a || block, anova() against the
# correlated term, confint() rows, AIC parameter count. Data seed
# 20261013 (dev/dblbar-data.R). Usage: Rscript ... base|lane
# Output: dev/dblbar-rev-log/misc-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- if (identical(arm, "base")) character(0) else
  "C:/Users/adf44/source/r/wt-dblbar-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(frmtmb))
sink(file.path(wt, sprintf("dev/dblbar-rev-log/misc-%s.txt", arm)),
     split = TRUE)
cat("lib:", find.package("frmtmb"), "\n")
dd <- dblbar_data()
tc <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
A <- frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd)
U <- frm(bf(y ~ f + (1 + f | g)), family = gaussian(), data = dd)
cat("lkj prior on || block:", substr(tc({
  frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd,
      prior = set_prior("lkj(2)", class = "cor", group = "g"))
  "accepted"
}), 1, 160), "\n")
cat("lkj class-wide prior on || fit:", substr(tc({
  frm(bf(y ~ f + (1 + f || g)), family = gaussian(), data = dd,
      prior = set_prior("lkj(2)", class = "cor"))
  "accepted"
}), 1, 160), "\n")
a <- tc(anova(A, U))
print(a)
cat("nobs-par: A", attr(logLik(A), "df"), "U", attr(logLik(U), "df"), "\n")
ci <- tc(suppressWarnings(confint(A)))
print(if (is.character(ci)) ci else rownames(ci))
cv <- tc(suppressWarnings(confint_varcorr(A)))
print(cv)
pr <- tc(suppressWarnings(confint(A, method = "profile",
                                  parm = "sd_g__fb")))
print(pr)
sink()
