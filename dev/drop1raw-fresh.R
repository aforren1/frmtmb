# A fit saved with saveRDS() by dev/drop1raw-survey.R, reloaded in a
# fresh session that has no `d`. Arg: base | lane.
# Log: dev/drop1raw-out/fresh-<arm>.txt
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drop1raw-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
res <- function(lab, expr) {
  r <- tryCatch({
    v <- q(expr)
    paste("ok", if (is.character(v)) paste(v, collapse = " | ") else "")
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 120)))
  cat(sprintf("%-46s %s\n", lab, r))
}
f <- readRDS(sprintf("dev/drop1raw-out/fits-%s.rds", arm))
cat("exists('d'):", exists("d"), "\n")
fit <- f$fit
res("drop1() log", rownames(drop1(fit)))
res("drop1() plain", rownames(drop1(f$fitp)))
res("update(. ~ . - x)", format(nobs(update(fit, . ~ . - x))))
res("update() plain", format(nobs(update(f$fitp, . ~ . - x))))
res("add1(~ . + I(x3^2))", rownames(add1(fit, ~ . + I(x3^2))))
res("step()", deparse1(formula(step(fit, trace = 0))))
res("influence() default data", format(nrow(influence(f$fitr,
                                                      groups = "g")$fixed)))
res("emmeans(~ x)", format(nrow(as.data.frame(emmeans::emmeans(fit, ~ x)))))
res("marginaleffects avg_slopes()", format(nrow(
  marginaleffects::avg_slopes(fit))))
res("conditional_effects()", paste(names(conditional_effects(fit)),
                                   collapse = ","))
res("predict()", format(length(predict(fit)[, 1])))
