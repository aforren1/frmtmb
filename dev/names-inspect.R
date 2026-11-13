# Lane names: inspect where a fit keeps the order its dpar formulas were
# written in (skew case of dev/names-cases.R). Output on stdout.
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("dev/names-cases.R")
cases <- names_cases(frmtmb::bf, list(
  categorical = categorical(), skew_normal = skew_normal(),
  gaussian = gaussian(), poisson = poisson()))
cs <- cases$skew
fit <- suppressWarnings(frm(cs$f, family = cs$fam, data = cs$d))
str(fit$spec$responses[[1]], max.level = 1)
print(names(fit$spec$responses[[1]]$dpars))
print(names(fit$spec$responses[[1]]$pforms))
print(names(fit$frame$linpreds))
for (bk in fit$frame$re_blocks) {
  cat(bk$term_label, "|", vapply(bk$components, `[[`, "", "lp_key"), "\n")
}
print(names(cs$f$pforms))
print(fit$spec$responses[[1]]$family$dpars)
