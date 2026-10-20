# marginaleffects::predictions() on a multinomial and a categorical fit,
# warnings let through, on the lane build. Output: dev/catmulti-out/me.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib", "C:/Users/adf44/source/r/rellib-r9", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
source("dev/catmulti-data.R")
dm <- catmulti_mult_data(12)
fm <- frm(bf(Yc | trials(N) ~ x), family = multinomial(), data = dm)
r <- tryCatch(withCallingHandlers(marginaleffects::predictions(fm),
  warning = function(w) invokeRestart("muffleWarning")),
  error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(r)) r else dim(r))
dc <- catmulti_cat_data(11)
fc <- frm(bf(f ~ x), family = categorical(), data = dc)
r <- tryCatch(marginaleffects::predictions(fc),
  error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(r)) r else head(as.data.frame(r)[, 1:5]))
r <- tryCatch(conditional_effects(fm), error = function(e) paste("ERROR:", conditionMessage(e)))
print(if (is.character(r)) r else lapply(r, function(d) c(dim(d), levels(d$cats__))))
r <- tryCatch(dharma_residuals(fm), error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)
