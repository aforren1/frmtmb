# brms 2.23.0 conditional_effects() on the saved multinomial fit, and
# insight::get_predicted() on a frmtmb categorical fit for comparison.
# Output: dev/catmulti-out/brms-ce.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(brms))
fm <- readRDS("dev/catmulti-out/brms-mult-bare.rds")
r <- tryCatch(conditional_effects(fm), error = function(e) conditionMessage(e))
print(if (is.character(r)) r else lapply(r, dim))
r <- tryCatch(conditional_effects(fm, categorical = TRUE),
              error = function(e) conditionMessage(e))
print(if (is.character(r)) r else lapply(r, function(d) c(dim(d), names(d))))
r <- tryCatch(insight::get_predicted(fm), error = function(e) conditionMessage(e))
print(if (is.character(r)) r else c(class(r), length(r), dim(r)))
fc <- readRDS("dev/catmulti-out/brms-cat-refb.rds")
r <- tryCatch(insight::get_predicted(fc), error = function(e) conditionMessage(e))
print(if (is.character(r)) r else c(class(r), length(r), dim(r)))
