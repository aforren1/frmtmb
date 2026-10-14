# Lane optstart, item R1: the units and the mi() placeholders of
# fit_imp2 (brms_missings), and the outer Hessian in natural units
# against the raw-scale one.
#   Rscript dev/optstart-r1-peek.R <lib>
.libPaths(c(commandArgs(TRUE)[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
data("nhanes", package = "mice")
f <- suppressWarnings(suppressMessages(frm(bf(bmi | mi() ~ age * mi(chl)) +
  bf(chl | mi() ~ age) + set_rescor(FALSE), data = nhanes)))
cat("par_units:\n"); print(setNames(f$par_units, names(f$opt$par)))
print(frmtmb:::predictor_units(f$frame))
lp <- f$frame$linpreds[[1]]
cat("colnames:", colnames(lp$X), "\n")
str(lp$mi, max.level = 2)
print(head(as.matrix(lp$X)))
cat("template beta names:", names(f$frame$par_template$beta), "\n")
for (k in names(f$frame$linpreds)) cat(k, "idx", f$frame$linpreds[[k]]$idx,
  "spread", frmtmb:::predictor_col_spread(f$frame$linpreds[[k]]), "\n")
cat("sd(chl, observed)", sd(nhanes$chl, na.rm = TRUE), "\n")
