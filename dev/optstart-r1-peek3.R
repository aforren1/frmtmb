# Lane optstart, item R1: the units of fit_imp2 after the mi() fix
.libPaths(c("C:/Users/adf44/source/r/wt-optstart-lib", "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
data("nhanes", package = "mice")
f <- frm(bf(bmi | mi() ~ age * mi(chl)) + bf(chl | mi() ~ age) +
  set_rescor(FALSE), data = nhanes)
print(setNames(f$par_units, names(f$opt$par)))
print(names(par_template(f)$beta))
for (k in names(f$frame$linpreds)) print(frmtmb:::predictor_col_spread(f$frame$linpreds[[k]], f$frame))
str(f$frame$y)
str(f$frame$mi_map)
