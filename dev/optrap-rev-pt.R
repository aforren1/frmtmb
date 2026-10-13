# Reviewer of lane optrap: does par_template() show the unit start of the
# group sds (NEWS)? sleepstudy in ms (unit 64).
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
pt <- par_template(Reaction ~ Days + (Days | Subject), data = sleepstudy)
print(pt$theta)
f <- frm(Reaction ~ Days + (Days | Subject), data = sleepstudy,
         dry_run = "objective")
print(f$obj$par)
