# Exact shapes behind the failures in brmssuite-spotcheck2.R, so the
# audit can quote what frmtmb returns and not only that it differed.
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
load("dev/brms-suite/brms/data/epilepsy.rda")
fit <- frm(count ~ zBase * Trt + (1 | patient), epilepsy, family = poisson())

cat("## variables(fit)\n"); print(variables(fit))
cat("\n## names of hypothesis() result\n")
h <- hypothesis(fit, "zBase > Trt1")
print(names(h)); print(class(h$hypothesis)); print(h$hypothesis)
cat("\n## posterior_summary(fit)\n")
print(tryCatch(posterior_summary(fit), error = function(e) conditionMessage(e)))
cat("\n## fitted(): are brms arguments swallowed?\n")
a <- fitted(fit)
b <- fitted(fit, re_formula = NA)
cc <- fitted(fit, scale = "linear")
d <- fitted(fit, dpar = "inv")
e <- fitted(fit, this_argument_does_not_exist = 42)
cat(sprintf("fitted(fit)[1]                  = %.6f\n", a[1]))
cat(sprintf("fitted(fit, re_formula = NA)[1] = %.6f  identical: %s\n",
            b[1], identical(a, b)))
cat(sprintf("fitted(fit, scale = 'linear')[1]= %.6f  identical: %s\n",
            cc[1], identical(a, cc)))
cat(sprintf("fitted(fit, dpar = 'inv')[1]    = %.6f  identical: %s\n",
            d[1], identical(a, d)))
cat(sprintf("fitted(fit, nonsense = 42)[1]   = %.6f  identical: %s\n",
            e[1], identical(a, e)))
cat("\n## names(summary(fit)) and names(VarCorr(fit))\n")
print(names(summary(fit))); print(names(VarCorr(fit)))
cat("\n## ngrps\n"); print(ngrps(fit))
cat("\n## class(family(fit))\n"); print(class(family(fit)))
