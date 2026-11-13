# univariate y | mi() against the complete-case fit: anova and df.residual
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
set.seed(9)
n <- 120
du <- data.frame(x = rnorm(n), z = rnorm(n))
du$y <- du$x + rnorm(n)
du$y[c(2, 4, 6, 8)] <- NA
u0 <- frm(bf(y | mi() ~ x), data = du)
cc <- frm(bf(y ~ x + z), data = du)
cat(arm, "df.residual mi fit", df.residual(u0), " complete-case fit",
    df.residual(cc), "\n")
cat(arm, "anova(mi fit, complete-case fit):",
    tryCatch({a <- anova(u0, cc); sprintf("ran, Chisq %.4f", a$Chisq[2])},
             error = function(e) substr(conditionMessage(e), 1, 120)), "\n")
