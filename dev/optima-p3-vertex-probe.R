# Lane optima c1 probe: a simplex at a vertex (weight at 1) on
# ls ~ mo(income) * age, the review's final-check construction.
.libPaths(c("C:/Users/adf44/source/r/wt-optima-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1234)
lev <- c("below_20", "20_to_40", "40_to_100", "greater_100")
income <- factor(sample(lev, 100, TRUE), levels = lev, ordered = TRUE)
d <- data.frame(income = income,
                ls = c(30, 60, 70, 75)[income] + rnorm(100, sd = 7))
d$age <- rnorm(100, mean = 40, sd = 10)
fit <- frm(ls ~ mo(income) * age, data = d)
tab <- frmtmb:::summary_mo_frame(fit, 0.95)
print(tab, digits = 10)
print(attr(tab, "face"))
