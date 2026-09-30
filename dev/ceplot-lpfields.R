# Lane ceplot: the fields of a linear predictor, to find where every
# variable a model reads is stored.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(200), fc = factor(sample(c("a", "b", "c"), 200, TRUE)),
                z = rnorm(200))
d$yo <- factor(sample(1:4, 200, TRUE), ordered = TRUE)
ff <- frm(bf(yo ~ x + cs(fc)) + sratio(), data = d)
lp <- ff$frame$linpreds[[1]]
print(names(lp))
str(lp$cs, max.level = 2)
print(names(ff$spec$responses[[1]]))
print(ff$spec$responses[[1]]$dpars)
