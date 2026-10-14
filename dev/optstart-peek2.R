# Lane optstart: the fields of a random-effect block and of Z
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(50), g = factor(rep(1:10, 5)), h = factor(rep(1:5, 10)))
d$y <- 1 + 2 * d$x + rnorm(10)[d$g] + rnorm(50)
f <- suppressWarnings(frm(y ~ x + (1 | g) + (1 + x | h), data = d))
fr <- f$frame
for (bk in fr$re_blocks) { cat("---\n"); str(bk[setdiff(names(bk), "components")], max.level = 1) }
str(fr$linpreds[[1]]$Z)
