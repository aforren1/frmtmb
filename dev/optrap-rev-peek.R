.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib", "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(101)
d3 <- data.frame(x = rnorm(300), g1 = factor("only"))
d3$y <- rnorm(300, 1 + 2 * d3$x, 1)
fr <- frm(bf(y ~ x + (1 | g1)), family = gaussian(), data = d3, dry_run = "frame")
bk <- fr$re_blocks[[1]]
str(bk[setdiff(names(bk), c("Z", "components"))], max.level = 1)
