# Lane postfix: frm_multiple() with a fixed-effect level missing from the
# first imputation, on one build.  Rscript dev/postfix-repro-levels.R [lib]
a <- commandArgs(TRUE)
.libPaths(c(if (length(a)) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
set.seed(8)
n <- 80
x <- rnorm(n)
y <- rnorm(n, 1 + 0.5 * x, 1)
imf <- lapply(1:3, function(i) data.frame(y = y, x = x))
imf[[1]]$f <- factor(rep(c("a", "b"), 40))
imf[[2]]$f <- factor(rep(c("a", "b", "c"), length.out = 80))
imf[[3]]$f <- imf[[2]]$f
r <- tryCatch(frm_multiple(bf(y ~ x + f), family = gaussian(), data = imf),
              error = function(e) paste("ERROR", conditionMessage(e)))
print(if (is.character(r)) r else class(r))
