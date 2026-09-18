source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(frmtmb))
set.seed(20260917); n <- 80
d <- data.frame(x = rnorm(n), g = factor(rep_len(paste0("g", 1:8), n)),
                h = factor(rep_len(paste0("h", 1:4), n)))
d$y <- rnorm(n, 1 + 0.5 * d$x, 1)
f2 <- frm(bf(y ~ x + (1 | g) + (1 | h)), family = gaussian(), data = d)
nd <- d[1:3, ]
a <- tryCatch(predict(f2, newdata = nd, re_formula = NULL),
              error = function(e) conditionMessage(e))
b <- tryCatch(predict(f2, newdata = nd, re_formula = ~ (1 | h)),
              error = function(e) conditionMessage(e))
cc <- tryCatch(predict(f2, newdata = nd, re_formula = NA),
               error = function(e) conditionMessage(e))
cat("NULL     :", if (is.character(a)) a else paste(round(a, 8), collapse=","), "\n")
cat("~(1|h)   :", if (is.character(b)) b else paste(round(b, 8), collapse=","), "\n")
cat("NA       :", if (is.character(cc)) cc else paste(round(cc, 8), collapse=","), "\n")
cat("identical(NULL, ~(1|h)):", identical(a, b), "\n")
