# Reviewer, claim 5 follow-up: cat() on categorical, cat(x = 4) in brms.
# Log: dev/aterms2-rev-log-05b-cat.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(505)
n <- 150
d <- data.frame(x = rnorm(n))
d$o <- factor(cut(d$x + rnorm(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE))
d$yi <- as.integer(d$o)
r <- function(e) tryCatch({suppressWarnings(e); "ok"},
                          error = function(e) conditionMessage(e))
cat("frm(yi ~ x, categorical):", r(frm(yi ~ x, data = d, family = categorical())), "\n")
cat("frm(o ~ x, categorical):", r(frm(o ~ x, data = d, family = categorical())), "\n")
cat("frm(o | cat(4) ~ x, categorical):", r(frm(o | cat(4) ~ x, data = d, family = categorical())), "\n")
cat("frm(o | thres(3) ~ x, categorical):", r(frm(o | thres(3) ~ x, data = d, family = categorical())), "\n")
cat("brms cat(x = 4):", r(brms::standata(brms::bf(yi | cat(x = 4) ~ x), d, family = brms::cumulative())), "\n")
