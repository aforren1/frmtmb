# Lane selimit2 fix round: what a linear predictor entry and the frame hold
.libPaths(c("C:/Users/adf44/source/r/wt-selimit2-lib", "C:/Users/adf44/source/r/rellib-r21", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(21)
d <- data.frame(f = factor(rep(c("a", "b"), each = 150)))
d$y <- ifelse(d$f == "a", 0L, rpois(300, 4))
fit <- suppressWarnings(frm(bf(y ~ 1, zi ~ f), family = zero_inflated_poisson(), data = d))
lp <- fit$frame$linpreds
print(names(lp)); str(lp[[2]], max.level = 1, give.attr = FALSE)
print(names(fit$frame))
str(fit$frame$y, max.level = 1); str(fit$frame$responses, max.level = 2)
print(fit$spec$responses[[1]]$family$family)
li <- fit$cache$se_limit; print(li$jw)
