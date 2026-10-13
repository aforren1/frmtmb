.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib", "C:/Users/adf44/source/r/rellib-r7", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1)
d <- data.frame(x = rnorm(100))
d$y <- 0.5 * d$x + ifelse(rbinom(100, 1, 0.4) == 1, rnorm(100, 4, 1), rnorm(100, 0, 1))
f <- frm(bf(y ~ x), family = mixture(gaussian, gaussian), data = d)
fam <- f$frame$spec$responses[[1]]$family
str(fam[setdiff(names(fam), c("loglik","dpars_fn"))], max.level = 1, give.attr = FALSE)
for (k in names(f$frame$linpreds)) { l <- f$frame$linpreds[[k]]; cat(k, l$dpar, l$resp, l$link$name, "\n") }
print(f$par_units); print(f$opt$par)
