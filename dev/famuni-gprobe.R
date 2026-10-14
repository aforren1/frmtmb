.libPaths(c("C:/Users/adf44/source/r/wt-famuni-lib", "C:/Users/adf44/source/r/rellib-r8", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(1); d <- data.frame(x = rnorm(80), g = gl(8, 10)); d$y <- 1 + d$x + (rexp(80)^(-0.1) - 1) / 0.1
fit <- frm(bf(y ~ x + (1|g), sigma ~ x), family = gen_extreme_value(), data = d)
print(names(fit$opt$par)); print(outer_par_names(fit)); print(fit$estimates$betad)
lp <- fit$frame$linpreds[["y.xi"]]; print(lp$idx); print(fit$frame$betad_fixed_idx)
print(drop(fit$obj$gr(fit$opt$par)))
