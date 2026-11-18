# Lane alapse: how a fixed quantile and weights appear in the frame, and
# whether eval_dpars() gives mu = X b.
#   Rscript dev/alapse-probe2.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/alapse-common.R")
set.seed(5)
d <- data.frame(x = stats::rnorm(100), w = rep(1:2, 50))
d$y <- 1 + 0.5 * d$x + stats::rnorm(100)
f <- suppressWarnings(frm(bf(y | weights(w) ~ x, quantile = 0.5),
                          family = asym_laplace(), data = d))
print(names(f$frame$linpreds))
str(f$frame$aterm_values, max.level = 2)
dp <- ns$eval_dpars(f)
str(dp)
lp <- f$frame$linpreds[["y.mu"]]
print(range(dp$y$mu - drop(lp$X %*% f$estimates$beta[lp$idx])))
print(names(f$opt$par))
print(f$frame$spec$responses$y$fixed_dpars)
print(names(f$spec$responses$y))
print(f$spec$responses$y$dpar_fixed)
print(f$frame[["sparse_x"]])
print(is.null(f$obj$env$random))
print(class(ns$fit_prior_entries(f)))
