# Lane alapse: what an asym_laplace fit carries (names, frame slots).
#   Rscript dev/alapse-probe1.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/alapse-common.R")
set.seed(5)
d <- data.frame(x = stats::rnorm(100), z = stats::rnorm(100))
d$y <- 1 + 0.5 * d$x + stats::rnorm(100)
r <- cap(frm(bf(y ~ x, sigma ~ z, quantile ~ z), family = asym_laplace(),
             data = d))
f <- r$v
print(r$wc)
print(f$opt[c("par", "objective", "convergence", "message")])
print(names(f$frame))
print(names(f$frame$linpreds))
str(f$frame$linpreds[[1]], max.level = 1)
print(names(f$frame$par_template))
print(ns$outer_par_names(f))
str(ns$outer_par_map(f))
print(names(f$frame$aterm_values[[1]]))
print(dim(f$frame$linpreds[[1]]$X))
print(fixef(f))
print(f$par_units)
H <- f$obj$he(f$opt$par)
print(round(H, 4))
