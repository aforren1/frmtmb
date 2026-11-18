# Lane alapse: why ald_cens_sep() does not fire on the all-censored fit.
lib_arg <- "lane"
source("dev/alapse-common.R")
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)),
                g4 = factor(rep(1:20, each = 4)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$ya <- 1 + 0.5 * d$x + rnorm(n)
d$ca <- "right"
d$yac <- rep(min(d$ya), n)
f <- suppressWarnings(frm(yac | cens(ca) ~ x, family = asym_laplace(),
                          data = d))
str(f$frame$aterm_values)
print(ns$fit_prior_density(f))
print(f$opt$objective)
print(str(ns$eval_dpars(f, b = NULL)))
print(ns$fit_outer_box(f))
print(f$frame$linpreds[["yac.sigma"]][c("par", "idx")])
print(f$frame$betad_fixed_idx)
print(names(f$spec$responses))
print(tryCatch(ns$ald_cens_sep(f), error = function(e) e))
