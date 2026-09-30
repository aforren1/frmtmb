.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib", "C:/Users/adf44/source/r/rellib-r3", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(1); n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n)); u <- rlogis(n) + 0.8 * d$x
d$y <- ifelse(runif(n) < 0.3, 0L, 1L + (u > -1) + (u > 0.3) + (u > 1.5))
f <- frm(bf(y ~ x, disc ~ 1 + z), family = hurdle_cumulative(), data = d,
         prior = set_prior("normal(0, 1)", class = "Intercept", dpar = "disc"))
lp <- f$frame$linpreds[[frmtmb:::linpred_key("y", "disc")]]
str(lp[c("par", "idx")]); print(colnames(lp$X))
ent <- frmtmb:::resolve_prior_input(list(frame = f$frame, spec = f$spec), f$prior)$entries
str(ent)
print(summary(f))
