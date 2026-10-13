# brms 2.23.0's update(prior = ): what prior the refit carries.
# Usage: Rscript dev/moparity-brms-update.R > dev/moparity-log/brms-update.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(11)
dd <- data.frame(x = rnorm(60), z = rnorm(60))
dd$y <- rnorm(60, 1 + 0.5 * dd$x, 1)
pr0 <- c(prior(normal(0, 0.05), class = b, coef = x),
         prior(normal(0, 10), class = Intercept))
fit <- brm(y ~ x + z, data = dd, prior = pr0, chains = 1, iter = 400,
           refresh = 0, seed = 1, backend = "rstan")
show <- function(lab, f) {
  cat("\n==", lab, "\n")
  ps <- prior_summary(f)
  print(ps[ps$source == "user", c("prior", "class", "coef", "source")])
  cat("b_x mean", round(mean(as_draws_df(f)$b_x), 4), "\n")
}
show("stored", fit)
f1 <- update(fit, prior = prior(exponential(2), class = sigma), refresh = 0,
             seed = 1)
show("update(prior = sigma prior)", f1)
f2 <- update(fit, prior = prior(normal(0, 1), class = b, coef = x),
             refresh = 0, seed = 1)
show("update(prior = new b_x prior)", f2)
f3 <- tryCatch(update(fit, prior = prior(normal(0, 1), class = b,
                                          coef = nosuch),
                      refresh = 0, seed = 1),
               error = function(e) e)
if (inherits(f3, "error")) cat("\n== nosuch coef: ERROR", conditionMessage(f3),
                               "\n") else show("update(prior = nosuch coef)", f3)
f4 <- update(fit, formula. = ~ . - x, prior = prior(exponential(2),
                                                     class = sigma),
             refresh = 0, seed = 1)
cat("\n== update(~ . - x, prior = sigma prior)\n")
ps <- prior_summary(f4)
print(ps[ps$source == "user", c("prior", "class", "coef", "source")])
