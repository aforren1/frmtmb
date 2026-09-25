# Edge cases of simulate(newdata = ): an NA covariate, and a grouping
# column absent from newdata while a same-named object exists in the
# formula environment (the global environment here).
source("dev/simnewdata-review/rv-prelude.R")
set.seed(7)
n <- 120
d <- data.frame(x = rnorm(n), g = factor(rep(1:12, 10)))
d$y <- 0.5 + 0.8 * d$x + rnorm(12)[d$g] + rnorm(n, 0, 0.5)
fit <- frm(bf(y ~ x + (1 | g)), data = d)
cat("NA covariate:", try_msg(print(simulate(fit, nsim = 2, seed = 1,
    newdata = data.frame(x = c(NA, 1), g = factor(c("1", "2")))))), "\n")
g <- factor(rep("3", 5))
r <- tryCatch(simulate(fit, nsim = 2, seed = 1, re_formula = NA,
                       newdata = data.frame(x = c(0, 1))),
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("group column absent, global g of length 5:\n"); print(r)
g <- factor(c("3", "4"))
r <- tryCatch(simulate(fit, nsim = 2000, seed = 1, re_formula = NA,
                       newdata = data.frame(x = c(0, 0))),
              error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(r)) cat("global g of length 2:", r, "\n") else
  cat("global g of length 2 (levels 3, 4): answers; cor of the two rows",
      round(cor(unlist(r[1, ]), unlist(r[2, ])), 3), "\n")
rm(g)
cat("group column absent, no global g:", try_msg(simulate(fit, nsim = 2,
    re_formula = NA, newdata = data.frame(x = c(0, 1)))), "\n")
cat("group column absent, NULL, allow_new_levels:", try_msg(simulate(fit,
    nsim = 2, newdata = data.frame(x = c(0, 1)), allow_new_levels = TRUE)),
    "\n")
