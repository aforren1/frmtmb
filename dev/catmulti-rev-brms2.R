# Reviewer: brms 2.23.0 on a two-category categorical with
# posterior_linpred(transform = TRUE), and brms's own dpar order for
# refcat = NA on levels "a b", "1x", "c". Log dev/catmulti-rev-out/brms2.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
lv <- c("a b", "1x", "c")
d <- data.frame(y = factor(rep(lv, 10), levels = lv), x = rnorm(30))
for (rc in list(NULL, NA)) {
  bt <- brms:::validate_formula(bf(y ~ x), data = d, family = categorical(refcat = rc))
  cat("refcat", format(rc), "family$dpars:", bt$family$dpars, "\n")
}
set.seed(3)
d2 <- data.frame(x = rnorm(200))
d2$y <- factor(ifelse(runif(200) < plogis(0.3 + d2$x), "b", "a"))
f <- brm(y ~ x, data = d2, family = categorical(), chains = 1, iter = 400,
         refresh = 0, backend = "rstan", seed = 2)
r <- tryCatch(dim(posterior_linpred(f, transform = TRUE)),
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat("2-cat posterior_linpred(transform = TRUE):", r, "\n")
r <- tryCatch(dim(posterior_linpred(f)), error = function(e) conditionMessage(e))
cat("2-cat posterior_linpred():", r, "\n")
cat("DONE\n")
