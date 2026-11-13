# brms 2.23.0 on the two fits of lane ppdraw (backlog 3269, 3277):
# posterior_predict() against posterior_epred(), and NA draws.
# Usage: Rscript dev/ppdraw-brms.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")

cat("\n== A. zero_inflated_binomial, no zero inflation in the data\n")
set.seed(7)
n <- 200
dz <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
dz$y <- rbinom(n, dz$N, plogis(0.3 + 0.5 * dz$x))
ba <- brm(y | trials(N) ~ x, family = zero_inflated_binomial(), data = dz,
          chains = 2, iter = 2000, seed = 11, refresh = 0)
zi <- as_draws_df(ba)$zi
cat("zi posterior quantiles (0.5, 0.9, 0.99, 1):",
    format(quantile(zi, c(0.5, 0.9, 0.99, 1)), digits = 4), "\n")
cat("prior on zi:\n"); print(prior_summary(ba))
set.seed(1)
pp <- posterior_predict(ba)
cat("posterior_predict mean", mean(pp), " posterior_epred mean",
    mean(posterior_epred(ba)), " NA share", mean(is.na(pp)), "\n")

cat("\n== B. mixture(cumulative, cumulative)\n")
set.seed(7)
do <- data.frame(x = rnorm(400))
ko <- rbinom(400, 1, 0.3)
do$yo <- cut(do$x + ifelse(ko == 1, 2, 0) + rlogis(400),
             c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
bb <- tryCatch(brm(bf(yo ~ x), family = mixture(cumulative(), cumulative()),
                   data = do, chains = 2, iter = 2000, seed = 11,
                   refresh = 0),
               error = function(e) {
                 cat("brm error:", conditionMessage(e), "\n")
                 NULL
               })
if (!is.null(bb)) {
  print(summary(bb))
  set.seed(13)
  pb <- posterior_predict(bb, newdata = data.frame(x = c(-1, 0, 1)))
  cat("NA share", mean(is.na(pb)), "\n")
  print(table(pb, useNA = "always"))
}
