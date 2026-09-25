# Rscript dev/fams-brms-pp-check.R > dev/fams-brms-pp-check.txt
#
# Is brms's posterior_predict_hurdle_negbinomial draw the zero-truncated
# NB? Compare its positive-part pmf (Monte Carlo) with the exact one.
set.seed(20260925)
mu <- 1.5
shape <- 0.8
N <- 2e6
t <- -log(1 - runif(N) * (1 - exp(-mu)))
d <- rnbinom(N, mu = mu - t, size = shape) + 1
k <- 1:6
emp <- tabulate(d, 6) / N
ex <- dnbinom(k, mu = mu, size = shape) /
  (1 - dnbinom(0, mu = mu, size = shape))
se <- sqrt(ex * (1 - ex) / N)
print(round(cbind(k, emp, ex, z = (emp - ex) / se), 5))
cat("mean draw", mean(d), "exact truncated mean",
    mu / (1 - dnbinom(0, mu = mu, size = shape)), "\n")
