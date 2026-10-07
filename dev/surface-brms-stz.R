# Lane surface, item 7: what brms 2.23.0 declares for a sum-to-zero
# threshold vector, alone and as a mixture component, and where its
# class "Intercept" prior lands. stancode() only; nothing compiles.
#
#   Rscript dev/surface-brms-stz.R > dev/surface-out/brms-stz.txt
source("dev/surface-env.R")
surface_env("base")
set.seed(3)
dp <- data.frame(x = rnorm(200))
dp$y <- sample(1:5, 200, TRUE)
pr <- brms::set_prior("normal(0, 3)", class = "Intercept")
sc1 <- brms::stancode(y ~ x, data = dp,
                      family = brms::cumulative(threshold = "sum_to_zero"),
                      prior = pr)
cat("## one family, sum_to_zero, prior normal(0, 3) on class Intercept\n")
l <- strsplit(sc1, "\n")[[1]]
cat(grep("Intercept|stz|sum", l, value = TRUE), sep = "\n")
cat("\n## mixture(cumulative(sum_to_zero), cumulative()), prior on mu1\n")
pr2 <- brms::set_prior("normal(0, 3)", class = "Intercept", dpar = "mu1")
sc2 <- brms::stancode(y ~ x, data = dp,
                      family = brms::mixture(
                        brms::cumulative(threshold = "sum_to_zero"),
                        brms::cumulative()), prior = pr2)
l <- strsplit(sc2, "\n")[[1]]
cat(grep("Intercept_mu1|stz|sum", l, value = TRUE), sep = "\n")
cat("\n## default_prior, one family\n")
print(brms::default_prior(y ~ x, data = dp,
                          family = brms::cumulative(threshold =
                                                      "sum_to_zero")))
